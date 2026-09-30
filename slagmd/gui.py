"""The window. Six tabs, one background monitor thread, nothing that blocks the interface.

Design rules, because a GUI that freezes is worse than a command line:
  * Every action that touches the disk, the scheduler or a trajectory runs on a worker
    thread. The only thing the main thread does is draw.
  * Workers never touch a widget. They put messages on a queue and the main thread drains
    it on a timer. Tkinter is not thread safe and breaking that rule produces crashes that
    look random.
  * Every number shown next to a result is shown with the tool's verdict on whether the
    number is allowed to exist. There is no place in this interface that prints a diffusion
    coefficient on its own.
"""
from __future__ import annotations

import queue
import threading
import traceback
from dataclasses import asdict
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from . import analysis as an
from . import scheduler as sched
from .composition import VERBEECK_2025_TABLE5, Composition
from .forcefield import ForceField
from .inputgen import EqBlock, Protocol
from .project import Project

NAVY = "#0F2A44"
TEAL = "#14607A"
BAND = "#E8F0F4"
GOOD = "#1B7F3B"
WARN = "#B2600A"
BAD = "#A32020"

STATUS_COLOUR = {
    "COMPLETED": GOOD, "RUNNING": TEAL, "PENDING": WARN, "NOT_SUBMITTED": "#555",
    "FAILED": BAD, "BLOCKED": BAD, "TIMEOUT": BAD, "OUT_OF_MEMORY": BAD,
    "CANCELLED": "#777", "UNKNOWN": WARN,
}


class App(tk.Tk):
    def __init__(self, project_root: Path | None = None):
        super().__init__()
        self.title("slagmd  molecular dynamics of molten oxide slags")
        self.geometry("1180x780")
        self.minsize(960, 640)
        self.msgq: queue.Queue = queue.Queue()
        self.project: Project | None = None
        self.monitor_on = tk.BooleanVar(value=False)
        self._monitor_thread: threading.Thread | None = None

        self._style()
        self._header()
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True, padx=8, pady=(0, 6))
        self._tab_project()
        self._tab_cluster()
        self._tab_runs()
        self._tab_queue()
        self._tab_analysis()
        self._tab_log()
        self._statusbar()

        self.after(200, self._drain)
        if project_root:
            self._open_project(Path(project_root))

    # ------------------------------------------------------------------ chrome
    def _style(self) -> None:
        s = ttk.Style(self)
        try:
            s.theme_use("clam")
        except tk.TclError:
            pass
        s.configure("TNotebook.Tab", padding=(14, 7))
        s.configure("Head.TLabel", foreground=NAVY,
                    font=("DejaVu Sans", 15, "bold"))
        s.configure("Sub.TLabel", foreground="#556674", font=("DejaVu Sans", 9))
        s.configure("Sect.TLabel", foreground=TEAL,
                    font=("DejaVu Sans", 10, "bold"))
        s.configure("Treeview", rowheight=23)
        s.configure("Treeview.Heading", font=("DejaVu Sans", 9, "bold"))

    def _header(self) -> None:
        bar = tk.Frame(self, bg=BAND)
        bar.pack(fill="x")
        inner = tk.Frame(bar, bg=BAND)
        inner.pack(fill="x", padx=12, pady=8)
        ttk.Label(inner, text="slagmd", style="Head.TLabel",
                  background=BAND).pack(anchor="w")
        ttk.Label(inner, background=BAND, style="Sub.TLabel",
                  text="build, submit, watch, repair and analyse molecular dynamics of "
                       "CaO-Al2O3-SiO2 melts.  Protocol and force fields from the "
                       "published Verbeeck and Bellemans scripts.").pack(anchor="w")
        self.proj_label = ttk.Label(inner, background=BAND, style="Sub.TLabel",
                                    text="no project open")
        self.proj_label.pack(anchor="w", pady=(4, 0))

    def _statusbar(self) -> None:
        bar = tk.Frame(self, bg="#EFEFEF")
        bar.pack(fill="x", side="bottom")
        self.status = ttk.Label(bar, text="ready", background="#EFEFEF",
                                style="Sub.TLabel")
        self.status.pack(side="left", padx=10, pady=3)
        self.busy = ttk.Progressbar(bar, mode="indeterminate", length=130)
        self.busy.pack(side="right", padx=10, pady=3)

    # ------------------------------------------------------------ tab: project
    def _tab_project(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="1  Project")
        top = ttk.Frame(f)
        top.pack(fill="x", padx=14, pady=12)
        ttk.Label(top, text="Project folder", style="Sect.TLabel").pack(anchor="w")
        ttk.Label(top, style="Sub.TLabel",
                  text="One folder holds every run, its inputs, its logs and the state "
                       "file. Closing this window loses nothing: the state is on disk "
                       "after every change, so reopening the folder resumes exactly "
                       "where it stopped.").pack(anchor="w", pady=(0, 6))
        row = ttk.Frame(top)
        row.pack(fill="x")
        self.root_var = tk.StringVar()
        ttk.Entry(row, textvariable=self.root_var, width=74).pack(side="left")
        ttk.Button(row, text="Browse", command=self._browse).pack(side="left", padx=6)
        ttk.Button(row, text="Open or create",
                   command=lambda: self._open_project(Path(self.root_var.get()))
                   ).pack(side="left")

        env = ttk.LabelFrame(f, text="  what this machine has  ")
        env.pack(fill="x", padx=14, pady=10)
        probe = sched.probe_cluster()
        txt = [f"scheduler detected     {probe['backend']}",
               f"LAMMPS on PATH         "
               f"{probe.get('lammps_exe', 'NO') if probe['lammps_on_path'] else 'NO'}",
               f"packmol on PATH        {'yes' if probe['packmol_on_path'] else 'no'}",
               f"cores here             {probe['cores_here']}"]
        if not probe["lammps_on_path"]:
            txt.append("")
            txt.append("LAMMPS is not installed here. That is fine for building inputs "
                       "and for analysing trajectories")
            txt.append("brought back from a cluster. To actually run, either set the "
                       "backend to slurm or pbs on a")
            txt.append("login node, or point the executable at a LAMMPS build.")
        ttk.Label(env, text="\n".join(txt), style="Sub.TLabel",
                  font=("DejaVu Sans Mono", 9), justify="left").pack(anchor="w",
                                                                    padx=10, pady=8)

        ff = ttk.LabelFrame(f, text="  force fields available, and what they were fitted to  ")
        ff.pack(fill="both", expand=True, padx=14, pady=(0, 12))
        box = tk.Text(ff, height=12, wrap="word", font=("DejaVu Sans", 9),
                      relief="flat", background="#FCFCFC")
        box.pack(fill="both", expand=True, padx=10, pady=8)
        for nm in ForceField.available():
            o = ForceField.load(nm)
            probs = o.check_self_consistency()
            box.insert("end", f"{nm}\n", "b")
            box.insert("end", f"  {o.label}\n"
                              f"  form {o.form},  charges "
                              f"{', '.join(f'{k} {v:+g}' for k, v in o.charges.items())}\n"
                              f"  source {o.source}\n"
                              f"  self consistency: "
                              f"{'PASSES' if not probs else '; '.join(probs)}\n\n")
        box.tag_config("b", font=("DejaVu Sans", 10, "bold"), foreground=NAVY)
        box.configure(state="disabled")

    def _browse(self) -> None:
        d = filedialog.askdirectory(title="choose or create a project folder")
        if d:
            self.root_var.set(d)

    def _open_project(self, root: Path) -> None:
        try:
            self.project = Project(root)
        except Exception as exc:                       # noqa: BLE001
            messagebox.showerror("could not open the project", str(exc))
            return
        self.root_var.set(str(root))
        self.proj_label.configure(
            text=f"project  {root}      {len(self.project.runs)} runs")
        self._cluster_to_form()
        self._refresh_queue()
        self._refresh_analysis_runs()
        self.log(f"opened project {root}, {len(self.project.runs)} runs")

    # ------------------------------------------------------------ tab: cluster
    def _tab_cluster(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="2  Cluster")
        ttk.Label(f, text="Where the jobs run", style="Sect.TLabel").pack(
            anchor="w", padx=14, pady=(12, 0))
        ttk.Label(f, style="Sub.TLabel", justify="left",
                  text="The defaults are the resources the published study used: two "
                       "nodes of 128 cores, about ten hours for the whole workflow.\n"
                       "Their own release contains no submission script at all, so this "
                       "tab is the part of the work that did not exist before.").pack(
            anchor="w", padx=14, pady=(0, 8))

        self.cvars: dict[str, tk.Variable] = {}
        grid = ttk.Frame(f)
        grid.pack(fill="x", padx=14)
        fields = [
            ("backend", "scheduler", ["local", "slurm", "pbs"]),
            ("partition", "partition or queue", None),
            ("account", "account or project", None),
            ("nodes", "nodes", None),
            ("tasks_per_node", "tasks per node", None),
            ("walltime", "walltime hh:mm:ss", None),
            ("memory_per_node", "memory per node, blank for default", None),
            ("lammps_exe", "LAMMPS executable", None),
            ("mpi_launcher", "launcher", ["srun", "mpirun", "mpiexec", ""]),
            ("email", "email for failure notices", None),
        ]
        for i, (key, label, choices) in enumerate(fields):
            r, c = divmod(i, 2)
            cell = ttk.Frame(grid)
            cell.grid(row=r, column=c, sticky="w", padx=(0, 26), pady=3)
            ttk.Label(cell, text=label, width=30, style="Sub.TLabel").pack(side="left")
            var = tk.StringVar()
            self.cvars[key] = var
            if choices is not None:
                ttk.Combobox(cell, textvariable=var, values=choices, width=18,
                             state="readonly").pack(side="left")
            else:
                ttk.Entry(cell, textvariable=var, width=21).pack(side="left")

        for key, label in [("modules", "modules to load, one per line"),
                           ("pre_commands", "extra shell lines before LAMMPS"),
                           ("extra_directives", "extra scheduler directives")]:
            ttk.Label(f, text=label, style="Sub.TLabel").pack(anchor="w", padx=14,
                                                              pady=(8, 0))
            t = tk.Text(f, height=3, width=100, font=("DejaVu Sans Mono", 9))
            t.pack(anchor="w", padx=14)
            self.cvars[key] = t                                     # type: ignore[assignment]

        self.cluster_warn = ttk.Label(f, text="", style="Sub.TLabel", foreground=BAD,
                                      wraplength=1000, justify="left")
        self.cluster_warn.pack(anchor="w", padx=14, pady=(8, 0))

        row = ttk.Frame(f)
        row.pack(anchor="w", padx=14, pady=10)
        ttk.Button(row, text="Save", command=self._cluster_from_form).pack(side="left")
        ttk.Button(row, text="Check for problems",
                   command=self._check_cluster).pack(side="left", padx=6)
        ttk.Button(row, text="Preview the submit script",
                   command=self._preview_submit).pack(side="left")

    def _cluster_to_form(self) -> None:
        if not self.project:
            return
        c = self.project.cluster
        for k, v in self.cvars.items():
            val = getattr(c, k, "")
            if isinstance(v, tk.Text):
                v.delete("1.0", "end")
                v.insert("1.0", "\n".join(val or []))
            else:
                v.set(str(val))

    def _cluster_from_form(self) -> None:
        if not self.project:
            return messagebox.showwarning("no project", "open a project first")
        c = self.project.cluster
        try:
            for k, v in self.cvars.items():
                if isinstance(v, tk.Text):
                    setattr(c, k, [l for l in v.get("1.0", "end").splitlines()
                                   if l.strip()])
                elif k in ("nodes", "tasks_per_node"):
                    setattr(c, k, int(v.get() or 1))
                else:
                    setattr(c, k, v.get())
        except ValueError as exc:
            return messagebox.showerror("bad value", str(exc))
        self.project.save()
        self._check_cluster()
        self.log("cluster settings saved")

    def _check_cluster(self) -> None:
        if not self.project:
            return
        probs = sched.config_problems(self.project.cluster)
        self.cluster_warn.configure(
            text=("\n".join("PROBLEM: " + p for p in probs) if probs
                  else "no problems found in these settings"),
            foreground=BAD if probs else GOOD)

    def _preview_submit(self) -> None:
        if not self.project:
            return
        txt = sched.submit_script(self.project.cluster, "example_01_hot",
                                 self.project.root / "runs" / "example",
                                 "in.01_hot", "log.01_hot")
        self._show_text("the submit script that will be generated", txt)

    # --------------------------------------------------------------- tab: runs
    def _tab_runs(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="3  Compositions")
        ttk.Label(f, text="Add a run", style="Sect.TLabel").pack(anchor="w", padx=14,
                                                                  pady=(12, 0))
        ttk.Label(f, style="Sub.TLabel", justify="left",
                  text="One run is one composition at one temperature with one force "
                       "field. Mole fractions are resolved to exact integer atom counts, "
                       "and oxygen is counted from the cations rather than rounded, so "
                       "the box is exactly charge neutral.").pack(anchor="w", padx=14,
                                                                  pady=(0, 8))
        g = ttk.Frame(f)
        g.pack(fill="x", padx=14)
        self.rvars: dict[str, tk.Variable] = {}
        spec = [("name", "run name", "S1_1600C", 18),
                ("comp", "composition, oxide then mole fraction",
                 "CaO 0.453 Al2O3 0.125 SiO2 0.422", 40),
                ("temp_C", "temperature in C", "1600", 8),
                ("atoms", "target atoms", "5000", 8),
                ("density", "initial density g/cm3", "2.6", 8),
                ("seed", "random seed", "458165", 10)]
        for i, (key, label, default, width) in enumerate(spec):
            cell = ttk.Frame(g)
            cell.grid(row=i // 2, column=i % 2, sticky="w", padx=(0, 26), pady=3)
            ttk.Label(cell, text=label, width=34, style="Sub.TLabel").pack(side="left")
            var = tk.StringVar(value=default)
            self.rvars[key] = var
            ttk.Entry(cell, textvariable=var, width=width).pack(side="left")

        row = ttk.Frame(f)
        row.pack(anchor="w", padx=14, pady=6)
        ttk.Label(row, text="force field", style="Sub.TLabel").pack(side="left")
        self.ff_var = tk.StringVar(value="bouhadja")
        ttk.Combobox(row, textvariable=self.ff_var, values=ForceField.available(),
                     width=12, state="readonly").pack(side="left", padx=(6, 20))
        self.visc_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(row, text="sample the stress tensor, so viscosity can be computed",
                        variable=self.visc_var).pack(side="left")
        self.nve_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row, text="add a short NVE energy conservation check",
                        variable=self.nve_var).pack(side="left", padx=(14, 0))

        row2 = ttk.Frame(f)
        row2.pack(anchor="w", padx=14, pady=8)
        ttk.Button(row2, text="Check this composition",
                   command=self._check_comp).pack(side="left")
        ttk.Button(row2, text="Add run", command=self._add_run).pack(side="left", padx=6)
        ttk.Button(row2, text="Add the ten published benchmark compositions",
                   command=self._add_benchmark).pack(side="left", padx=6)
        ttk.Button(row2, text="Prepare every unprepared run",
                   command=self._prepare_all).pack(side="left", padx=6)

        self.comp_out = tk.Text(f, height=20, wrap="word",
                                font=("DejaVu Sans Mono", 9), background="#FCFCFC")
        self.comp_out.pack(fill="both", expand=True, padx=14, pady=(4, 12))

    def _parse_run_form(self):
        from .composition import parse_composition_line
        mf = parse_composition_line(self.rvars["comp"].get())
        return (self.rvars["name"].get().strip(), mf,
                float(self.rvars["temp_C"].get()) + 273.15,
                int(self.rvars["atoms"].get()), float(self.rvars["density"].get()),
                int(self.rvars["seed"].get()))

    def _check_comp(self) -> None:
        self.comp_out.delete("1.0", "end")
        try:
            _n, mf, T, atoms, dens, _s = self._parse_run_form()
            c = Composition(mf, atoms)
            ff = ForceField.load(self.ff_var.get())
            ff.check_ready(c.species)
            out = [c.summary(), ""]
            ok, scale = c.charges_are_scaled_formal(ff.charges)
            out.append(f"net box charge under {ff.name}: "
                       f"{c.net_charge(ff.charges):+.3e}")
            out.append(f"charges are a uniform {scale:.4f} of the formal charges: {ok}")
            out.append(f"cubic box at {dens} g/cm3: {c.box_length_A(dens):.2f} A")
            out.append(f"molar mass of the box: {c.molar_mass_g():.1f} g")
            proto = Protocol(target_T=T)
            if self.visc_var.get():
                proto.stress_every_steps = 5
            out.append("")
            out.append(f"protocol: {proto.total_ns():.2f} ns total, "
                       f"{proto.equilibrate_ns():.0f} ns equilibration in "
                       f"{len(proto.eq_blocks)} blocks, {proto.produce_ns:.0f} ns "
                       f"production")
            out.append(f"estimated cost: {proto.est_core_hours(c.n_atoms):.0f} core hours")
            for w in proto.warnings():
                out.append(f"WARNING: {w}")
            self.comp_out.insert("end", "\n".join(out))
        except Exception as exc:                        # noqa: BLE001
            self.comp_out.insert("end", f"REFUSED\n\n{exc}")

    def _add_run(self) -> None:
        if not self.project:
            return messagebox.showwarning("no project", "open a project first")
        try:
            name, mf, T, atoms, dens, seed = self._parse_run_form()
            self.project.add_run(name, mf, T, self.ff_var.get(), atoms, dens, seed,
                                 nve_check=self.nve_var.get(),
                                 stress_for_viscosity=self.visc_var.get())
        except Exception as exc:                        # noqa: BLE001
            return messagebox.showerror("could not add the run", str(exc))
        self.log(f"added run {name}")
        self._refresh_queue()

    def _add_benchmark(self) -> None:
        if not self.project:
            return messagebox.showwarning("no project", "open a project first")
        added = self.project.add_verbeeck_benchmark(
            self.ff_var.get(), (1600,), int(self.rvars["atoms"].get()))
        self.comp_out.delete("1.0", "end")
        self.comp_out.insert("end",
            "Queued the ten compositions of the published benchmark table at 1600 C, "
            "where the published temperature range allows it.\n\n"
            "This exists so the tool can be checked against a study whose answers are in "
            "print. If it cannot reproduce their self diffusion coefficients, the tool is "
            "wrong, and that is worth finding out before using it on anything new.\n\n"
            + "\n".join(f"  {a}" for a in added)
            + f"\n\n{len(added)} runs added.")
        self._refresh_queue()

    def _prepare_all(self) -> None:
        if not self.project:
            return
        todo = [r for r in self.project.runs.values() if not r.prepared]
        if not todo:
            return messagebox.showinfo("nothing to do", "every run is already prepared")
        self._work(self._prepare_worker, todo)

    def _prepare_worker(self, todo) -> None:
        for r in todo:
            self.msgq.put(("status", f"preparing {r.name}"))
            try:
                notes = self.project.prepare(r)                 # type: ignore[union-attr]
                self.msgq.put(("comp", f"--- {r.name}\n" + "\n".join(notes) + "\n"))
            except Exception as exc:                            # noqa: BLE001
                self.msgq.put(("comp", f"--- {r.name}\nFAILED: {exc}\n"))
        self.msgq.put(("refresh", None))

    # -------------------------------------------------------------- tab: queue
    def _tab_queue(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="4  Queue")
        bar = ttk.Frame(f)
        bar.pack(fill="x", padx=14, pady=10)
        ttk.Checkbutton(bar, text="watch and repair automatically",
                        variable=self.monitor_on,
                        command=self._toggle_monitor).pack(side="left")
        ttk.Label(bar, text="  jobs at once", style="Sub.TLabel").pack(side="left",
                                                                       padx=(16, 4))
        self.conc = tk.IntVar(value=4)
        ttk.Spinbox(bar, from_=1, to=40, textvariable=self.conc,
                    width=4).pack(side="left")
        ttk.Label(bar, text="  poll every", style="Sub.TLabel").pack(side="left",
                                                                     padx=(16, 4))
        self.poll_s = tk.IntVar(value=60)
        ttk.Spinbox(bar, from_=5, to=1800, increment=5, textvariable=self.poll_s,
                    width=5).pack(side="left")
        ttk.Label(bar, text="seconds", style="Sub.TLabel").pack(side="left", padx=4)
        ttk.Button(bar, text="Poll now", command=self._poll_once).pack(side="left",
                                                                       padx=(18, 4))
        ttk.Button(bar, text="Submit next", command=self._fill_once).pack(side="left")
        ttk.Button(bar, text="Write report", command=self._write_report).pack(
            side="left", padx=4)

        cols = ("run", "composition", "T_C", "ff", "progress", "stage", "status",
                "job", "attempts", "note")
        self.tree = ttk.Treeview(f, columns=cols, show="headings", height=16)
        widths = {"run": 150, "composition": 230, "T_C": 46, "ff": 78, "progress": 62,
                  "stage": 100, "status": 96, "job": 90, "attempts": 62, "note": 300}
        for c in cols:
            self.tree.heading(c, text=c.replace("_", " "))
            self.tree.column(c, width=widths[c], anchor="w")
        self.tree.pack(fill="both", expand=True, padx=14)
        for st, col in STATUS_COLOUR.items():
            self.tree.tag_configure(st, foreground=col)
        self.tree.bind("<Double-1>", lambda e: self._run_detail())

        row = ttk.Frame(f)
        row.pack(anchor="w", padx=14, pady=8)
        for label, cmd in [("Details and history", self._run_detail),
                           ("Read the diagnosis", self._show_diagnosis),
                           ("Unblock and retry", self._reset_selected),
                           ("Cancel", self._cancel_selected)]:
            ttk.Button(row, text=label, command=cmd).pack(side="left", padx=(0, 6))

    def _selected_run(self):
        sel = self.tree.selection()
        if not sel or not self.project:
            messagebox.showinfo("nothing selected", "select a run in the table first")
            return None
        return self.project.runs.get(self.tree.item(sel[0], "values")[0])

    def _refresh_queue(self) -> None:
        if not self.project:
            return
        self.tree.delete(*self.tree.get_children())
        for row in self.project.summary_rows():
            self.tree.insert("", "end", values=tuple(row[c] for c in
                                                     self.tree["columns"]),
                             tags=(row["status"],))
        n = self.project.n_active()
        blocked = sum(1 for r in self.project.runs.values() if r.is_blocked())
        done = sum(1 for r in self.project.runs.values() if r.is_done())
        self.proj_label.configure(
            text=f"project  {self.project.root}      {len(self.project.runs)} runs, "
                 f"{done} finished, {n} in flight, {blocked} blocked")

    def _toggle_monitor(self) -> None:
        if not self.project:
            self.monitor_on.set(False)
            return messagebox.showwarning("no project", "open a project first")
        if self.monitor_on.get():
            self._monitor_thread = threading.Thread(target=self._monitor_loop,
                                                    daemon=True)
            self._monitor_thread.start()
            self.log("monitor started: it will submit, watch, diagnose and resubmit")
        else:
            self.log("monitor stopped. Nothing already submitted is affected.")

    def _monitor_loop(self) -> None:
        import time
        while self.monitor_on.get() and self.project:
            try:
                for m in self.project.fill_queue(self.conc.get()):
                    self.msgq.put(("log", m))
                for m in self.project.poll_all():
                    self.msgq.put(("log", m))
                self.project.write_report()
                self.msgq.put(("refresh", None))
            except Exception:                            # noqa: BLE001
                self.msgq.put(("log", "monitor error:\n" + traceback.format_exc()))
            for _ in range(max(1, self.poll_s.get())):
                if not self.monitor_on.get():
                    return
                time.sleep(1)

    def _poll_once(self) -> None:
        if self.project:
            self._work(lambda: [self.msgq.put(("log", m))
                                for m in self.project.poll_all()])   # type: ignore

    def _fill_once(self) -> None:
        if self.project:
            self._work(lambda: [self.msgq.put(("log", m))
                                for m in self.project.fill_queue(self.conc.get())])

    def _write_report(self) -> None:
        if self.project:
            p = self.project.write_report()
            self._show_text(f"report  {p}", p.read_text())

    def _run_detail(self) -> None:
        r = self._selected_run()
        if not r or not self.project:
            return
        out = [f"run           {r.name}",
               f"composition   " + "  ".join(f"{k} {v}" for k, v in
                                             sorted(r.mole_fractions.items())),
               f"temperature   {r.temperature_K:.2f} K  "
               f"({r.temperature_K-273.15:.0f} C)",
               f"force field   {r.force_field}",
               f"atoms         {r.target_atoms} target",
               f"folder        {self.project.run_dir(r.name)}",
               f"prepared      {r.prepared}", ""]
        if r.protocol_overrides:
            out += ["protocol overrides applied, usually by an automatic repair:",
                    f"  {r.protocol_overrides}", ""]
        if r.cluster_overrides:
            out += ["cluster overrides applied:", f"  {r.cluster_overrides}", ""]
        for s in r.stage_list():
            st = r.stages[s]
            out.append(f"--- {s}   {st.status}   job {st.job_id or '-'}   "
                       f"attempts {st.attempts}")
            out += [f"      {h}" for h in st.history]
        if r.notes:
            out += ["", "notes from preparation:"] + [f"  {n}" for n in r.notes]
        self._show_text(f"{r.name}", "\n".join(out))

    def _show_diagnosis(self) -> None:
        r = self._selected_run()
        if not r or not self.project:
            return
        d = sorted(self.project.run_dir(r.name).glob("diagnosis.*.txt"))
        if not d:
            return messagebox.showinfo("no diagnosis", "nothing has failed in this run")
        self._show_text(d[-1].name, d[-1].read_text())

    def _reset_selected(self) -> None:
        r = self._selected_run()
        if not r or not self.project:
            return
        s = r.current_stage()
        if not s:
            return messagebox.showinfo("nothing to do", "this run has finished")
        self.project.reset_stage(r, s)
        self.log(f"{r.name}/{s} reset, it will be submitted again")
        self._refresh_queue()

    def _cancel_selected(self) -> None:
        r = self._selected_run()
        if not r or not self.project:
            return
        for m in self.project.cancel_run(r):
            self.log(m)
        self._refresh_queue()

    # ----------------------------------------------------------- tab: analysis
    def _tab_analysis(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="5  Analysis")
        bar = ttk.Frame(f)
        bar.pack(fill="x", padx=14, pady=10)
        ttk.Label(bar, text="run", style="Sub.TLabel").pack(side="left")
        self.an_run = tk.StringVar()
        self.an_combo = ttk.Combobox(bar, textvariable=self.an_run, width=26,
                                     state="readonly")
        self.an_combo.pack(side="left", padx=(6, 16))
        ttk.Button(bar, text="Refresh list",
                   command=self._refresh_analysis_runs).pack(side="left")
        ttk.Label(bar, text="  or a trajectory file", style="Sub.TLabel").pack(
            side="left", padx=(16, 4))
        self.an_file = tk.StringVar()
        ttk.Entry(bar, textvariable=self.an_file, width=34).pack(side="left")
        ttk.Button(bar, text="Browse",
                   command=lambda: self.an_file.set(
                       filedialog.askopenfilename() or "")).pack(side="left", padx=4)

        bar2 = ttk.Frame(f)
        bar2.pack(fill="x", padx=14, pady=(0, 8))
        for label, cmd in [
            ("Everything", lambda: self._analyse("all")),
            ("Structure: RDF and coordination", lambda: self._analyse("rdf")),
            ("Speciation: Qn and oxygen types", lambda: self._analyse("spec")),
            ("Diffusion with error bars", lambda: self._analyse("msd")),
            ("Conductivity, both routes", lambda: self._analyse("cond")),
            ("Viscosity", lambda: self._analyse("visc")),
            ("Was it equilibrated?", lambda: self._analyse("equil")),
        ]:
            ttk.Button(bar2, text=label, command=cmd).pack(side="left", padx=(0, 5))

        self.an_out = tk.Text(f, wrap="word", font=("DejaVu Sans Mono", 9),
                              background="#FCFCFC")
        self.an_out.pack(fill="both", expand=True, padx=14, pady=(0, 6))
        self.an_out.tag_config("h", font=("DejaVu Sans Mono", 10, "bold"),
                               foreground=NAVY)
        self.an_out.tag_config("good", foreground=GOOD)
        self.an_out.tag_config("bad", foreground=BAD)
        row = ttk.Frame(f)
        row.pack(anchor="w", padx=14, pady=(0, 10))
        ttk.Button(row, text="Save these results beside the run",
                   command=self._save_analysis).pack(side="left")
        ttk.Button(row, text="Make the figures",
                   command=lambda: self._analyse("figs")).pack(side="left", padx=6)

    def _refresh_analysis_runs(self) -> None:
        if not self.project:
            return
        names = [r.name for r in self.project.runs.values()]
        self.an_combo.configure(values=names)
        if names and not self.an_run.get():
            self.an_run.set(names[0])

    def _traj_path(self) -> Path | None:
        if self.an_file.get():
            return Path(self.an_file.get())
        if not self.project or not self.an_run.get():
            return None
        p = self.project.run_dir(self.an_run.get()) / "prod.lammpstrj"
        return p if p.exists() else None

    def _analyse(self, what: str) -> None:
        p = self._traj_path()
        if p is None:
            return messagebox.showinfo(
                "no trajectory",
                "There is no prod.lammpstrj for that run yet. Either the production "
                "stage has not finished, or browse to a trajectory file brought back "
                "from a cluster.")
        self._work(self._analyse_worker, what, p)

    def _analyse_worker(self, what: str, path: Path) -> None:
        put = self.msgq.put
        try:
            put(("an_clear", None))
            put(("status", f"reading {path.name}"))
            T = 1873.0
            charges = None
            if self.project and self.an_run.get() in (self.project.runs or {}):
                r = self.project.runs[self.an_run.get()]
                T = r.temperature_K
                charges = ForceField.load(r.force_field).charges
            traj = an.read_lammps_dump(path)
            put(("an", (f"trajectory {path.name}\n", "h")))
            put(("an", (f"  {traj.n_frames} frames, {traj.n_atoms} atoms, "
                        f"{traj.times_ps[-1]:.1f} ps, species "
                        f"{', '.join(traj.species_list())}\n"
                        f"  coordinates are "
                        f"{'unwrapped, correct for MSD' if traj.unwrapped else 'WRAPPED, they will be unwrapped here'}\n"
                        f"  mean box {traj.mean_box()[0]:.2f} A, "
                        f"volume {traj.volume_A3():.0f} A^3\n"
                        f"  temperature taken as {T:.0f} K\n\n", None)))

            if what in ("all", "rdf", "figs"):
                put(("status", "radial distribution functions"))
                rdf = an.compute_rdf(traj, rmax=10.0, nbins=400)
                put(("an", ("structure\n", "h")))
                put(("an", ("  pair      first peak A   first min A   coordination\n",
                            None)))
                for k in sorted(rdf.g):
                    put(("an", (f"  {k:9s} {rdf.first_peak[k]:11.3f} "
                                f"{rdf.first_min[k]:13.3f} "
                                f"{rdf.cn_at_first_min[k]:14.2f}\n", None)))
                for w in rdf.warnings:
                    put(("an", (f"  WARNING: {w}\n", "bad")))
                put(("an", ("\n", None)))

            if what in ("all", "spec"):
                put(("status", "oxygen speciation and Qn"))
                sp = an.compute_speciation(traj)
                put(("an", ("speciation\n", "h")))
                put(("an", (f"  cutoffs measured from this trajectory: "
                            f"{', '.join(f'{k} {v:.2f} A' for k, v in sp.cutoffs.items())}\n",
                            None)))
                put(("an", (f"  bridging      {100*sp.frac_bridging:6.2f} %\n"
                            f"  non bridging  {100*sp.frac_non_bridging:6.2f} %\n"
                            f"  free          {100*sp.frac_free:6.2f} %\n"
                            f"  triclustered  {100*sp.frac_tricluster:6.2f} %\n"
                            f"  NBO per network former  {sp.nbo_per_T:.3f}\n", None)))
                for f_, dist in sp.qn_distribution.items():
                    put(("an", (f"  Q(n) for {f_}: "
                                + "  ".join(f"Q{k} {100*v:.1f}%"
                                            for k, v in sorted(dist.items()))
                                + f"   mean {sp.mean_qn[f_]:.2f}\n", None)))
                for f_, c in sp.coordination.items():
                    put(("an", (f"  coordination of {f_} by oxygen: {c:.3f}\n", None)))
                for n in sp.notes:
                    put(("an", (f"  NOTE: {n}\n", "bad")))
                put(("an", ("\n", None)))

            if what in ("all", "msd", "cond", "figs"):
                put(("status", "diffusion"))
                put(("an", ("diffusion\n", "h")))
                Ds: dict[str, float] = {}
                for sp_name in traj.species_list():
                    try:
                        dwe = an.diffusion_with_error(traj, sp_name, n_blocks=5)
                    except ValueError as exc:
                        put(("an", (f"  {sp_name}: {exc}\n", "bad")))
                        continue
                    Ds[sp_name] = dwe.D_mean_m2_s
                    tag = "good" if dwe.whole.diffusive else "bad"
                    put(("an", (f"  {sp_name}\n", None)))
                    put(("an", (f"      D = {dwe.D_m2_s:.4e} m2/s "
                                f"= {dwe.D_cm2_s if hasattr(dwe,'D_cm2_s') else dwe.D_m2_s*1e4:.4e} cm2/s\n",
                                None)))
                    put(("an", (f"      blocks  {dwe.D_mean_m2_s:.4e} "
                                f"+/- {dwe.D_stderr_m2_s:.2e} m2/s over "
                                f"{dwe.n_blocks} blocks\n", None)))
                    put(("an", (f"      log-log slope {dwe.whole.loglog_slope:.3f}, "
                                f"R squared {dwe.whole.r_squared:.5f}\n", None)))
                    put(("an", (f"      {dwe.whole.verdict}\n", tag)))
                    put(("an", (f"      {dwe.verdict}\n", None)))
                put(("an", ("\n", None)))

                if what in ("all", "cond") and charges:
                    put(("status", "conductivity"))
                    counts = {s: int((traj.elements == s).sum())
                              for s in traj.species_list()}
                    put(("an", ("electrical conductivity\n", "h")))
                    try:
                        ne = an.conductivity_nernst_einstein(
                            counts, charges, Ds, traj.volume_A3(), T)
                        put(("an", (f"  Nernst-Einstein, assumes ions move independently:"
                                    f"  {ne:.3f} S/m\n", None)))
                        ex, info = an.conductivity_exact_einstein(traj, charges, T)
                        put(("an", (f"  exact Einstein, includes cross correlations:   "
                                    f"  {ex:.3f} S/m\n", None)))
                        put(("an", (f"  block error on the exact value: "
                                    f"{info['relative_stderr_percent']:.0f} percent\n",
                                    None)))
                        if ne:
                            ratio = ex / ne
                            put(("an", (f"  Haven ratio, exact over Nernst-Einstein: "
                                        f"{ratio:.3f}\n", None)))
                            # The interpretation depends on the number, so print the one
                            # that applies. Until 29 Sep 2026 the "below 1" sentence was
                            # printed whatever the ratio was, including when it was 1.96.
                            if ratio < 0.9:
                                put(("an", ("  Below 1: the ions do NOT move independently. "
                                            "The shortfall is exactly the correlation that "
                                            "Nernst-Einstein throws away.\n", None)))
                            elif ratio > 1.1:
                                put(("an", ("  Above 1. In a real melt that is unusual and "
                                            "normally means the exact value is not "
                                            "converged, because it rests on one collective "
                                            "coordinate. Check the block error before "
                                            "reading anything into it.\n", "bad")))
                            else:
                                put(("an", ("  Close to 1: the ions are moving more or less "
                                            "independently, so the approximation is fine "
                                            "here.\n", None)))
                        put(("an", (f"  {info['verdict']}\n",
                                    "bad" if info["relative_stderr_percent"] > 30
                                    else "good")))
                    except Exception as exc:                 # noqa: BLE001
                        put(("an", (f"  could not compute: {exc}\n", "bad")))
                    put(("an", ("\n", None)))

            if what in ("all", "visc", "figs"):
                put(("status", "viscosity"))
                put(("an", ("viscosity\n", "h")))
                sf = path.parent / "stress.txt"
                if not sf.exists():
                    put(("an", ("  There is no stress.txt, so viscosity cannot be "
                                "computed from this run.\n"
                                "  Tick the stress sampling box when adding the run. It "
                                "costs almost nothing on a\n"
                                "  production run that is happening anyway, and the "
                                "published scripts for this system\n"
                                "  never compute viscosity at all, which is why it is "
                                "worth having.\n", "bad")))
                else:
                    for fn, nm in ((an.viscosity_green_kubo, "Green-Kubo"),
                                   (an.viscosity_einstein_helfand,
                                    "Einstein-Helfand")):
                        try:
                            v = fn(sf, T, traj.volume_A3())
                            put(("an", (f"  {nm}: "
                                        f"{'could not be determined' if v.eta_Pa_s is None else f'{v.eta_Pa_s:.4e} Pa s'}\n",
                                        None)))
                            put(("an", (f"      {v.verdict}\n",
                                        "bad" if "NOT" in v.verdict or "not" in
                                        v.verdict[:40] else "good")))
                        except Exception as exc:              # noqa: BLE001
                            put(("an", (f"  {nm}: {exc}\n", "bad")))
                put(("an", ("\n", None)))

            if what in ("all", "equil"):
                put(("status", "equilibration"))
                put(("an", ("was it equilibrated?\n", "h")))
                logs = sorted(path.parent.glob("log.*"))
                if not logs:
                    put(("an", ("  no LAMMPS log found beside the trajectory\n", "bad")))
                for lg in logs:
                    if lg.suffix in (".err", ".out"):
                        continue
                    try:
                        d = an.read_lammps_log(lg)
                        if not d:
                            continue
                        rep = an.equilibration_report(d)
                        put(("an", (f"  {lg.name}\n", None)))
                        for k, v in rep.items():
                            if k == "_verdicts":
                                continue
                            put(("an", (f"      {k:9s} mean {v['mean']:12.4f}  "
                                        f"sd {v['std']:10.4f}  drift "
                                        f"{v['drift_percent_over_tail']:+7.2f} %\n",
                                        None)))
                        for vd in rep["_verdicts"]:
                            put(("an", (f"      {vd}\n",
                                        "good" if "settled" in vd else "bad")))
                    except Exception as exc:                  # noqa: BLE001
                        put(("an", (f"  {lg.name}: {exc}\n", "bad")))

            if what == "figs":
                put(("status", "drawing figures"))
                out = self._figures(traj, path.parent, T)
                put(("an", (f"\nfigures written to {out}\n", "good")))

            put(("status", "analysis finished"))
        except Exception:                                     # noqa: BLE001
            put(("an", ("\nANALYSIS FAILED\n" + traceback.format_exc(), "bad")))

    def _figures(self, traj, outdir: Path, T: float) -> Path:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        outdir = Path(outdir) / "figures"
        outdir.mkdir(exist_ok=True)
        rdf = an.compute_rdf(traj, rmax=10.0, nbins=400)
        fig, ax = plt.subplots(figsize=(7, 4.4))
        for k in sorted(rdf.g):
            if "O" in k:
                ax.plot(rdf.r, rdf.g[k], lw=1.4, label=k)
        ax.set_xlabel("r  (Angstrom)")
        ax.set_ylabel("g(r)")
        ax.set_title("partial radial distribution functions")
        ax.legend(fontsize=8, ncol=2)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / "rdf.png", dpi=140)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 4.4))
        for sp_name in traj.species_list():
            t, m = an.compute_msd(traj, sp_name)
            ax.loglog(t[1:], m[1:], lw=1.4, label=sp_name)
        ax.plot([], [], " ")
        lo = max(1e-3, traj.times_ps[1])
        ax.loglog([lo, traj.times_ps[-1] / 2],
                  [1.0, 1.0 * (traj.times_ps[-1] / 2) / lo],
                  "k--", lw=0.9, label="slope 1, diffusive")
        ax.set_xlabel("time  (ps)")
        ax.set_ylabel("MSD  (Angstrom squared)")
        ax.set_title("mean squared displacement, log-log. A slope of 1 is diffusive.")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25, which="both")
        fig.tight_layout()
        fig.savefig(outdir / "msd_loglog.png", dpi=140)
        plt.close(fig)
        return outdir

    def _save_analysis(self) -> None:
        if not self.project or not self.an_run.get():
            return
        d = self.project.run_dir(self.an_run.get())
        p = d / "ANALYSIS.txt"
        p.write_text(self.an_out.get("1.0", "end"))
        self.log(f"analysis saved to {p}")
        messagebox.showinfo("saved", str(p))

    # ---------------------------------------------------------------- tab: log
    def _tab_log(self) -> None:
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="6  Log")
        self.logbox = tk.Text(f, wrap="word", font=("DejaVu Sans Mono", 9),
                              background="#FCFCFC")
        self.logbox.pack(fill="both", expand=True, padx=14, pady=12)

    def log(self, msg: str) -> None:
        self.logbox.insert("end", msg.rstrip() + "\n")
        self.logbox.see("end")

    # ------------------------------------------------------------- plumbing
    def _work(self, fn, *args) -> None:
        self.busy.start(12)

        def runner():
            try:
                fn(*args)
            except Exception:                                # noqa: BLE001
                self.msgq.put(("log", traceback.format_exc()))
            finally:
                self.msgq.put(("done", None))
        threading.Thread(target=runner, daemon=True).start()

    def _drain(self) -> None:
        try:
            while True:
                kind, payload = self.msgq.get_nowait()
                if kind == "log":
                    self.log(str(payload))
                elif kind == "status":
                    self.status.configure(text=str(payload))
                elif kind == "refresh":
                    self._refresh_queue()
                    self._refresh_analysis_runs()
                elif kind == "comp":
                    self.comp_out.insert("end", str(payload) + "\n")
                    self.comp_out.see("end")
                elif kind == "an_clear":
                    self.an_out.delete("1.0", "end")
                elif kind == "an":
                    text, tag = payload
                    self.an_out.insert("end", text, tag or ())
                    self.an_out.see("end")
                elif kind == "done":
                    self.busy.stop()
        except queue.Empty:
            pass
        self.after(200, self._drain)

    def _show_text(self, title: str, body: str) -> None:
        w = tk.Toplevel(self)
        w.title(title)
        w.geometry("980x640")
        t = tk.Text(w, wrap="none", font=("DejaVu Sans Mono", 9))
        ys = ttk.Scrollbar(w, orient="vertical", command=t.yview)
        xs = ttk.Scrollbar(w, orient="horizontal", command=t.xview)
        t.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        t.grid(row=0, column=0, sticky="nsew")
        ys.grid(row=0, column=1, sticky="ns")
        xs.grid(row=1, column=0, sticky="ew")
        w.rowconfigure(0, weight=1)
        w.columnconfigure(0, weight=1)
        t.insert("1.0", body)


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="slagmd graphical workflow manager")
    ap.add_argument("project", nargs="?", help="project folder to open")
    a = ap.parse_args(argv)
    App(Path(a.project) if a.project else None).mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
