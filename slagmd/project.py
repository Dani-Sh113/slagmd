"""The state machine. A project on disk, a queue of runs, and the loop that keeps it moving.

Everything is written to disk the moment it changes. There is no state that lives only in
memory, and no batch that is held up to be written at the end. If the laptop's internet
drops, or the session dies, or the GUI is closed, the next thing that opens this project
reads project.json and carries on from the same place. That is the whole design.

One run is one composition at one temperature with one force field. It moves through the
five stages in order. A stage is only considered done when the input script's own marker
is in its log, never because the scheduler said COMPLETED.
"""
from __future__ import annotations

import json
import os
import shutil
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import diagnose as diag
from . import scheduler as sched
from .builder import build_random, build_with_packmol, check_structure, write_data_file
from .composition import Composition
from .forcefield import ForceField
from .inputgen import EqBlock, Protocol, stage_order, write_inputs

SCHEMA = 2


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _proto_to_dict(p: Protocol) -> dict:
    d = {k: v for k, v in asdict(p).items() if k != "eq_blocks"}
    d["eq_blocks"] = [asdict(b) for b in p.eq_blocks]
    return d


def _proto_from_dict(d: dict) -> Protocol:
    d = dict(d)
    blocks = d.pop("eq_blocks", None)
    p = Protocol(**d)
    if blocks is not None:
        p.eq_blocks = [EqBlock(**b) if isinstance(b, dict) else b for b in blocks]
    return p


def atomic_write(path: Path, text: str) -> None:
    """Write via a temp file then rename, so a crash mid-write cannot corrupt the state."""
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


# ============================================================================ run
@dataclass
class StageState:
    name: str
    status: str = "NOT_SUBMITTED"
    job_id: str = ""
    attempts: int = 0
    submitted_at: str = ""
    finished_at: str = ""
    last_diagnosis: str = ""
    history: list[str] = field(default_factory=list)

    def log(self, msg: str) -> None:
        self.history.append(f"{_now()}  {msg}")


@dataclass
class Run:
    name: str
    mole_fractions: dict[str, float]
    temperature_K: float
    force_field: str
    target_atoms: int = 5000
    initial_density: float = 2.6
    seed: int = 12345
    protocol_overrides: dict = field(default_factory=dict)
    cluster_overrides: dict = field(default_factory=dict)
    stages: dict[str, StageState] = field(default_factory=dict)
    stage_sequence: list[str] = field(default_factory=list)
    prepared: bool = False
    created_at: str = field(default_factory=_now)
    notes: list[str] = field(default_factory=list)

    def stage_list(self) -> list[str]:
        """Stages in run order. The order is recorded when the run is prepared, because
        the number of equilibration blocks is a protocol choice and not a constant."""
        if self.stage_sequence:
            return [s for s in self.stage_sequence if s in self.stages]
        return sorted(self.stages)

    def current_stage(self) -> str | None:
        """The first stage that is not COMPLETED."""
        for s in self.stage_list():
            if self.stages[s].status != sched.COMPLETED:
                return s
        return None

    def is_done(self) -> bool:
        return self.current_stage() is None and self.prepared

    def is_blocked(self) -> bool:
        st = self.current_stage()
        return bool(st and self.stages[st].status == "BLOCKED")

    def progress(self) -> str:
        done = sum(1 for s in self.stage_list()
                   if self.stages[s].status == sched.COMPLETED)
        return f"{done}/{len(self.stage_list())}"


# ======================================================================== project
class Project:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.state_file = self.root / "project.json"
        self.cluster = sched.ClusterConfig()
        self.protocol = Protocol()
        self.runs: dict[str, Run] = {}
        self.event_log = self.root / "events.log"
        self.auto_diagnose = True
        self.max_attempts = diag.MAX_ATTEMPTS
        if self.state_file.exists():
            self.load()
        else:
            probe = sched.probe_cluster()
            self.cluster.backend = probe["backend"]
            if "lammps_exe" in probe:
                self.cluster.lammps_exe = probe["lammps_exe"]
            if probe["backend"] == "local":
                self.cluster.nodes = 1
                self.cluster.tasks_per_node = max(1, min(8, probe["cores_here"]))
                self.cluster.mpi_launcher = "mpirun" if shutil.which("mpirun") else ""
                self.cluster.modules = []
            self.save()

    # ------------------------------------------------------------------ persistence
    def save(self) -> None:
        payload = {
            "schema": SCHEMA,
            "saved_at": _now(),
            "cluster": asdict(self.cluster),
            "protocol_defaults": _proto_to_dict(self.protocol),
            "auto_diagnose": self.auto_diagnose,
            "max_attempts": self.max_attempts,
            "runs": {
                name: {**{k: v for k, v in asdict(r).items() if k != "stages"},
                       "stages": {s: asdict(st) for s, st in r.stages.items()}}
                for name, r in self.runs.items()},
        }
        atomic_write(self.state_file, json.dumps(payload, indent=2))

    def load(self) -> None:
        d = json.loads(self.state_file.read_text())
        self.cluster = sched.ClusterConfig(**d.get("cluster", {}))
        self.protocol = _proto_from_dict(d.get("protocol_defaults", {}))
        self.auto_diagnose = d.get("auto_diagnose", True)
        self.max_attempts = d.get("max_attempts", diag.MAX_ATTEMPTS)
        self.runs = {}
        for name, rd in d.get("runs", {}).items():
            stages = {s: StageState(**sd) for s, sd in rd.pop("stages", {}).items()}
            self.runs[name] = Run(stages=stages, **rd)

    def event(self, msg: str) -> None:
        with self.event_log.open("a") as fh:
            fh.write(f"{_now()}  {msg}\n")

    # ------------------------------------------------------------------- run setup
    def run_dir(self, name: str) -> Path:
        return self.root / "runs" / name

    def add_run(self, name: str, mole_fractions: dict[str, float],
                temperature_K: float, force_field: str = "bouhadja",
                target_atoms: int = 5000, initial_density: float = 2.6,
                seed: int = 12345, protocol_overrides: dict | None = None,
                nve_check: bool = False, stress_for_viscosity: bool = False) -> Run:
        if name in self.runs:
            raise ValueError(f"a run called {name!r} already exists")
        po = dict(protocol_overrides or {})
        po["target_T"] = temperature_K
        if nve_check:
            po.setdefault("nve_check_ns", 0.1)
        if stress_for_viscosity:
            po.setdefault("stress_every_steps", 5)
        r = Run(name=name, mole_fractions=dict(mole_fractions),
                temperature_K=temperature_K, force_field=force_field,
                target_atoms=target_atoms, initial_density=initial_density,
                seed=seed, protocol_overrides=po)
        self.runs[name] = r
        self.save()
        self.event(f"run {name} added, {mole_fractions} at {temperature_K} K, {force_field}")
        return r

    def protocol_for(self, run: Run) -> Protocol:
        p = _proto_from_dict(_proto_to_dict(self.protocol))
        for k, v in run.protocol_overrides.items():
            if k == "eq_blocks":
                p.eq_blocks = [EqBlock(**b) if isinstance(b, dict) else b for b in v]
            elif k == "_scale_eq_pdamp":
                # a remedy asking for a gentler barostat across every equilibration block
                for b in p.eq_blocks:
                    b.pdamp_fs *= float(v)
            elif k == "_scale_eq_tdamp":
                for b in p.eq_blocks:
                    b.tdamp_fs *= float(v)
            elif hasattr(p, k):
                setattr(p, k, v)
        p.target_T = run.temperature_K
        return p

    def cluster_for(self, run: Run) -> sched.ClusterConfig:
        c = sched.ClusterConfig(**asdict(self.cluster))
        for k, v in run.cluster_overrides.items():
            if hasattr(c, k):
                setattr(c, k, v)
        return c

    def prepare(self, run: Run, use_packmol: bool = False,
                progress=None) -> list[str]:
        """Build the structure, write the deck, write the submit scripts. Returns notes."""
        rundir = self.run_dir(run.name)
        rundir.mkdir(parents=True, exist_ok=True)
        notes: list[str] = []

        def say(m: str) -> None:
            notes.append(m)
            if progress:
                progress(m)

        comp = Composition(run.mole_fractions, run.target_atoms)
        ff = ForceField.load(run.force_field)
        ff.check_ready(comp.species)
        say(f"force field {ff.label} passed its consistency identities")

        ok, scale = comp.charges_are_scaled_formal(ff.charges)
        if ok:
            say(f"charges are a uniform {scale:.4f} of the formal charges, so any "
                f"stoichiometric composition is exactly neutral")
        else:
            say("WARNING: the charges are not a uniform scaling of the formal charges. "
                "Neutrality was checked for this box only.")
        q = comp.net_charge(ff.charges)
        if abs(q) > 1e-6:
            raise ValueError(f"net box charge {q:+.6f}, refusing to build")

        say(comp.summary())

        if use_packmol:
            try:
                struct = build_with_packmol(comp, rundir / "pack",
                                           run.initial_density, run.seed)
                say("initial configuration built by Packmol, as in Verbeeck 2025")
            except Exception as exc:                  # noqa: BLE001
                say(f"packmol unavailable or failed ({exc}), using the internal packer")
                struct = build_random(comp, run.initial_density, run.seed)
        else:
            struct = build_random(comp, run.initial_density, run.seed)
            say("initial configuration built by the internal random packer")

        problems = check_structure(struct, comp, ff.charges)
        if problems:
            raise ValueError("the starting structure failed its checks:\n  "
                             + "\n  ".join(problems))
        say(f"structure checked: {struct.n_atoms} atoms, box {struct.box:.2f} A, "
            f"closest pair {struct.min_distance():.2f} A, "
            f"density {struct.density_g_cm3():.3f} g/cm3")

        write_data_file(struct, ff.charges, rundir / "melt.data",
                        title=f"slagmd {run.name}")

        proto = self.protocol_for(run)
        for w in proto.warnings():
            say(f"PROTOCOL WARNING: {w}")
        paths = write_inputs(rundir, run.name, comp, ff, proto,
                             data_file="melt.data", date=_now())
        say(f"wrote {len(paths)} stage scripts, total simulated time "
            f"{proto.total_ns():.2f} ns")

        cfg = self.cluster_for(run)
        for prob in sched.config_problems(cfg):
            say(f"CLUSTER WARNING: {prob}")
        for stage in paths:
            script = sched.submit_script(cfg, f"{run.name}_{stage}", rundir,
                                         f"in.{stage}", f"log.{stage}")
            sp = rundir / f"submit.{stage}.sh"
            sp.write_text(script)
            sp.chmod(0o755)

        run.stage_sequence = stage_order(proto)
        run.stages = {s: StageState(name=s) for s in run.stage_sequence if s in paths}
        run.prepared = True
        run.notes = notes
        self.save()
        self.event(f"run {run.name} prepared, {len(paths)} stages")
        return notes

    # -------------------------------------------------------------------- the loop
    def submit_next(self, run: Run) -> str:
        """Submit the first stage that is not complete. Returns a message."""
        if not run.prepared:
            return f"{run.name}: not prepared yet"
        stage = run.current_stage()
        if stage is None:
            return f"{run.name}: all stages complete"
        st = run.stages[stage]
        if st.status in (sched.PENDING, sched.RUNNING):
            return f"{run.name}/{stage}: already {st.status}"
        if st.status == "BLOCKED":
            return (f"{run.name}/{stage}: BLOCKED, needs a human. "
                    f"{st.last_diagnosis}")
        if st.attempts >= self.max_attempts:
            st.status = "BLOCKED"
            st.log(f"gave up after {st.attempts} attempts")
            self.save()
            return f"{run.name}/{stage}: {st.attempts} attempts, giving up"

        # never resubmit a stage whose predecessor did not really finish
        order = run.stage_list()
        i = order.index(stage)
        if i > 0:
            prev = order[i - 1]
            if not diag.stage_completed(self.run_dir(run.name), prev):
                run.stages[prev].status = "NOT_SUBMITTED"
                self.save()
                return (f"{run.name}: {prev} has no completion marker, so {stage} cannot "
                        f"start. Rolled {prev} back to be resubmitted.")

        cfg = self.cluster_for(run)
        try:
            job_id, reply = sched.submit(cfg, self.run_dir(run.name),
                                         f"submit.{stage}.sh")
        except sched.SubmitError as exc:
            st.status = "BLOCKED"
            st.last_diagnosis = str(exc)
            st.log(f"submission refused: {exc}")
            self.save()
            return f"{run.name}/{stage}: submission failed, {exc}"
        st.job_id = job_id
        st.status = sched.PENDING
        st.attempts += 1
        st.submitted_at = _now()
        st.log(f"submitted as {job_id} (attempt {st.attempts})")
        self.save()                      # written BEFORE anything else happens
        self.event(f"{run.name}/{stage} submitted as {job_id}")
        return f"{run.name}/{stage}: submitted as {job_id}"

    def poll(self, run: Run) -> list[str]:
        """Update one run's state from the scheduler and the logs. Returns messages."""
        msgs: list[str] = []
        if not run.prepared:
            return msgs
        rundir = self.run_dir(run.name)
        cfg = self.cluster_for(run)
        for stage in run.stage_list():
            st = run.stages[stage]
            if st.status in (sched.COMPLETED, "BLOCKED", "NOT_SUBMITTED"):
                continue
            if not st.job_id:
                continue
            state, detail = sched.status(cfg, st.job_id)

            # the marker is the only real authority on whether the stage finished
            if diag.stage_completed(rundir, stage):
                if st.status != sched.COMPLETED:
                    st.status = sched.COMPLETED
                    st.finished_at = _now()
                    st.log(f"completion marker found in the log ({detail})")
                    msgs.append(f"{run.name}/{stage}: COMPLETE")
                    self.event(f"{run.name}/{stage} completed")
                continue

            if state in (sched.PENDING, sched.RUNNING):
                if st.status != state:
                    st.status = state
                    st.log(f"scheduler says {state}: {detail}")
                continue

            # it left the queue without the marker. Something went wrong.
            st.status = state if state in diag.__dict__ or state else sched.FAILED
            st.finished_at = _now()
            d = diag.diagnose_job(rundir, stage)
            rem = d.remedy
            if rem is None:
                st.status = sched.COMPLETED if "normal completion" in d.evidence \
                    else sched.FAILED
                st.log(f"left the queue as {state}, no remedy needed. {d.evidence}")
                msgs.append(f"{run.name}/{stage}: ended as {state}, {d.evidence}")
                self.save()
                continue

            st.last_diagnosis = f"{rem.code}: {rem.plain}"
            st.log(f"scheduler state {state}. diagnosis {rem.code}. {rem.what_happened} "
                   f"WHY: {rem.why} ACTION: {rem.action}")
            (rundir / f"diagnosis.{stage}.{st.attempts}.txt").write_text(
                f"stage      {stage}\nattempt    {st.attempts}\n"
                f"scheduler  {state}  {detail}\n\ncode       {rem.code}\n\n"
                f"WHAT HAPPENED\n{rem.what_happened}\n\nWHY\n{rem.why}\n\n"
                f"ACTION\n{rem.action}\n\nEVIDENCE\n{d.evidence}\n\n"
                f"LOG TAIL\n{d.log_tail}\n")
            msgs.append(f"{run.name}/{stage}: FAILED, {rem.code}, {rem.plain}")
            self.event(f"{run.name}/{stage} failed {rem.code}")

            if not (self.auto_diagnose and rem.auto_retry):
                st.status = "BLOCKED"
                msgs.append(f"{run.name}/{stage}: BLOCKED for a human. {rem.action}")
                self.save()
                continue
            if st.attempts >= self.max_attempts:
                st.status = "BLOCKED"
                st.log(f"{st.attempts} attempts reached, not retrying")
                msgs.append(f"{run.name}/{stage}: BLOCKED, {st.attempts} attempts")
                self.save()
                continue

            msgs += self._apply_remedy(run, stage, rem)
        self.save()
        return msgs

    def _apply_remedy(self, run: Run, stage: str, rem: diag.Remedy) -> list[str]:
        msgs: list[str] = []
        if rem.protocol_changes:
            run.protocol_overrides.update(rem.protocol_changes)
            msgs.append(f"{run.name}: protocol changed {rem.protocol_changes}")
        if rem.cluster_changes:
            newc = diag.apply_cluster_change(
                asdict(self.cluster_for(run)), rem.cluster_changes)
            for k in rem.cluster_changes:
                run.cluster_overrides[k] = newc[k]
            msgs.append(f"{run.name}: cluster changed to "
                        f"{ {k: newc[k] for k in rem.cluster_changes} }")

        if rem.rebuild_structure or rem.protocol_changes:
            try:
                keep = {s: run.stages[s] for s in run.stage_list()}
                notes = self.prepare(run)
                for s, old in keep.items():
                    if s in run.stages:
                        run.stages[s].attempts = old.attempts
                        run.stages[s].history = old.history
                        if not rem.rebuild_structure and s != stage \
                                and old.status == sched.COMPLETED:
                            run.stages[s].status = sched.COMPLETED
                msgs.append(f"{run.name}: rebuilt with the new settings")
            except Exception as exc:                  # noqa: BLE001
                run.stages[stage].status = "BLOCKED"
                run.stages[stage].last_diagnosis = f"rebuild failed: {exc}"
                msgs.append(f"{run.name}: rebuild FAILED, {exc}")
                self.save()
                return msgs

        target = rem.resume_from_stage or stage
        order = run.stage_list()
        if target in order:
            for s in order[order.index(target):]:
                run.stages[s].status = "NOT_SUBMITTED"
                run.stages[s].job_id = ""
        self.save()
        msgs.append(self.submit_next(run))
        return msgs

    # ---------------------------------------------------------------- bulk actions
    def poll_all(self) -> list[str]:
        msgs: list[str] = []
        for run in list(self.runs.values()):
            msgs += self.poll(run)
        return msgs

    def n_active(self) -> int:
        return sum(1 for r in self.runs.values() for s in r.stage_list()
                   if r.stages[s].status in (sched.PENDING, sched.RUNNING))

    def fill_queue(self, max_concurrent: int = 4) -> list[str]:
        """Keep up to max_concurrent jobs in flight, submitting the next stage of
        whichever runs are idle. This is the throughput knob."""
        msgs: list[str] = []
        for run in self.runs.values():
            if self.n_active() >= max_concurrent:
                break
            if not run.prepared or run.is_done() or run.is_blocked():
                continue
            stage = run.current_stage()
            if stage and run.stages[stage].status == "NOT_SUBMITTED":
                msgs.append(self.submit_next(run))
        return msgs

    def cancel_run(self, run: Run) -> list[str]:
        out = []
        cfg = self.cluster_for(run)
        for s in run.stage_list():
            st = run.stages[s]
            if st.status in (sched.PENDING, sched.RUNNING) and st.job_id:
                out.append(sched.cancel(cfg, st.job_id))
                st.status = sched.CANCELLED
                st.log("cancelled by the user")
        self.save()
        return out

    def reset_stage(self, run: Run, stage: str) -> None:
        """Unblock a stage after a human has fixed whatever it was."""
        st = run.stages[stage]
        st.status = "NOT_SUBMITTED"
        st.attempts = 0
        st.job_id = ""
        st.last_diagnosis = ""
        st.log("reset by the user")
        self.save()

    # -------------------------------------------------------------- import a study
    def add_verbeeck_benchmark(self, force_field: str = "bouhadja",
                              temperatures_C: tuple[int, ...] = (1600,),
                              target_atoms: int = 5000) -> list[str]:
        """Queue the ten benchmark compositions of Verbeeck 2025 Table 5.

        This exists so the tool can be validated against a published study rather than
        against nothing. If it cannot reproduce their self-diffusion coefficients, the
        tool is wrong, and that is worth knowing before using it on anything new.
        """
        from .composition import VERBEECK_2025_TABLE5
        added = []
        for label, mf, (tlo, thi) in VERBEECK_2025_TABLE5:
            for tc in temperatures_C:
                if not (tlo <= tc <= thi):
                    continue
                name = f"{label}_{tc}C_{force_field}"
                if name in self.runs:
                    continue
                self.add_run(name, mf, tc + 273.15, force_field,
                             target_atoms=target_atoms)
                added.append(name)
        return added

    # ------------------------------------------------------------------- reporting
    def summary_rows(self) -> list[dict]:
        rows = []
        for r in self.runs.values():
            stage = r.current_stage()
            st = r.stages.get(stage) if stage else None
            rows.append({
                "run": r.name,
                "composition": "  ".join(f"{k} {v:.3f}"
                                         for k, v in sorted(r.mole_fractions.items())),
                "T_C": f"{r.temperature_K - 273.15:.0f}",
                "ff": r.force_field,
                "progress": r.progress(),
                "stage": stage or "done",
                "status": st.status if st else "COMPLETE",
                "job": st.job_id if st else "",
                "attempts": st.attempts if st else 0,
                "note": st.last_diagnosis if st else "",
            })
        return rows

    def write_report(self) -> Path:
        """A plain text report of the whole project, for reading on a phone."""
        lines = [f"slagmd project report   {_now()}",
                 f"root  {self.root}", "",
                 f"backend {self.cluster.backend}, {self.cluster.nodes} nodes x "
                 f"{self.cluster.tasks_per_node} tasks, walltime "
                 f"{self.cluster.walltime}", ""]
        for row in self.summary_rows():
            lines.append(f"{row['run']:28s} {row['T_C']:>5s} C  {row['ff']:9s} "
                         f"{row['progress']:>5s}  {row['stage']:14s} {row['status']:12s} "
                         f"{row['note']}")
        lines.append("")
        blocked = [r for r in self.runs.values() if r.is_blocked()]
        if blocked:
            lines.append("NEEDS A HUMAN:")
            for r in blocked:
                s = r.current_stage()
                lines.append(f"  {r.name}/{s}: {r.stages[s].last_diagnosis}")
        p = self.root / "REPORT.txt"
        atomic_write(p, "\n".join(lines) + "\n")
        return p
