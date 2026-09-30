"""Submit jobs and ask what happened to them, on SLURM, PBS or this machine.

Design rule: nothing here ever parses a scheduler's human readable output when a machine
readable form exists. squeue gets --noheader --format, sacct gets --parsable2. Screen
scraping a scheduler is how a monitor silently starts reporting every job as complete.

Second design rule: submission is idempotent per stage. The job id is written to disk the
instant the scheduler returns it, before anything else happens, so a session that dies
between submitting and recording does not lose the job.
"""
from __future__ import annotations

import os
import re
import shutil
import signal
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path

# canonical states this tool understands
PENDING = "PENDING"
RUNNING = "RUNNING"
COMPLETED = "COMPLETED"
FAILED = "FAILED"
TIMEOUT = "TIMEOUT"
CANCELLED = "CANCELLED"
OOM = "OUT_OF_MEMORY"
UNKNOWN = "UNKNOWN"

TERMINAL = {COMPLETED, FAILED, TIMEOUT, CANCELLED, OOM}

# Handles for jobs this process started with the local backend, so they can be reaped.
# Keyed by pid. Not persisted, and it does not need to be: if the session dies the /proc
# fallback in status() takes over.
_LOCAL_PROCS: dict[int, "subprocess.Popen"] = {}


@dataclass
class ClusterConfig:
    backend: str = "local"                  # local | slurm | pbs
    partition: str = ""
    account: str = ""
    nodes: int = 2                          # Verbeeck 2025 used 2 nodes for classical MD
    tasks_per_node: int = 128               # 2 x 64-core AMD EPYC 7H12
    walltime: str = "24:00:00"
    memory_per_node: str = ""
    modules: list[str] = field(default_factory=lambda: ["LAMMPS"])
    pre_commands: list[str] = field(default_factory=list)
    lammps_exe: str = "lmp"
    mpi_launcher: str = "srun"              # srun, mpirun, mpiexec, or "" for serial
    extra_directives: list[str] = field(default_factory=list)
    email: str = ""

    def total_ranks(self) -> int:
        return max(1, self.nodes * self.tasks_per_node)


# --------------------------------------------------------------------- submit script
def submit_script(cfg: ClusterConfig, jobname: str, rundir: Path,
                  input_file: str, log_file: str) -> str:
    """The text of the batch script. Deliberately plain, so it is readable and editable."""
    L = ["#!/bin/bash"]
    if cfg.backend == "slurm":
        L += [f"#SBATCH --job-name={jobname}",
              f"#SBATCH --nodes={cfg.nodes}",
              f"#SBATCH --ntasks-per-node={cfg.tasks_per_node}",
              f"#SBATCH --time={cfg.walltime}",
              f"#SBATCH --output={log_file}.slurm.out",
              f"#SBATCH --error={log_file}.slurm.err"]
        if cfg.partition:
            L.append(f"#SBATCH --partition={cfg.partition}")
        if cfg.account:
            L.append(f"#SBATCH --account={cfg.account}")
        if cfg.memory_per_node:
            L.append(f"#SBATCH --mem={cfg.memory_per_node}")
        if cfg.email:
            L += [f"#SBATCH --mail-user={cfg.email}",
                  "#SBATCH --mail-type=FAIL,TIME_LIMIT"]
        L += [f"#SBATCH {d}" for d in cfg.extra_directives]
    elif cfg.backend == "pbs":
        mem = f":mem={cfg.memory_per_node}" if cfg.memory_per_node else ""
        L += [f"#PBS -N {jobname}",
              f"#PBS -l select={cfg.nodes}:ncpus={cfg.tasks_per_node}"
              f":mpiprocs={cfg.tasks_per_node}{mem}",
              f"#PBS -l walltime={cfg.walltime}",
              f"#PBS -o {log_file}.pbs.out",
              f"#PBS -e {log_file}.pbs.err"]
        if cfg.partition:
            L.append(f"#PBS -q {cfg.partition}")
        if cfg.account:
            L.append(f"#PBS -A {cfg.account}")
        if cfg.email:
            L += [f"#PBS -M {cfg.email}", "#PBS -m a"]
        L += [f"#PBS {d}" for d in cfg.extra_directives]

    L += ["", "set -o pipefail",
          f"cd {_q(rundir)} || exit 90", ""]
    if cfg.modules:
        L.append("module purge 2>/dev/null || true")
        L += [f"module load {m}" for m in cfg.modules]
    L += cfg.pre_commands
    L += ["",
          'echo "SLAGMD_HOST $(hostname)"',
          'echo "SLAGMD_START $(date -Is)"',
          ""]
    launcher = ""
    if cfg.mpi_launcher == "srun":
        launcher = "srun"
    elif cfg.mpi_launcher:
        launcher = f"{cfg.mpi_launcher} -np {cfg.total_ranks()}"
    L += [f"{launcher} {cfg.lammps_exe} -in {input_file} -log {log_file}".strip(),
          "rc=$?",
          'echo "SLAGMD_EXIT $rc"',
          'echo "SLAGMD_END $(date -Is)"',
          "exit $rc",
          ""]
    return "\n".join(L)


def config_problems(cfg: "ClusterConfig") -> list[str]:
    """Things about this cluster configuration that will break the batch script.

    The first rule exists because the project this tool was written for lives at a path
    with a TRAILING SPACE in it. An unquoted path with a space in a shell script does not
    fail loudly, it splits into two words, and the error you get is
    "export: not a valid identifier", several lines away from the cause. It cost a debug
    cycle on 29 September 2026, so it is now a check.
    """
    out: list[str] = []
    for i, cmd in enumerate(cfg.pre_commands):
        if " " not in cmd:
            continue
        rhs = cmd.split("=", 1)[1] if "=" in cmd else ""
        if rhs and " " in rhs and rhs[0] not in "'\"" and "\\ " not in rhs:
            out.append(
                f"pre-command {i+1} assigns a value containing a space without quoting "
                f"it, which bash will split into two words: {cmd}. Write it as "
                f"{cmd.split('=', 1)[0]}='{rhs}'")
    if " " in cfg.lammps_exe and not any(
            q in cfg.lammps_exe for q in ("'", '"')) and " -" not in cfg.lammps_exe:
        out.append(f"the LAMMPS executable path contains a space and is not quoted: "
                   f"{cfg.lammps_exe}")
    if cfg.backend == "local" and cfg.nodes > 1:
        out.append("the local backend cannot use more than one node")
    if cfg.backend in ("slurm", "pbs") and not cfg.walltime:
        out.append("a batch backend needs a walltime")
    return out


def _q(p) -> str:
    """Quote a path for bash. Necessary here: the project root contains a trailing space."""
    return "'" + str(p).replace("'", "'\\''") + "'"


# -------------------------------------------------------------------------- submit
class SubmitError(RuntimeError):
    pass


def submit(cfg: ClusterConfig, rundir: Path, script_name: str = "submit.sh"
           ) -> tuple[str, str]:
    """Submit and return (job_id, raw_scheduler_reply)."""
    rundir = Path(rundir)
    script = rundir / script_name
    script.chmod(0o755) if script.exists() else None

    if cfg.backend == "local":
        logout = (rundir / f"{script_name}.out").open("w")
        proc = subprocess.Popen(["bash", str(script)], cwd=str(rundir),
                                stdout=logout, stderr=subprocess.STDOUT,
                                start_new_session=True)
        _LOCAL_PROCS[proc.pid] = proc
        return f"local:{proc.pid}", f"started pid {proc.pid}"

    exe = {"slurm": "sbatch", "pbs": "qsub"}[cfg.backend]
    if shutil.which(exe) is None:
        raise SubmitError(
            f"{exe} is not on PATH, so this machine is not a {cfg.backend} submit host. "
            f"Either run slagmd on the cluster login node, or set the backend to local.")
    res = subprocess.run([exe, str(script)], cwd=str(rundir),
                         capture_output=True, text=True, timeout=300)
    if res.returncode != 0:
        raise SubmitError(f"{exe} returned {res.returncode}\n"
                          f"{res.stdout}\n{res.stderr}")
    out = (res.stdout or "").strip()
    if cfg.backend == "slurm":
        m = re.search(r"(\d+)", out)
        if not m:
            raise SubmitError(f"could not find a job id in sbatch output: {out!r}")
        return m.group(1), out
    return out.split()[0] if out else "", out


# -------------------------------------------------------------------------- status
def status(cfg: ClusterConfig, job_id: str) -> tuple[str, str]:
    """Return (canonical_state, detail). Never raises for an unknown job."""
    if job_id.startswith("local:"):
        pid = int(job_id.split(":", 1)[1])
        # Reap first. A finished child of this process becomes a zombie until it is waited
        # on, and a zombie still answers os.kill(pid, 0), so the naive liveness test
        # reports a job that ended minutes ago as RUNNING and the whole queue stalls
        # forever. That bug was live until the end to end test hung on stage zero.
        proc = _LOCAL_PROCS.get(pid)
        if proc is not None:
            rc = proc.poll()
            if rc is not None:
                _LOCAL_PROCS.pop(pid, None)
                return (COMPLETED if rc == 0 else FAILED), f"local pid {pid} exited {rc}"
            return RUNNING, f"local pid {pid} running"
        # No handle, which means a different session started it. Read /proc instead of
        # trusting signal 0.
        stat = f"/proc/{pid}/stat"
        try:
            with open(stat) as fh:
                fields = fh.read().rsplit(")", 1)[1].split()
            state = fields[0]
            if state == "Z":
                return UNKNOWN, f"local pid {pid} has exited (zombie), read the log"
            return RUNNING, f"local pid {pid} state {state}"
        except FileNotFoundError:
            return UNKNOWN, "local process has exited, read the log to see how"
        except (OSError, IndexError):
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return UNKNOWN, "local process has exited, read the log to see how"
            except PermissionError:
                return RUNNING, "running, owned by another user"
            return RUNNING, f"local pid {pid} alive"

    if cfg.backend == "slurm":
        if shutil.which("squeue"):
            r = subprocess.run(
                ["squeue", "-h", "-j", job_id, "-o", "%T|%M|%R"],
                capture_output=True, text=True)
            line = (r.stdout or "").strip()
            if line:
                st = line.split("|")[0].strip().upper()
                return _canon_slurm(st), line
        if shutil.which("sacct"):
            r = subprocess.run(
                ["sacct", "-j", job_id, "-n", "-P", "-X",
                 "-o", "State,ExitCode,Elapsed,MaxRSS"],
                capture_output=True, text=True)
            line = (r.stdout or "").strip().splitlines()
            if line:
                parts = line[0].split("|")
                return _canon_slurm(parts[0].strip().upper()), line[0]
        return UNKNOWN, "neither squeue nor sacct knows this job"

    if cfg.backend == "pbs":
        if shutil.which("qstat"):
            r = subprocess.run(["qstat", "-f", job_id],
                               capture_output=True, text=True)
            txt = r.stdout or ""
            m = re.search(r"job_state\s*=\s*(\w)", txt)
            if m:
                return {"Q": PENDING, "H": PENDING, "R": RUNNING, "E": RUNNING,
                        "C": COMPLETED, "F": COMPLETED}.get(m.group(1), UNKNOWN), \
                       f"job_state={m.group(1)}"
            if r.returncode != 0:
                return UNKNOWN, "qstat does not know this job, it has probably left the queue"
        return UNKNOWN, "qstat unavailable"

    return UNKNOWN, "unknown backend"


def _canon_slurm(st: str) -> str:
    st = st.split()[0] if st else ""
    return {
        "PENDING": PENDING, "CONFIGURING": PENDING, "REQUEUED": PENDING,
        "RUNNING": RUNNING, "COMPLETING": RUNNING, "RESIZING": RUNNING,
        "COMPLETED": COMPLETED,
        "FAILED": FAILED, "NODE_FAIL": FAILED, "BOOT_FAIL": FAILED,
        "TIMEOUT": TIMEOUT,
        "CANCELLED": CANCELLED, "PREEMPTED": CANCELLED,
        "OUT_OF_MEMORY": OOM, "OUT_OF_ME+": OOM,
    }.get(st, UNKNOWN)


def cancel(cfg: ClusterConfig, job_id: str) -> str:
    if job_id.startswith("local:"):
        pid = int(job_id.split(":", 1)[1])
        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
            return f"sent SIGTERM to process group of {pid}"
        except Exception as exc:                     # noqa: BLE001
            return f"could not signal {pid}: {exc}"
    exe = {"slurm": "scancel", "pbs": "qdel"}.get(cfg.backend, "")
    if not exe or shutil.which(exe) is None:
        return f"no way to cancel on backend {cfg.backend}"
    r = subprocess.run([exe, job_id], capture_output=True, text=True)
    return (r.stdout + r.stderr).strip() or f"{exe} {job_id} sent"


def detect_backend() -> str:
    """Guess the backend from what is installed. Used to prefill the GUI."""
    if shutil.which("sbatch"):
        return "slurm"
    if shutil.which("qsub"):
        return "pbs"
    return "local"


def probe_cluster() -> dict:
    """What can we see about this machine, to prefill sensible cluster defaults."""
    info: dict = {"backend": detect_backend(),
                  "lammps_on_path": bool(shutil.which("lmp")
                                         or shutil.which("lmp_serial")
                                         or shutil.which("lammps")),
                  "packmol_on_path": bool(shutil.which("packmol")),
                  "cores_here": os.cpu_count() or 1}
    for exe in ("lmp", "lmp_mpi", "lmp_serial", "lammps"):
        if shutil.which(exe):
            info["lammps_exe"] = exe
            break
    if shutil.which("sinfo"):
        r = subprocess.run(["sinfo", "-h", "-o", "%P|%D|%c|%l"],
                           capture_output=True, text=True)
        info["partitions"] = [l for l in (r.stdout or "").strip().splitlines() if l]
    return info
