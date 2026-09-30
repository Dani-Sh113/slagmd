"""Read a dead job's log and say what killed it, then prescribe the change that fixes it.

Each rule is a real LAMMPS failure with a real cause. The prescription is a concrete edit
to the protocol or the cluster config, not advice. The project driver applies it and
resubmits, and records both, so the history shows what was changed and why.

The important part is the last rule. If nothing matches, the tool says so and does NOT
resubmit. A resubmission loop that does not understand the failure burns the allocation
and produces the same crash, and that is worse than stopping.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Remedy:
    code: str
    what_happened: str            # plain English, for Danish
    why: str                      # the physics or the mechanics
    action: str                   # what the tool will do
    protocol_changes: dict = field(default_factory=dict)
    cluster_changes: dict = field(default_factory=dict)
    rebuild_structure: bool = False
    resume_from_stage: str | None = None   # None means resume the same stage
    auto_retry: bool = True
    plain: str = ""               # one line for the GUI queue table


@dataclass
class Diagnosis:
    matched: bool
    remedy: Remedy | None
    evidence: str                 # the lines from the log that triggered it
    log_tail: str


# ------------------------------------------------------------------- the rule table
# (regex, builder). Order matters: the most specific patterns come first.
def _rules():
    return [
        (r"ERROR.*Lost atoms|Lost atoms:",
         lambda m: Remedy(
             code="LOST_ATOMS",
             what_happened="LAMMPS lost atoms. The count of atoms it can find no longer "
                           "matches the count it started with.",
             why="Two atoms got close enough that the repulsive part of the potential "
                 "threw one of them across the box in a single step. In an oxide melt "
                 "this is almost always the first dynamics step after a bad starting "
                 "structure, or a timestep too long for the temperature. At 4000 K the "
                 "atoms are moving fast and 2 fs is already near the limit.",
             action="Halve the timestep for the hot stage, minimise harder, and rebuild "
                    "the starting structure with a larger minimum separation.",
             protocol_changes={"timestep_fs": 1.0},
             rebuild_structure=True,
             resume_from_stage="00_minimise",
             plain="atoms flew apart, smaller timestep and a cleaner start")),

        (r"ERROR.*Non-numeric (atom coords|pressure)|nan|NaN",
         lambda m: Remedy(
             code="NAN",
             what_happened="A coordinate or a pressure became not-a-number.",
             why="An energy or a force overflowed. Same root cause as lost atoms: a pair "
                 "far inside the repulsive wall. It can also mean a pair_coeff with the "
                 "wrong units, which makes the forces enormous from step zero.",
             action="Rebuild the structure, minimise, and verify the pair coefficients "
                    "against the force field table before resubmitting.",
             protocol_changes={"timestep_fs": 1.0},
             rebuild_structure=True,
             resume_from_stage="00_minimise",
             plain="energies overflowed, rebuild and check the coefficients")),

        (r"ERROR.*Out of range atoms.*comm|Out of range atoms - cannot compute",
         lambda m: Remedy(
             code="OUT_OF_RANGE",
             what_happened="Atoms moved further in one step than the communication "
                           "cutoff, so a processor could not find its neighbours.",
             why="Either the system is genuinely too hot for the timestep, or the "
                 "neighbour list is being rebuilt too rarely for how fast atoms move.",
             action="Rebuild the neighbour list every step, widen the skin, and reduce "
                    "the timestep.",
             protocol_changes={"timestep_fs": 1.0},
             plain="atoms outstripped the neighbour list")),

        (r"ERROR.*Neighbor list overflow|boost neigh_modify one",
         lambda m: Remedy(
             code="NEIGH_OVERFLOW",
             what_happened="The neighbour list ran out of room.",
             why="With an 11 Angstrom cutoff in a dense oxide melt each atom has several "
                 "hundred neighbours, and LAMMPS's default page size is not always "
                 "enough.",
             action="Raise neigh_modify one and page, and resubmit the stage.",
             # These two were missing until 29 Sep 2026. The action said it would raise
             # them and the remedy carried no protocol change, so the tool resubmitted the
             # stage unchanged and the same failure repeated to the attempt cap. Caught
             # while writing the manual, which had to document the broken promise.
             protocol_changes={"neigh_one": 8000, "neigh_page": 200000},
             plain="neighbour list too small, raising one and page")),

        (r"ERROR.*Domain too large for neighbor bins|Bad matrix inversion|"
         r"Non-numeric box dimensions",
         lambda m: Remedy(
             code="BOX_BLEW_UP",
             what_happened="The simulation box became an impossible size.",
             why="The barostat lost control. In NPT a barostat damping constant that is "
                 "too small fights the thermostat and the volume oscillates until it "
                 "diverges. It is the classic failure of an aggressive NPT start.",
             action="Increase the barostat damping so it responds more slowly, and "
                    "restart the stage.",
             protocol_changes={"hot_pdamp_fs": 800.0, "prod_pdamp_fs": 4000.0,
                               "_scale_eq_pdamp": 4.0},
             plain="barostat lost control, slowing it down")),

        (r"ERROR.*Ewald.*accuracy|KSpace accuracy must be|"
         r"Cannot use Ewald with 2d",
         lambda m: Remedy(
             code="KSPACE",
             what_happened="The Ewald sum could not reach the requested accuracy.",
             why="A 1e-6 accuracy with an 11 Angstrom real space cutoff needs a lot of "
                 "reciprocal vectors. On a small box the two requests can be "
                 "incompatible.",
             action="Loosen the Ewald accuracy one decade, to 1e-5, which is still tight "
                    "for a melt, and note it in the run record.",
             protocol_changes={"ewald_accuracy": 1.0e-5},
             plain="Ewald could not hit 1e-6, loosening to 1e-5")),

        (r"ERROR.*Cannot open|No such file or directory.*\.data|"
         r"Cannot open restart file",
         lambda m: Remedy(
             code="MISSING_INPUT",
             what_happened="A file the stage needed was not there.",
             why="Usually the previous stage never wrote its restart file, because it "
                 "died. Occasionally the run directory is on a filesystem the compute "
                 "node cannot see.",
             action="Go back and check the previous stage, do not resubmit this one.",
             auto_retry=False,
             plain="a required input file is missing, previous stage probably failed")),

        (r"ERROR.*(Incorrect args for pair coefficients|"
         r"All pair coeffs are not set|Pair style .* requires)",
         lambda m: Remedy(
             code="BAD_PAIR_COEFF",
             what_happened="LAMMPS rejected the pair coefficients.",
             why="The generated deck and the pair style disagree. born/coul/long takes "
                 "A rho sigma C D, buck/coul/long takes A rho C. A force field written "
                 "in the wrong form gives exactly this.",
             action="This is a bug in the generator, not a physics problem. Stop and fix "
                    "the force field form, do not resubmit.",
             auto_retry=False,
             plain="generated pair_coeff is malformed, this is a tool bug")),

        (r"slurmstepd.*(exceeded memory|oom-kill|Out Of Memory)|"
         r"Killed.*out of memory|MemoryError|std::bad_alloc",
         lambda m: Remedy(
             code="OOM",
             what_happened="The job was killed for using too much memory.",
             why="Ewald plus a long neighbour list plus per-atom dump buffers on a large "
                 "box adds up. Spreading the same atoms over more nodes gives each rank "
                 "less to hold.",
             action="Double the node count and resubmit the same stage.",
             cluster_changes={"nodes": "x2"},
             plain="ran out of memory, doubling the nodes")),

        (r"DUE TO TIME LIMIT|CANCELLED AT.*DUE TO TIME|"
         r"job killed.*walltime|=>> PBS: job killed: walltime",
         lambda m: Remedy(
             code="WALLTIME",
             what_happened="The job hit the wall clock limit and was killed part way "
                           "through.",
             why="This is the commonest failure and the least serious, because a restart "
                 "file means nothing is lost. A 50 ns production run of 5000 atoms is "
                 "not a 24 hour job on two nodes.",
             action="Double the walltime and resume from the last restart file. If the "
                    "queue has a maximum walltime, split the stage instead.",
             cluster_changes={"walltime": "x2"},
             plain="hit the walltime, doubling it and resuming")),

        (r"srun: error|mpirun noticed that process rank|"
         r"Signal: Segmentation fault|signal 11|BAD TERMINATION OF ONE OF YOUR",
         lambda m: Remedy(
             code="MPI_CRASH",
             what_happened="One MPI rank died and took the job with it.",
             why="Could be a node fault, which is not your problem and a resubmission "
                 "fixes it, or a real segfault in LAMMPS, which a resubmission will "
                 "reproduce. Try once, then stop.",
             action="Resubmit once unchanged. If it fails the same way twice, stop and "
                    "report the node list.",
             plain="an MPI rank died, one retry")),

        (r"ERROR.*Unknown (command|pair style|fix|compute)|"
         r"Invalid (fix|compute) style",
         lambda m: Remedy(
             code="MISSING_PACKAGE",
             what_happened="LAMMPS does not recognise something in the deck.",
             why="The LAMMPS build on this cluster was compiled without the package that "
                 "provides it. born/coul/long needs the KSPACE package.",
             action="Stop. Load a LAMMPS module built with KSPACE, or ask the cluster "
                    "support. Resubmitting the same deck will fail identically.",
             auto_retry=False,
             plain="this LAMMPS build lacks a needed package")),

        (r"command not found|module.*not found|Unable to locate a modulefile",
         lambda m: Remedy(
             code="ENVIRONMENT",
             what_happened="The batch script could not find LAMMPS or a module.",
             why="Module names differ on every cluster. Nothing about the physics is "
                 "wrong.",
             action="Stop and fix the module list and executable name in the cluster "
                    "settings.",
             auto_retry=False,
             plain="module or executable name is wrong for this cluster")),
    ]


TAIL_CHARS = 4000


def diagnose_text(text: str) -> Diagnosis:
    """Diagnose from a blob of log text."""
    tail = text[-TAIL_CHARS:]
    for pattern, build in _rules():
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            # keep the three lines around the match as evidence
            lines = text.splitlines()
            idx = 0
            for i, line in enumerate(lines):
                if re.search(pattern, line, re.IGNORECASE):
                    idx = i
                    break
            ev = "\n".join(lines[max(0, idx - 2): idx + 3])
            return Diagnosis(matched=True, remedy=build(m), evidence=ev, log_tail=tail)

    # an unrecognised end. Distinguish "died silently" from "actually finished".
    if "SLAGMD_STAGE_DONE" in text or "Total wall time" in text:
        return Diagnosis(matched=False, remedy=None,
                         evidence="the log contains a normal completion marker",
                         log_tail=tail)
    return Diagnosis(
        matched=False,
        remedy=Remedy(
            code="UNKNOWN",
            what_happened="The job stopped and the log has no recognised error and no "
                          "completion marker.",
            why="Something outside LAMMPS ended it, a node reboot, a filesystem outage, "
                "or an error message this tool has not been taught yet.",
            action="STOPPING. Read the last lines of the log yourself. The tool will not "
                   "resubmit a failure it does not understand, because that burns the "
                   "allocation and reproduces the crash.",
            auto_retry=False,
            plain="unrecognised failure, stopping rather than guessing"),
        evidence="no known error pattern and no completion marker",
        log_tail=tail)


def diagnose_job(rundir: Path, stage: str) -> Diagnosis:
    """Read every log a stage could have written and diagnose from all of them together."""
    rundir = Path(rundir)
    blobs: list[str] = []
    candidates = [
        rundir / f"log.{stage}", rundir / f"log.{stage}.slurm.err",
        rundir / f"log.{stage}.slurm.out", rundir / f"log.{stage}.pbs.err",
        rundir / f"log.{stage}.pbs.out",
        rundir / f"submit.{stage}.sh.out",      # where shell and module errors land
        rundir / "submit.sh.out", rundir / "log.lammps",
    ]
    for c in candidates:
        if c.exists():
            try:
                blobs.append(f"----- {c.name} -----\n" + c.read_text(errors="replace"))
            except OSError:
                pass
    if not blobs:
        return Diagnosis(
            matched=False,
            remedy=Remedy(
                code="NO_LOG",
                what_happened="No log file exists at all.",
                why="The job never started. The batch script may not have been submitted, "
                    "or it failed before LAMMPS was reached.",
                action="Check the scheduler's own output and the submit script.",
                auto_retry=False,
                plain="no log written, the job never ran"),
            evidence="no log files found in the run directory",
            log_tail="")
    return diagnose_text("\n\n".join(blobs))


def stage_completed(rundir: Path, stage: str) -> bool:
    """Did this stage actually finish? Checked by its own marker, not by the queue.

    The scheduler saying COMPLETED means the batch script exited zero. LAMMPS can exit
    zero having done almost nothing. The marker is printed by the input script itself as
    its last action, so it is the only trustworthy signal.
    """
    rundir = Path(rundir)
    for name in (f"log.{stage}", "log.lammps"):
        p = rundir / name
        if p.exists():
            try:
                if f"SLAGMD_STAGE_DONE {stage}" in p.read_text(errors="replace"):
                    return True
            except OSError:
                pass
    return False


def apply_cluster_change(cfg_dict: dict, changes: dict) -> dict:
    """Apply a remedy's cluster changes. Supports the 'x2' doubling form."""
    out = dict(cfg_dict)
    for key, val in changes.items():
        if val == "x2":
            cur = out.get(key)
            if key == "walltime" and isinstance(cur, str) and ":" in cur:
                parts = [int(x) for x in cur.split(":")]
                while len(parts) < 3:
                    parts.append(0)
                secs = parts[0] * 3600 + parts[1] * 60 + parts[2]
                secs = min(secs * 2, 7 * 24 * 3600)
                out[key] = f"{secs//3600}:{(secs%3600)//60:02d}:{secs%60:02d}"
            elif isinstance(cur, int):
                out[key] = cur * 2
        else:
            out[key] = val
    return out


MAX_ATTEMPTS = 3
