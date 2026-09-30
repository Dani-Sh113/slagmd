"""A fake LAMMPS, so the whole pipeline can be tested on a laptop with no LAMMPS.

It is not a simulator. It reads the generated input script well enough to know which stage
it is, writes a log in LAMMPS's own format with the completion marker, writes a restart
file, and for the production stage writes a trajectory of Brownian walkers with a diffusion
coefficient we choose. Because we choose it, the analysis can be tested against a known
answer instead of against a plausible looking plot.

It also has a FAIL_MODE environment variable that makes it emit each real LAMMPS error
message verbatim, which is how the diagnosis rules are tested. A diagnosis rule that has
never seen its own error message is a guess.

Usage, standing in for the real binary:
    python3 -m slagmd.mocklammps -in in.04_produce -log log.04_produce
    SLAGMD_MOCK_FAIL=LOST_ATOMS python3 -m slagmd.mocklammps -in in.01_randomise ...
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np

FAIL_MESSAGES = {
    "LOST_ATOMS": ("ERROR: Lost atoms: original 5002 current 4998 "
                   "(src/thermo.cpp:488)\nLast command: run             25000000"),
    "NAN": ("ERROR on proc 3: Non-numeric atom coords - simulation unstable "
            "(src/domain.cpp:210)"),
    "OUT_OF_RANGE": ("ERROR: Out of range atoms - cannot compute PPPM "
                     "(src/KSPACE/pppm.cpp:1927)"),
    "NEIGH_OVERFLOW": ("ERROR on proc 0: Neighbor list overflow, boost neigh_modify one "
                       "(src/npair_half_bin_newton.cpp:130)"),
    "BOX_BLEW_UP": ("ERROR: Domain too large for neighbor bins "
                    "(src/nbin_standard.cpp:97)"),
    "KSPACE": ("ERROR: KSpace accuracy must be > 0 (src/ewald.cpp:152)"),
    "MISSING_INPUT": ("ERROR: Cannot open restart file restart.02_cool: "
                      "No such file or directory (src/read_restart.cpp:135)"),
    "BAD_PAIR_COEFF": ("ERROR: Incorrect args for pair coefficients "
                       "(src/pair_born_coul_long.cpp:401)"),
    "OOM": ("slurmstepd: error: Detected 1 oom-kill event(s) in StepId=8812345.batch. "
            "Some of your processes may have been killed by the cgroup out-of-memory "
            "handler."),
    "WALLTIME": ("slurmstepd: error: *** JOB 8812345 ON node042 CANCELLED AT "
                 "2026-09-29T14:02:11 DUE TO TIME LIMIT ***"),
    "MPI_CRASH": ("srun: error: node071: task 44: Segmentation fault\n"
                  "srun: error: Terminating job step 8812345.0"),
    "MISSING_PACKAGE": ("ERROR: Unknown pair style born/coul/long "
                        "(src/force.cpp:271)\nPair style born/coul/long requires the "
                        "KSPACE package"),
    "ENVIRONMENT": ("Lmod has detected the following error: Unable to locate a modulefile "
                    "for 'LAMMPS'\nlmp: command not found"),
    "SILENT": "",
}


def _read_data(path: Path):
    """Pull atom count, box and per atom type out of a data file we wrote ourselves."""
    txt = path.read_text().splitlines()
    n = 0
    box = 40.0
    masses: dict[int, str] = {}
    types: list[int] = []
    section = ""
    for line in txt:
        s = line.strip()
        if s.endswith("atoms") and s.split()[0].isdigit():
            n = int(s.split()[0])
        if "xlo xhi" in s:
            lo, hi = float(s.split()[0]), float(s.split()[1])
            box = hi - lo
        if s.startswith(("Masses", "Atoms")):
            section = s.split()[0]
            continue
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if section == "Masses" and len(parts) >= 2 and parts[0].isdigit():
            comment = line.split("#")[1].strip() if "#" in line else f"T{parts[0]}"
            masses[int(parts[0])] = comment
        elif section == "Atoms" and len(parts) >= 6:
            types.append(int(parts[1]))
    return n, box, masses, types


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-in", dest="infile", required=True)
    ap.add_argument("-log", dest="logfile", default="log.lammps")
    ap.add_argument("--D", type=float, default=None,
                    help="diffusion coefficient in m^2/s for the fake production run")
    args, _unknown = ap.parse_known_args(argv)

    infile = Path(args.infile)
    logf = Path(args.logfile)
    stage = infile.name.replace("in.", "")
    text = infile.read_text()

    fail = os.environ.get("SLAGMD_MOCK_FAIL", "")
    fail_stage = os.environ.get("SLAGMD_MOCK_FAIL_STAGE", "")
    do_fail = bool(fail) and (not fail_stage or fail_stage == stage)
    # SLAGMD_MOCK_FAIL_ONCE makes the fault transient, so the automatic recovery path can
    # be tested for real: fail on the first attempt, succeed on the resubmission. A
    # permanent fault only ever tests the attempt cap.
    if do_fail and os.environ.get("SLAGMD_MOCK_FAIL_ONCE"):
        sentinel = infile.parent / f".mock_failed_{stage}"
        if sentinel.exists():
            do_fail = False
        else:
            sentinel.write_text("failed once\n")

    lines = ["LAMMPS (27 Jun 2024)  [SLAGMD MOCK, NOT A REAL SIMULATION]",
             f"reading input {infile.name}"]

    # find the run length and the dump settings from the deck
    m = re.search(r"^run\s+(\d+)", text, re.M)
    nsteps = int(m.group(1)) if m else 1000
    m = re.search(r"^timestep\s+([\d.]+)", text, re.M)
    dt_fs = float(m.group(1)) if m else 2.0
    m = re.search(r"^dump\s+\w+\s+all\s+custom\s+(\d+)\s+(\S+)", text, re.M)
    dump_every, dump_name = (int(m.group(1)), m.group(2)) if m else (0, "")
    m = re.search(r"dump_modify\s+\w+\s+element\s+(.+)", text)
    elements = m.group(1).split() if m else []
    m = re.search(r"^read_data\s+(\S+)", text, re.M)
    data_name = m.group(1) if m else "melt.data"
    m = re.search(r"npt\s+temp\s+([\d.]+)\s+([\d.]+)", text)
    T0, T1 = (float(m.group(1)), float(m.group(2))) if m else (1873.0, 1873.0)

    workdir = infile.parent
    datafile = workdir / data_name
    if datafile.exists():
        natoms, box, masses, types = _read_data(datafile)
    else:
        # a later stage reads a restart; recover the shape from the stashed json
        import json
        meta = workdir / ".mock_state.json"
        if meta.exists():
            d = json.loads(meta.read_text())
            natoms, box, types = d["natoms"], d["box"], d["types"]
            masses = {int(k): v for k, v in d["masses"].items()}
        else:
            natoms, box, masses, types = 500, 25.0, {1: "O"}, [1] * 500

    import json
    (workdir / ".mock_state.json").write_text(json.dumps(
        {"natoms": natoms, "box": box, "types": types,
         "masses": {str(k): v for k, v in masses.items()}}))

    if do_fail:
        lines.append("Setting up run ...")
        lines.append(FAIL_MESSAGES.get(fail, "ERROR: something unrecognised happened"))
        logf.write_text("\n".join(lines) + "\n")
        sys.stderr.write(FAIL_MESSAGES.get(fail, "") + "\n")
        return 1

    # ------------------------------------------------- a plausible thermo table
    rng = np.random.default_rng(abs(hash(stage)) % (2**31))
    nthermo = 60
    lines += ["", "Step Time Temp Press Volume Density PotEng KinEng TotEng Enthalpy Lx"]
    vol0 = box ** 3
    for i in range(nthermo):
        f = i / max(1, nthermo - 1)
        T = T0 + (T1 - T0) * f + rng.normal(0, 8)
        vol = vol0 * (1.0 + 0.002 * (1 - f) + rng.normal(0, 0.0004))
        dens = 2.65 * vol0 / vol
        pe = -1.35e5 + rng.normal(0, 60)
        ke = 1.5 * natoms * 0.0019872 * T
        lines.append(
            f"{int(nsteps*f)} {nsteps*f*dt_fs/1000:.3f} {T:.2f} {rng.normal(1,300):.2f} "
            f"{vol:.2f} {dens:.5f} {pe:.2f} {ke:.2f} {pe+ke:.2f} {pe+ke+1:.2f} "
            f"{vol**(1/3):.4f}")
    lines += ["", f"Loop time of {nsteps*1e-4:.4f} on 256 procs for {nsteps} steps",
              "", f"Total wall time: 0:{int(nsteps*1e-4)//60:02d}:{int(nsteps*1e-4)%60:02d}"]

    # ------------------------------------------- a trajectory with a KNOWN D and RDF
    if dump_every and dump_name:
        D_m2_s = args.D if args.D is not None else 1.0e-9     # a typical melt value
        D_A2_ps = D_m2_s * 1.0e8                              # m^2/s -> A^2/ps
        dt_ps = dump_every * dt_fs / 1000.0
        nframes = max(2, min(400, nsteps // dump_every + 1))
        sigma = np.sqrt(2.0 * D_A2_ps * dt_ps)                # per axis, per frame

        el = [elements[t - 1] if 0 < t <= len(elements)
              else masses.get(t, f"T{t}") for t in types]
        pos = rng.random((natoms, 3)) * box
        out = []
        for fr in range(nframes):
            if fr:
                pos = pos + rng.normal(0.0, sigma, size=(natoms, 3))
            out += ["ITEM: TIMESTEP", str(fr * dump_every),
                    "ITEM: NUMBER OF ATOMS", str(natoms),
                    "ITEM: BOX BOUNDS pp pp pp",
                    f"0.0 {box:.6f}", f"0.0 {box:.6f}", f"0.0 {box:.6f}",
                    "ITEM: ATOMS id type element xu yu zu"]
            for i in range(natoms):
                out.append(f"{i+1} {types[i]} {el[i]} "
                           f"{pos[i,0]:.4f} {pos[i,1]:.4f} {pos[i,2]:.4f}")
        (workdir / dump_name).write_text("\n".join(out) + "\n")
        lines.append(f"# MOCK: wrote {dump_name}, {nframes} frames, "
                     f"true D = {D_m2_s:.4e} m^2/s per species")

    # stress series for the viscosity path, white noise so the true answer is known
    if "stressout" in text:
        nst = 200000
        w = rng.normal(0.0, 2000.0, size=(nst, 3))            # atm
        rows = ["# step pxy pxz pyz pxx pyy pzz (atm)"]
        for i in range(0, nst, 1):
            rows.append(f"{i*5} {w[i,0]:.3f} {w[i,1]:.3f} {w[i,2]:.3f} 0 0 0")
        (workdir / "stress.txt").write_text("\n".join(rows) + "\n")
        lines.append("# MOCK: wrote stress.txt as white noise, so Green-Kubo must give "
                     "a viscosity that does NOT plateau. That is the point.")

    for rst in re.findall(r"^write_restart\s+(\S+)", text, re.M):
        (workdir / rst).write_text("mock restart\n")
    for mark in re.findall(r'print\s+"SLAGMD_STAGE_DONE (\S+)"', text):
        lines.append(f"SLAGMD_STAGE_DONE {mark}")

    logf.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
