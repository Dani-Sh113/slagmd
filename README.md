# slagmd

**Build, submit, watch, repair and analyse molecular dynamics of molten oxide slags.**

`slagmd` takes a slag composition in mole fractions and carries it all the way to diffusion
coefficients, viscosity and structural speciation: it resolves the composition to an exactly
charge-neutral set of atom counts, packs a starting configuration, writes a multi-stage LAMMPS
deck, generates and submits the scheduler script, watches the job, diagnoses it when it dies,
repairs and resubmits it, and then analyses the trajectory.

It is built around one idea: **a simulation tool should refuse to give you a number it cannot
stand behind.** Most of the design effort in this repository went into the refusals, which are
listed in full below.

---

## Why it exists

Two published studies of CaO-Al<sub>2</sub>O<sub>3</sub>-SiO<sub>2</sub> melts released their
complete LAMMPS input files, 84 files in total. Those files were read end to end, and they
contain **no job submission layer at all**: no `#SBATCH`, no `#PBS`, no `module load`, no `srun`,
no `qsub`, no job arrays, no Makefile and no workflow manager. The whole launch layer is a command
line typed by hand, for a study that amounts to roughly 754 chained jobs. Building the starting
configuration is a six-step manual procedure in a graphical program, including a trial-and-error
loop.

**That gap is what this fills.** The physics is theirs and is credited below. The orchestration,
the validation, the failure diagnosis and the analysis are this repository.

---

## Install

```bash
git clone https://github.com/Dani-Sh113/slagmd.git
cd slagmd
pip install -r requirements.txt
```

Needs **Python 3.10 or later**, plus `numpy`, `matplotlib` and `tkinter`.

On Debian and Ubuntu, `tkinter` is a separate system package:

```bash
sudo apt install python3-tk
```

**LAMMPS is not required** to build inputs, to run the test suite, or to analyse a trajectory.
It is only needed to actually run a simulation.

---

## Run it

```bash
python3 run_slagmd.py                  # open the window
python3 run_slagmd.py /path/to/project # and open a project straight away
python3 tests/test_slagmd.py           # the self-audit, 221 checks
```

---

## The six tabs

| tab | what it does |
|:--|:--|
| **Project** | opens or creates a project folder, and reports what this machine has |
| **Cluster** | scheduler, nodes, walltime, modules and executable; previews the batch script it will generate and checks the settings for problems |
| **Compositions** | mole fractions to exact charge-neutral atom counts, cost in core hours, and queueing. One button queues ten published benchmark compositions |
| **Queue** | watches every job, and when one dies it diagnoses and repairs it |
| **Analysis** | radial distribution functions, Q<sup>n</sup> and oxygen speciation, diffusion, ionic conductivity and viscosity from a finished trajectory, with the figures |
| **Log** | everything that happened |

---

## What it refuses to do

**This is the important part, and it is deliberate.** Each of these is a silent-wrong-answer
failure mode that the tool turns into a loud refusal.

- It **refuses to build a deck from unverified force-field parameters.** Every parameter set
  carries its printed source, and one marked `verified: false` is rejected outright.
- It **refuses a species the force field has no parameters for**, rather than defaulting to
  zero. LAMMPS will happily run a melt held together by nothing but Coulomb and tell you nothing
  is wrong.
- It **refuses a box that is not more than twice the Coulomb cutoff**, where an atom would
  interact with its own periodic image, and warns between two and two and a half times.
- It **refuses to call a fitted slope a diffusion coefficient** when the log-log slope of the
  mean squared displacement is not within 0.15 of 1. A slope that is not 1 is not diffusion.
- It **refuses to quote a Green-Kubo viscosity** when the running integral has not plateaued,
  and reports the scatter instead.
- It **refuses to resubmit a failure it does not recognise.** A retry loop that does not
  understand the failure burns the allocation and faithfully reproduces the crash.
- It **never rounds oxygen.** Oxygen is counted from the cations, so the box is exactly neutral
  rather than nearly neutral.

---

## Crash proofing

State is written to `project.json` after every change, through a temporary file and an atomic
rename, so a crash cannot corrupt it. The scheduler's job number is recorded **before** anything
else happens, so a session dying between submitting and recording cannot orphan a job. Reopening
the project folder resumes exactly where it stopped.

**A stage counts as finished only when the input script's own marker appears in its log**, never
because the scheduler reported COMPLETED. LAMMPS can exit zero having done almost nothing.

---

## How it is verified

`python3 tests/test_slagmd.py` runs **221 checks** and needs no network, no cluster and no LAMMPS.
It ships a mock LAMMPS that reproduces every failure mode the diagnoser claims to recognise.

**The analysis is checked against known answers, not against itself:**

| what is checked | against what |
|:--|:--|
| diffusion coefficient | a synthetic Brownian walk whose true `D` is known exactly |
| radial distribution function and coordination number | a simple cubic lattice, whose first peak must be at 3.0 Å with coordination 6 |
| mean squared displacement by FFT | the direct double loop, agreeing to 5 parts in 10<sup>14</sup> |
| log-log diffusive-regime test | a ballistic trajectory, whose slope must be 2, and a diffusive one, whose slope must be 1 |
| viscosity | a stress signal with an analytically known integral |
| charge neutrality | every composition, exactly, with no rounding of oxygen |

**Two independent reproductions of published claims are also in the test set.** A paper's stated
result that one potential's cation-oxygen wells are "about 66 % and 57 %" deeper than another's
reproduces from the raw parameters as 65.8 %, 56.6 % and 56.6 %; and the parameters transcribed
here from the printed tables match the authors' own released `pair_coeff` lines on all forty
numbers to better than one part in ten thousand.

---

## Honest limitations

**Read this before trusting any number it produces.**

1. **It has never been run against a real LAMMPS binary.** Every test uses the bundled mock. The
   acceptance test it still needs is one benchmark composition run for real and compared against
   published diffusion coefficients.
2. **The damping schedule between stages is reconstructed from the authors' released input files,
   not from either paper.** Both papers state that thermostat and barostat damping were changed
   between stages and place the schedule only in supplementary material.
3. **Viscosity is computed nowhere in the released source material.** The Green-Kubo and
   Einstein-Helfand implementations here are checked against synthetic signals with known answers,
   but they have no published slag value to be validated against.
4. **The force fields cover Ca, Al, Si and O only.** Iron and titanium are not parameterised, and
   the tool will refuse a composition containing them rather than guess.
5. **`[MINE]` settings are defaults, not published values.** Every line a generated deck writes is
   tagged `[THEIRS]`, `[PAPER]` or `[MINE]` precisely so this can never be confused.

---

## Layout

```
slagmd/
  forcefields.json   verified parameters, with provenance and self-consistency identities
  composition.py     mole fractions to exact charge-neutral atom counts
  forcefield.py      loading, validation, LAMMPS pair coefficients
  builder.py         random packing with a minimum separation, and the data file
  inputgen.py        the multi-stage deck, with the reconstructed damping cascade
  scheduler.py       SLURM, PBS and local submission and status
  diagnose.py        14 real failure patterns, each with a cause and a remedy
  analysis.py        RDF, speciation, diffusion, conductivity, viscosity
  project.py         the resumable state machine
  mocklammps.py      a fake LAMMPS for testing, including every failure mode
  gui.py             the window
tests/test_slagmd.py the self-audit
docs/MANUAL.md       the full user manual, A to Z
```

---

## Defects found in the published input files

Recorded here so that nobody has to rediscover them. All are handled by this tool.

1. `production.file` ends with `unfix` calls for fixes it never defines, so LAMMPS aborts after
   the physics has been written. Cosmetic for a single run, fatal for a dependency chain. This
   tool defines every fix it unfixes.
2. `MSCD.py` hard-codes one charge set while offering another in its header.
3. The Packmol example is not charge neutral: 295 oxygens where the stoichiometry needs 290.
4. Four stale run-length comments. The step counts themselves are correct.

---

## Provenance and citation

**The physics is not mine.** The interatomic potentials, the simulation protocol and the benchmark
compositions come from published work, and every number in `forcefields.json` records the printed
table it was transcribed from.

If you use this tool, **cite the original science**:

- M. Matsui, *Mineralogical Magazine* **58A** (1994) 571–572.
- B. Guillot and N. Sator, *Geochimica et Cosmochimica Acta* **71** (2007) 1249–1265.
- M. Bouhadja, N. Jakse and A. Pasturel, *Journal of Chemical Physics* **138** (2013).
- H. Verbeeck, I. Bellemans et al., "Molecular dynamics studies of the structural and transport
  properties of CaO-Al<sub>2</sub>O<sub>3</sub>-SiO<sub>2</sub> melts: comparing the accuracy of
  the empirical force fields", *Computational Materials Science* **259** (2025) 114150.
- H. Verbeeck et al., *Materials & Design* **264** (2026) 115709.

The authors' own input files are published as open data at Mendeley Data,
`10.17632/v7m3k6yksw.1` and `10.17632/gkcjjk744s.1`.

**Those files are not redistributed in this repository.** What is here is the set of numerical
parameter values, transcribed from published tables and cross-checked against the released
`pair_coeff` lines, together with a citation for each. Download the originals from Mendeley Data
if you need them, and observe whatever licence they carry.

---

## Licence

The code in this repository is released under the **MIT Licence**. See [LICENSE](LICENSE).

**The MIT licence covers this code only.** It does not and cannot cover the published potential
parameters, the protocol or the experimental data they were fitted to, which belong to their
authors and are cited above.
