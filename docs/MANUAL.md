---
title: "slagmd: user manual"
subtitle: "Build, submit, watch, repair and analyse molecular dynamics of molten oxide slags"
author: "Danish Ilyas"
date: "29 September 2026, version 1.0.0"
---

# Contents

| Section | What is in it |
|---|---|
| How to read this manual | Entry points, and the two conventions used throughout |
| 1. What this tool is for | The gap in the published release, and what slagmd replaces |
| 2. Install and first run | Requirements, the audit suite, what a project folder is |
| 3. The five minute tour | A complete worked example on a laptop, with no cluster and no LAMMPS |
| 4. Every tab, field by field | Project, Cluster, Compositions, Queue, Analysis, Log |
| 5. The concepts the tool is built on | Force fields, oxygen counting, box size, the protocol, the damping cascade |
| 6. Every refusal, and what to do about it | Every error and warning the tool can raise, with its fix |
| 7. The failure diagnosis system | All the failure patterns, and why it will not guess |
| 8. The analysis, and how to read the output | Structure, speciation, diffusion, conductivity, viscosity, equilibration |
| 9. Running on a real cluster | A SLURM example end to end, and what changes elsewhere |
| 10. Extending it | New force fields, new oxides, new failure patterns |
| 11. Troubleshooting | Everything that is neither a refusal nor a job failure |
| 12. Glossary | Every term used |
| Provenance | What is theirs, what is ours, and the known defects on both sides |

# How to read this manual

You do not need to read it all. Pick your entry point.

| If you are | Read |
|---|---|
| Trying it for the first time on a laptop | Sections 2 and 3, in that order, then stop |
| Setting it up on a cluster | Sections 2, 4.2 and 9 |
| Wondering why it refused to do something | Section 6, then 5 |
| Looking at a number it printed and unsure whether to trust it | Section 8 |
| A job died and you want to know what happened | Section 7 |
| Adding a force field, an oxide, or a new error pattern | Section 10 |
| Lost on a word | Section 12, the glossary |

Two conventions run through the whole document.

**Provenance labels.** Everything the tool generates is labelled `[THEIRS]`, `[PAPER]` or
`[MINE]`. `[THEIRS]` means the setting was copied from an input file that Verbeeck and
Bellemans published. `[PAPER]` means it is stated in the printed paper. `[MINE]` means it
is slagmd's own addition and must never be quoted as theirs. The same distinction is
carried through this manual: where behaviour comes from the published study it is said so
and the source is named.

**Refusals.** slagmd stops rather than guessing more often than most tools do. Every one of
those stops is listed in section 6 with its cause and its fix. A refusal is not a bug.

# 1. What this tool is for

## 1.1 The problem

Héléna Verbeeck, Inge Bellemans and co-workers published two papers on molecular dynamics of
CaO-Al2O3-SiO2 melts, and with each paper they released their simulation scripts on Mendeley
Data. That release is unusually generous. It contains the LAMMPS input files for all three
force fields, the Packmol inputs, the post-processing scripts, and the VASP inputs for their
ab initio comparison. Eighty four files in total.

What it does not contain is any way to run them at scale.

Every one of those 84 files was searched for `#SBATCH`, `#PBS`, `#$ -`, `sbatch`, `qsub`,
`squeue`, `module load`, `srun`, `mpirun`, `mpiexec`, `walltime`, `partition`, job arrays,
Makefiles and workflow managers such as Snakemake, Nextflow, ASE or AiiDA. Zero hits. There
is no scheduler script, no job array, no dependency logic, no retry, no module list, and no
record of how many cores any published run used. The entire launch layer in the public
record is one command line, quoted here exactly as their own README gives it:

```
mympirun lmp -in equilibration_step1.file -var temp 1873 -var input_data cool_1873.out
```

That single line has to be typed once per stage, per temperature, per composition. Their
high throughput study covers 40 compositions at 5 temperatures, trimmed by a liquidus
screen to 119 usable pairs. Cooling is shared across the five temperatures, so the job count
is 40 cooling jobs plus 119 chains of 6 stages each, which is **754 chained jobs**. Each
chain is strictly sequential: equilibration step 2 cannot start until step 1 has written its
output file.

Building the starting structure is worse, because it is not scripted at all. Their own
Packmol README describes it as a sequence of manual actions in a graphical program:

> Convert the CAT+Ox.xyz into a POSCAR file. Reshape the cell to obtain a cubic system.
> Convert the POSCAR file into a LAMMPS data file. Add 2 Angström to each lattice side to
> avoid overlapping atoms. Check the rdf (in Ovito) to confirm that there are no overlapping
> atoms. Adjust the lattice size in the LAMMPS datafile if needed untul there is no overlap
> anymore.

Six steps in OVITO, the last one a trial and error loop, standing between a composition and
a file LAMMPS can read. No script exists for any of it. The data file headers confirm the
tool: `# LAMMPS data file written by OVITO Basic 3.8.5`.

There are three further gaps that cost real time:

1. **Diffusion coefficients are transcribed by hand.** Their analysis script prints the
   numbers to the screen and draws a plot. Their README then asks the user to store the
   values in a text file themselves.
2. **Nothing checks that a stage finished.** There is no exit code handling, no check for a
   half written output file, and no `read_restart` anywhere in the release, so resuming a
   killed run is a manual act.
3. **Viscosity is never computed from the simulations at all.** Searching all 84 files for
   `viscos`, `green.?kubo`, `stress` and `ave/correlate` gives no match. Where viscosity
   appears in their papers it is a FactSage screening number, not an MD result.

## 1.2 What slagmd replaces

| Step in their workflow | How it is done in the release | How slagmd does it |
|---|---|---|
| Choose atom counts for a composition | By hand, and the released Packmol example is not charge neutral, 295 oxygens where the stoichiometry needs 290 | `composition.py` resolves mole fractions to exact integer counts, and counts oxygen from the cations so the box is exactly neutral |
| Build the starting structure | Packmol in two passes, then six manual OVITO steps with a trial and error loop | `builder.py` packs it directly with a minimum separation criterion, or calls Packmol when it is installed, and checks the result before anything is submitted |
| Write the LAMMPS deck | Seven hand edited files per composition and force field, each re-declaring the force field | `inputgen.py` generates every stage from one protocol object, with the provenance of each line marked |
| Submit | One command line typed by hand, 754 times | `scheduler.py` writes a batch script and submits it, on SLURM, PBS or the local machine |
| Know whether a stage worked | Nothing | The input script prints its own completion marker, and only that marker counts as finished |
| Handle a failure | Nothing | `diagnose.py` reads the log, names the failure, applies a specific fix and resubmits, or stops if it does not recognise it |
| Get the numbers out | Printed to the screen, then typed into a text file | `analysis.py` computes them and reports, beside every number, whether that number is allowed to exist |
| Viscosity | Not computed | Green-Kubo and Einstein-Helfand, from four extra lines in a production run that is happening anyway |

The protocol and the force fields are **theirs**. The three force fields were transcribed
from the printed parameter tables and then checked, number by number, against the
`pair_coeff` lines in their released input files: 40 numbers, agreement better than one part
in ten thousand. Their five stage protocol, including the thermostat and barostat damping
cascade that appears in neither paper, was recovered from the released scripts and is
reproduced exactly.

The submission layer, the structure builder, the failure diagnosis, the error bars, the
regime tests and the viscosity are **ours**.

## 1.3 What it is not

- It is not a molecular dynamics engine. LAMMPS does the physics.
- It is not a force field fitting tool. It refuses to use a parameter set that does not
  carry a printed source.
- It does not know about liquidus temperatures. Their 200 to 119 reduction came from
  FactSage, which is commercial software, and no script for it was released. slagmd warns
  when a target temperature looks low and leaves the judgement to you.
- It has never been run against a real LAMMPS. Everything in it is verified against a mock
  LAMMPS and against analytic ground truths where the right answer is known exactly. The
  first real use should be one of the ten published benchmark compositions, so the answer
  can be compared with a published one.

# 2. Install and first run

## 2.1 What you need

| Requirement | Why | What happens without it |
|---|---|---|
| Python 3.10 or later | The code uses modern type syntax | It will not import |
| numpy | Every array operation, the FFT based displacement routine | It will not import |
| matplotlib | The two figures on the Analysis tab | Everything works except the `Make the figures` button |
| tkinter | The window itself | The package still imports and can be driven from Python, but `run_slagmd.py` fails |
| LAMMPS | Running the simulations | **Everything else still works.** You can build inputs, check them, generate batch scripts, and analyse trajectories brought back from a cluster. You just cannot run |
| Packmol | Optional. Their own way of building the starting structure | slagmd uses its own packer instead and says so |

tkinter is part of the Python standard library but some Linux distributions package it
separately. On Debian or Ubuntu it is `python3-tk`. On Anaconda it is already there.

## 2.2 Install

There is nothing to install. The tool is a folder of plain Python files with no build step,
no `pip install`, and no virtual environment of its own.

```bash
cd ~/slagmd
python3 -c "import numpy, matplotlib, tkinter; print('dependencies ok')"
```

Quote the path. This folder lives under a directory whose name ends in a space, and an
unquoted path with a space in it does not fail loudly, it splits into two words several
lines away from the cause. That trap cost a debugging cycle on 29 September 2026 and it is
now one of the things the tool checks for you, in `scheduler.config_problems`.

## 2.3 Run the audit first

```bash
python3 tests/test_slagmd.py
```

This is the whole confidence in the tool, so run it before you trust any number it gives
you. It takes a couple of minutes and prints one line per check. The last three lines look
like this:

```
========================================================================
221 checks passed, 0 failed
========================================================================
```

The exit code is 0 when everything passed and 1 otherwise, so it can be chained.

### What "221 checks passed" actually means

It is not 209 unit tests of the "does this function return a number" kind. The checks are
grouped by what they protect against, and most of them compare against something whose
answer is known independently.

| Group | What it checks | Examples of the ground truth used |
|---|---|---|
| A, force field tables | That the transcribed parameters are not mistyped | Matsui's own functional form implies A equals 4.184 times rho for all ten pairs. Guillot's cation-oxygen A must equal Matsui's A times exp(sigma over rho), which is algebra, not an assumption. The paper's claim that Bouhadja's Ca-O well is about 66 percent deeper is recomputed from the raw parameters. All 40 Bouhadja numbers are compared against the `pair_coeff` lines in their released `cooling.file` |
| B, composition | Charge neutrality and exact stoichiometry | The ten published benchmark compositions, each resolved to atom counts, each neutral to better than 1e-9, each within 0.2 percent of 5000 atoms and 0.001 of the requested mole fractions |
| C, structure | That a starting configuration is fit to run | A deliberate 0.3 Angström overlap must be caught. A box too small for the Coulomb cutoff must be refused. A net charge must be caught |
| D, the deck | That the generated LAMMPS input reproduces their protocol | Their arithmetic: 2,500,000 steps is 5 ns, 1,063,500 steps is 2.127 ns at exactly 1 K/ps, 25,000,000 steps is 50 ns, five equilibration blocks give 50 ns, and the cost estimate reproduces their own 2560 core hours |
| E, analysis | Every property against an exact answer | A Brownian walk with a chosen D, because the mean squared displacement of a free walk is exactly 6Dt. A simple cubic lattice, whose first RDF peak is at the lattice spacing and whose coordination number is 6. An ideal gas, whose g(r) is 1 everywhere. White noise and an Ornstein-Uhlenbeck process, whose Green-Kubo viscosities are known in closed form. A caged oscillator, which must be **refused** as a diffusion coefficient |
| F, diagnosis | Every failure rule against the real error text | Each rule is tested against the verbatim LAMMPS or SLURM message, and an unrecognised failure must refuse to retry |
| G, the pipeline | The whole thing end to end, twice with a failure | A run is taken through all its stages under a mock LAMMPS, the diffusion coefficient is recovered from the written trajectory, then a transient walltime failure is injected and the tool must diagnose it, double the walltime, resubmit and finish, with no human involved |
| H, the window | The interface driven as a person would | A run added through the form, prepared, executed, analysed, and the figures checked to be real images rather than empty files |

If any check fails, the failing names are listed at the end. Do not use the tool until they
are fixed.

## 2.4 Open the window

```bash
python3 run_slagmd.py                       # open with no project
python3 run_slagmd.py /path/to/project      # open a project folder straight away
```

Those are the only two forms. There is no other command line flag.

## 2.5 What a project folder is

One folder holds everything about one study.

```
myproject/
  project.json          all the state: cluster settings, protocol, every run and stage
  events.log            an append only record of what happened and when
  REPORT.txt            a plain text summary, rewritten on every poll
  runs/
    S1_1600C/
      melt.data         the starting structure, LAMMPS data file
      in.00_minimise    one input script per stage
      in.01_hot
      ...
      submit.00_minimise.sh    one batch script per stage
      ...
      protocol.json     the full manifest: composition, force field, charges, cost estimate
      log.08_produce    LAMMPS logs
      prod.lammpstrj    the production trajectory
      stress.txt        the stress series, if viscosity sampling was on
      diagnosis.*.txt   one file per failure, written for a human to read
      ANALYSIS.txt      whatever the Analysis tab last printed, if you saved it
      figures/          rdf.png and msd_loglog.png
```

`project.json` is written after every single change, through a temporary file and a rename,
so a crash in the middle of a write cannot corrupt it. Closing the window loses nothing.
Reopening the same folder resumes at exactly the same place. There is no state that lives
only in memory.

# 3. The five minute tour

This section takes you from an empty folder to a diffusion coefficient with an error bar, on
a laptop, with no cluster and no LAMMPS. It works because slagmd ships with a mock LAMMPS: a
small program that reads the generated input script, writes a log in LAMMPS's own format, and
for the production stage writes a trajectory of random walkers whose diffusion coefficient is
**chosen in advance**. Because it is chosen, you can check that the analysis gets it back.

The mock's default is D = 1.0e-9 m^2/s for every species. Remember that number. It is the
answer you are looking for at the end.

Follow this exactly. It takes about five minutes of waiting in total, most of it in the last
step.

## Step 1. Open the window

```bash
cd ~/slagmd
python3 run_slagmd.py
```

You get a window with six numbered tabs along the top and a header that says what the tool
is. The status bar at the bottom says `ready`.

## Step 2. Make a project (tab 1, Project)

1. Type a folder path into the **Project folder** box, for example `/tmp/slagdemo`. It does
   not have to exist. Or press **Browse** and pick one.
2. Press **Open or create**.

The header line now reads `project /tmp/slagdemo   0 runs`.

Below that, the box titled **what this machine has** already told you, before you did
anything, what it found: the scheduler it detected, whether LAMMPS is on your PATH, whether
Packmol is, and how many cores this machine has. If LAMMPS is not there it says so in plain
words and explains that building inputs and analysing trajectories still work.

The bottom box lists the three force fields, each with its label, its charge set, its printed
source, and the result of its own consistency check. All three should say `PASSES`.

## Step 3. Point it at the mock LAMMPS (tab 2, Cluster)

Set these six fields. Leave everything else alone.

| Field | Value |
|---|---|
| scheduler | `local` |
| nodes | `1` |
| tasks per node | `1` |
| LAMMPS executable | `python3 -m slagmd.mocklammps` |
| launcher | leave **blank**, pick the empty entry in the dropdown |
| modules to load, one per line | leave **empty** |

Then in the box labelled **extra shell lines before LAMMPS**, put one line, with the quotes:

```
export PYTHONPATH=~/slagmd
```

The quotes are not optional. That path contains a space, and without them bash splits it into
two words.

Press **Save**. The warning line below should read `no problems found in these settings` in
green.

Now press **Preview the submit script**. A window opens showing exactly the shell script that
will be written. Read it. It is short, and every line in it is one you can understand. This
is the layer their release does not have.

## Step 4. Add a composition (tab 3, Compositions)

The fields come prefilled with the first of the ten published benchmark compositions. Change
one thing only:

| Field | Value |
|---|---|
| run name | `S1_1600C` |
| composition, oxide then mole fraction | `CaO 0.453 Al2O3 0.125 SiO2 0.422` |
| temperature in C | `1600` |
| target atoms | `1400`  (change this, the default is 5000) |
| initial density g/cm3 | `2.6` |
| random seed | `458165` |
| force field | `bouhadja` |
| sample the stress tensor | leave **ticked** |
| add a short NVE energy conservation check | leave unticked |

1400 atoms instead of 5000 is purely so the tour finishes quickly. For real work use 5000,
which is what the published study used.

Press **Check this composition**. Nothing is submitted, nothing is written. You get:

```
requested  Al2O3 0.1250  CaO 0.4530  SiO2 0.4220
achieved   Al2O3 0.1257  CaO 0.4531  SiO2 0.4212   (max error 0.00084)
formula units  Al2O3x63  CaOx227  SiO2x211
atoms          Si 211  Al 126  Ca 227  O 838   total 1402
ideal NBO/T from stoichiometry  0.973

net box charge under bouhadja: +5.684e-14
charges are a uniform 0.6000 of the formal charges: True
cubic box at 2.6 g/cm3: 27.29 A
molar mass of the box: 31830.5 g

protocol: 107.13 ns total, 50 ns equilibration in 5 blocks, 50 ns production
estimated cost: 718 core hours
```

Read that carefully, because it is the whole philosophy of the tool in one panel.

- You asked for 1400 atoms and got **1402**, because the mole fractions were resolved to whole
  formula units and oxygen was then counted from the cations. It was never rounded.
- The net box charge is 5.7e-14, which is zero to the last bit of floating point
  precision on a sum of 1402 terms. Not "small", not "close enough".
- The charges are reported as a uniform 0.6 of the formal charges, which is a stronger
  statement than checking this one box: it means every stoichiometric mixture with this force
  field is exactly neutral.
- It tells you, before you spend anything, that this will cost about 718 core hours.

Press **Add run**. Then press **Prepare every unprepared run**. The panel fills with the
preparation record:

```
force field Bouhadja, Jakse and Pasturel (2013), BMH, fitted to AIMD of low-silica
CaO-Al2O3-SiO2 melts passed its consistency identities
charges are a uniform 0.6000 of the formal charges, so any stoichiometric composition is
exactly neutral
requested  Al2O3 0.1250  CaO 0.4530  SiO2 0.4220
achieved   Al2O3 0.1257  CaO 0.4531  SiO2 0.4212   (max error 0.00084)
formula units  Al2O3x63  CaOx227  SiO2x211
atoms          Si 211  Al 126  Ca 227  O 838   total 1402
ideal NBO/T from stoichiometry  0.973
initial configuration built by the internal random packer
structure checked: 1402 atoms, box 27.29 A, closest pair 1.60 A, density 2.600 g/cm3
wrote 9 stage scripts, total simulated time 107.13 ns
```

Nine stages: minimise, hot, cool, five equilibration blocks, production.

## Step 5. Run it (tab 4, Queue)

The table now shows one row, `S1_1600C`, status `NOT_SUBMITTED`, progress `0/9`.

Tick **watch and repair automatically**. Set **jobs at once** to 1 and **poll every** to 5
seconds, because the mock finishes in under a second and there is no point waiting a minute
between polls.

Watch the table. The status column walks through `PENDING`, `RUNNING`, `COMPLETED` for each
stage in turn, and the progress column climbs `1/9`, `2/9`, and so on. On a laptop the whole
nine stage chain takes about 15 seconds under the mock. When it reaches `9/9` the stage column
reads `done`.

Untick **watch and repair automatically**.

If it stops short, select the row and press **Details and history**. Every stage carries its
own timestamped history, including what the scheduler said and what the tool decided.

## Step 6. Read a diffusion coefficient (tab 5, Analysis)

Pick `S1_1600C` in the **run** dropdown. Press **Everything**.

This is the slow step. On a 400 frame, 1400 atom trajectory it takes a few minutes, because
the radial distribution function compares every pair in every frame. The status bar tells you
which property it is working on.

It starts by telling you what it read:

```
trajectory prod.lammpstrj
  400 frames, 1402 atoms, 3990.0 ps, species Ca, Si, Al, O
  coordinates are unwrapped, correct for MSD
  mean box 27.29 A, volume 20329 A^3
  temperature taken as 1873 K
```

Then, after the structure and speciation blocks, the part you came for:

```
diffusion
  Ca
      D = 1.0786e-09 m2/s = 1.0786e-05 cm2/s
      blocks  1.0150e-09 +/- 8.32e-12 m2/s over 5 blocks
      log-log slope 1.050, R squared 0.99989
      diffusive, log-log slope 1.05, R squared 0.99989. The Einstein relation is valid
      over this window.
      error bar 0.8 percent over 5 independent blocks. Well determined.
  Si
      D = 1.0205e-09 m2/s = 1.0205e-05 cm2/s
      blocks  1.0068e-09 +/- 1.98e-11 m2/s over 5 blocks
      log-log slope 1.015, R squared 0.99999
      diffusive, log-log slope 1.02, R squared 0.99999. The Einstein relation is valid
      over this window.
      error bar 2.0 percent over 5 independent blocks. Well determined.
  O
      D = 1.0036e-09 m2/s = 1.0036e-05 cm2/s
      blocks  1.0045e-09 +/- 1.28e-11 m2/s over 5 blocks
      log-log slope 1.002, R squared 1.00000
      diffusive, log-log slope 1.00, R squared 1.00000. The Einstein relation is valid
      over this window.
      error bar 1.3 percent over 5 independent blocks. Well determined.
```

**That is the point of the tour.** The mock was told to make every species diffuse at
1.0e-9 m^2/s. The analysis got back 1.00e-9 to 1.08e-9, with a log-log slope of 1 and an
error bar of one or two percent. The chain from input generation, through the scheduler,
through the written trajectory, through the analysis, is intact.

Your numbers will differ in the last digits because the mock uses a random seed drawn from
the stage name, but they should sit within a few percent of 1.0e-9 and every log-log slope
should be within 0.05 of 1.

The rest of the panel prints conductivity by both routes, viscosity by both routes, and an
equilibration verdict on every stage log. Section 8 explains how to read each one.

## Step 7. Draw the figures

Press **Make the figures**. Two PNG files appear in `runs/S1_1600C/figures/`:

| File | What it shows |
|---|---|
| `rdf.png` | Every partial radial distribution function involving oxygen |
| `msd_loglog.png` | Mean squared displacement of every species on log axes, with a dashed reference line of slope 1 |

Press **Save these results beside the run** to write the whole panel to
`runs/S1_1600C/ANALYSIS.txt`.

## What the mock can and cannot tell you

The mock is honest about being a fake, and you should be too.

| Property | Meaningful under the mock? | Why |
|---|---|---|
| Diffusion | **Yes, exactly** | The walkers have a chosen D, so the answer is known |
| Conductivity | Partly | The walkers are uncorrelated by construction, so the exact Einstein route should agree with Nernst-Einstein on average. On one trajectory it will not, and the tool will tell you the block error is too large to trust. That is the correct behaviour |
| Viscosity | As arithmetic only | The mock writes white noise, whose Green-Kubo viscosity is known in closed form, so it checks the integration and the unit conversion. It is not a melt viscosity |
| Equilibration report | Yes | The mock writes a real thermo table, and the report reads it exactly as it would a real one |
| **RDF, coordination, speciation, Qn** | **No** | This is important. The mock's atoms are free walkers that wander far outside the box, so the minimum image convention folds their separations into nonsense. You will see a first peak near 0.04 Angström, coordination numbers of zero, and 100 percent free oxygen, **and the tool will print two warnings saying exactly that**. Those are artefacts of the mock, not of the analysis. Structure numbers are only meaningful on a real trajectory |

# 4. Every tab, field by field

The window has six tabs, numbered in the order you normally use them. A background thread
does all the work, so nothing freezes. The progress bar at the bottom right spins while a
worker is busy.

## 4.1 Tab 1, Project

| Control | What it does |
|---|---|
| Project folder | The path to the folder that holds this study. It does not have to exist |
| Browse | A folder picker |
| Open or create | Opens the folder. If `project.json` is there it is loaded and you resume. If not, a new project is created and the cluster settings are prefilled from what this machine has |

Below that are two read only boxes.

**what this machine has** reports the scheduler detected on this machine, whether a LAMMPS
executable is on your PATH and what it is called, whether Packmol is on your PATH, and the
core count. When LAMMPS is missing it says so and explains what still works.

**force fields available** lists every force field marked `verified: true`, with its label,
its form, its charge set, its printed source, and the outcome of its own consistency check.
If any of them says anything other than `PASSES`, stop and read section 10.1, because a
parameter has been mistyped.

## 4.2 Tab 2, Cluster

### First, what a scheduler is

On your laptop you start a program by typing its name, and it runs immediately. A
supercomputer cannot work that way, because a thousand people share it. Instead you write a
small shell script that says what you want to run and how much of the machine you need, and
you hand that script to a **scheduler**. The scheduler puts it in a queue, waits until the
resources you asked for are free, then runs it on your behalf and saves the output to a file.
You are not there when it runs.

The two schedulers you will meet are **SLURM** and **PBS**. They do the same job with
different words. slagmd speaks both, plus a third option, `local`, which just runs the thing
on this machine with no queue at all.

| Word | What it means, plainly |
|---|---|
| **Scheduler** | The program that owns the queue and decides when your job runs. SLURM or PBS |
| **Batch script** | The shell script you hand to the scheduler. The lines starting `#SBATCH` (SLURM) or `#PBS` are the resource request. Everything after them is what actually runs |
| **Partition** or **queue** | A named slice of the machine with its own rules. Typically one partition for short test jobs, one for long production jobs, one for machines with big memory, one for nodes with GPUs. You must pick one that exists on your cluster, and only your cluster's documentation knows the names |
| **Account** or **project** | Who gets billed for the core hours. Many clusters refuse a job with no account |
| **Node** | One physical computer in the cluster. It has its own memory. A job asking for 2 nodes gets two whole computers that talk to each other over the network |
| **Task** or **rank** | One process of your parallel program. LAMMPS splits the simulation box into as many pieces as there are tasks. `nodes 2` and `tasks per node 128` gives 256 tasks |
| **Walltime** | The clock time limit. Not processor time. If your job is still running when the walltime expires, the scheduler kills it, however much work is left. Format `hh:mm:ss` |
| **Module** | Clusters keep many versions of every program installed side by side and none of them are active by default. `module load LAMMPS/2Aug2023-foss-2022a` is what puts one specific build on your PATH. The names are different on every cluster, which is why this field cannot have a useful default |
| **Launcher** | The command that actually starts the parallel processes. `srun` on SLURM, `mpirun` or `mpiexec` elsewhere. Some Flemish clusters use a wrapper called `mympirun`, which is what their own README uses |

### The fields

| Field | Meaning | A sensible value |
|---|---|---|
| scheduler | `local`, `slurm` or `pbs` | Whatever was detected. `local` on a laptop |
| partition or queue | The named slice of the cluster | Blank on a laptop. On a cluster, a name from its own documentation. There is no safe default |
| account or project | Who is billed | Blank on a laptop. Your project code on a cluster |
| nodes | Whole computers | `2`. That is what the published study used for classical MD |
| tasks per node | Parallel processes per node | `128`. Their nodes were two 64 core AMD EPYC 7H12 chips at 2.6 GHz |
| walltime hh:mm:ss | Clock limit | `24:00:00` is the default. For a real 50 ns production run of 5000 atoms on 2 nodes, ask for more. See section 9 |
| memory per node | Leave blank unless the cluster demands it | Blank |
| LAMMPS executable | The command that runs LAMMPS | `lmp`. Prefilled if one was found on your PATH |
| launcher | `srun`, `mpirun`, `mpiexec` or blank | `srun` on SLURM. Blank means run it directly with no MPI, which is right for `local` |
| email for failure notices | Adds a mail directive on failure and on hitting the time limit | Your address, or blank |
| modules to load, one per line | Each becomes a `module load` line, preceded by `module purge` | Your cluster's LAMMPS module. Empty on a laptop |
| extra shell lines before LAMMPS | Arbitrary shell, inserted after the module loads | Blank, unless you need to set a variable. Quote anything containing a space |
| extra scheduler directives | Raw directives, one per line, without the `#SBATCH` prefix | `--exclusive` if you want the whole node |

### The three buttons

**Save** writes the settings into `project.json` and immediately runs the check.

**Check for problems** runs `scheduler.config_problems`, which looks for four specific things
that will break the generated script:

1. A pre-command that assigns a value containing a space without quoting it. This is the
   trailing space trap. Bash splits it into two words and the error you get is
   `export: not a valid identifier`, several lines from the cause.
2. A LAMMPS executable path containing an unquoted space.
3. The `local` backend asked for more than one node, which it cannot do.
4. A batch backend with no walltime.

**Preview the submit script** shows the exact text that will be written for a stage, without
writing anything. Use it before you submit 754 jobs.

## 4.3 Tab 3, Compositions

One **run** is one composition, at one temperature, with one force field. It moves through
nine stages in order.

| Field | Meaning | Sensible value |
|---|---|---|
| run name | The folder name under `runs/` and the job name prefix. Must be unique in the project | `S1_1600C`, or `S1_1600C_bouhadja` if you will compare force fields |
| composition, oxide then mole fraction | Oxide names and mole fractions. Accepts `CaO 0.453 Al2O3 0.125 SiO2 0.422` or `CaO=0.453,Al2O3=0.125,SiO2=0.422`. If they do not sum to 1 they are normalised and you are told | Any of the supported oxides: CaO, Al2O3, SiO2, MgO, FeO, Fe2O3, TiO2, Na2O, K2O |
| temperature in C | Production temperature. Converted to kelvin internally | `1600`. Their five temperatures were 1400 to 1600 C in 50 C steps |
| target atoms | How big the box should be. Not exact, because formula units are whole numbers | `5000`, which is what they used. Below about 1000 the box gets close to twice the Coulomb cutoff and the tool starts warning |
| initial density g/cm3 | Only a starting guess. The protocol runs NPT at 1 atm from 4000 K down, so the box finds its own volume | `2.6`, a reasonable middle for these liquids near 1600 C |
| random seed | Feeds both the structure packing and the LAMMPS velocity creation | Vary it per replicate. Their scripts hard code `458165` in every copy, which means a repeated run reproduces the same trajectory rather than giving an independent one |
| force field | `bouhadja`, `matsui` or `guillot` | `bouhadja`. Their own force field comparison paper concludes it reproduces transport in CaO-Al2O3-SiO2 better than the other two |
| sample the stress tensor, so viscosity can be computed | Adds four lines to the production deck that write the off diagonal pressure every 5 steps | Leave it ticked. It costs almost nothing on a run that is happening anyway, and their 84 released files contain nothing that computes viscosity |
| add a short NVE energy conservation check | Appends a 0.1 ns constant energy run after production | Tick it once per new system. It is the cheapest check in molecular dynamics and almost nobody does it |

### The four buttons

| Button | What it does | Does it write anything? |
|---|---|---|
| Check this composition | Resolves the atom counts, checks the force field, checks neutrality, computes the box edge, estimates the cost, and lists any protocol warnings | No |
| Add run | Records the run in `project.json` | Yes, the state file |
| Add the ten published benchmark compositions | Queues the ten compositions of their Table 5, at 1600 C, for those whose published temperature range includes it | Yes |
| Prepare every unprepared run | For each run: builds the structure, writes the data file, writes all nine input scripts, writes all nine batch scripts, writes `protocol.json` | Yes, everything |

`Add the ten published benchmark compositions` exists for one reason. It lets you check the
tool against a study whose answers are in print. If it cannot reproduce their self diffusion
coefficients, the tool is wrong, and that is worth finding out before using it on anything
new.

## 4.4 Tab 4, Queue

The table has one row per run.

| Column | Meaning |
|---|---|
| run | The run name |
| composition | Mole fractions as given |
| T_C | Temperature in Celsius |
| ff | Force field |
| progress | Stages completed out of total, for example `4/9` |
| stage | The first stage that is not complete, or `done` |
| status | That stage's state, colour coded |
| job | The scheduler's job number, or `local:<pid>` |
| attempts | How many times this stage has been submitted. The cap is 3 |
| note | The last diagnosis, one line |

| Status | Colour | Meaning |
|---|---|---|
| NOT_SUBMITTED | grey | Waiting its turn |
| PENDING | amber | In the queue, not yet started |
| RUNNING | teal | Running |
| COMPLETED | green | Finished, **and its own completion marker was found in the log** |
| FAILED, TIMEOUT, OUT_OF_MEMORY, CANCELLED | red or grey | Left the queue without finishing |
| BLOCKED | red | Needs a human. Either the failure is one the tool refuses to retry, or three attempts were used |
| UNKNOWN | amber | The scheduler does not recognise the job any more. Read the log |

### The controls

| Control | What it does | Sensible value |
|---|---|---|
| watch and repair automatically | Starts the background monitor: it submits, polls, diagnoses failures, applies remedies, resubmits, and rewrites `REPORT.txt` every cycle | On, once you trust the settings |
| jobs at once | How many stages may be in flight across the whole project. This is the throughput knob | `4` on a cluster. `1` with the mock |
| poll every | Seconds between cycles | `60` on a cluster, because `squeue` is a shared service and hammering it is antisocial. `5` with the mock |
| Poll now | One cycle immediately | |
| Submit next | Submit the next idle stage, without starting the monitor | |
| Write report | Writes and shows `REPORT.txt`, a plain text summary that reads well on a phone | |
| Details and history | Everything about the selected run, including a timestamped history per stage and any protocol or cluster overrides an automatic repair applied | |
| Read the diagnosis | Shows the most recent `diagnosis.*.txt` for the selected run: what happened, why, what the tool did, the log lines that triggered it, and the log tail | |
| Unblock and retry | Resets the current stage to not submitted and clears its attempt count. Use after **you** have fixed whatever it was | |
| Cancel | Cancels every pending or running job of the selected run | |

Turning the monitor off does not affect anything already submitted. Those jobs keep running
on the cluster. Closing the window does not affect them either.

## 4.5 Tab 5, Analysis

| Control | What it does |
|---|---|
| run | Pick a run in this project. It looks for `runs/<name>/prod.lammpstrj` |
| or a trajectory file, Browse | Analyse any LAMMPS dump, for example one copied back from a cluster by hand |
| Everything | Runs every analysis in order |
| Structure: RDF and coordination | Partial radial distribution functions, first peak, first minimum, coordination number |
| Speciation: Qn and oxygen types | Bridging, non bridging, free and triclustered oxygen, plus Qn per network former |
| Diffusion with error bars | Mean squared displacement, the Einstein fit, the regime test, and block error bars |
| Conductivity, both routes | Nernst-Einstein and exact Einstein, and their ratio |
| Viscosity | Green-Kubo and Einstein-Helfand, from `stress.txt` |
| Was it equilibrated? | Reads every `log.*` beside the trajectory and reports drift in density, temperature, volume, pressure and total energy |
| Save these results beside the run | Writes the panel to `ANALYSIS.txt` |
| Make the figures | Writes `figures/rdf.png` and `figures/msd_loglog.png` |

Two things to know. The temperature is taken from the run record, not from the trajectory, and
the panel says which value it used. The charges are taken from the run's force field, which is
why conductivity only appears when you analyse a run rather than a loose file.

There is no place in this panel that prints a diffusion coefficient on its own. Every number
comes with the tool's verdict on whether it is allowed to exist.

## 4.6 Tab 6, Log

Everything the session did, in order, including full Python tracebacks if a worker thread
threw. It is the first place to look when a button appeared to do nothing.

This is the session log, held in the window. The project also keeps `events.log` on disk,
which survives the window closing and records run creation, preparation, submission,
completion and failure.

# 5. The concepts the tool is built on

## 5.1 What a force field is

A molecular dynamics program only knows how to do one thing: work out the force on every
atom, move each atom a little, and repeat. The force comes from a formula that says how the
energy of a pair of atoms depends on how far apart they are. That formula, plus its numbers,
is the **force field**. It is not measured during the simulation. It is decided in advance,
and everything the simulation says afterwards is a consequence of it.

For an oxide melt the formula has three parts.

| Part | Formula | What it does physically |
|---|---|---|
| Coulomb | q_i q_j / r | Ions attract or repel because they carry charge. This is the long ranged part, and it is why the Ewald sum exists |
| Repulsion | A exp((sigma - r) / rho), or A exp(-r / rho) | Two atoms cannot occupy the same place. The electron clouds push back, very steeply |
| Dispersion | minus C / r^6 | The weak attraction that exists between any two atoms |

Two functional forms appear in this work. **Born-Mayer-Huggins**, which slagmd calls `bmh`
and which LAMMPS calls `born/coul/long`, puts a distance `sigma` inside the exponential.
**Buckingham**, `buck`, does not. They are algebraically the same thing with the prefactor
absorbed, and that fact is used later as a consistency check.

The charges are not the formal charges. In a real silicate the electrons are partly shared,
so a fully ionic model is too strongly bound. Every one of these three force fields therefore
uses **effective charges**, a fixed fraction of the formal ones.

| Force field | Form | Ca | Al | Si | O | Fraction of formal | Fitted to |
|---|---|---|---|---|---|---|---|
| Matsui (1994) | bmh | 0.945 | 1.4175 | 1.89 | -0.945 | 0.4725 | 26 CMAS crystals |
| Guillot and Sator (2007) | buck | 0.945 | 1.4175 | 1.89 | -0.945 | 0.4725 | Silicate melt density and structure |
| Bouhadja, Jakse and Pasturel (2013) | bmh | 1.2 | 1.8 | 2.4 | -1.2 | 0.6 | Ab initio MD of low silica CaO-Al2O3-SiO2 melts |

Bouhadja's larger charges make the Coulomb attraction stronger, which is why its cation-oxygen
energy wells are deeper: about 66 percent deeper for Ca-O and about 57 percent for Al-O and
Si-O. Those percentages are a claim in their force field comparison paper, and the audit suite
recomputes them from the raw parameters rather than taking them on trust.

Guillot has one further feature worth knowing about: it has **no short range terms at all**
between cation pairs. Those interactions are purely Coulombic. That is stated explicitly in
the paper, and the tool treats the absence as deliberate rather than as missing data. When it
writes the deck it emits an explicit zero line with the comment `# Al-Ca Coulomb only`, so
nobody reading the deck later has to wonder.

## 5.2 Why the tool refuses unverified parameters

Every force field in `forcefields.json` carries a `verified` flag and a `source` string naming
the printed table it came from. `ForceField.load` refuses outright to load one whose flag is
false:

```
force field 'X' is marked verified=false. Refusing to build a deck from unverified
parameters. Add a printed source first.
```

The reason is that a wrong parameter does not produce an error. It produces a plausible
number. A force field recalled from memory, or interpolated, or estimated, gives you a melt
that runs perfectly happily for 50 ns and then hands you a diffusion coefficient that is
wrong by a factor you cannot see. There is no downstream check that catches it. The only
place to catch it is before the run.

On top of the flag there are **algebraic identities**, run every time a force field is loaded
for a real deck. These are the checks that catch a typing mistake in a transcribed table,
which is by far the most likely way a wrong number gets in.

| Force field | The identity | Why it must hold |
|---|---|---|
| Matsui | A_ij equals 4.184 times rho_ij, for all ten pairs | Matsui's own functional form has the prefactor f0 times (B_i + B_j), with f0 = 4.184 kJ per mol per Angström and rho_ij = B_i + B_j. So A and rho are not independent numbers, and if one is mistyped the identity breaks |
| Guillot | For the three cation-oxygen pairs, A_Guillot equals A_Matsui times exp(sigma_Matsui / rho) | This is what rewriting Born-Mayer-Huggins as Buckingham at fixed rho gives you. It is algebra |

The Guillot identity is tested at 0.5 percent, not 0.1 percent, and the reason is written into
the code. Matsui's A is printed to four decimal places and the exponential factor is about
four million, so a rounding in the last printed digit moves the answer by a couple of tenths
of a percent. The measured gap on 29 September 2026 was 0.21 to 0.22 percent for all three
pairs, the same size and the same sign, which is the signature of rounding. A genuine typing
mistake is a transposed digit, which is tens of percent or more, and still fails the test.

The Guillot O-O pair deliberately does **not** satisfy the identity, because Guillot's O-O rho
is 0.265 Angström against Matsui's 0.276. That pair is a real refit, not a rewriting, so the
check skips it rather than pretending otherwise. This is recorded in the JSON so nobody
rediscovers it as a bug.

## 5.3 Why oxygen is counted rather than rounded

Suppose you want 5000 atoms of CaO 0.453, Al2O3 0.125, SiO2 0.422. The obvious approach is to
multiply out how many of each oxide you need, round each to a whole number, and count up the
atoms.

That is wrong, and the failure is quiet.

Once each oxide count is rounded independently, the oxygen count no longer follows exactly
from the cations. Stoichiometry is broken. Break stoichiometry and you break charge
neutrality. And a LAMMPS box with a net charge and an Ewald sum does not fail. The Ewald
method assumes a neutral cell, so a charged one is handled by adding a uniform neutralising
background, silently. The simulation runs. The numbers are wrong.

Their own released Packmol example has exactly this defect. It asks for 40 Al, 80 Ca, 75 Si,
which needs 3(20) + 80 + 2(75) = 290 oxygen for stoichiometry, and requests **295**. With
Bouhadja charges that is a net minus 6 e in the cell. Their README does say the file is a
template whose numbers must be changed, so it never entered a production run, but it is a
clean illustration of how easy the mistake is.

slagmd does it the other way round.

1. Work out how many whole **formula units** of each oxide gets closest to the target atom
   count, using largest remainder rounding so the mole fractions stay as close as whole
   numbers allow.
2. Count the cations from the formula units.
3. Count the oxygen **from the cations**. It is never rounded, never adjusted, never chosen.

The atom total then comes out near the target rather than exactly on it, 1402 instead of 1400
in the tour, and that is the right trade. It also reports the achieved mole fractions and the
largest error, so you can see what the integer arithmetic cost you. For the ten published
benchmark compositions at 5000 atoms that error is under 0.001 in every case.

There is one further check, and it is stronger than testing one box. `charges_are_scaled_formal`
asks whether every charge in the force field is the **same** multiple of its formal charge.
If it is, neutrality holds for any stoichiometric mixture whatever, not just the one in front
of you. All three force fields pass: Matsui and Guillot at 0.4725, Bouhadja at 0.6.

## 5.4 Why the box must be more than twice the Coulomb cutoff

The simulation box is small and periodic. The box you see is surrounded by identical copies of
itself in every direction, which is how a few thousand atoms stand in for a bulk liquid.

To work out the force between two atoms in a periodic box you use the **minimum image
convention**: of all the copies of atom B, interact with the nearest one. That is
unambiguous only as long as the interaction range is less than half the box edge. If the
cutoff is larger than L over 2, an atom can reach more than one copy of the same neighbour,
and the question "which one is nearest" has more than one answer.

Their Coulomb real space cutoff is 11 Angström, which they took from Guillot and Sator. So the
hard requirement is L greater than 22 Angström.

| Box edge | What the tool does | Why |
|---|---|---|
| L less than or equal to 2 times the cutoff | **Refuses.** It tells you the minimum atom count that would fix it | The answer is not noisy, it is wrong. An atom interacts with its own periodic image |
| L between 2 and 2.5 times the cutoff | **Warns.** The run is valid and goes ahead | Finite size effects on diffusion and viscosity are measurable at this size. Their study used about 5000 atoms for this reason |
| L greater than 2.5 times the cutoff | Silent | Fine |

At 2.6 g per cm3, 5000 atoms of a CaO-Al2O3-SiO2 melt gives a box of about 41 Angström, well
clear. The refusal only fires on small test boxes, and it fires with a number: it tells you
how many atoms you need.

The same rule appears again in the analysis. `compute_rdf` will not compute g(r) beyond half
the smallest box edge. If you ask for more it cuts the range and says so, because beyond L
over 2 the minimum image convention double counts and g(r) stops meaning anything.

## 5.5 The five stage protocol

Their papers name five stages. The released scripts spread those five over seven LAMMPS
invocations, and slagmd generates nine stage scripts. The mapping is this.

| Their stage | Their file | slagmd stage | Length | What it is for |
|---|---|---|---|---|
| Initialization | `cooling.file`, first part | `00_minimise` | none | Remove any residual overlap, create velocities, minimise the energy |
| Randomization | `cooling.file` | `01_hot` | 5 ns at 4000 K, NPT | Forget the starting configuration completely. Whatever the packer produced must leave no trace |
| Cooling | `cooling.file` | `02_cool` | 2.127 ns, linear ramp 4000 K to 1873 K at 1 K per ps | Bring it down slowly enough that the structure follows |
| Equilibration | `equilibration_step1` to `step5` | `03_equil1` to `07_equil5` | 5 blocks of 10 ns, 50 ns total | Let the density settle, and hand over to production with the thermostat barely touching the dynamics |
| Production | `production.file` | `08_produce` | 50 ns | Measure everything |
| (not theirs) | | `09_nve_check` | 0.1 ns, optional | A constant energy run, to prove the timestep and the Ewald accuracy are not leaking energy |

Total 107.13 ns for a 1600 C run, which at 2 fs per step is 53.6 million steps.

Every stage reads the previous stage's restart file and writes its own. Every stage ends by
printing its own completion marker, `SLAGMD_STAGE_DONE <stage>`. That marker, and nothing
else, is what counts as finished. The scheduler saying `COMPLETED` only means the batch script
exited zero, and LAMMPS can exit zero having done almost nothing.

Three numbers in that table are worth pausing on, because they show why reading the released
scripts mattered more than reading the papers.

- `run 2500000` at 2 fs is 5 ns, and both papers say 5 ns. Their own comment on that line says
  `#2000 ps = 4 ns`, which is wrong twice over. There are four such stale comments in the
  release. The step counts are right and the comments are not, so slagmd generates the comment
  from the step count and the two can never drift apart.
- `run 1063500` at 2 fs is 2127 ps for a 2127 K drop, which is exactly 1 K per ps as the paper
  states.
- Five equilibration blocks of 5,000,000 steps is 5 times 10 ns, which is the 50 ns the paper
  quotes without saying it was five blocks.

## 5.6 Why the damping cascade exists

This is the part of their protocol that appears in neither paper, only in the released
scripts, and it is the most useful thing in the release.

Both the thermostat and the barostat are feedback controllers. The thermostat watches the
temperature and nudges the atoms to bring it to the target. The barostat watches the pressure
and changes the box volume. Each has a **damping time**: roughly, how long it takes to react.

There is a tension built into that.

- A **short** damping time means a tight, aggressive controller. The temperature and the
  density reach their targets fast. But the controller is also interfering with the motion of
  the atoms all the time, and the motion of the atoms is exactly what you are trying to
  measure. A tight thermostat distorts diffusion.
- A **long** damping time means a loose controller that barely interferes. The dynamics are
  nearly free. But it takes a long time to bring the density anywhere.

Their solution is to do both, in sequence. Start tight, to get the density where it belongs
quickly, then loosen it in steps until, by the time production starts, the controllers are
almost not acting. Recovered from their scripts:

| Stage | Thermostat damping | Barostat damping | Length | Their comment in slagmd |
|---|---|---|---|---|
| Hot, 4000 K | 100 fs | 200 fs | 5 ns | |
| Cooling | 100 fs | 200 fs | 2.127 ns | |
| Equilibration 1 | 50 fs | 100 fs | 10 ns | tight thermostat, let the density find itself |
| Equilibration 2 | 100 fs | 100 fs | 10 ns | loosen the thermostat |
| Equilibration 3 | 200 fs | 200 fs | 10 ns | loosen both |
| Equilibration 4 | 200 fs | 800 fs | 10 ns | loosen the barostat further |
| Equilibration 5 | 500 fs | 1000 fs | 10 ns | nearly free dynamics before production |
| Production | 1000 fs | 2000 fs | 50 ns | |

Both columns increase monotonically, by a factor of 20 on the thermostat and 20 on the
barostat from the tightest block to production. The audit suite checks that monotonicity, so
an edit that breaks the cascade fails a test.

The paper describes this in words and says the schedule is in the Supplementary Material. The
released scripts are the only place in the public record where the six pairs are actually
given. In the scripts they are written as multiples of the timestep, for example
`$(25*dt)`, which is where the 50 fs comes from.

slagmd also adds one thing of its own here, marked `[MINE]`: each equilibration block writes a
running average of density and volume to `equil<N>_running.txt`, so the analysis can **prove**
the density stopped drifting rather than assume it.

# 6. Every refusal, and what to do about it

A **refusal** stops the work and shows you a message. A **warning** lets the work go ahead and
tells you what it is worried about. In the window a refusal appears as `REFUSED` at the top of
the panel, or as an error dialog, and is written to the Log tab. Warnings appear inline,
prefixed `WARNING`, `PROTOCOL WARNING`, `CLUSTER WARNING` or `NOTE`.

None of these are bugs. Every one exists because the alternative is a plausible wrong number.

## 6.1 Refusals about the composition

From `composition.py`. All of these raise `CompositionError` and appear as `REFUSED` when you
press **Check this composition** or **Add run**.

| Message you see | What caused it | Fix |
|---|---|---|
| `oxide 'XYZ' is not in the supported list: [...]` | An oxide name that is not one of the nine known ones, or a typo such as `AL2O3` | Use one of CaO, Al2O3, SiO2, MgO, FeO, Fe2O3, TiO2, Na2O, K2O. Names are case sensitive. To add a tenth, see section 10.2 |
| `mole fractions sum to zero` | Every fraction was zero, or negative values cancelled out | Give positive fractions |
| `target_atoms below 50 is not a melt, it is a cluster` | A target atom count under 50 | Use at least a few hundred, and in practice at least 1000 so the box clears twice the Coulomb cutoff |
| `cannot pair up oxide and fraction in '...'` | An odd number of tokens in the composition line, usually a missing number or a stray word | Write it as name, number, name, number. Commas, equals signs and colons are all accepted as separators |
| `'abc' is not a number` | A fraction that will not parse as a float, often a comma used as a decimal point | Use a full stop for the decimal point |
| `the force field has no charge for ['Fe2']. Do not guess one.` | The composition contains a species the chosen force field was never parameterised for | Choose a different force field, or add a verified parameter set. Do not invent a charge |
| `density must be positive` | A zero or negative initial density | Use something near 2.6 g per cm3 |

There is one thing this module does **without** refusing. If your mole fractions do not sum to
1, they are normalised and the original total is recorded rather than rejected. The panel then
shows both the requested and the achieved fractions, so you can see what happened.

## 6.2 Refusals about the force field

From `forcefield.py`, raising `ForceFieldError`.

| Message you see | What caused it | Fix |
|---|---|---|
| `no force field 'X'. Available: [...]` | A name that is not a key in `forcefields.json`, or a key starting with an underscore, which is metadata | Use a name from the list the message prints |
| `force field 'X' is marked verified=false. Refusing to build a deck from unverified parameters. Add a printed source first.` | The parameter set has no confirmed printed source | Find the table it came from, check every number against it, record the source and the identity it should satisfy, then set the flag. Section 10.1 |
| `<label> has no effective charge for ['Fe2']. This force field was not parameterised for those species. Do not invent a charge, choose a different force field or add a verified parameter set.` | An iron or titanium bearing composition with Matsui, Guillot or Bouhadja, none of which cover them | Either extend the force field with verified parameters, or drop the species |
| `<label> has no parameters for the pairs ['Fe2-O', 'Fe2-Si']. LAMMPS would run with zero short-range repulsion for them and the melt would be physically meaningless. Refusing.` | A pair is missing from the table | Same as above. This is the important refusal: LAMMPS would happily run with A equal to zero, and you would get a liquid held together by nothing but Coulomb attraction, which collapses |
| `the stored parameters failed their own consistency identity:` followed by lines such as `matsui Al-O: A=0.7916 but 4.184*rho=0.7196, so either A or rho is mistyped` | A number in `forcefields.json` has been edited or mistyped | Go back to the printed table. Do not adjust the tolerance to make the test pass |

The Guillot cation-cation pairs are absent **on purpose**, and the tool knows the difference.
`missing_pairs` skips a cation-cation pair when the force field carries a `cation_cation` note,
which Guillot's entry does. A missing pair with no such note is a refusal.

## 6.3 Refusals about the starting structure

From `builder.py`. `build_random` and `build_with_packmol` raise `BuildError`.
`check_structure` returns a list of problems instead, and `Project.prepare` turns any non
empty list into a `ValueError` that stops preparation.

| Message you see | What caused it | Fix |
|---|---|---|
| `box edge 20.6 A is not more than twice the 11 A Coulomb cutoff. With 602 atoms an atom would interact with its own periodic image. Increase target_atoms to at least 735.` | Too few atoms for the cutoff | Use the number the message gives you, or more. 5000 is the published value |
| `could not place a Si atom after 4000 tries at box 27.57 A. The box is too small for the separation criteria. Lower the initial density or reduce min_sep.` | The packer ran out of room. Almost always an initial density far too high | Lower the initial density. Remember the protocol runs NPT from 4000 K down, so the starting density is only a guess and the box finds its own volume |
| `packmol is not on PATH. Use build_random instead.` | You asked for Packmol and it is not installed | Nothing to do. When called through `Project.prepare` this is caught, reported in the notes, and the internal packer is used instead. It never fails silently into something different |
| `packmol produced no output.` followed by Packmol's own output | Packmol ran and failed, usually because the tolerance cannot be met | Read the Packmol output shown. Loosen the tolerance or lower the density |
| `packmol placed 4988 atoms, expected 5000` | Packmol did not place everything | Same |
| `the starting structure failed its checks:` then `atom count 1400 does not match composition 1402` | Internal inconsistency between the packer and the composition | This is a tool bug. Report it |
| ... `species Ca: structure has 226, needs 227` | Same | Same |
| ... `net box charge is +0.945000, not zero. Ewald would be wrong.` | The charge set and the atom counts disagree | Almost always a force field whose charges are not a uniform scaling of the formal charges. The panel warns about that separately |
| ... `closest pair is 0.412 A. Below 1 A the repulsive term explodes and LAMMPS will lose atoms on the first step.` | Two atoms were placed on top of each other | Change the seed. If it repeats, the minimum separations are too small for that density |
| ... `box edge 21.0 A is not more than twice the 11 A cutoff, so the minimum image convention is not single valued` | Same as the first row, reached by a different path | Same |
| `net box charge +0.945000, refusing to build` | Raised by `Project.prepare` before the structure is even built | Same |
| `a run called 'S1_1600C' already exists` | Duplicate run name | Pick another name. The run name is the folder name |

And one warning, not a refusal:

| Warning | Meaning |
|---|---|
| `box edge 22.6 A is valid but under 28 A, which is 2.5 times the cutoff. Finite size effects on diffusion and viscosity are measurable at this size. Verbeeck 2025 used about 5000 atoms for this reason.` | The run is valid and goes ahead. If you are publishing a transport number, use a bigger box |

## 6.4 Protocol warnings

From `Protocol.warnings()` in `inputgen.py`. None of these stop anything. All of them are
printed by **Check this composition** and again during **Prepare**, and all of them are stored
in `protocol.json` so they travel with the run.

| Condition | What it tells you |
|---|---|
| timestep over 2.5 fs | Longer than the 2 fs they used. The Si-O stretch needs about ten steps per period |
| target temperature under 1200 K | Says what that is in Celsius and tells you to check the liquidus. Below it you have a supercooled liquid and the diffusion coefficients will not converge |
| production under 5 ns | Short against their 50 ns. Silicon is the slowest species and needs tens of ns at 1600 C to leave the cage regime |
| total equilibration under 5 ns | The density will still be drifting during production |
| fewer than 2 equilibration blocks | They used five with the damping loosened at each step. One block equilibrates more slowly and leaves more thermostat interference in the dynamics |
| randomisation temperature under 3000 K | May not erase the memory of the random insertion |
| trajectory sampling less often than every 20 ps | Few points for the mean squared displacement fit. They used 10 ps |
| stress sampled less often than every 20 steps | The stress autocorrelation decays in well under a picosecond. Sampling it less often than about every 10 steps aliases it away, and the Green-Kubo viscosity will be **wrong rather than noisy** |

Their own protocol, unmodified, produces no warnings. The audit checks that.

## 6.5 Refusals in the analysis

From `analysis.py`. Some raise, some return a verdict instead of a number. The distinction
matters: a raise means it could not compute anything, a verdict means it computed something
and is telling you not to quote it.

| Message | Raise or verdict | What caused it | Fix |
|---|---|---|---|
| `no frames found in <path>` | raise | The file is not a LAMMPS dump, or it is empty, or truncated before the first frame | Check the file |
| `no atoms of species 'Mg' in the trajectory` | raise | Asked for a species the dump does not contain | Check the species list the panel printed |
| `fitting window has fewer than 4 points, the run is too short` | raise | A trajectory with only a handful of frames | Sample more often, or run longer |
| `no oxygen atoms were classified, check the cutoffs` | raise | Speciation found no oxygen at all, meaning the dump has no species called `O` | Make sure the dump has an `element` column, or that the type to element mapping is right |
| `only 50 stress samples in stress.txt. Green-Kubo needs tens of thousands. The stress must be sampled every few steps, not every picosecond: the stress autocorrelation decays in well under 1 ps and a 1 ps sampling interval aliases it away completely.` | raise | The stress file is far too short or too sparse | Set the stress sampling interval to 5 steps, which is what the tick box does |
| `too few stress samples for Einstein-Helfand` | raise | Same | Same |
| `rmax 10.0 A exceeds half the smallest box edge, cut to 13.64 A. Beyond L/2 the minimum image convention double counts and g(r) is meaningless.` | warning | You asked for a longer range than the box allows | Nothing. The range was cut for you |
| `the MSD slope is not positive. Nothing is diffusing. Either the run is far too short or the melt has frozen at this temperature.` | verdict | A flat or falling mean squared displacement | Check the temperature against the liquidus. Run longer |
| `log-log slope is 0.42, not 1. The atoms are still rattling in their cages over this window, so this is NOT a diffusion coefficient. Run longer, or fit later in the trajectory. This is the single commonest way a wrong D gets published.` | verdict | The fit window is in the cage regime | Run longer. Silicon needs the most |
| `log-log slope is 1.31, above 1. Either the fit window is still in the ballistic regime at very short times, or the MSD is dominated by a drifting centre of mass.` | verdict | Fitting too early, or centre of mass drift | Use `com yes` on the LAMMPS compute, which the generated deck does, or subtract the drift |
| `the regime is diffusive ... but the linear fit R squared is only 0.9820. The statistics are poor. Report an error bar from block averages, do not quote a single figure.` | verdict | A noisy mean squared displacement | Use the block error bar, which is printed on the next line anyway |
| `the error bar is 41 percent of the value. This D is not measured, it is estimated. Run longer or average independent trajectories before quoting it.` | verdict | Block to block scatter over 25 percent | Run longer, or run replicates with different seeds |
| `too few frames to split into independent blocks, so there is no error bar. Treat the value as indicative only.` | verdict | Fewer than about 60 frames | Sample more often |
| `the block to block scatter gives an error of 79 percent. This is NOT a measured conductivity. The exact Einstein relation needs either a very long trajectory or many independent runs, because it depends on one collective coordinate. Quote the Nernst-Einstein value and state that the exact value could not be resolved.` | verdict | Inherent to the method on one trajectory | Run replicates. See section 8.5 |
| `NOT converged. Across 2.0 to 10.0 ps the running integral scatters by 34 percent. Green-Kubo has not plateaued. Do not quote this number.` | verdict | The stress autocorrelation has a long tail that has not died | Run much longer, average independent trajectories, or use the Einstein-Helfand value |
| `the requested plateau window contains fewer than 5 points, so no viscosity is reported. Widen the window or sample more often.` | verdict, and the viscosity is returned as nothing at all | The plateau window falls outside the data | Widen it |
| `over the fit window the Helfand mean square displacement grows as t^0.41, not t^1, so this is not yet a viscosity. If the exponent is below 1 the window is too LATE ... If it is above 1 the window is too early` | verdict | The Einstein-Helfand fit window is in the wrong place | Move it. The default is deliberately early, 0.2 to 2 percent of the series. See section 8.6 |
| `density is still drifting by +0.84 percent over the second half of the run. It is not equilibrated. Every diffusion coefficient and conductivity measured here is at the wrong density.` | verdict | Equilibration was too short | Add equilibration time. This is the check almost nobody performs |
| `temperature fluctuates by 310 K around 2395 K, which is over 5 percent. Check the thermostat damping.` | verdict | Either a genuinely loose thermostat, or you pointed the report at a cooling stage log, where the temperature is **meant** to be changing | Read which log it was reporting on |
| `47.2 percent of oxygens are triclustered. That is high. It is real in alumina rich melts, but it can also mean the cutoff is too generous.` | note | Suspicious speciation | Check the cutoff against the radial distribution function first minimum |
| `100.0 percent free oxygen. In a silicate melt free O is usually a few percent at most, so a large number suggests the cutoff is too tight or the composition is very basic.` | note | Same | Same. On a mock trajectory this note fires because the mock is not a melt |

## 6.6 Refusals about submission

From `scheduler.py`, raising `SubmitError`, which puts the stage into `BLOCKED`.

| Message | What caused it | Fix |
|---|---|---|
| `sbatch is not on PATH, so this machine is not a slurm submit host. Either run slagmd on the cluster login node, or set the backend to local.` | You set the backend to SLURM on a machine that has no SLURM | Run slagmd on the login node, or use `local` |
| `sbatch returned 1` followed by the scheduler's own output | The scheduler rejected the script: a partition that does not exist, an account with no allocation, a walltime over the queue limit | Read the scheduler's message. It is quoted in full |
| `could not find a job id in sbatch output: '...'` | `sbatch` succeeded but printed something unexpected | Unusual. Report it with the output quoted |

## 6.7 The attempt cap

Any stage that has been submitted three times is set to `BLOCKED` and not submitted again.
Three is `diagnose.MAX_ATTEMPTS`. The Queue tab shows the attempt count, and `Unblock and
retry` resets both the status and the count, which you should only press after changing
something.

# 7. The failure diagnosis system

## 7.1 How it works

When a stage leaves the queue **without** its completion marker, the tool does this, in order:

1. Reads every log the stage could have written: the LAMMPS log, the scheduler's standard
   output and error files under both SLURM and PBS names, the batch script's own output, and
   `log.lammps` as a fallback. All of them are concatenated and searched together, because a
   module error lands in a different file from a LAMMPS error.
2. Matches the text against a table of patterns, most specific first.
3. If one matches, it builds a **remedy**: a plain English description of what happened, an
   explanation of why, the action it will take, a set of changes to the protocol or the cluster
   settings, and a flag saying whether it may retry.
4. Writes `diagnosis.<stage>.<attempt>.txt` in the run folder, containing all of that plus the
   log lines that triggered the match and the last 4000 characters of the log.
5. Applies the changes, rebuilds the deck if the changes affect it, rolls the affected stages
   back to not submitted, and resubmits, **if** the remedy allows a retry and the attempt count
   is under three.

The Queue tab shows a one line version in the `note` column, and `Read the diagnosis` opens
the full file.

## 7.2 The patterns

The table below has 15 rows. Thirteen of them are regular expressions matched against the log.
The fourteenth, `UNKNOWN`, is the fall through for a log that ends without a recognised error
and without a completion marker. The fifteenth, `NO_LOG`, is reached before any matching is
attempted, when there is no log file at all. The tool's own summary counts 14 failure patterns,
which is the thirteen plus `UNKNOWN`.

Every one of the thirteen is tested in the audit suite against the **verbatim** error text that
LAMMPS or SLURM actually prints, taken from `mocklammps.py`. A diagnosis rule that has never
seen its own error message is a guess.

| Code | What the log says | What it means, plainly |
|---|---|---|
| `LOST_ATOMS` | `ERROR: Lost atoms: original 5002 current 4998 (src/thermo.cpp:488)` | LAMMPS can no longer find all the atoms it started with. Two atoms got close enough that the repulsive part of the potential threw one of them across the box in a single step. In an oxide melt this is almost always the first dynamics step after a bad starting structure, or a timestep too long for 4000 K |
| `NAN` | `ERROR on proc 3: Non-numeric atom coords - simulation unstable (src/domain.cpp:210)` | A coordinate or a pressure became not a number. An energy or a force overflowed. Same root cause as lost atoms, a pair far inside the repulsive wall. It can also mean a pair coefficient with the wrong units, which makes the forces enormous from step zero |
| `OUT_OF_RANGE` | `ERROR: Out of range atoms - cannot compute PPPM (src/KSPACE/pppm.cpp:1927)` | Atoms moved further in one step than the distance processors exchange information over, so a processor could not find its neighbours. Either too hot for the timestep, or the neighbour list is rebuilt too rarely |
| `NEIGH_OVERFLOW` | `ERROR on proc 0: Neighbor list overflow, boost neigh_modify one (src/npair_half_bin_newton.cpp:130)` | The list of nearby atoms ran out of room. With an 11 Angström cutoff in a dense oxide melt each atom has several hundred neighbours, and LAMMPS's default page size is not always enough |
| `BOX_BLEW_UP` | `ERROR: Domain too large for neighbor bins (src/nbin_standard.cpp:97)` | The simulation box became an impossible size. The barostat lost control: a barostat damping time too short fights the thermostat and the volume oscillates until it diverges. The classic failure of an aggressive constant pressure start |
| `KSPACE` | `ERROR: KSpace accuracy must be > 0 (src/ewald.cpp:152)` | The Ewald sum could not reach the accuracy it was asked for. A 1e-6 accuracy with an 11 Angström real space cutoff needs a lot of reciprocal vectors, and on a small box the two requests can be incompatible |
| `MISSING_INPUT` | `ERROR: Cannot open restart file restart.02_cool: No such file or directory (src/read_restart.cpp:135)` | A file this stage needed was not there. Usually the previous stage never wrote its restart file because it died. Occasionally the run directory is on a filesystem the compute node cannot see |
| `BAD_PAIR_COEFF` | `ERROR: Incorrect args for pair coefficients (src/pair_born_coul_long.cpp:401)` | LAMMPS rejected the pair coefficients. The deck and the pair style disagree: `born/coul/long` takes A rho sigma C D, `buck/coul/long` takes A rho C. A force field written in the wrong form gives exactly this |
| `OOM` | `slurmstepd: error: Detected 1 oom-kill event(s) in StepId=8812345.batch.` | The job was killed for using too much memory. Ewald plus a long neighbour list plus per atom dump buffers on a large box adds up. Spreading the same atoms over more nodes gives each process less to hold |
| `WALLTIME` | `slurmstepd: error: *** JOB 8812345 ON node042 CANCELLED AT 2026-09-29T14:02:11 DUE TO TIME LIMIT ***` | The job hit the clock limit and was killed part way through. The commonest failure and the least serious, because a restart file means nothing is lost. A 50 ns production run of 5000 atoms is not a 24 hour job on two nodes |
| `MPI_CRASH` | `srun: error: node071: task 44: Segmentation fault` | One parallel process died and took the job with it. This could be a faulty node, which is not your problem and a resubmission fixes, or a real crash inside LAMMPS, which a resubmission will reproduce |
| `MISSING_PACKAGE` | `ERROR: Unknown pair style born/coul/long (src/force.cpp:271)` and `Pair style born/coul/long requires the KSPACE package` | LAMMPS does not recognise something in the deck, because the build on this cluster was compiled without the package that provides it. `born/coul/long` needs KSPACE |
| `ENVIRONMENT` | `Lmod has detected the following error: Unable to locate a modulefile for 'LAMMPS'` or `lmp: command not found` | The batch script could not find LAMMPS or a module. Module names differ on every cluster. Nothing about the physics is wrong |
| `UNKNOWN` | Nothing recognised, and no completion marker | Something outside LAMMPS ended it: a node reboot, a filesystem outage, or an error message this tool has not been taught yet |
| `NO_LOG` | There is no log file at all | The job never started. The batch script may not have been submitted, or it failed before LAMMPS was reached |

And what it does about each one:

| Code | What the tool does | Retry? |
|---|---|---|
| `LOST_ATOMS` | Halves the timestep to 1 fs, rebuilds the starting structure, and restarts from `00_minimise` | Yes |
| `NAN` | Halves the timestep, rebuilds the structure, restarts from `00_minimise`, and tells you to verify the pair coefficients | Yes |
| `OUT_OF_RANGE` | Halves the timestep | Yes |
| `NEIGH_OVERFLOW` | Sets `neigh_one` to 8000 and `neigh_page` to 200000 in the protocol, which makes the generated deck emit `neigh_modify one 8000 page 200000`, then resubmits the stage | Yes |
| `BOX_BLEW_UP` | Increases barostat damping: hot to 800 fs, production to 4000 fs, and multiplies every equilibration block's barostat damping by 4. Rebuilds and resubmits | Yes |
| `KSPACE` | Loosens the Ewald accuracy one decade to 1e-5, which is still tight for a melt, and records the change in the run history | Yes |
| `MISSING_INPUT` | **Stops.** Tells you to check the previous stage | **No** |
| `BAD_PAIR_COEFF` | **Stops.** This is a bug in the generator, not a physics problem | **No** |
| `OOM` | Doubles the node count and resubmits the same stage | Yes |
| `WALLTIME` | Doubles the walltime, capped at seven days, and resumes. If your queue has a maximum walltime below that, split the stage instead | Yes |
| `MPI_CRASH` | Resubmits once unchanged. The attempt cap of three then stops the loop. If it fails the same way twice, stop and report the node list | Yes, but bounded |
| `MISSING_PACKAGE` | **Stops.** Load a LAMMPS module built with KSPACE, or ask cluster support. Resubmitting the same deck will fail identically | **No** |
| `ENVIRONMENT` | **Stops.** Fix the module list and the executable name on the Cluster tab | **No** |
| `UNKNOWN` | **Stops**, and says so: `STOPPING. Read the last lines of the log yourself.` | **No** |
| `NO_LOG` | **Stops.** Check the scheduler's own output and the submit script | **No** |

One case is not a failure at all. If the log contains `SLAGMD_STAGE_DONE` or LAMMPS's own
`Total wall time`, the tool reports that a normal completion marker is present and produces no
remedy. This is what happens when a job is killed a moment after finishing.

## 7.3 Why it refuses to retry a failure it does not recognise

This is the single most important design decision in `diagnose.py`, and it is worth stating
plainly.

A retry loop that does not understand the failure has two possible outcomes. Either the failure
was transient, in which case the retry works, or it was not, in which case the retry reproduces
the crash. There is no way to tell in advance. So the loop will resubmit the same broken job
until it hits the attempt cap, and then, in the worst case, do it again for the next 118
compositions.

The cost of guessing is real. Their production stage is a 50 ns run of 5000 atoms on 2 nodes of
128 cores. A job that reaches the walltime before crashing burns its whole allocation. Three
attempts across 119 chains is a serious amount of somebody's compute budget spent producing the
same error message.

The cost of stopping is that you read a log. One log, once, and then you teach the tool the
pattern so it never needs you again. Section 10.3 says how.

So the rule is: **if slagmd cannot name the failure, it stops and tells you it is stopping.** It
does not retry hopefully, and it does not mark the run as finished.

Four recognised failures also refuse to retry, and they are the four where a retry is certainly
useless: a missing input file, a malformed pair coefficient, a LAMMPS build lacking a package,
and a wrong module name. Nothing about the job will be different the second time.

## 7.4 A worked example, from the audit suite

This exact sequence runs in group G of the audit, with the fault injected deliberately.

1. A run is prepared. Its walltime is `24:00:00`.
2. Stage `01_hot` is submitted. The mock LAMMPS is told to fail once with SLURM's real walltime
   message, and to succeed on the second attempt.
3. The stage leaves the queue with no completion marker. The tool reads the logs, matches
   `DUE TO TIME LIMIT`, and builds the `WALLTIME` remedy.
4. It writes `diagnosis.01_hot.1.txt`.
5. It applies `walltime: x2`, which gives `48:00:00`, stores it as a cluster override on that
   run, rewrites the batch scripts, rolls `01_hot` back to not submitted, and resubmits.
6. The second attempt succeeds. The chain continues through all remaining stages and finishes.

The audit checks that the walltime ended at exactly `48:00:00`, that the diagnosis file exists,
that the run finished with no human involved, and that the failed stage took exactly two
attempts and not three.

The same test is then run with a **permanent** fault. With `MISSING_PACKAGE`, the run must stop
after exactly **one** attempt, not three. With `MPI_CRASH`, which is allowed one retry, it must
stop at the attempt cap. Both are checked.

# 8. The analysis, and how to read the output

A note on units before anything else. Inside the analysis module, distance is in Angström,
time in picoseconds and temperature in kelvin. Results come out in SI: diffusion in m^2 per s
and also cm^2 per s, conductivity in S per m, viscosity in Pa s.

A note on where the numbers come from. Both papers' own diffusion coefficients are printed in
units of 1e-6 cm^2 per s, which is 1e-10 m^2 per s. Their high throughput dataset's first
composition, at 1673 K, reports D_Ca = 0.4664, D_Al = 0.0265, D_Si = 0.0142, D_O = 0.0256 in
those units. In slagmd's units that is 4.7e-11, 2.7e-12, 1.4e-12 and 2.6e-12 m^2 per s. Use
that as an order of magnitude anchor, and note the ordering: calcium is faster than everything
else by more than a factor of ten, and silicon is the slowest. If your run comes out with
silicon faster than calcium, something is wrong before you look at any other number.

## 8.1 Radial distribution functions and coordination numbers

### What it is

The radial distribution function, g(r), answers one question: if I stand on an atom of type A,
how many atoms of type B are at distance r, compared with how many there would be if the liquid
were a structureless gas at the same density? A value of 1 means "no more and no fewer than
random". A peak at 1.6 Angström means "there is a shell of neighbours at 1.6 Angström".

The **first peak** is the bond length. The **first minimum** after it is the edge of the first
shell. The **coordination number** is how many B atoms are inside that shell, obtained by
integrating 4 pi rho g(r) r^2 out to the first minimum, which is the definition their force
field comparison paper uses.

### What the tool prints

```
structure
  pair      first peak A   first min A   coordination
  Si-O            1.612         2.284           4.02
  Al-O            1.760         2.480           4.31
  Ca-O            2.380         3.180           6.40
  O-O             2.640         3.520          10.5
```

One row per pair, for every pair of species present. The numbers above are illustrative
of a real calcium aluminosilicate melt, not output from a particular run.

### How to tell a good number from a bad one

The tool does **not** check these against expected values, so this is on you. These are the
values silicate melt chemistry expects, not numbers taken from their tables:

| Pair | Expected first peak | Expected coordination | If it is different |
|---|---|---|---|
| Si-O | about 1.6 Angström | 4, tightly | A coordination well away from 4 means the cutoff is wrong or the melt has not equilibrated. Silicon in a silicate is tetrahedral and stubbornly so |
| Al-O | about 1.75 Angström | 4 to 5 | Aluminium is the interesting one. It sits at 4 when there is enough calcium to charge balance it, and drifts up when there is not. A value near 6 in a calcium rich melt is suspicious |
| Ca-O | about 2.4 Angström | 6 to 7, broadly | Calcium's shell is genuinely ill defined, so a spread is normal |
| O-O | about 2.6 Angström | 9 to 12 | |

Three failure signatures to watch for:

- **A first peak below 1 Angström.** Physically impossible. It means the coordinates you gave
  it are not what it thinks they are, most often unwrapped coordinates being folded by the
  minimum image convention against a box they have long since left. This is what the mock
  trajectory produces. **Since 29 September 2026 the tool warns about it**, naming the pairs
  and telling you not to use the coordination numbers. It also warns separately when a
  coordination number comes out near zero, which cannot happen in a condensed melt.
- **A warning that rmax was cut.** Harmless, but it tells you the box is smaller than you
  thought.
- **A coordination number of exactly 0.00.** The first minimum was found in the wrong place,
  so the integral covers nothing.

One piece of history worth knowing. The peak finder used to return the **second** coordination
shell, because it began its search one bin past the first non empty bin. On a sharply peaked
g(r), such as a crystal or a very cold melt, the first non empty bin **is** the first peak. The
bug was caught by testing against a simple cubic lattice, whose answers are 3.0 Angström and 6,
and where the tool was reporting 4.245 Angström and 18. That test is now check number one of
group E.

## 8.2 Oxygen speciation and Qn

### What it is

In a silicate melt the oxygen atoms do the structural work, and they come in four kinds.

| Class | Definition | What it means |
|---|---|---|
| **Bridging** | Bonded to exactly two network formers | The glue. A Si-O-Si or Si-O-Al link holding the network together |
| **Non bridging** | Bonded to exactly one network former | A broken link. The network ends here, and the charge is balanced by a nearby calcium |
| **Free** | Bonded to no network former | Surrounded only by modifier cations. Rare in a silicate melt |
| **Triclustered** | Bonded to three or more | Real in high alumina melts, not an artefact. Their group counts them, so slagmd does too |

A network former is Si, Al, Fe(III) or Ti. A modifier is Ca, Mg, Na, K or Fe(II).

**Qn** is the same information seen from the cation's side: for each network former, how many
**bridging** oxygens it has. Q4 is a fully connected tetrahedron, in the middle of the network.
Q0 is an isolated tetrahedron floating free. The Qn distribution is the standard way to describe
how depolymerised a melt is, and depolymerisation is what controls viscosity and diffusion.

Note the definition carefully. Qn counts bridging oxygens, not simply neighbours. A silicon with
four oxygens of which two are non bridging is Q2, not Q4.

**NBO/T** is the number of non bridging oxygens per network former, a single number summarising
the same thing.

### The cutoff is the whole argument

Deciding whether an oxygen is "bonded" to a silicon means picking a distance. Get it wrong and
every number above changes. So by default slagmd does not use a literature value. It computes
the former-oxygen partial radial distribution function from your own trajectory and takes the
**first minimum** as the cutoff, then prints which cutoffs it used.

```
speciation
  cutoffs measured from this trajectory: Si 2.28 A, Al 2.48 A
  bridging       58.31 %
  non bridging   38.92 %
  free            1.44 %
  triclustered    1.33 %
  NBO per network former  0.982
  Q(n) for Si: Q0 1.2%  Q1 9.4%  Q2 28.1%  Q3 39.0%  Q4 22.3%   mean 2.72
  Q(n) for Al: Q0 2.1%  Q1 12.0%  Q2 30.4%  Q3 37.2%  Q4 18.3%   mean 2.58
  coordination of Si by oxygen: 4.021
  coordination of Al by oxygen: 4.312
```

This is what a plausible panel looks like for a calcium aluminosilicate near 1600 C. Their own
hand picked cutoffs, for comparison, are Si-O 2.0, Al-O 2.5, Ca-O 3.0 and O to Si or Al 2.5
Angström, and their own script comment next to them reads `# --> cutoff should be adjusted`.
You can pass explicit cutoffs instead of measuring them, from Python, if you want to match
theirs exactly.

### How to tell a good number from a bad one

| Check | What it should be | What it means if it is not |
|---|---|---|
| The four fractions sum to 1 | Exactly | A tool bug. The audit checks this to 1e-9 |
| Free oxygen | A few percent at most in a silicate melt | Over 10 percent triggers a note. Either the cutoff is too tight, so real bonds are being missed, or the composition really is very basic |
| Triclustered oxygen | A few percent, more in alumina rich melts | Over 5 percent triggers a note. Check the cutoff against the RDF first minimum |
| Measured NBO/T against the ideal | Should be close, but not identical | The panel gives the measured value and the composition panel gives the ideal value from stoichiometry alone. They differ because aluminium can act as a former or as a modifier depending on how much charge balancing calcium is around, which stoichiometry cannot know. That gap is a real result, not an error |
| Qn distributions sum to 1 | Exactly, per species | A tool bug |
| Mean Qn | Between 0 and 4, and falling as you add CaO | A melt with more calcium is more broken up, so mean Qn should drop |

## 8.3 Diffusion

### What it is

Watch one atom. Over time t it wanders away from where it started. Square that distance and
average over every atom and every possible starting point, and you get the **mean squared
displacement**. In a liquid, for long enough times, it grows in a straight line:

MSD = 6 D t

That is the Einstein relation, and D is the self diffusion coefficient. The 6 is because there
are three dimensions and the factor of 2 from the random walk. slagmd fits a straight line to
the middle of the mean squared displacement, from 30 percent to 80 percent of the available
lag, and divides the slope by 6, exactly as their equation 4 does.

Two implementation points. The mean squared displacement uses **every** available starting
point, not just the first frame, which is what makes a long run give a smooth curve. And it
uses unwrapped coordinates, `xu yu zu`, which the generated deck dumps. Wrapped coordinates give
a wrong answer the moment an atom crosses the box boundary, and in constant pressure dynamics
the box length is changing too, so unwrapping afterwards from image flags is awkward and easy to
get wrong. If you hand it a wrapped dump anyway it unwraps by following each atom and undoing
any jump larger than half the box, and it says in the panel which kind of coordinates it found.

### The log-log slope test, in plain words

The Einstein relation is only true at long times. At short times an atom is rattling inside a
cage of its neighbours, and the mean squared displacement flattens off instead of rising. If you
fit a straight line through that flat part, you get a slope, you divide it by 6, and you get a
number with the units of a diffusion coefficient that is not one. Nothing about the arithmetic
complains.

There is a simple test. Plot the mean squared displacement against time on axes where both are
logarithmic. A relationship of the form MSD proportional to t^n becomes a straight line of
slope n. So:

| Log-log slope | Regime | What it means |
|---|---|---|
| about 2 | Ballistic | Very short times. Atoms have not yet hit a neighbour |
| well below 1 | Caged | Atoms are rattling in place. **Not diffusion** |
| **1** | **Diffusive** | The Einstein relation holds. This is the only regime where D exists |
| above 1 | Suspect | Either you are still in the ballistic regime, or the whole box is drifting and you are measuring the drift |

slagmd measures that slope over the same window it fits, and if it is not within **0.15** of 1
it refuses to call the answer a diffusion coefficient. Their own analysis does the same thing in
a different way: their script searches forward through the trajectory until the log-log slope is
within 0.01 of unity, then applies a 2 percent deviation criterion and a safety buffer.

### The block error bar, in plain words

A single fitted line has no error bar worth having. You can compute the scatter of the points
around the fit, but it is useless as an uncertainty, because consecutive points of a mean
squared displacement are not independent measurements. Every point shares starting positions
with its neighbours. So the curve can be beautifully smooth, with an R squared of 0.99999, and
still sit on top of a badly determined slope.

The cheap honest answer is to cut the trajectory into pieces that do not share starting points,
fit each piece separately, and look at how much the answers differ. That is a **block average**.
slagmd uses five blocks by default, and reports the mean of the block values and the standard
error of that mean.

The audit suite tests this the right way round. It builds a random walk with a known D, and
checks not that the answer is within some arbitrary tolerance, but that the **true value lies
inside the error bar the tool itself reported**. If it did not, the error bar would be the bug.

Their own scripts quote a standard deviation over sliding fit windows instead. Blocks are
stricter, because sliding windows on one trajectory are not independent of each other.

### What the tool prints, and how to read it

```
  Si
      D = 1.0205e-09 m2/s = 1.0205e-05 cm2/s
      blocks  1.0068e-09 +/- 1.98e-11 m2/s over 5 blocks
      log-log slope 1.015, R squared 0.99999
      diffusive, log-log slope 1.02, R squared 0.99999. The Einstein relation is valid
      over this window.
      error bar 2.0 percent over 5 independent blocks. Well determined.
```

| Line | Read it as |
|---|---|
| `D = ...` | From the whole trajectory. This is the number to quote, **if** the two verdicts below are good |
| `blocks ... +/- ...` | The independent estimate and its uncertainty. It should agree with the line above. A large disagreement means the trajectory is not stationary |
| `log-log slope`, `R squared` | The evidence. Slope within 0.15 of 1, R squared above 0.995 |
| First verdict | Whether this is a diffusion coefficient at all |
| Second verdict | Whether it is well enough determined to quote |

| Relative error bar | Verdict | What to do |
|---|---|---|
| under 10 percent | Well determined | Quote it |
| 10 to 25 percent | Usable | Quote it **with** the error bar, never as a bare number |
| over 25 percent | Not measured, estimated | Run longer, or average independent trajectories with different seeds |
| no blocks possible | No error bar | Fewer than about 60 frames. Treat the value as indicative only |

## 8.4 Electrical conductivity, two routes

Their force field comparison paper defines conductivity two ways, and slagmd computes both,
because the difference between them is itself a result.

### Nernst-Einstein, their equation 7

Take the diffusion coefficient of each species, weight it by the square of that species' charge
and by how many of them there are, add it up, divide by volume and temperature. This is the
approximate route, and its assumption is stated in its own name in the panel: it **assumes the
ions move independently**.

It is cheap, because it needs nothing that was not already computed for diffusion.

### Exact Einstein, their equation 8

Instead of adding up independent mean squared displacements, track the mean squared displacement
of the **total charge**: add up the position of every ion weighted by its charge, and watch how
that single quantity wanders. Cross correlations between species are then included
automatically, because they are inside the square.

### The ratio is the physics

The panel prints both and their ratio:

```
electrical conductivity
  Nernst-Einstein, assumes ions move independently:  153.769 S/m
  exact Einstein, includes cross correlations:       301.666 S/m
  block error on the exact value: 24 percent
  Haven ratio, exact over Nernst-Einstein: 1.962
```

The ratio of the exact value to the approximate one is the **Haven ratio**. A value below 1
means the ions do not move independently: when a calcium goes one way, its neighbours tend to go
with it, and that correlated motion carries no net charge. The gap between the two numbers is
exactly the correlation the approximation throws away.

**Read the block error before you read the ratio.** The exact Einstein route has a problem built
into it. The collective charge displacement is a single coordinate, not one per atom, so however
many thousand ions you have, one trajectory gives only a handful of independent samples of it.
Measured on 120 synthetic replicas of a 600 ion box over 400 ps, the scatter from one trajectory
to the next was 79 percent while the average over replicas was correct to 5 percent.

So a single trajectory exact Einstein conductivity is close to meaningless, and the tool says so
rather than printing four significant figures without comment. Above 30 percent block error the
verdict tells you to quote the Nernst-Einstein value and state that the exact one could not be
resolved. The way to get it properly is replicates: several runs with different seeds, averaged.

One honest caution about this panel. The value shown above is from the mock, whose walkers are
uncorrelated by construction, which is why its ratio is nowhere near 1: that is the 24 percent
block error being honest.

The interpretation printed under the Haven ratio now depends on the number. Below 0.9 it says
the ions are not moving independently and that the shortfall is the correlation Nernst-Einstein
discards. Above 1.1 it warns that in a real melt this normally means the exact value has not
converged, and tells you to check the block error first. In between it says the approximation is
fine. Until 29 September 2026 the "below 1" sentence was printed whatever the ratio was.

For a real slag at 1600 C, a plausible conductivity is somewhere between 0.1 and 500 S per m,
which is the band the audit suite uses as a sanity check.

### One trap from their own scripts

Their conductivity script hard codes the Matsui and Guillot charge set while its own header
offers Bouhadja as an option. Run it on a Bouhadja trajectory and every charge is wrong by the
ratio 1.2698, so the conductivity is wrong by about 1.61 times. slagmd takes the charges from
the force field the run actually used, never from a literal.

## 8.5 Viscosity, two routes

This is the part their released work does not do at all. Searching all 84 files for `viscos`,
`green.?kubo`, `stress` and `ave/correlate` gives nothing: no pressure tensor compute, no stress
output, no integration. Viscosity appears in their papers only as a FactSage screening number.

It costs four lines in a production run that is happening anyway, which is what the stress
sampling tick box adds.

### Green-Kubo

Shear viscosity is the time integral of the autocorrelation of the off diagonal parts of the
pressure tensor, multiplied by volume over kB T. In words: when a fluid is sheared, the stress
takes some time to relax, and the longer it takes, the more viscous the fluid is. The integral
of "how long does a stress fluctuation remember itself" **is** the viscosity.

### The plateau test, in plain words

The integral has to be taken to infinite time, and you do not have infinite time. What you do in
practice is compute the **running integral**, plot it against the upper limit, and look for the
part where it stops changing. That flat part is the answer. If there is no flat part, there is no
answer.

Green-Kubo viscosity in a viscous melt is notorious for this. The autocorrelation has a long,
slowly decaying tail, the noise in the tail is integrated along with the signal, and the running
integral wanders instead of settling. Anyone who quotes a single Green-Kubo viscosity without
showing that it plateaued is quoting noise.

So slagmd returns the whole running integral, not one number, measures the scatter across the
window you asked for, and **refuses to call it converged if that scatter is over 10 percent**.

```
  Green-Kubo: 3.1797e-04 Pa s
      plateau found, scatter 1.5 percent across 2.0 to 10.0 ps. Still report it with an
      error bar from independent runs, not from this one window.
```

versus

```
  Green-Kubo: 8.4e-02 Pa s
      NOT converged. Across 1.0 to 10.0 ps the running integral scatters by 34 percent.
      Green-Kubo has not plateaued. Do not quote this number. Either run much longer,
      average over many independent trajectories, or use the Einstein-Helfand form which
      integrates the noise instead of differentiating it.
```

Even when it does plateau, the verdict tells you the error bar must come from independent runs
and not from that one window. That is correct and you should follow it.

### Einstein-Helfand

The same physics, integrated once. Define the running integral of the off diagonal pressure, and
its mean squared displacement grows linearly in time with a slope proportional to the viscosity.
This is to Green-Kubo what the Einstein relation is to a velocity autocorrelation.

It behaves better in practice, because integrating suppresses noise while differentiating
amplifies it. That is precisely why it is worth computing both: **they fail in different
directions, so agreement between them is real evidence.**

### The fit window, and a bug worth knowing about

The Einstein-Helfand fit window default is 0.2 percent to 2 percent of the series, which looks
absurdly early until you understand why.

The linear regime starts once the lag is much longer than the stress relaxation time, which in a
melt is well under a picosecond. It **ends** long before the series does, because at a lag equal
to a third of the trajectory there are hardly any independent starting points left and the
Helfand displacement is itself a random walk. Fitting at 20 to 60 percent of the series, which
is the obvious choice and was this function's first default, gave an answer **five times too
small** on a signal whose exact viscosity is known in closed form.

Worse, the log-log slope test could not see it, because the slope was being measured over the
whole series rather than over the fit window, so it returned the same number whatever window you
chose. It is now measured on the window, which is the one thing it is for.

So: always read the reported exponent.

| Reported exponent | Meaning | What to do |
|---|---|---|
| within 0.2 of 1 | Linear regime confirmed | Compare against Green-Kubo |
| below 1 | The window is too **late**, too few independent starting points | Move it earlier |
| above 1 | The window is too **early**, inside the stress relaxation time | Move it later |

### How to tell a good viscosity from a bad one

| Check | Good |
|---|---|
| Green-Kubo plateau scatter | under 10 percent |
| Einstein-Helfand exponent | within 0.2 of 1 |
| The two values against each other | Same order of magnitude. If they differ by a factor of two, neither is converged |
| Absolute value for a CaO-Al2O3-SiO2 melt at 1600 C | Of order 0.1 to 10 Pa s. Their own high throughput screen kept only compositions below 10 Pa s. A value of 1e-4 Pa s is water-like and means you are analysing noise, which is exactly what the mock produces |
| Number of stress samples | Tens of thousands. The tool refuses below 100 and explains the aliasing |

## 8.6 The equilibration report

### What it is

A simple question that almost nobody asks: was the box actually equilibrated when the
measurement was taken?

It matters because density is not a detail. Every transport coefficient scales with it. If the
density is still drifting during production, every diffusion coefficient and conductivity you
measured is at the wrong density, and no amount of averaging fixes that.

The report reads a LAMMPS log, pulls out every thermodynamic table, and for each of density,
temperature, volume, pressure and total energy takes the **second half** of the series, fits a
straight line, and reports the drift across that tail as a percentage of the mean.

### What it prints

```
  log.08_produce
      Density   mean       2.6485  sd     0.0014  drift   +0.14 %
      Temp      mean    1873.9743  sd     7.7007  drift   +0.17 %
      Volume    mean   20340.4977  sd    10.8883  drift   -0.14 %
      Press     mean     -57.6213  sd   299.4708  drift -516.89 %
      TotEng    mean -127178.5580  sd    64.6886  drift   +0.05 %
      density and temperature look settled over the tail
```

### How to read it

| Quantity | What to look at | Good |
|---|---|---|
| Density | The drift | Under 0.5 percent. Above that it says the box is not equilibrated, in those words |
| Temp | The standard deviation against the mean | Under 5 percent. In a 1402 atom box at 1873 K, 8 K of scatter is normal and healthy |
| Volume | Should mirror density | |
| Press | **Ignore the drift column** | Pressure in a small box fluctuates by hundreds of atmospheres around a mean near 1, so a drift expressed as a percentage of that mean is a meaningless number like minus 517 percent. Look at the mean, which should be near your target pressure, and at the standard deviation, which will be large |
| TotEng | The drift | Under about 0.1 percent. This is where a bad timestep shows up |

Two things to watch for. The report runs on **every** `log.*` beside the trajectory, including
the cooling stage, where the temperature is deliberately ramping from 4000 K to the target. It
will correctly report a large temperature spread there and tell you to check the thermostat
damping. That is not a problem with your run, it is the report being applied to a stage whose
temperature is meant to change. Read the log name.

And if you ticked the NVE check, look at the total energy drift in `log.09_nve_check`. A drift
larger than about 1e-4 per atom per ns means the timestep is too long or the Ewald accuracy too
loose, and every transport number above it is suspect.

# 9. Running on a real cluster

## 9.1 Where slagmd has to live

slagmd must run on a machine that can submit jobs, which means the cluster's **login node**. It
checks: if the backend is SLURM and `sbatch` is not on the PATH, it refuses to submit and tells
you why rather than failing later.

That has a practical consequence. The window is a graphical program, so you need either X
forwarding (`ssh -X`), or a remote desktop session, or you drive the project from Python instead
of the window. Everything the buttons do is a method on `Project`, so a headless session is a
few lines:

```python
import sys
sys.path.insert(0, "/path/to/AUTOMATION")
from slagmd.project import Project
p = Project("/scratch/myslag")
p.cluster.backend = "slurm"
p.cluster.partition = "cpu_rome"
p.cluster.account = "2026_045"
p.cluster.nodes = 2
p.cluster.tasks_per_node = 128
p.cluster.walltime = "72:00:00"
p.cluster.modules = ["cluster/dodrio/cpu_rome",
                     "LAMMPS/2Aug2023-foss-2022a-kokkos"]
p.cluster.lammps_exe = "lmp"
p.cluster.mpi_launcher = "srun"
p.save()
p.add_verbeeck_benchmark("bouhadja", temperatures_C=(1600,), target_atoms=5000)
for run in p.runs.values():
    p.prepare(run)
p.fill_queue(max_concurrent=4)
```

Then poll it from a loop, or reopen the window later and turn on the monitor. Reopening the
project folder picks up exactly where it stopped, because all the state is in `project.json`.

## 9.2 A realistic SLURM setup, end to end

The example below uses the sort of names a Flemish cluster uses, because that is the system the
published study ran on: the Flemish Supercomputer Center, two nodes of 128 cores each, two 64
core AMD EPYC 7H12 chips at 2.6 GHz per node. **Every one of these names is a placeholder. Get
the real ones from your own cluster's documentation.**

### Step 1. Find out what your cluster calls things

```bash
sinfo -o "%P %D %c %l"          # partitions, node count, cores per node, time limit
sacctmgr show assoc where user=$USER format=account,partition
module avail 2>&1 | grep -i lammps
```

slagmd runs the first of these for you and puts the result in the project's environment probe,
so the Project tab may already show your partitions.

### Step 2. Fill in the Cluster tab

| Field | Value in this example |
|---|---|
| scheduler | `slurm` |
| partition or queue | `cpu_rome` |
| account or project | `2026_045` |
| nodes | `2` |
| tasks per node | `128` |
| walltime hh:mm:ss | `72:00:00` |
| memory per node | blank |
| LAMMPS executable | `lmp` |
| launcher | `srun` |
| email for failure notices | your address |
| modules | `cluster/dodrio/cpu_rome` then `LAMMPS/2Aug2023-foss-2022a-kokkos`, one per line |
| extra scheduler directives | `--exclusive` |

Press **Check for problems**, then **Preview the submit script**.

### Step 3. Read the generated batch script

This is exactly what those settings produce for the production stage:

```bash
#!/bin/bash
#SBATCH --job-name=S1_1600C_08_produce
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=72:00:00
#SBATCH --output=log.08_produce.slurm.out
#SBATCH --error=log.08_produce.slurm.err
#SBATCH --partition=cpu_rome
#SBATCH --account=2026_045
#SBATCH --mail-user=you@example.com
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --exclusive

set -o pipefail
cd '/scratch/gent/vo/001/slag/runs/S1_1600C' || exit 90

module purge 2>/dev/null || true
module load cluster/dodrio/cpu_rome
module load LAMMPS/2Aug2023-foss-2022a-kokkos

echo "SLAGMD_HOST $(hostname)"
echo "SLAGMD_START $(date -Is)"

srun lmp -in in.08_produce -log log.08_produce
rc=$?
echo "SLAGMD_EXIT $rc"
echo "SLAGMD_END $(date -Is)"
exit $rc
```

Nine things in that script are deliberate.

| Line | Why it is there |
|---|---|
| `#SBATCH --output=...slurm.out` and `--error=...slurm.err` | Separate from the LAMMPS log, because a module error or an out of memory kill never reaches the LAMMPS log. The diagnosis reads both |
| `--mail-type=FAIL,TIME_LIMIT` | Only failures. Nobody wants 754 success emails |
| `set -o pipefail` | So a failure in a pipeline is not hidden by a later success |
| `cd '...' \|\| exit 90` | **Quoted**, because this project's own path contains a space. Exit code 90 is distinctive, so a failure to change directory is unmistakable |
| `module purge` before the loads | Clusters preload things. A stale module is a common cause of a mysterious crash |
| `2>/dev/null \|\| true` on the purge | Some sites have no module system at all, and the script must not die on that |
| `SLAGMD_HOST` and `SLAGMD_START` | So a failure can be traced to a node. When `MPI_CRASH` fires twice on the same node name, that is a hardware ticket, not your bug |
| `rc=$?` then `exit $rc` | The script exits with LAMMPS's own exit code, so the scheduler records the truth |
| `SLAGMD_EXIT $rc` echoed | Because the exit code alone is not enough. LAMMPS can exit zero having done nothing, which is why the real test is the completion marker inside the LAMMPS log |

### Step 4. Prepare and launch

Tab 3, **Add the ten published benchmark compositions**, then **Prepare every unprepared run**.
Check the notes for `CLUSTER WARNING` or `PROTOCOL WARNING` lines. Then tab 4: set **jobs at
once**, set **poll every** to 60 seconds, and tick **watch and repair automatically**.

Leave it. Close the window if you like, and reopen it later.

### Step 5. Choose the walltime honestly

Their own figure is the anchor: the whole multi stage classical workflow, about 5000 atoms, took
about 10 hours on 2 nodes of 128 cores, which is 2560 core hours for roughly 107 ns. That is what
`est_core_hours` scales from, linearly in atoms and in simulated time.

Split across the stages, production is about half the total, so plan for something like five
hours of wall time for production on that hardware, and ask for a good deal more than you plan
for. A walltime that is too short is the cheapest failure in the table, because the restart file
means nothing is lost and the tool doubles it for you. A walltime that is too long only costs you
queue priority.

If your queue has a maximum walltime below what a stage needs, doubling will not help. Split the
stage instead: shorten `produce_ns`, or shorten each equilibration block and add more of them.

## 9.3 What changes at a different site

| Thing to change | SLURM | PBS |
|---|---|---|
| Backend | `slurm` | `pbs` |
| Partition field becomes | `--partition` | `-q`, the queue name |
| Account field becomes | `--account` | `-A` |
| Nodes and tasks become | `--nodes` and `--ntasks-per-node` | one `-l select=N:ncpus=M:mpiprocs=M` line |
| Memory becomes | `--mem` | `:mem=` inside the same select line |
| Launcher | usually `srun` | usually `mpirun`, and slagmd adds `-np <total ranks>` for you |
| Status is read with | `squeue -h -j <id>` with a bar separated format string, falling back to `sacct -j <id> -n -P -X` | `qstat -f <id>`, parsing `job_state` |
| Cancel is | `scancel` | `qdel` |

The PBS version of the same script looks like this:

```bash
#!/bin/bash
#PBS -N S1_1600C_08_produce
#PBS -l select=2:ncpus=64:mpiprocs=64:mem=240gb
#PBS -l walltime=48:00:00
#PBS -o log.08_produce.pbs.out
#PBS -e log.08_produce.pbs.err
#PBS -q q_long
#PBS -A PROJ123

set -o pipefail
cd '/home/u/slag/runs/S1_1600C' || exit 90

module purge 2>/dev/null || true
module load lammps/2024.06

echo "SLAGMD_HOST $(hostname)"
echo "SLAGMD_START $(date -Is)"

mpirun -np 128 lmp -in in.08_produce -log log.08_produce
rc=$?
echo "SLAGMD_EXIT $rc"
echo "SLAGMD_END $(date -Is)"
exit $rc
```

Three site specific things that are **not** covered and that you may need to add by hand in
**extra shell lines before LAMMPS** or **extra scheduler directives**:

1. **A launcher that is not srun, mpirun or mpiexec.** Their own README uses `mympirun`, a
   Flemish cluster wrapper that derives the rank count from the scheduler environment. Type it
   into the launcher field and slagmd will append `-np <total ranks>`, which `mympirun` does not
   need but tolerates. If that is a problem, put the whole command in the LAMMPS executable field
   and leave the launcher blank.
2. **A container.** Put the `singularity exec ...` or `apptainer exec ...` prefix into the LAMMPS
   executable field. If it contains a space and no quotes, the checker will warn you, and the
   warning is worth heeding.
3. **Scratch filesystems.** If your cluster requires the working directory to be on a fast
   parallel filesystem, put the whole project there. The batch script changes directory to the
   run folder, so the project root is where the input and output live. A common failure mode is
   a login node home directory that compute nodes cannot see, and that shows up as
   `MISSING_INPUT`.

## 9.4 Which LAMMPS build you need

The generated deck needs the **KSPACE** package, which provides `ewald` and the
`born/coul/long` and `buck/coul/long` pair styles. That requirement is written into the header of
every generated input file. Their own scripts carry the same note.

If the module you loaded was built without it, the job dies with `MISSING_PACKAGE` and the tool
stops instead of retrying, because retrying the same deck against the same build will fail
identically.

Their runs used LAMMPS's 27 June 2024 version. Any reasonably recent version should behave the
same for these pair styles, but the version is worth recording, and `protocol.json` records the
force field, the charges and the full protocol beside every run.

# 10. Extending it

## 10.1 Adding a force field

Force fields live in `slagmd/forcefields.json`. Nothing is compiled, so adding one is editing a
file.

### The entry

```json
"newff": {
  "label": "Author (year), form, what it was fitted to",
  "form": "bmh",
  "verified": true,
  "source": "the printed table this came from, named exactly",
  "charges": {"Ca": 1.2, "Al": 1.8, "Si": 2.4, "O": -1.2},
  "charge_source": "the sentence or table the charges came from",
  "identity_check": {
    "rule": "the algebraic identity these numbers must satisfy, in words",
    "tolerance": 0.001,
    "result_2026_09_29": "PASSES for all 10 pairs"
  },
  "pairs": {
    "Al-Al": {"A": 0.2798, "sigma": 1.5704, "rho": 0.068, "C": 1355.62},
    "Al-O":  {"A": 0.7237, "sigma": 2.6067, "rho": 0.164, "C": 3336.01}
  }
}
```

| Field | Rules |
|---|---|
| The key | Lower case. A key starting with an underscore is metadata and is never offered as a force field |
| `label` | Shown in the interface. Say what it was fitted to, because that is what decides whether it suits your melt |
| `form` | `bmh` for A exp((sigma - r) / rho) minus C over r^6, which LAMMPS calls `born/coul/long`. `buck` for A exp(minus r over rho) minus C over r^6, `buck/coul/long`. There is no third option |
| `verified` | **Set this to `true` only when every number has been checked against a printed table.** With `false`, loading it for a real deck is refused outright |
| `source` | A real citation. This string is printed into the header of every generated input file, so a deck can always be traced back |
| `charges` | One per species. A species with no charge here is refused rather than given a default |
| `pairs` | Key is the two species joined by a hyphen. The order does not matter, the loader sorts it, so `O-Si` and `Si-O` are the same key |
| `cation_cation` | **Presence of this key** tells the tool that cation-cation short range terms are deliberately absent, as in Guillot. Its value is a note for a human. Without it, a missing pair is a refusal |

### Units

| Quantity | Unit in the JSON |
|---|---|
| `A` | kJ per mol, exactly as printed in the source table |
| `C` | Angström^6 kJ per mol |
| `sigma`, `rho` | Angström |

The tool converts A and C to kcal per mol for LAMMPS's `units real` by dividing by 4.184, which
is exact by the definition of the thermochemical calorie. Do not pre-convert. Put the printed
values in.

### The consistency identity, and one thing to be careful about

The point of an identity is that it catches a typing mistake without needing a second copy of
the table. Matsui's is that A equals 4.184 times rho, because his own functional form makes them
the same number. Guillot's is that its cation-oxygen A equals Matsui's A times the exponential of
sigma over rho, because Buckingham is Born-Mayer-Huggins rewritten.

**A new force field cannot slip through unchecked, but you still have to give it a check.**

`ForceField.check_self_consistency` knows which force fields have an implemented check. If you
add a name it does not know, it reports one of two problems rather than passing silently:

- if your JSON **declares** an `identity_check` block but nothing implements it, it says so and
  tells you to implement it or remove the block;
- if there is **no** check at all, it says the numbers have been transcribed but never
  cross-checked, and to treat them as unconfirmed.

Either way the message comes back from `check_self_consistency`, and `check_ready` refuses to
build a deck. Until 29 September 2026 a new force field simply passed, while its JSON appeared
to declare a check, which quietly undercut the whole verified-parameters guarantee.

**The strongest check is not algebra, it is the authors' own file.** If the people who published
the force field also released a LAMMPS input, put their `pair_coeff` values into a
`released_pair_coeff_kcal` block and the tool will compare against them directly. That is what
`bouhadja` does, and it catches a single wrong digit in any of the forty numbers, including in
parameters no algebraic identity can constrain. Prefer it whenever a released file exists.

If you only have an algebraic identity, add a branch:

```python
if self.name == "newff":
    for k, p in self.pairs.items():
        expect = <whatever the identity says>
        if abs(p["A"] - expect) > tol:
            problems.append(f"newff {k}: A={p['A']} but identity gives {expect:.4f}")
```

and then add a check to group A of `tests/test_slagmd.py` that restates it independently, so the
identity is asserted twice by two different pieces of code.

If your source has no such identity, the next best thing is what was done for Bouhadja: compare
every number against a **second** independent copy, in that case the `pair_coeff` lines in their
released input files. That comparison is 40 numbers and it lives in group D of the audit.

### Choosing a tolerance honestly

Set the tolerance from the number of digits your source prints, not from what makes the test
pass. The Guillot identity is tested at 0.5 percent because Matsui's A is printed to four decimal
places and the exponential factor is about four million, so the last printed digit moves the
answer by a couple of tenths of a percent. The measured gap was 0.21 to 0.22 percent for all
three pairs: same size, same sign, which is rounding. A real typing mistake is a transposed digit,
tens of percent or more, and still fails at 0.5 percent. Write that reasoning into a comment, as
the existing code does, so the next person does not loosen it further without thinking.

## 10.2 Adding an oxide

Four places in `slagmd/composition.py`, and one more in `slagmd/analysis.py`.

### 1. `OXIDES`

```python
OXIDES: dict[str, tuple[str, int, int]] = {
    "BaO": ("Ba", 1, 1),     # cation symbol, cations per formula unit, oxygens per unit
}
```

For an oxide where the same element appears in two oxidation states, use a distinct cation
symbol, as iron does with `Fe2` for Fe(II) and `Fe3` for Fe(III). The charge differs, so they
cannot share a species.

### 2. `ATOMIC_MASS`

```python
ATOMIC_MASS = {"Ba": 137.327, ...}
```

IUPAC 2021 standard atomic weights, in g per mol. This dictionary does double duty: its **key
order** sets the order species appear in, which sets the LAMMPS atom type numbers. Adding a
species in the middle renumbers the types of everything after it, which will not match any
trajectory you have already produced. Add at the end.

### 3. `FORMAL_CHARGE`

```python
FORMAL_CHARGE = {"Ba": 2.0, ...}
```

Used only to check that a force field's partial charges are a uniform scaling of the formal
charges. If you leave it out, that check returns false and you lose the strong neutrality
statement.

### 4. `NETWORK_FORMERS`, if it is one

```python
NETWORK_FORMERS = ("Si", "Al", "Fe3", "Ti")
```

Only add a cation here if it genuinely carries bridging oxygen. This tuple decides the oxygen
speciation classes, the Qn counting, and the ideal NBO per tetrahedron.

### 5. Nothing else

`NETWORK_FORMERS` is defined once, in `composition.py`, and `analysis.py` imports it. Editing
the one tuple changes the composition panel and the speciation analysis together. Until 29
September 2026 there were two copies and they could drift apart, which is why this section
used to tell you to edit both.

### Then

The force field must have a charge and every pair for the new species, or the tool refuses with
the message in section 6.2. That is the point: adding an oxide to the list does not make it
runnable, it makes it *askable*, and the force field then has to be able to answer.

## 10.3 Adding a failure pattern

Two files, and the second one is what makes the first trustworthy.

### 1. Add a rule to `diagnose.py`

The table in `_rules()` is a list of `(regex, builder)` pairs, matched **in order**, so put a
specific pattern above a general one.

```python
(r"ERROR.*Bond atoms .* missing|Bond atom missing",
 lambda m: Remedy(
     code="BOND_MISSING",
     what_happened="LAMMPS could not find an atom that a bond refers to.",
     why="Plain English. What actually went wrong physically or mechanically. "
         "This is what somebody reads at 2am, so write it for them.",
     action="What the tool is about to do, or what the human must do.",
     protocol_changes={"timestep_fs": 1.0},
     cluster_changes={"nodes": "x2"},
     rebuild_structure=False,
     resume_from_stage=None,
     auto_retry=True,
     plain="one line for the queue table")),
```

| Field | What it controls |
|---|---|
| `code` | The short name. Appears in the queue note, the history and the diagnosis file name |
| `what_happened` | One or two sentences. No jargon |
| `why` | The cause. The audit requires this to be more than 40 characters, on purpose |
| `action` | What will be done. The audit requires more than 20 characters |
| `protocol_changes` | A dict merged into the run's protocol overrides. Any field of `Protocol` works. Two special keys exist: `_scale_eq_pdamp` and `_scale_eq_tdamp` multiply the damping of **every** equilibration block by a factor |
| `cluster_changes` | A dict merged into the run's cluster overrides. The string `"x2"` doubles the current value, and it understands `hh:mm:ss` walltimes as well as integers. Walltime doubling is capped at seven days |
| `rebuild_structure` | True rebuilds the starting structure from scratch. Use it when the failure suggests the structure itself was bad |
| `resume_from_stage` | The stage to go back to. `None` means retry the same stage. Setting it rolls that stage and everything after it back to not submitted |
| `auto_retry` | **Set this to False whenever a retry cannot possibly help.** A missing package, a wrong module name, a malformed pair coefficient, a missing input file |
| `plain` | One short line, lower case, for the queue table |

If `protocol_changes` or `rebuild_structure` is set, the whole deck is rewritten before the
resubmission, and stages that had already completed keep their completed status unless the
structure was rebuilt.

### 2. Add the verbatim message to `mocklammps.py`

```python
FAIL_MESSAGES = {
    "BOND_MISSING": ("ERROR: Bond atoms 4012 4013 missing on proc 7 at step 120000 "
                     "(src/ntopo_bond_all.cpp:59)"),
}
```

This is not optional and it is not decoration. The audit suite loops over every entry in
`FAIL_MESSAGES`, feeds it to the diagnoser, and requires that the right code comes back with a
cause and an action. So adding your message here **automatically** adds a test for your rule, and
the key must match the `code` exactly.

Copy the message from a real log. Do not write it from memory. A diagnosis rule that has never
seen its own error message is a guess, and the whole value of this table is that each pattern has
been matched against text LAMMPS or SLURM actually emits.

### 3. Run the audit

```bash
python3 tests/test_slagmd.py
```

The count goes up by two, to 211, and both new checks should pass. If your regex is too greedy it
will shadow a later rule and one of the existing checks will fail instead, which is exactly what
you want to find out.

### If a failure should stop rather than retry

Set `auto_retry=False` and add the code to the loop in group F of the audit that asserts the tool
treats it as needing a human:

```python
for code in ("BAD_PAIR_COEFF", "MISSING_PACKAGE", "ENVIRONMENT", "MISSING_INPUT"):
```

## 10.4 Changing the protocol

The protocol is the `Protocol` dataclass in `inputgen.py`. Its defaults **are** their published
protocol, and the fields are labelled in the source with `[THEIRS]` or `[MINE]`. Change them per
run rather than globally, by passing `protocol_overrides` to `add_run`, so their protocol stays
the reference.

```python
p.add_run("test_short", {"CaO": 0.5, "SiO2": 0.5}, 1873.15, "bouhadja",
          protocol_overrides={"produce_ns": 5.0, "hot_ns": 1.0})
```

The equilibration cascade is a list of `EqBlock(tdamp_fs, pdamp_fs, ns, label)`. Replacing it
with fewer than two blocks triggers a warning, and replacing it with one block of fixed damping
triggers a different one, both of which say why in the panel. If you edit `THEIR_EQ_CASCADE`
itself, edit the table in the module docstring in the same commit, because the audit checks the
arithmetic and a person checks the documentation.

# 11. Troubleshooting

Refusals are in section 6 and job failures are in section 7. This section is for the things that
are neither.

| Symptom | Likely cause | What to do |
|---|---|---|
| `python3 run_slagmd.py` fails with `ModuleNotFoundError: No module named 'tkinter'` | tkinter is packaged separately on some Linux distributions | Install `python3-tk`, or drive the project from Python as in section 9.1 |
| `_tkinter.TclError: no display name and no $DISPLAY environment variable` | No graphical display, typical over plain `ssh` | Use `ssh -X`, or a remote desktop, or drive it from Python |
| A button appears to do nothing | The work is on a background thread and something threw | Look at tab 6, Log. Full tracebacks land there |
| The window is responsive but the progress bar never stops | A worker is genuinely still running. A full analysis of a 400 frame, 1400 atom trajectory takes several minutes | Watch the status bar at the bottom left, which names the property being computed |
| `no project open` warnings on every button | You opened the window but never pressed **Open or create** | Go to tab 1 and open a project |
| A local job shows `RUNNING` forever | This was a real bug and is fixed. A finished child process becomes a zombie until it is waited on, and a zombie still answers a liveness probe, so the naive test reported a job that ended minutes ago as running and the whole queue stalled. The status routine now reaps the process first, and falls back to reading `/proc` when a different session started it | If you still see it, the process really is alive |
| `local process has exited, read the log to see how` | A local job ended and this session has no handle on it, because a different session started it | Read the log. The status is honestly reported as unknown rather than guessed |
| `<stage> has no completion marker, so <next stage> cannot start. Rolled <stage> back to be resubmitted.` | A stage was marked complete by the scheduler but its own marker is not in the log | Nothing. This is the tool protecting you. It will rerun the stage rather than build on a stage that did not really finish |
| A stage sits at `PENDING` for hours | Normal on a busy cluster | `squeue -j <id>` will say why, often `Priority` or `Resources` |
| Everything is `BLOCKED` with `submission failed, sbatch is not on PATH` | slagmd is not running on a submit host | Move to the login node, or set the backend to `local` |
| The analysis says there is no trajectory | Production has not finished, or you are looking at the wrong run | Check the Queue tab. Or browse to a trajectory file brought back from a cluster by hand |
| Viscosity says there is no `stress.txt` | The stress sampling box was not ticked when the run was added | Ticking it now will not help that run, because the production stage has already run. Add a new run with it ticked |
| `Make the figures` does nothing and the Log shows a matplotlib import error | matplotlib is not installed | Install it. Everything else works without it |
| The RDF first peak is absurdly small, and coordination numbers are zero | You are analysing a mock trajectory, or a dump whose coordinates do not belong to the box in its header. The tool now says so in a warning | Section 3, the last table. The mock's structure numbers are meaningless by construction |
| Numbers differ slightly between two runs of the audit | Several checks use random data and the mock draws its seed from the stage name | Expected. The tolerances are set with that in mind. A **failure**, rather than a different digit, is real |
| A generated shell script fails with `export: not a valid identifier` | An unquoted path containing a space, several lines above the error | This is the trailing space trap. **Check for problems** on the Cluster tab catches it. Quote the value |
| A path fails in bash with `No such file or directory` and the path looks right | Same cause. This folder's own path ends in a space | Quote every path, always |
| `project.json` looks stale | It should not. It is written after every change, through a temporary file and a rename | If it really is stale, check `events.log`, which is append only and independent |
| Two windows open on the same project | Both write `project.json` and the last writer wins | Do not. There is no lock |
| Analysis is very slow | The radial distribution function compares every pair in every frame. That cost grows with the square of the atom count and linearly with frames | Analyse fewer frames. `read_lammps_dump` takes `stride` and `max_frames` from Python. There is no control for this in the window |
| The equilibration report complains about temperature on a cooling stage | It is reporting on `log.02_cool`, where the temperature is meant to be ramping from 4000 K to the target | Read the log name above the numbers |
| The pressure drift is reported as minus 500 percent | Pressure in a small box fluctuates by hundreds of atmospheres around a mean near 1 atm, so a drift as a percentage of that mean is meaningless | Ignore that one column. Section 8.6 |

## 11.1 If a number looks wrong

Work through this in order.

1. **Was it equilibrated?** Press `Was it equilibrated?` and read the density drift for the
   production log. If the density was moving, nothing else matters.
2. **Was it diffusive?** Read the log-log slope. Not within 0.15 of 1 means there is no diffusion
   coefficient to argue about.
3. **What is the error bar?** Over 25 percent means the number is an estimate.
4. **Do the two independent routes agree?** The generated deck computes the mean squared
   displacement inside LAMMPS as well, into `msd_lammps.txt`, specifically so it can be compared
   against the Python. If they disagree the Python is wrong, and that is worth knowing. For
   viscosity, compare Green-Kubo against Einstein-Helfand, which fail in different directions.
5. **Is the structure sane?** Si-O at 1.6 Angström with a coordination number of 4. If that is
   wrong the melt is wrong, and no transport number from it means anything.
6. **Run the audit.** `python3 tests/test_slagmd.py`. If the tool's own arithmetic has broken,
   it will say so.
7. **Compare against a published benchmark.** Queue one of the ten benchmark compositions with
   `bouhadja` at 1600 C and compare the self diffusion coefficients against the published ones.
   That is the acceptance test, and it has not yet been run against a real LAMMPS.

## 11.2 After a crash or a dropped connection

Nothing special is required.

1. Reopen the same project folder.
2. Press **Poll now**. Every stage's real state is read back from the scheduler and from the
   logs.
3. Turn the monitor back on.

Jobs already on the cluster were never affected by the window closing. Stages that finished while
the window was shut are detected by their completion markers. Stages that failed while the window
was shut are diagnosed on the first poll.

# 12. Glossary

Every term used in this manual, in alphabetical order.

| Term | Meaning |
|---|---|
| **Account** | The project code a cluster bills core hours to. Many clusters refuse a job without one |
| **Angström** | 1e-10 m. The natural unit of atomic distance. A silicon to oxygen bond is about 1.6 of them |
| **Attempt cap** | The limit on how many times slagmd will resubmit one stage. Three |
| **Autocorrelation** | How much a fluctuating quantity resembles itself a time t later. Its integral is what Green-Kubo turns into a transport coefficient |
| **Ballistic regime** | Very short times, before an atom has hit a neighbour. The mean squared displacement grows as t squared |
| **Barostat** | The part of the algorithm that holds the pressure at its target by changing the box volume |
| **Batch script** | The shell script handed to a scheduler. The `#SBATCH` or `#PBS` lines are the resource request |
| **BMH, Born-Mayer-Huggins** | A pair potential of the form A exp((sigma minus r) over rho) minus C over r^6, plus Coulomb. LAMMPS calls it `born/coul/long` |
| **Block average** | Cutting a trajectory into pieces that share no starting points, analysing each, and using the spread as an error bar |
| **Bridging oxygen** | An oxygen bonded to exactly two network formers. The links that hold a silicate network together |
| **Buckingham** | A pair potential of the form A exp(minus r over rho) minus C over r^6, plus Coulomb. The same thing as BMH with the prefactor absorbed. LAMMPS calls it `buck/coul/long` |
| **Cage regime** | Times short enough that an atom is still rattling inside a shell of neighbours. The mean squared displacement flattens. **Not** diffusion |
| **CMAS** | Shorthand for the CaO-MgO-Al2O3-SiO2 system, the chemistry of most slags |
| **Coordination number** | How many neighbours of a given type are inside the first shell. Four for silicon in a silicate |
| **Core hour** | One processor core running for one hour. The currency of a cluster allocation |
| **Coulomb cutoff** | The distance beyond which the real space part of the electrostatic sum is handed to the reciprocal space part. 11 Angström here |
| **Damping time** | How quickly a thermostat or barostat reacts. Short means tight and intrusive, long means loose and slow |
| **Dispersion** | The weak attraction between any two atoms, the minus C over r^6 term |
| **Effective charge** | A partial charge smaller than the formal one, used because electrons in a real silicate are partly shared. Here 0.4725 or 0.6 of formal |
| **Einstein relation** | Mean squared displacement equals 6 D t at long times. The definition of the self diffusion coefficient |
| **Einstein-Helfand** | Viscosity from the mean squared displacement of the time integral of the stress. The integrated form of Green-Kubo, and better behaved because integrating suppresses noise |
| **Equilibration** | Running the system until its properties stop changing, before measuring anything |
| **Ewald summation** | The standard way to compute electrostatics in a periodic box: a short ranged part in real space and a long ranged part in reciprocal space |
| **Formal charge** | The textbook oxidation state. Plus 4 for silicon, minus 2 for oxygen |
| **Force field** | The formula and the numbers that give the energy of a pair of atoms as a function of their separation. Decided before the simulation, not during it |
| **Free oxygen** | An oxygen bonded to no network former, surrounded only by modifier cations. Rare in a silicate melt |
| **g(r)** | See radial distribution function |
| **Green-Kubo** | Getting a transport coefficient from the time integral of an autocorrelation function. Used here for viscosity |
| **Haven ratio** | The exact conductivity divided by the Nernst-Einstein estimate. Below 1 means the ions do not move independently |
| **Liquidus** | The temperature above which a composition is fully liquid. Below it you have crystals, or a supercooled liquid, and transport coefficients stop converging |
| **LAMMPS** | The molecular dynamics program that does the physics. slagmd generates its input and reads its output |
| **Login node** | The machine you land on when you connect to a cluster. Where you submit jobs from, and not where you run them |
| **Minimum image convention** | In a periodic box, interact with the nearest copy of each neighbour. Unambiguous only while the interaction range is under half the box edge |
| **Modifier** | A cation that breaks the network rather than forming it: Ca, Mg, Na, K, Fe(II) |
| **Module** | A cluster's mechanism for putting one specific version of a program on your PATH. `module load LAMMPS/...` |
| **MSD, mean squared displacement** | The average squared distance an atom has moved in a time t, averaged over atoms and over starting points |
| **NBO/T** | Non bridging oxygens per network former. A one number summary of how broken up a melt is |
| **Nernst-Einstein** | Electrical conductivity estimated from the individual diffusion coefficients, assuming the ions move independently |
| **Network former** | A cation that carries bridging oxygen and builds the network: Si, Al, Fe(III), Ti |
| **Node** | One physical computer in a cluster, with its own memory |
| **Non bridging oxygen** | An oxygen bonded to exactly one network former. A broken link in the network |
| **NPT** | Constant number of atoms, pressure and temperature. The ensemble used throughout here, at 1 atm |
| **NVE** | Constant number of atoms, volume and energy. Used here only as an optional check: if energy drifts, the timestep or the Ewald accuracy is wrong |
| **OVITO** | The visualisation program their workflow uses for six manual steps of structure preparation. slagmd replaces those with code |
| **Packmol** | A program that packs molecules into a box without overlaps. Their way of building the starting structure |
| **Partition**, **queue** | A named slice of a cluster with its own rules and limits |
| **PBS** | One of the two schedulers slagmd speaks. Directives start `#PBS` |
| **Periodic boundary conditions** | The box is surrounded by copies of itself, so a few thousand atoms stand in for a bulk liquid |
| **Qn** | For a network former, the number of bridging oxygens it carries. Q4 is fully connected, Q0 is isolated |
| **Radial distribution function** | g(r). How many neighbours are at distance r, relative to a structureless liquid of the same density. Its first peak is a bond length |
| **Rank** | One process of a parallel program. See task |
| **Remedy** | In slagmd, a diagnosed failure's prescription: what happened, why, what to change, and whether to retry |
| **Restart file** | A binary snapshot LAMMPS writes so a later stage can pick up exactly where an earlier one stopped |
| **Scheduler** | The program that owns a cluster's queue and decides when your job runs. SLURM or PBS |
| **Self diffusion coefficient** | D. How fast one species wanders through the melt, in m^2 per s |
| **Slag** | A molten oxide mixture. Here CaO-Al2O3-SiO2, the chemistry of steelmaking and of metal recovery |
| **SLURM** | One of the two schedulers slagmd speaks. Directives start `#SBATCH` |
| **Speciation** | The classification of every oxygen into bridging, non bridging, free or triclustered |
| **Stage** | One LAMMPS invocation in the chain. slagmd generates nine, from minimise to production |
| **Stoichiometry** | The whole number ratios a formula demands. CaO needs one oxygen per calcium, Al2O3 needs three per two aluminium |
| **Stress tensor** | The pressure written as a matrix. Its off diagonal parts are what a shear viscosity is computed from |
| **Task** | One process of a parallel program. LAMMPS splits the box into as many pieces as there are tasks |
| **Thermostat** | The part of the algorithm that holds the temperature at its target |
| **Timestep** | How far the simulation advances in one step. 2 fs here, which is about a tenth of a silicon to oxygen vibration period |
| **Triclustered oxygen** | An oxygen bonded to three or more network formers. Real in alumina rich melts |
| **Trajectory** | The file of atomic positions written every so many steps. Here `prod.lammpstrj`, every 10 ps |
| **Unwrapped coordinates** | Positions that keep counting up as an atom crosses the box boundary, `xu yu zu` in a LAMMPS dump. What a mean squared displacement must use |
| **Walltime** | The clock time limit on a job. Exceed it and the scheduler kills the job, however much work is left |
| **Wrapped coordinates** | Positions folded back into the box, `x y z`. Correct for structure, wrong for displacement |

# Provenance

**Theirs.** The protocol, the five stages, the damping cascade, all three force fields with their
charges, the Ewald settings, the 11 Angström cutoff, the forced splitting parameter, the tail
correction, the overlap deletion, the minimisation tolerances, the 200 bin radial distribution
function, the fixed 2.5 Angström oxygen coordination cutoff, the ten benchmark compositions, the
group naming, the velocity seed, and the hardware figures the cost estimate scales from. Sources:

- H. Verbeeck, I. Bellemans et al., *Molecular dynamics studies of the structural and transport
  properties of CaO-Al2O3-SiO2 melts: comparing the accuracy of the empirical force fields*,
  Computational Materials Science **259** (2025) 114150. Scripts: Mendeley Data
  `10.17632/v7m3k6yksw.1`.
- H. Verbeeck et al., *Self-diffusion coefficients in CaO-Al2O3-SiO2 melts: a high-throughput
  molecular dynamics approach*, Materials & Design **264** (2026) 115709. Scripts: Mendeley Data
  `10.17632/gkcjjk744s.1`.
- Original force fields: M. Matsui, Mineralogical Magazine **58A** (1994) 571; B. Guillot and
  N. Sator, Geochimica et Cosmochimica Acta **71** (2007) 1249; M. Bouhadja, N. Jakse and
  A. Pasturel, Journal of Chemical Physics **138** (2013).

Both datasets were downloaded on 29 September 2026, 84 of 84 files, verified against the
published manifest sizes. The detailed reading of every script is in
`../research/A5_verbeeck_scripts.md`.

**Ours.** The job submission layer, the scheduler abstraction, the composition resolver, the
structure builder, the structure checks, the box against cutoff rule, the failure diagnosis
table, the completion marker convention, the resumable state machine, the block error bars, the
log-log regime test, the Green-Kubo plateau test, the Einstein-Helfand route, the equilibration
report, the stress sampling, the LAMMPS side mean squared displacement cross check, the optional
NVE check, the mock LAMMPS, the audit suite, and the window.

**Their defects, handled and recorded so nobody rediscovers them.**

| Defect | Effect | What slagmd does |
|---|---|---|
| `production.file` ends with `unfix` calls for fixes it never defines | LAMMPS aborts after the physics is written, so the job exits non zero and any dependency chain breaks. Their published results are unaffected | Generates the teardown from the fixes it actually created, so the class of bug cannot recur. The audit checks that every `unfix` has a matching `fix` |
| Their conductivity script hard codes the Matsui charge set while offering Bouhadja in its header | Conductivity out by about 1.61 times if used as labelled | Takes charges from the force field the run actually used |
| The Packmol example is not charge neutral, 295 oxygens where stoichiometry needs 290 | A net minus 6 e cell, which Ewald would smear into a neutralising background | Counts oxygen from the cations, never from a literal |
| Four stale run length comments, for example `run 2500000 #2000 ps = 4 ns` where the truth is 5 ns | Misleading to a reader. The step counts are right | Generates every comment from the step count, so the two cannot drift |
| The Guillot script's citation comment points at Bouhadja's paper | A copy and paste slip in a comment | Each force field carries its own verified source string |
| The velocity seed 458165 is identical in all four copies | Repeated runs of one composition are the same trajectory, not independent replicates | Exposes the seed per run |

**Five bugs the audit found in slagmd itself, all fixed.** Listed because each of them is now a
test, and because a tool that claims no bugs has not been tested.

1. The radial distribution function peak finder returned the **second** coordination shell.
   Caught by a simple cubic lattice, where the answers are 3.0 Angström and 6 and the tool was
   reporting 4.245 Angström and 18.
2. The Einstein-Helfand viscosity was five times too small, because the fit window was too late
   in the series, and the log-log test was measured over the whole series so it could not see it.
3. A finished local job showed as `RUNNING` forever, because a zombie process still answers a
   liveness probe.
4. The displacement routine grew with the square of the frame count and was unusable on a stress
   series. Replaced with the standard FFT algorithm, checked against the direct sum.
5. The trailing space in this project's own path broke the generated shell script. Now a check.

**Status.** Complete and audited on 29 September 2026: 221 checks, 0 failures. Never yet run
against a real LAMMPS. The first real use should be one of the ten benchmark compositions with
`bouhadja` at 1600 C, so the self diffusion coefficients can be compared with the published ones.
That is the acceptance test.
