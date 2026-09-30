"""Generate the LAMMPS input deck for the melt protocol.

This version is NOT a reconstruction from the papers. The real scripts were downloaded
from the two Mendeley Data releases, so every setting below is taken from a file the
authors published, and the provenance is marked in the generated deck:

  [THEIRS]  copied from their released input files
  [PAPER]   stated in the printed paper
  [MINE]    slagmd's own addition, not theirs, never to be cited as theirs

Datasets:
  10.17632/v7m3k6yksw.1  force field comparison scripts   (Comput. Mater. Sci. 259, 114150)
  10.17632/gkcjjk744s.1  high throughput self diffusion   (Mater. Des. 264, 115709)

Their protocol, recovered from cooling.file, equilibration_step1 to 5 and production.file,
with dt = 2 fs throughout, damping written in the scripts as multiples of dt:

  stage                 Tdamp     Pdamp     steps        ns
  hot 4000 K NPT        100 fs    200 fs     2,500,000     5.0
  cool to target        100 fs    200 fs     1,063,500     2.127   (1 K/ps)
  equilibrate 1          50 fs    100 fs     5,000,000    10.0
  equilibrate 2         100 fs    100 fs     5,000,000    10.0
  equilibrate 3         200 fs    200 fs     5,000,000    10.0
  equilibrate 4         200 fs    800 fs     5,000,000    10.0
  equilibrate 5         500 fs   1000 fs     5,000,000    10.0
  production           1000 fs   2000 fs    25,000,000    50.0

Five equilibration blocks of 10 ns, which is the 50 ns the paper quotes, with the
thermostat and barostat loosened at every step. That schedule is the thing the papers push
to Supplementary Material, so the release is its only public source.

Three defects in their released scripts are worked around here, each marked in the code:
  1. production.file calls unfix on RDF and coordination fixes it never defined, so every
     production job exits with an error even though the physics finished.
  2. MSCD.py hard codes the Matsui charge set while its header offers Bouhadja, which
     would put the conductivity out by roughly the square of 2.4/1.89.
  3. The Packmol example is not charge neutral, 295 oxygens where the stoichiometry needs
     290, and its cell matches no production run.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

from .composition import Composition
from .forcefield import ForceField

DS_FF = "Mendeley Data 10.17632/v7m3k6yksw.1"
DS_HT = "Mendeley Data 10.17632/gkcjjk744s.1"
PAPER_FF = "Verbeeck et al., Comput. Mater. Sci. 259 (2025) 114150"
PAPER_HT = "Verbeeck et al., Mater. Des. 264 (2026) 115709"


@dataclass
class EqBlock:
    """One equilibration block: temperature held, damping in fs, length in ns."""
    tdamp_fs: float
    pdamp_fs: float
    ns: float
    label: str = ""


# Their cascade, verbatim. Do not edit these without editing the note above.
THEIR_EQ_CASCADE = [
    EqBlock(50.0, 100.0, 10.0, "tight thermostat, let the density find itself"),
    EqBlock(100.0, 100.0, 10.0, "loosen the thermostat"),
    EqBlock(200.0, 200.0, 10.0, "loosen both"),
    EqBlock(200.0, 800.0, 10.0, "loosen the barostat further"),
    EqBlock(500.0, 1000.0, 10.0, "nearly free dynamics before production"),
]


@dataclass
class Protocol:
    target_T: float = 1873.0             # 1600 C
    melt_T: float = 4000.0               # THEIRS
    pressure_atm: float = 1.0            # THEIRS
    timestep_fs: float = 2.0             # THEIRS
    cool_rate_K_per_ps: float = 1.0      # THEIRS
    hot_ns: float = 5.0                  # THEIRS, 2,500,000 steps
    hot_tdamp_fs: float = 100.0          # THEIRS
    hot_pdamp_fs: float = 200.0          # THEIRS
    eq_blocks: list[EqBlock] = field(
        default_factory=lambda: [EqBlock(**asdict(b)) for b in THEIR_EQ_CASCADE])
    produce_ns: float = 50.0             # THEIRS, 25,000,000 steps
    prod_tdamp_fs: float = 1000.0        # THEIRS
    prod_pdamp_fs: float = 2000.0        # THEIRS
    sample_every_ps: float = 10.0        # THEIRS, dump every 5000 steps
    coulomb_cutoff_A: float = 11.0       # THEIRS
    ewald_accuracy: float = 1.0e-6       # THEIRS
    gewald: float = 0.15                 # THEIRS, kspace_modify gewald 0.15
    pair_modify_tail: bool = True        # THEIRS, pair_modify tail yes
    delete_overlap_A: float = 0.5        # THEIRS, delete_atoms overlap 0.5 metals Ox
    coord_cutoff_A: float = 2.5          # THEIRS, coord/atom cutoff 2.5 group Si_Al
    thermo_every_ps: float = 2.0
    # neighbour list sizing. LAMMPS defaults are usually fine, but an 11 Angstrom cutoff in
    # a dense oxide melt gives each atom several hundred neighbours and the default page can
    # overflow. The NEIGH_OVERFLOW remedy raises these.
    neigh_one: int = 0                   # 0 leaves the LAMMPS default alone
    neigh_page: int = 0
    seed: int = 458165                   # THEIRS, the seed in their velocity create
    # MINE, and off by default so a run reproduces theirs exactly unless asked
    stress_every_steps: int = 0          # 5 enables Green-Kubo viscosity sampling
    nve_check_ns: float = 0.0            # a short NVE energy conservation test

    # ------------------------------------------------------------------ arithmetic
    def steps(self, ns: float) -> int:
        return int(round(ns * 1.0e6 / self.timestep_fs))

    def cool_ns(self) -> float:
        return (self.melt_T - self.target_T) / self.cool_rate_K_per_ps / 1000.0

    def equilibrate_ns(self) -> float:
        return sum(b.ns for b in self.eq_blocks)

    def total_ns(self) -> float:
        return (self.hot_ns + self.cool_ns() + self.equilibrate_ns()
                + self.produce_ns + self.nve_check_ns)

    def est_core_hours(self, n_atoms: int) -> float:
        """Rough cost, anchored on their own reported figure.

        Their whole multi-stage classical workflow took about 10 hours on 2 nodes of 128
        cores for about 5000 atoms, which is 2560 core hours for 117 ns. Scaling is linear
        in atoms and in simulated time for this kind of run.
        """
        their_core_hours = 2.0 * 128.0 * 10.0
        their_ns = 5.0 + 2.127 + 50.0 + 50.0
        their_atoms = 5000.0
        return their_core_hours * (self.total_ns() / their_ns) * (n_atoms / their_atoms)

    def warnings(self) -> list[str]:
        w: list[str] = []
        if self.timestep_fs > 2.5:
            w.append(f"a {self.timestep_fs} fs timestep is longer than the 2 fs they used. "
                     "The Si-O stretch needs about ten steps per period.")
        if self.target_T < 1200:
            w.append(f"{self.target_T:.0f} K is {self.target_T-273.15:.0f} C. Check the "
                     "liquidus. Below it this is a supercooled liquid and the diffusion "
                     "coefficients will not converge.")
        if self.produce_ns < 5:
            w.append(f"a {self.produce_ns:g} ns production run is short against their "
                     "50 ns. Silicon is the slowest species and needs tens of ns at "
                     "1600 C to leave the cage regime.")
        if self.equilibrate_ns() < 5:
            w.append(f"{self.equilibrate_ns():g} ns of equilibration against their 50 ns. "
                     "The density will still be drifting during production.")
        if len(self.eq_blocks) < 2:
            w.append("they used five equilibration blocks with the damping loosened at "
                     "each step. One block with fixed damping equilibrates more slowly "
                     "and leaves more thermostat interference in the dynamics.")
        if self.melt_T < 3000:
            w.append(f"randomising at {self.melt_T:.0f} K instead of their 4000 K may not "
                     "erase the memory of the random insertion.")
        if self.sample_every_ps > 20:
            w.append(f"sampling every {self.sample_every_ps:g} ps gives few points for the "
                     "MSD fit. They used 10 ps.")
        if self.stress_every_steps and self.stress_every_steps > 20:
            w.append("the stress autocorrelation decays in well under a picosecond. "
                     "Sampling it less often than about every 10 steps aliases it away "
                     "and the Green-Kubo viscosity will be wrong rather than noisy.")
        return w


HEADER = """# ---------------------------------------------------------------------------
# {name}   stage {stage}
# generated by slagmd on {date}
#
# Composition  {comp_line}
# Atoms        {natoms}
# Force field  {ff_label}
#              {ff_source}
#
# PROVENANCE OF EVERY SETTING BELOW
#   [THEIRS] copied from the authors' released input files, {ds}
#   [PAPER]  stated in {paper}
#   [MINE]   slagmd's own addition. Never cite a [MINE] line as theirs.
#
# PACKAGES REQUIRED IN LAMMPS BUILD: KSPACE
# ---------------------------------------------------------------------------
"""


def _preamble(proto: Protocol, comp: Composition, ff: ForceField,
              read: str, name: str, stage: str, date: str) -> list[str]:
    species = comp.species
    comp_line = "  ".join(
        f"{ox} {v:.4f}" for ox, v in sorted(comp.mole_fractions.items()))
    L = [HEADER.format(name=name, stage=stage, date=date or "unknown date",
                       comp_line=comp_line, natoms=comp.n_atoms, ff_label=ff.label,
                       ff_source=ff.source, ds=DS_FF, paper=PAPER_FF),
         "units           real                    # [THEIRS]",
         "atom_style      charge                  # [THEIRS]",
         "boundary        p p p                   # [THEIRS]",
         "",
         read,
         ""]
    # their group definitions, kept because their own analysis scripts expect them
    for i, sp in enumerate(species, start=1):
        gname = "Ox" if sp == "O" else sp
        L.append(f"group           {gname:<8s} type {i}        # [THEIRS] naming")
    cations = [s for s in species if s != "O"]
    L.append(f"group           metals   type "
             f"{' '.join(str(species.index(s)+1) for s in cations)}     # [THEIRS]")
    formers = [s for s in cations if s in ("Si", "Al", "Fe3", "Ti")]
    if formers:
        L.append(f"group           Si_Al    union {' '.join(formers)}"
                 f"        # [THEIRS] network formers, used for oxygen coordination")
    L += ["",
          "neighbor        2.0 bin                 # [THEIRS]",
          "neigh_modify    delay 0 every 1 check yes    # [THEIRS]"]
    if proto.neigh_one or proto.neigh_page:
        bits = []
        if proto.neigh_one:
            bits.append(f"one {proto.neigh_one}")
        if proto.neigh_page:
            bits.append(f"page {proto.neigh_page}")
        L.append("neigh_modify    " + " ".join(bits)
                 + "    # [MINE] raised after a neighbour list overflow")
    L += [
          "",
          f"{ff.pair_style(proto.coulomb_cutoff_A)}"
          f"     # [THEIRS] {proto.coulomb_cutoff_A:g} A real space cutoff"]
    if proto.pair_modify_tail:
        L.append("pair_modify     tail yes                # [THEIRS] long range dispersion "
                 "correction to energy and pressure")
    L.append("")
    L.append("# potential parameters converted to kcal/mol for units real, by dividing the")
    L.append("# published kJ/mol values by 4.184. Verified to agree with their own")
    L.append("# pair_coeff lines to every printed digit.")
    if ff.form == "bmh":
        L.append("pair_coeff      * * 0 1 0 0 0             # [THEIRS] zero everything first")
    else:
        L.append("pair_coeff      * * 0 1 0                 # [THEIRS] zero everything first")
    L += ff.pair_coeff_lines(species)
    L.append("")
    for i, sp in enumerate(species, start=1):
        L.append(f"set             type {i} charge {ff.charges[sp]:g}"
                 f"        # [THEIRS] {sp}")
    L += ["",
          f"kspace_style    ewald {proto.ewald_accuracy:g}              # [THEIRS]",
          f"kspace_modify   gewald {proto.gewald:g}             # [THEIRS]",
          "",
          f"timestep        {proto.timestep_fs:g}                     # [THEIRS]",
          f"variable        p equal {proto.pressure_atm:g}",
          "thermo_style    custom step time temp press density vol pe ke etotal",
          f"thermo          {max(1, proto.steps(proto.thermo_every_ps/1000.0))}",
          "thermo_modify   flush yes                # [MINE] so a dead job's log is "
          "readable",
          ""]
    return L


def _npt(tdamp: float, pdamp: float, t0: float, t1: float, tag: str) -> str:
    return (f"fix             eq all npt temp {t0:g} {t1:g} {tdamp:g} "
            f"iso ${{p}} ${{p}} {pdamp:g}    # {tag}")


def write_inputs(outdir: Path, name: str, comp: Composition, ff: ForceField,
                 proto: Protocol, data_file: str = "melt.data",
                 date: str = "") -> dict[str, Path]:
    """Write one input script per stage. Returns {stage_name: path} in run order."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ff.check_ready(comp.species)
    species = comp.species
    paths: dict[str, Path] = {}

    def emit(stage: str, read: str, body: list[str]) -> None:
        lines = _preamble(proto, comp, ff, read, name, stage, date) + body + [
            f"write_restart   restart.{stage}",
            f'print           "SLAGMD_STAGE_DONE {stage}"',
            ""]
        p = outdir / f"in.{stage}"
        p.write_text("\n".join(lines))
        paths[stage] = p

    # ------------------------------------------------------- 00 minimise
    emit("00_minimise", f"read_data       {data_file}", [
        "# Stage 0. [THEIRS] delete any residual overlap, then minimise.",
        "# Their delete_atoms overlap is the neat trick here: rather than fight the packing",
        "# algorithm, they simply remove a cation that landed on top of an oxygen. slagmd's",
        "# own packer does not create overlaps, so this normally deletes nothing, and the",
        "# atom count is checked afterwards so a silent deletion cannot pass unnoticed.",
        f"variable        n_before equal count(all)",
        f"delete_atoms    overlap {proto.delete_overlap_A:g} metals Ox    # [THEIRS]",
        "variable        n_after equal count(all)",
        'print           "SLAGMD_ATOMS_BEFORE ${n_before} AFTER ${n_after}"',
        "",
        "fix             initial all nve           # [THEIRS]",
        "run             0                         # [THEIRS]",
        "unfix           initial",
        f"velocity        all create {proto.melt_T:g} {proto.seed} rot yes dist gaussian"
        f"   # [THEIRS]",
        "minimize        1.0e-6 1.0e-8 10000 100000     # [THEIRS]",
        "reset_timestep  0",
        ""])

    # ------------------------------------------------------- 01 hot randomise
    emit("01_hot", "read_restart    restart.00_minimise", [
        f"# Stage 1. [THEIRS] {proto.hot_ns:g} ns at {proto.melt_T:g} K in NPT, "
        f"{proto.steps(proto.hot_ns)} steps.",
        "# The job of this stage is to forget the starting configuration completely.",
        _npt(proto.hot_tdamp_fs, proto.hot_pdamp_fs, proto.melt_T, proto.melt_T,
             "[THEIRS] Tdamp 100 fs, Pdamp 200 fs"),
        f"dump            DUMP all custom {proto.steps(10.0/1000.0)} hot.trj "
        f"id type xu yu zu",
        "dump_modify     DUMP sort id element " + " ".join(species),
        f"run             {proto.steps(proto.hot_ns)}",
        "unfix           eq",
        "undump          DUMP",
        "reset_timestep  0",
        ""])

    # ------------------------------------------------------- 02 cool
    emit("02_cool", "read_restart    restart.01_hot", [
        f"# Stage 2. [THEIRS] linear ramp {proto.melt_T:g} K to {proto.target_T:g} K at "
        f"{proto.cool_rate_K_per_ps:g} K/ps,",
        f"# which is {proto.cool_ns():.3f} ns and {proto.steps(proto.cool_ns())} steps.",
        _npt(proto.hot_tdamp_fs, proto.hot_pdamp_fs, proto.melt_T, proto.target_T,
             "[THEIRS]"),
        f"run             {proto.steps(proto.cool_ns())}",
        "unfix           eq",
        "reset_timestep  0",
        ""])

    # ------------------------------------------------------- 03..07 equilibrate
    prev = "02_cool"
    for i, blk in enumerate(proto.eq_blocks, start=1):
        stage = f"{2+i:02d}_equil{i}"
        emit(stage, f"read_restart    restart.{prev}", [
            f"# Equilibration block {i} of {len(proto.eq_blocks)}. [THEIRS] {blk.ns:g} ns, "
            f"Tdamp {blk.tdamp_fs:g} fs, Pdamp {blk.pdamp_fs:g} fs.",
            f"# {blk.label}",
            "# The cascade exists because a thermostat tight enough to pull the density to",
            "# equilibrium quickly is also tight enough to distort the dynamics you want",
            "# to measure. So they start tight and loosen it in five steps, ending nearly",
            "# free. This schedule appears in no paper, only in the released scripts.",
            _npt(blk.tdamp_fs, blk.pdamp_fs, proto.target_T, proto.target_T, "[THEIRS]"),
            "variable        dens equal density",
            "variable        vol equal vol",
            f"fix             avg all ave/time 100 "
            f"{max(1, proto.steps(blk.ns/50.0)//100)} {max(1, proto.steps(blk.ns/50.0))} "
            f"v_dens v_vol file equil{i}_running.txt      # [MINE] so the analysis can "
            f"PROVE the density stopped drifting",
            f"run             {proto.steps(blk.ns)}",
            "unfix           avg",
            "unfix           eq",
            "reset_timestep  0",
            ""])
        prev = stage

    # ------------------------------------------------------- 08 production
    dump_every = max(1, proto.steps(proto.sample_every_ps / 1000.0))
    n_sp = len(species)
    # their RDF list: like-pairs and every cation-oxygen pair, no cation-cation cross terms
    o_type = species.index("O") + 1
    rdf_pairs: list[str] = []
    for i in range(1, n_sp + 1):
        rdf_pairs.append(f"{i} {i}")
        if i != o_type:
            rdf_pairs.append(f"{i} {o_type}")
    body = [
        f"# Production. [THEIRS] {proto.produce_ns:g} ns, "
        f"{proto.steps(proto.produce_ns)} steps, Tdamp {proto.prod_tdamp_fs:g} fs, "
        f"Pdamp {proto.prod_pdamp_fs:g} fs.",
        "#",
        "# NOTE ON A DEFECT IN THEIR RELEASED production.file: it ends with unfix calls for",
        "# an RDF fix and four coordination fixes that the file never defines, so LAMMPS",
        "# aborts with 'Could not find unfix ID'. The physics is already written by then so",
        "# their published results are unaffected, but the job exits non-zero, which breaks",
        "# any dependency chain. Every fix unfixed below is defined above it.",
        "",
        _npt(proto.prod_tdamp_fs, proto.prod_pdamp_fs, proto.target_T, proto.target_T,
             "[THEIRS]"),
        "",
        "# [THEIRS] both wrapped and unwrapped coordinates are dumped. The unwrapped ones",
        "# are what the MSD must use: wrapped coordinates give a wrong MSD the moment an",
        "# atom crosses the boundary, and in NPT the box length is changing too, so",
        "# unwrapping afterwards from image flags is awkward and error prone.",
        f"dump            DUMP all custom {dump_every} prod.lammpstrj "
        f"id type element xu yu zu x y z",
        f"dump_modify     DUMP element {' '.join(species)} sort id",
        "",
        f"compute         RDF all rdf 200 {' '.join(rdf_pairs)}     # [THEIRS] 200 bins",
        f"fix             rdfout all ave/time 100 "
        f"{max(1, proto.steps(proto.produce_ns)//100//100)} "
        f"{proto.steps(proto.produce_ns)} c_RDF[*] file rdf_lammps.txt mode vector",
        "",
        f"compute         CordOx Ox coord/atom cutoff {proto.coord_cutoff_A:g} group Si_Al"
        f"   # [THEIRS] fixed 2.5 A",
        f"fix             CoordHist Ox ave/histo 5000 100 500000 0 16 17 c_CordOx "
        f"file Cord_N_Ox.txt mode vector      # [THEIRS]",
        "",
        "# [MINE] LAMMPS-side MSD per species, as an independent check on the python",
        "# analysis. If the two disagree the python is wrong, and that is worth knowing.",
    ]
    for i, s in enumerate(species, start=1):
        gname = "Ox" if s == "O" else s
        body.append(f"compute         msd{s} {gname} msd com yes")
    body += [
        f"fix             msdout all ave/time {dump_every} 1 {dump_every} "
        + " ".join(f"c_msd{s}[4]" for s in species)
        + " file msd_lammps.txt title3 \"# time_fs " + " ".join(species) + "\"",
        ""]
    if proto.stress_every_steps > 0:
        body += [
            "# [MINE] AND THIS IS THE GAP. Their 84 released files contain no compute of",
            "# the pressure tensor, no stress autocorrelation, no Green-Kubo, no",
            "# ave/correlate. Viscosity is never computed from MD anywhere in their work;",
            "# it enters only as a FactSage screening number. These four lines are the",
            "# whole cost of adding it to a production run that is happening anyway.",
            "variable        pxy equal pxy",
            "variable        pxz equal pxz",
            "variable        pyz equal pyz",
            f"fix             stressout all ave/time {proto.stress_every_steps} 1 "
            f"{proto.stress_every_steps} v_pxy v_pxz v_pyz file stress.txt "
            f"title3 \"# step pxy pxz pyz (atm)\"",
            ""]
    body += [
        f"run             {proto.steps(proto.produce_ns)}",
        "unfix           msdout",
        "unfix           CoordHist",
        "unfix           rdfout"]
    if proto.stress_every_steps > 0:
        body.append("unfix           stressout")
    body += ["unfix           eq", "undump          DUMP", ""]
    prod_stage = f"{3+len(proto.eq_blocks):02d}_produce"
    emit(prod_stage, f"read_restart    restart.{prev}", body)

    # ------------------------------------------------- optional NVE energy check
    if proto.nve_check_ns > 0:
        stage = f"{4+len(proto.eq_blocks):02d}_nve_check"
        emit(stage, f"read_restart    restart.{prod_stage}", [
            "# [MINE] A short NVE run. If total energy drifts by more than about 1e-4 per",
            "# atom per ns, the timestep is too long or the Ewald accuracy too loose, and",
            "# every transport number above is suspect. It is the cheapest check in",
            "# molecular dynamics and almost nobody does it.",
            "fix             nve all nve",
            f"run             {proto.steps(proto.nve_check_ns)}",
            "unfix           nve",
            ""])

    # ------------------------------------------------------------------ manifest
    import json
    (outdir / "protocol.json").write_text(json.dumps({
        "name": name,
        "composition_requested": comp.mole_fractions,
        "composition_achieved": comp.achieved_mole_fractions,
        "atom_counts": comp.counts,
        "species_to_lammps_type": {s: i + 1 for i, s in enumerate(species)},
        "n_atoms": comp.n_atoms,
        "force_field": ff.name,
        "force_field_label": ff.label,
        "force_field_source": ff.source,
        "charges": ff.charges,
        "protocol": {**{k: v for k, v in asdict(proto).items() if k != "eq_blocks"},
                     "eq_blocks": [asdict(b) for b in proto.eq_blocks]},
        "total_ns": proto.total_ns(),
        "estimated_core_hours": round(proto.est_core_hours(comp.n_atoms), 1),
        "stages": list(paths),
        "protocol_warnings": proto.warnings(),
        "provenance": {"force_field_scripts": DS_FF, "high_throughput_scripts": DS_HT,
                       "papers": [PAPER_FF, PAPER_HT]},
    }, indent=2))
    return paths


def stage_order(proto: Protocol) -> list[str]:
    """The stage names this protocol will generate, in run order."""
    order = ["00_minimise", "01_hot", "02_cool"]
    order += [f"{2+i:02d}_equil{i}" for i in range(1, len(proto.eq_blocks) + 1)]
    order.append(f"{3+len(proto.eq_blocks):02d}_produce")
    if proto.nve_check_ns > 0:
        order.append(f"{4+len(proto.eq_blocks):02d}_nve_check")
    return order
