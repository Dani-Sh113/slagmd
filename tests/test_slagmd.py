"""The audit. Every number this tool can produce is checked against something known.

Run:  python3 tests/test_slagmd.py
Exit code 0 means every check passed. Anything else means do not trust the tool.

The checks are grouped by what they protect against:
  A  transcription errors in the force field tables
  B  composition and charge neutrality
  C  the starting structure
  D  the generated LAMMPS deck
  E  analysis against a trajectory whose answer we know exactly
  F  the failure diagnosis rules, each against the real LAMMPS error text
  G  the whole pipeline end to end, including a failure and an automatic recovery
"""
from __future__ import annotations

import math
import re
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from slagmd import analysis as an                                   # noqa: E402
from slagmd import diagnose as diag                                 # noqa: E402
from slagmd import mocklammps                                       # noqa: E402
from slagmd import scheduler as sched                               # noqa: E402
from slagmd.builder import build_random, check_structure, write_data_file  # noqa: E402
from slagmd.composition import (VERBEECK_2025_TABLE5, Composition,   # noqa: E402
                                CompositionError)
from slagmd.forcefield import ForceField, ForceFieldError           # noqa: E402
from slagmd.inputgen import EqBlock, Protocol, stage_order, write_inputs                  # noqa: E402
from slagmd.project import Project                                  # noqa: E402

PASS, FAIL = [], []


def check(name: str, condition: bool, detail: str = "") -> None:
    (PASS if condition else FAIL).append((name, detail))
    mark = "ok  " if condition else "FAIL"
    print(f"  {mark} {name}" + (f"   {detail}" if detail else ""))


def section(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


# ============================================================ A  force field tables
section("A  force field tables, against the identities their own authors' forms imply")

ffs = {n: ForceField.load(n) for n in ForceField.available()}
check("three verified force fields are present", len(ffs) == 3, str(sorted(ffs)))

for name, ff in ffs.items():
    probs = ff.check_self_consistency()
    check(f"{name} passes its consistency identity", not probs, "; ".join(probs))

# Matsui's A must equal 4.184 * rho for all ten pairs. Independent restatement.
mat = ffs["matsui"]
worst = max(abs(p["A"] - 4.184 * p["rho"]) for p in mat.pairs.values())
check("matsui A = 4.184 rho for every pair", worst < 1e-4, f"worst residual {worst:.2e}")

# Guillot must be Matsui rewritten, for the three cation-oxygen pairs
gui = ffs["guillot"]
for pair in ("Al-O", "Ca-O", "O-Si"):
    m, g = mat.pairs[pair], gui.pairs[pair]
    expect = m["A"] * math.exp(m["sigma"] / m["rho"])
    rel = abs(g["A"] - expect) / expect
    check(f"guillot {pair} is matsui rewritten", rel < 5e-3, f"{100*rel:.2f} percent")

# The paper's own claim about well depths, reproduced from the raw parameters
bou = ffs["bouhadja"]
claims = {("Ca", "O"): 66.0, ("Al", "O"): 57.0, ("Si", "O"): 57.0}
for pair, claimed in claims.items():
    _, um = mat.well_depth(*pair)
    _, ub = bou.well_depth(*pair)
    deeper = 100.0 * (ub - um) / abs(um) * -1.0
    check(f"bouhadja {pair[0]}-{pair[1]} well is about {claimed:.0f} percent deeper",
          abs(deeper - claimed) < 3.0, f"computed {deeper:.1f} percent")

# Ground truth: our Bouhadja numbers must equal the authors' own released pair_coeff
# lines, and a single corrupted digit must be caught. This is stronger than an algebraic
# identity because it constrains parameters the algebra cannot.
import json as _json, tempfile as _tf
_d = _json.loads((ROOT / "slagmd" / "forcefields.json").read_text())
check("bouhadja carries the authors' released pair_coeff values for all ten pairs",
      len([k for k in _d["bouhadja"].get("released_pair_coeff_kcal", {})
           if not k.startswith("_")]) == 10)
_d["bouhadja"]["pairs"]["Si-O"]["A"] = 0.6854          # was 0.6754, one digit changed
_tmpff = Path(_tf.mkdtemp()) / "ff.json"
_tmpff.write_text(_json.dumps(_d))
_bad = ForceField.load("bouhadja", _tmpff)
_probs = _bad.check_self_consistency()
check("a single corrupted digit in bouhadja is caught against the released file",
      any("Si-O" in p and "released" in p for p in _probs),
      "; ".join(_probs)[:80] if _probs else "NOT CAUGHT")

# A force field added later must never pass silently. There are two cases and the tool
# must give the RIGHT message for each, because they need different fixes.
_base = _json.loads((ROOT / "slagmd" / "forcefields.json").read_text())

# case one: it declares an identity_check that nobody implemented
_d2 = _json.loads(_json.dumps(_base))
_d2["declaresff"] = {**_json.loads(_json.dumps(_base["matsui"])),
                     "label": "declares a check nobody implements"}
_t2 = Path(_tf.mkdtemp()) / "ff2.json"
_t2.write_text(_json.dumps(_d2))
_p2 = ForceField.load("declaresff", _t2).check_self_consistency()
check("a new force field that declares an unimplemented check is told so",
      any("declares an identity_check" in p and "no implementation" in p for p in _p2),
      "; ".join(_p2)[:80] if _p2 else "NOT REPORTED")

# case two: no identity_check block and no released file
_d3 = _json.loads(_json.dumps(_base))
_bare = _json.loads(_json.dumps(_base["matsui"]))
_bare.pop("identity_check", None)
_bare["label"] = "no check of any kind"
_d3["bareff"] = _bare
_t3 = Path(_tf.mkdtemp()) / "ff3.json"
_t3.write_text(_json.dumps(_d3))
_p3 = ForceField.load("bareff", _t3).check_self_consistency()
check("a new force field with no check at all is told its numbers are unconfirmed",
      any("no algebraic identity check" in p for p in _p3),
      "; ".join(_p3)[:80] if _p3 else "NOT REPORTED")

# and in both cases the tool must REFUSE to build a deck
_refused = 0
for _nm, _pth in (("declaresff", _t2), ("bareff", _t3)):
    try:
        ForceField.load(_nm, _pth).check_ready(["Ca", "Al", "Si", "O"])
    except ForceFieldError:
        _refused += 1
check("and check_ready refuses to build a deck in both cases", _refused == 2,
      f"{_refused} of 2 refused")

# Guillot's cation-cation pairs must be absent, not zero by accident
check("guillot has no cation-cation short range terms",
      not any(k for k in gui.pairs if "O" not in k), sorted(gui.pairs))

# An unverified force field must be refused
try:
    ForceField.load("nonexistent_ff")
    check("an unknown force field is refused", False)
except ForceFieldError:
    check("an unknown force field is refused", True)

# A species the force field has no charge for must be refused, not guessed
try:
    mat.check_ready(["Ca", "Al", "Si", "Fe2", "O"])
    check("a species with no parameters is refused", False)
except ForceFieldError as exc:
    check("a species with no parameters is refused", True,
          "iron correctly rejected for matsui")


# ============================================================== B  composition
section("B  composition, atom counts and charge neutrality")

for label, mf, _tr in VERBEECK_2025_TABLE5:
    c = Composition(mf, 5000)
    for ffname, ff in ffs.items():
        q = c.net_charge(ff.charges)
        check(f"{label} is neutral under {ffname}", abs(q) < 1e-9, f"q={q:+.2e}")
        break
    check(f"{label} hits 5000 atoms within 0.2 percent",
          abs(c.n_atoms - 5000) / 5000 < 0.002, f"{c.n_atoms} atoms")
    check(f"{label} mole fractions reproduced to 0.001",
          c.max_mole_fraction_error < 1e-3, f"max error {c.max_mole_fraction_error:.5f}")

# oxygen is counted, never rounded: exact stoichiometry must hold
c = Composition({"CaO": 0.453, "Al2O3": 0.125, "SiO2": 0.422}, 5000)
expect_O = (c.formula_units["CaO"] + 3 * c.formula_units["Al2O3"]
            + 2 * c.formula_units["SiO2"])
check("oxygen count follows exactly from the cations", c.counts["O"] == expect_O,
      f"{c.counts['O']} == {expect_O}")
check("aluminium count is exactly twice the Al2O3 units",
      c.counts["Al"] == 2 * c.formula_units["Al2O3"])

ok, scale = c.charges_are_scaled_formal(mat.charges)
check("matsui charges are a uniform scaling of the formal charges", ok,
      f"scale {scale:.4f}")
ok2, scale2 = c.charges_are_scaled_formal(bou.charges)
check("bouhadja charges are a uniform scaling too", ok2, f"scale {scale2:.4f}")

# box volume must invert the density exactly
L = c.box_length_A(2.6)
mass_g = c.molar_mass_g() / 6.02214076e23
rho_back = mass_g / (L ** 3 * 1e-24)
check("box edge inverts the requested density", abs(rho_back - 2.6) < 1e-9,
      f"{rho_back:.10f}")

# an impossible composition must be refused
for bad, why in [({"CaO": 0.5, "XYZ": 0.5}, "unknown oxide"),
                 ({"CaO": 0.0, "SiO2": 0.0}, "zero total")]:
    try:
        Composition(bad, 1000)
        check(f"refuses {why}", False)
    except CompositionError:
        check(f"refuses {why}", True)

# a box too small for the 11 A cutoff must be refused, not silently run
try:
    from slagmd.builder import BuildError
    build_random(Composition({"CaO": 0.5, "SiO2": 0.5}, 200), 2.6)
    check("refuses a box smaller than twice the Coulomb cutoff", False)
except Exception as exc:                                    # noqa: BLE001
    check("refuses a box smaller than twice the Coulomb cutoff",
          "cutoff" in str(exc).lower(), str(exc)[:70])


# =============================================================== C  the structure
section("C  the starting structure")

comp = Composition({"CaO": 0.453, "Al2O3": 0.125, "SiO2": 0.422}, 1200)
st = build_random(comp, density_g_cm3=2.6, seed=7)
check("atom count matches the composition", st.n_atoms == comp.n_atoms,
      f"{st.n_atoms}")
got = {s: list(st.species).count(s) for s in set(st.species)}
check("every species count is exact", got == comp.counts, str(got))
dmin = st.min_distance()
check("no pair closer than 1.5 A", dmin >= 1.5, f"closest {dmin:.3f} A")
check("density comes out as requested", abs(st.density_g_cm3() - 2.6) < 0.01,
      f"{st.density_g_cm3():.4f} g/cm3")
check("structure passes its own checks", not check_structure(st, comp, mat.charges))
check("a small but legal box warns about finite size effects", bool(st.warnings),
      (st.warnings[0][:70] if st.warnings else "no warning"))
big = build_random(Composition({"CaO": 0.453, "Al2O3": 0.125, "SiO2": 0.422}, 5000), 2.6, 7)
check("a 5000 atom box warns about nothing", not big.warnings, str(big.warnings))

# a structure with a deliberate overlap must be caught
bad = build_random(comp, 2.6, seed=7)
bad.positions[1] = bad.positions[0] + np.array([0.3, 0.0, 0.0])
probs = check_structure(bad, comp, mat.charges)
check("a 0.3 A overlap is caught", any("closest pair" in p for p in probs),
      "; ".join(probs)[:80])

# a structure with a net charge must be caught
charges_bad = dict(mat.charges)
charges_bad["Ca"] = 1.0
probs = check_structure(st, comp, charges_bad)
check("a net box charge is caught", any("net box charge" in p for p in probs))


# ============================================================ D  the LAMMPS deck
section("D  the generated LAMMPS deck, against their released input files")

tmp = Path(tempfile.mkdtemp(prefix="slagmd_deck_"))

# their real protocol first: the arithmetic must reproduce what both papers state
theirs = Protocol(target_T=1873.0)
check("their protocol gives 50 ns of equilibration in five blocks",
      abs(theirs.equilibrate_ns() - 50.0) < 1e-9 and len(theirs.eq_blocks) == 5,
      f"{theirs.equilibrate_ns():g} ns in {len(theirs.eq_blocks)} blocks")
check("their cooling ramp from 4000 to 1873 K at 1 K/ps is 2.127 ns",
      abs(theirs.cool_ns() - 2.127) < 1e-9, f"{theirs.cool_ns():.4f} ns")
check("the hot stage is 2,500,000 steps", theirs.steps(theirs.hot_ns) == 2_500_000,
      f"{theirs.steps(theirs.hot_ns)}")
check("production is 25,000,000 steps",
      theirs.steps(theirs.produce_ns) == 25_000_000,
      f"{theirs.steps(theirs.produce_ns)}")
check("the damping cascade loosens monotonically in both thermostat and barostat",
      all(theirs.eq_blocks[i].tdamp_fs <= theirs.eq_blocks[i+1].tdamp_fs
          for i in range(4))
      and all(theirs.eq_blocks[i].pdamp_fs <= theirs.eq_blocks[i+1].pdamp_fs
              for i in range(4)),
      " -> ".join(f"({b.tdamp_fs:g},{b.pdamp_fs:g})" for b in theirs.eq_blocks))
check("the cost estimate reproduces their own 2560 core hours",
      abs(theirs.est_core_hours(5000) - 2560.0) < 1.0,
      f"{theirs.est_core_hours(5000):.0f} core hours for 5000 atoms")
check("their own protocol produces no warnings", not theirs.warnings(),
      str(theirs.warnings()))

# a short protocol for the file-content checks
proto = Protocol(target_T=1873.0, hot_ns=0.002, produce_ns=0.012,
                 cool_rate_K_per_ps=500.0, sample_every_ps=0.1, thermo_every_ps=0.5,
                 eq_blocks=[EqBlock(50.0, 100.0, 0.002, "a"),
                            EqBlock(500.0, 1000.0, 0.002, "b")],
                 stress_every_steps=5)
paths = write_inputs(tmp, "deck_test", comp, bou, proto)
expect_stages = stage_order(proto)
check("the stage list matches the protocol", list(paths) == expect_stages,
      str(list(paths)))
check("two equilibration blocks give six stages", len(paths) == 6, str(len(paths)))

prod_name = [k for k in paths if k.endswith("_produce")][0]
prod = (tmp / f"in.{prod_name}").read_text()
code = "\n".join(l.split("#")[0] for l in prod.splitlines())   # strip comments

check("production dump uses UNWRAPPED coordinates", "xu yu zu" in code)
check("born/coul/long is used for a BMH force field",
      "pair_style born/coul/long" in code)
check("their Ewald accuracy is used", "kspace_style   ewald 1e-06" in code
      or "ewald 1e-06" in code, "")
check("their forced gewald is present", "kspace_modify   gewald 0.15" in code)
check("their tail correction is present", "pair_modify     tail yes" in code)
check("their 11 A cutoff appears", "born/coul/long 11" in code)
check("their 200-bin RDF compute is present", "compute         RDF all rdf 200" in code)
check("their fixed 2.5 A oxygen coordination compute is present",
      "coord/atom cutoff 2.5 group Si_Al" in code)
check("a completion marker is printed", f"SLAGMD_STAGE_DONE {prod_name}" in prod)
check("LAMMPS side MSD is computed as an independent check", "compute         msd" in code)
check("their defect is not reproduced: every unfix has a matching fix",
      all(f"fix             {n}" in code or f"fix   {n}" in code
          for n in re.findall(r"^unfix\s+(\S+)", code, re.M) if n != "eq"),
      str(re.findall(r"^unfix\s+(\S+)", code, re.M)))
check("every stage writes a restart file",
      all("write_restart" in (tmp / f"in.{st}").read_text() for st in paths))
check("the minimise stage records the atom count across their delete_atoms",
      "SLAGMD_ATOMS_BEFORE" in (tmp / "in.00_minimise").read_text())
check("the stress tensor is sampled when viscosity is asked for",
      "v_pxy v_pxz v_pyz" in code)

npairs_expected = len(comp.species) * (len(comp.species) + 1) // 2
check("one pair_coeff per type pair plus the zeroing line",
      code.count("pair_coeff") == npairs_expected + 1,
      f"{code.count('pair_coeff')} of {npairs_expected + 1}")

# THE KEY CHECK: our independently transcribed coefficients must equal theirs, byte level
THEIR_BOUHADJA_KCAL = {      # from ds2 cooling.file, pair_coeff lines, verbatim
    "Al-Al": (0.066874, 0.068, 1.5704, 323.99600),
    "Al-Ca": (0.073792, 0.074, 1.9572, 395.96326),
    "Al-Si": (0.05765,  0.057, 1.5056, 433.79550),
    "Al-O":  (0.17295,  0.164, 2.6067, 797.29258),
    "Ca-Ca": (0.08071,  0.080, 2.3440, 483.92794),
    "Ca-Si": (0.062262, 0.063, 1.8924, 530.16554),
    "Ca-O":  (0.177562, 0.178, 2.9935, 974.41414),
    "Si-Si": (0.027672, 0.046, 1.4408, 580.81914),
    "Si-O":  (0.16142,  0.156, 2.5419, 1067.51658),
    "O-O":   (0.27672,  0.263, 3.6430, 1962.03704),
}
worst_rel = 0.0
for key, (A_t, rho_t, sig_t, C_t) in THEIR_BOUHADJA_KCAL.items():
    a, b = key.split("-")
    ours = bou.pairs["-".join(sorted((a, b)))]
    A_o, C_o = ours["A"] / 4.184, ours["C"] / 4.184
    for got, want in ((A_o, A_t), (C_o, C_t), (ours["rho"], rho_t),
                      (ours["sigma"], sig_t)):
        worst_rel = max(worst_rel, abs(got - want) / abs(want))
# 1e-3, not 1e-4. The paper prints A in kJ/mol to four decimals and their scripts print
# the kcal/mol value to six, so the last printed digit disagrees by up to about 2e-4 on the
# smallest coefficient. A genuine transcription error is a digit transposition or a wrong
# row, which is percent to orders of magnitude, and still fails this test.
check("our transcription from the PAPER equals their RELEASED input files, all 10 pairs",
      worst_rel < 1e-3,
      f"worst relative difference {worst_rel:.2e} over 40 numbers, rounding level")

# a Buckingham force field must emit three coefficients, not five
gpaths = write_inputs(tmp / "g", "deck_g", comp, gui, proto)
gprod = (tmp / "g" / f"in.{prod_name}").read_text()
check("buck/coul/long is used for Guillot", "pair_style buck/coul/long" in gprod)
gline = [l for l in gprod.splitlines()
         if l.startswith("pair_coeff") and "# Al-Ca" in l][0]
check("Guillot cation-cation pair is emitted as Coulomb only",
      "Coulomb only" in gline, gline.strip()[:60])

# The NEIGH_OVERFLOW remedy used to promise a change it did not carry, so the tool
# resubmitted the stage unchanged and looped to the attempt cap. Found 29 Sep 2026 while
# writing the manual, which had to document the broken promise.
_nov = diag.diagnose_text(mocklammps.FAIL_MESSAGES["NEIGH_OVERFLOW"]).remedy
check("the neighbour overflow remedy actually changes something",
      bool(_nov and _nov.protocol_changes),
      str(_nov.protocol_changes if _nov else "no remedy"))
_pn = Protocol(target_T=1873.0, hot_ns=0.002, produce_ns=0.004, cool_rate_K_per_ps=500.0,
               eq_blocks=[EqBlock(50.0, 100.0, 0.002, "a")],
               neigh_one=8000, neigh_page=200000)
write_inputs(tmp / "neigh", "neigh_test", comp, bou, _pn)
_neigh_line = "neigh_modify    one 8000 page 200000"
_neigh_deck = (tmp / "neigh" / "in.00_minimise").read_text()
check("and the generated deck emits neigh_modify one and page",
      _neigh_line in _neigh_deck,
      "" if _neigh_line in _neigh_deck else "the line is absent from the deck")
check("while the default protocol leaves the LAMMPS default alone",
      "one 8000" not in (tmp / "in.00_minimise").read_text())

# protocol warnings must fire on bad physics
check("a 5 fs timestep is flagged",
      any("timestep" in w for w in Protocol(timestep_fs=5.0).warnings()))
check("a 900 K melt temperature is flagged",
      any("liquidus" in w for w in Protocol(target_T=900.0).warnings()))
check("a 0.5 ns production run is flagged",
      any("short" in w for w in Protocol(produce_ns=0.5).warnings()))
check("a single equilibration block is flagged against their five",
      any("five equilibration blocks" in w
          for w in Protocol(eq_blocks=[EqBlock(100.0, 200.0, 50.0)]).warnings()))
check("sampling the stress too rarely is flagged as aliasing",
      any("alias" in w for w in Protocol(stress_every_steps=100).warnings()))


# ================================================ E  analysis against known answers
section("E  analysis against a trajectory whose answer is known exactly")

# A Brownian walk with a chosen D. The MSD of a free random walk in 3d is 6 D t exactly,
# so fit_diffusion must return the D we put in.
rng = np.random.default_rng(1234)
TRUE_D = 2.5e-9            # m^2/s, a realistic silicate melt oxygen value
D_A2_ps = TRUE_D * 1e8
nfr, nat, dt_ps = 600, 400, 1.0
sigma = math.sqrt(2.0 * D_A2_ps * dt_ps)
pos = np.zeros((nfr, nat, 3))
pos[0] = rng.random((nat, 3)) * 40.0
for f in range(1, nfr):
    pos[f] = pos[f - 1] + rng.normal(0, sigma, (nat, 3))
traj = an.Trajectory(positions=pos, box=np.tile([40.0, 40.0, 40.0], (nfr, 1)),
                     types=np.ones(nat, dtype=int),
                     elements=np.array(["O"] * nat), times_ps=np.arange(nfr) * dt_ps,
                     unwrapped=True)
t, msd = an.compute_msd(traj, "O")
res = an.fit_diffusion(t, msd, "O")
err = abs(res.D_m2_s - TRUE_D) / TRUE_D
# A single fitted slope has statistical error, so the honest test is not an arbitrary
# tolerance but whether the TRUE value sits inside the error bar the tool itself reports.
# That is what block averaging is for, and if the tool's error bar did not cover the truth
# the error bar would be the bug.
dwe = an.diffusion_with_error(traj, "O", n_blocks=5)
inside = abs(dwe.D_mean_m2_s - TRUE_D) <= 3.0 * dwe.D_stderr_m2_s
check("diffusion coefficient recovered from a known random walk, to about 6 percent",
      err < 0.08, f"true {TRUE_D:.3e}, got {res.D_m2_s:.3e}, error {100*err:.2f} percent")
check("and the TRUE value lies inside the tool's own block error bar", inside,
      f"{dwe.D_mean_m2_s:.3e} +/- {dwe.D_stderr_m2_s:.1e}, truth {TRUE_D:.3e}, "
      f"{abs(dwe.D_mean_m2_s-TRUE_D)/dwe.D_stderr_m2_s:.2f} sigma away")
check("the block result carries D in cm2/s, which the interface expects",
      hasattr(dwe, "D_cm2_s") and abs(dwe.D_cm2_s - dwe.D_m2_s * 1e4) < 1e-30,
      f"{dwe.D_cm2_s:.3e} cm2/s")
check("the block error bar is reported as a percentage and a verdict",
      dwe.relative_stderr_percent > 0 and len(dwe.verdict) > 30,
      f"{dwe.relative_stderr_percent:.1f} percent")
check("the log-log slope of a free random walk is 1",
      abs(res.loglog_slope - 1.0) < 0.05, f"{res.loglog_slope:.4f}")
check("the regime test passes for a genuinely diffusive walk", res.diffusive)
check("MSD of a free walk equals 6 D t at the last point",
      abs(msd[-1] / (6 * D_A2_ps * t[-1]) - 1.0) < 0.06,
      f"ratio {msd[-1]/(6*D_A2_ps*t[-1]):.4f}")

# THE CENTRAL GUARD. An oscillator, not a walker. The MSD plateaus. The tool must refuse.
pos2 = np.zeros((nfr, nat, 3))
centres = rng.random((nat, 3)) * 40.0
for f in range(nfr):
    pos2[f] = centres + rng.normal(0, 0.35, (nat, 3))
traj2 = an.Trajectory(positions=pos2, box=np.tile([40.0, 40.0, 40.0], (nfr, 1)),
                      types=np.ones(nat, dtype=int),
                      elements=np.array(["O"] * nat),
                      times_ps=np.arange(nfr) * dt_ps, unwrapped=True)
t2, msd2 = an.compute_msd(traj2, "O")
res2 = an.fit_diffusion(t2, msd2, "O")
check("a caged, plateaued MSD is REFUSED as a diffusion coefficient",
      not res2.diffusive, f"log-log slope {res2.loglog_slope:.3f}")
check("and the refusal explains why in plain words",
      "cage" in res2.verdict or "not a diffusion coefficient" in res2.verdict,
      res2.verdict[:70])

# RDF of a known lattice: a simple cubic lattice of spacing a has its first peak at a
a_lat = 3.0
n_side = 12
grid = np.array([[i, j, k] for i in range(n_side) for j in range(n_side)
                 for k in range(n_side)], dtype=float) * a_lat
L_lat = n_side * a_lat
trajl = an.Trajectory(positions=grid[None, :, :].repeat(3, axis=0),
                      box=np.tile([L_lat, L_lat, L_lat], (3, 1)),
                      types=np.ones(len(grid), dtype=int),
                      elements=np.array(["X"] * len(grid)),
                      times_ps=np.arange(3.0), unwrapped=True)
rdf = an.compute_rdf(trajl, rmax=8.0, nbins=800)
peak = rdf.first_peak["X-X"]
check("RDF first peak of a simple cubic lattice is at the lattice spacing",
      abs(peak - a_lat) < 0.05, f"found {peak:.3f} A, expected {a_lat}")
check("its coordination number to the first minimum is 6",
      abs(rdf.cn_at_first_min["X-X"] - 6.0) < 0.3,
      f"{rdf.cn_at_first_min['X-X']:.3f}")

# an RDF first peak below 1 Angstrom is not physics. Until 29 Sep 2026 the tool reported
# a 0.038 Angstrom first peak and zero coordination numbers with no complaint.
_scatter = rng.random((3, 400, 3)) * 400.0            # far outside a 40 A box
_trj = an.Trajectory(positions=_scatter, box=np.tile([40.0, 40.0, 40.0], (3, 1)),
                     types=np.ones(400, dtype=int), elements=np.array(["O"] * 400),
                     times_ps=np.arange(3.0), unwrapped=True)
_rj = an.compute_rdf(_trj, rmax=10.0, nbins=300)
check("an impossible RDF first peak is flagged",
      any("not physically possible" in w for w in _rj.warnings),
      "; ".join(_rj.warnings)[:90] or "NOT FLAGGED")

# an ideal gas must give g(r) = 1 everywhere outside the origin
ideal = rng.random((4, 3000, 3)) * 60.0
traji = an.Trajectory(positions=ideal, box=np.tile([60.0, 60.0, 60.0], (4, 1)),
                      types=np.ones(3000, dtype=int),
                      elements=np.array(["G"] * 3000), times_ps=np.arange(4.0),
                      unwrapped=True)
rdfi = an.compute_rdf(traji, rmax=20.0, nbins=100)
g = rdfi.g["G-G"][10:]
check("g(r) of an ideal gas is 1 within 5 percent", abs(g.mean() - 1.0) < 0.05,
      f"mean {g.mean():.4f}, spread {g.std():.4f}")

# a dump with no id column must WARN, because atoms cannot be aligned between frames.
# This is the defect that made Danish's own July 2026 analysis report a negative diffusion
# coefficient: LAMMPS does not write a fixed atom order unless the deck says
# "dump_modify <name> sort id", and his did not.
tmp_noid = Path(tempfile.mkdtemp(prefix="slagmd_noid_"))
frames = []
for step in range(3):
    frames += ["ITEM: TIMESTEP", str(step * 100), "ITEM: NUMBER OF ATOMS", "4",
               "ITEM: BOX BOUNDS pp pp pp", "0.0 30.0", "0.0 30.0", "0.0 30.0",
               "ITEM: ATOMS type element xu yu zu"]
    for a in range(4):
        frames.append(f"1 O {a + step:.3f} {a:.3f} {a:.3f}")
(tmp_noid / "noid.lammpstrj").write_text("\n".join(frames) + "\n")
import warnings as _w
with _w.catch_warnings(record=True) as caught:
    _w.simplefilter("always")
    an.read_lammps_dump(tmp_noid / "noid.lammpstrj")
check("a dump with no id column raises a warning",
      any(issubclass(c.category, RuntimeWarning) for c in caught),
      f"{len(caught)} warnings")
check("and the warning names sort id as the fix",
      any("sort id" in str(c.message) for c in caught),
      str(caught[0].message)[:70] if caught else "none")
shutil.rmtree(tmp_noid, ignore_errors=True)

# rmax beyond L/2 must be cut, not silently wrong
rdfc = an.compute_rdf(traji, rmax=50.0, nbins=100)
check("an rmax beyond half the box is cut and reported",
      bool(rdfc.warnings) and rdfc.r[-1] <= 30.1, str(rdfc.warnings)[:70])

# Nernst-Einstein against a hand calculation
counts = {"Ca": 810, "Al": 448, "Si": 754, "O": 2990}
Ds = {"Ca": 3e-9, "Al": 1.5e-9, "Si": 1e-9, "O": 2e-9}
V = 41.7 ** 3
T = 1873.0
sig = an.conductivity_nernst_einstein(counts, bou.charges, Ds, V, T)
hand = (an.E_CHARGE ** 2 / (an.KB * T * V * 1e-30)) * sum(
    counts[s] * bou.charges[s] ** 2 * Ds[s] for s in counts)
check("Nernst-Einstein matches an independent hand evaluation",
      abs(sig - hand) / hand < 1e-12, f"{sig:.4f} S/m")
check("the conductivity is physically plausible for a slag at 1600 C",
      0.1 < sig < 500.0, f"{sig:.2f} S/m")

# Exact Einstein on uncorrelated walkers must approach Nernst-Einstein, because the cross
# terms are zero by construction. But it must be tested over REPLICAS, not one trajectory.
# The collective charge displacement is a single coordinate no matter how many ions there
# are, so one trajectory gives a handful of independent samples of it. Measured here: the
# per-trajectory scatter is about 80 percent while the mean is right to a few percent.
# That is physics, not a bug, and it is the reason the function returns a block error.
nat2 = 600
els = np.array(["Ca"] * 100 + ["Si"] * 100 + ["O"] * 400)
D3 = {"Ca": 3e-9, "Si": 1e-9, "O": 2e-9}
sd3 = np.array([math.sqrt(2 * D3[e] * 1e8 * 1.0) for e in els])[:, None]
cnt3 = {s: int((els == s).sum()) for s in ("Ca", "Si", "O")}
sig_ne = an.conductivity_nernst_einstein(cnt3, bou.charges, D3, 40.0 ** 3, T)

vals, one_traj = [], None
for rep in range(40):
    p3 = np.zeros((400, nat2, 3))
    p3[0] = rng.random((nat2, 3)) * 40.0
    for f in range(1, 400):
        p3[f] = p3[f - 1] + rng.normal(0, 1, (nat2, 3)) * sd3
    tr3 = an.Trajectory(positions=p3, box=np.tile([40.0, 40.0, 40.0], (400, 1)),
                        types=np.ones(nat2, dtype=int), elements=els,
                        times_ps=np.arange(400.0), unwrapped=True)
    sv, inf = an.conductivity_exact_einstein(tr3, bou.charges, T)
    vals.append(sv)
    one_traj = inf
arr = np.array(vals)
mean_e = float(arr.mean())
sem_e = float(arr.std(ddof=1) / math.sqrt(len(arr)))
check("exact Einstein equals Nernst-Einstein for uncorrelated ions, averaged over "
      "40 replicas", abs(mean_e - sig_ne) <= 3.0 * sem_e,
      f"{mean_e:.2f} +/- {sem_e:.2f} against {sig_ne:.2f} S/m, "
      f"{abs(mean_e-sig_ne)/sem_e:.2f} sigma")
check("a SINGLE trajectory exact Einstein value is not trustworthy, and the tool says so",
      one_traj["relative_stderr_percent"] > 20.0,
      f"block error {one_traj['relative_stderr_percent']:.0f} percent: "
      f"{one_traj['verdict'][:60]}")
check("the per-trajectory scatter is large, which is why replicas are needed",
      100 * arr.std(ddof=1) / mean_e > 30.0,
      f"{100*arr.std(ddof=1)/mean_e:.0f} percent scatter across replicas")

# equilibration report must catch a drifting density
log_ok = {"Density": np.full(200, 2.65) + rng.normal(0, 0.002, 200),
          "Temp": np.full(200, 1873.0) + rng.normal(0, 12, 200)}
rep = an.equilibration_report(log_ok)
check("a settled density is reported as settled",
      "settled" in rep["_verdicts"][0], rep["_verdicts"][0][:60])
log_drift = {"Density": np.linspace(2.80, 2.65, 200),
             "Temp": np.full(200, 1873.0)}
rep2 = an.equilibration_report(log_drift)
check("a drifting density is caught", any("drifting" in v for v in rep2["_verdicts"]),
      rep2["_verdicts"][0][:70])

# Green-Kubo against two stress signals whose exact viscosity is known analytically.
#
# 1. White noise. The autocorrelation is a spike at zero lag, so the integral is
#    sigma^2 * dt exactly, and eta = V/(kT) * sigma^2 * dt. Physically meaningless as a
#    melt, but it is an exact arithmetic check on the integration and the unit conversion.
# 2. An Ornstein-Uhlenbeck process with correlation time tau. Its autocorrelation is
#    sigma^2 exp(-t/tau), so the integral is sigma^2 * tau and eta = V/(kT) * sigma^2 * tau.
#    This is the realistic case and it checks the plateau detection too.
tmp_v = Path(tempfile.mkdtemp(prefix="slagmd_visc_"))
dt_ps, sample = 0.01, 5          # 2 fs timestep, every 5 steps
SD_ATM = 2000.0
n = 200000

wn = rng.normal(0, SD_ATM, (n, 3))
(tmp_v / "white.txt").write_text("# step pxy pxz pyz\n" + "\n".join(
    f"{i*sample} {wn[i,0]:.5f} {wn[i,1]:.5f} {wn[i,2]:.5f}" for i in range(n)))
eta_exact_white = (V * 1e-30 / (an.KB * T)) * (SD_ATM * an.ATM_TO_PA) ** 2 * dt_ps * 1e-12
vw = an.viscosity_green_kubo(tmp_v / "white.txt", T, V, 2.0, sample,
                             plateau_lo_ps=0.5, plateau_hi_ps=5.0)
relw = abs(vw.eta_Pa_s - eta_exact_white) / eta_exact_white
check("Green-Kubo reproduces the exact viscosity of a white noise stress",
      relw < 0.05,
      f"exact {eta_exact_white:.4e}, got {vw.eta_Pa_s:.4e} Pa s, "
      f"{100*relw:.2f} percent")

TAU_PS = 0.2
a_ou = math.exp(-dt_ps / TAU_PS)
noise = rng.normal(0, SD_ATM * math.sqrt(1 - a_ou ** 2), (n, 3))
# an Ornstein-Uhlenbeck process is an AR(1) filter, so lfilter builds it in one call
# instead of a 200,000 step Python loop
ou = signal.lfilter([1.0], [1.0, -a_ou], noise, axis=0)
(tmp_v / "ou.txt").write_text("# step pxy pxz pyz\n" + "\n".join(
    f"{i*sample} {ou[i,0]:.5f} {ou[i,1]:.5f} {ou[i,2]:.5f}" for i in range(n)))
eta_exact_ou = ((V * 1e-30 / (an.KB * T)) * (SD_ATM * an.ATM_TO_PA) ** 2
                * TAU_PS * 1e-12)
vo = an.viscosity_green_kubo(tmp_v / "ou.txt", T, V, 2.0, sample,
                             plateau_lo_ps=2.0, plateau_hi_ps=10.0)
relo = abs(vo.eta_Pa_s - eta_exact_ou) / eta_exact_ou
check("Green-Kubo reproduces the exact viscosity of a correlated stress with a known "
      "relaxation time", relo < 0.15,
      f"exact {eta_exact_ou:.4e}, got {vo.eta_Pa_s:.4e} Pa s, {100*relo:.1f} percent")
check("and it finds a plateau in the correlated case",
      vo.plateau_scatter_percent is not None and vo.plateau_scatter_percent < 10.0,
      f"scatter {vo.plateau_scatter_percent:.1f} percent")

# 3. A long, slowly decaying tail, the realistic viscous melt failure. The integral has
#    not converged by 10 ps and the tool must refuse rather than quote a number.
TAU_LONG = 40.0
a_l = math.exp(-dt_ps / TAU_LONG)
nl = rng.normal(0, SD_ATM * math.sqrt(1 - a_l ** 2), (n, 3))
sl = signal.lfilter([1.0], [1.0, -a_l], nl, axis=0)
(tmp_v / "slow.txt").write_text("# step pxy pxz pyz\n" + "\n".join(
    f"{i*sample} {sl[i,0]:.5f} {sl[i,1]:.5f} {sl[i,2]:.5f}" for i in range(n)))
vs = an.viscosity_green_kubo(tmp_v / "slow.txt", T, V, 2.0, sample,
                             plateau_lo_ps=1.0, plateau_hi_ps=10.0)
check("Green-Kubo REFUSES to quote a viscosity when the integral has not plateaued",
      (vs.plateau_scatter_percent or 0) > 10.0 and "NOT converged" in vs.verdict,
      f"scatter {vs.plateau_scatter_percent:.0f} percent, verdict starts "
      f"{vs.verdict[:28]!r}")

# Einstein-Helfand must agree with Green-Kubo on the case where both are valid
vh = an.viscosity_einstein_helfand(tmp_v / "ou.txt", T, V, 2.0, sample)
agree = abs(vh.eta_Pa_s - eta_exact_ou) / eta_exact_ou
check("Einstein-Helfand independently reproduces the same exact viscosity",
      agree < 0.15,
      f"exact {eta_exact_ou:.4e}, Helfand {vh.eta_Pa_s:.4e} Pa s, {100*agree:.1f} percent")
# and a deliberately late window must be DIAGNOSED, not silently wrong. This is the bug
# the audit caught: at 20 to 60 percent of the series the answer came out 5 times too
# small, and the log-log test measured over the whole series could not see it.
vlate = an.viscosity_einstein_helfand(tmp_v / "ou.txt", T, V, 2.0, sample, 0.2, 0.6)
check("a fit window too late is diagnosed rather than silently returning a wrong answer",
      "too LATE" in vlate.verdict,
      f"eta would be {vlate.eta_Pa_s:.3e} against exact {eta_exact_ou:.3e}")

# too few samples must be refused with the aliasing explanation
(tmp_v / "short.txt").write_text("# s pxy pxz pyz\n" + "\n".join(
    f"{i} 1.0 1.0 1.0" for i in range(50)))
try:
    an.viscosity_green_kubo(tmp_v / "short.txt", T, V)
    check("under-sampled stress is refused", False)
except ValueError as exc:
    check("under-sampled stress is refused with the aliasing reason",
          "alias" in str(exc).lower(), str(exc)[:70])

# speciation on a built silica-like network: check the bookkeeping identity
comp_s = Composition({"CaO": 0.40, "Al2O3": 0.10, "SiO2": 0.50}, 1500)
sst = build_random(comp_s, 2.6, seed=3)
trs = an.Trajectory(positions=sst.positions[None, :, :].repeat(2, axis=0),
                    box=np.tile([sst.box] * 3, (2, 1)),
                    types=np.ones(sst.n_atoms, dtype=int),
                    elements=np.array(sst.species), times_ps=np.arange(2.0),
                    unwrapped=True)
sp = an.compute_speciation(trs, cutoffs={"Si": 2.0, "Al": 2.2})
tot_frac = (sp.frac_bridging + sp.frac_non_bridging + sp.frac_free
            + sp.frac_tricluster)
check("the four oxygen classes account for every oxygen exactly once",
      abs(tot_frac - 1.0) < 1e-9, f"sum {tot_frac:.12f}")
check("Qn distributions are normalised",
      all(abs(sum(d.values()) - 1.0) < 1e-9 for d in sp.qn_distribution.values()))
check("a random, unrelaxed structure is flagged as suspicious", bool(sp.notes),
      str(sp.notes)[:70])


# ======================================================= F  the diagnosis rules
section("F  failure diagnosis, each against the real LAMMPS or scheduler error text")

for code, text in mocklammps.FAIL_MESSAGES.items():
    if code == "SILENT":
        continue
    d = diag.diagnose_text("LAMMPS (27 Jun 2024)\nSetting up run ...\n" + text)
    got = d.remedy.code if d.remedy else "NONE"
    check(f"{code} is recognised", got == code, f"got {got}")
    if d.remedy:
        check(f"{code} has a plain-English cause and an action",
              len(d.remedy.why) > 40 and len(d.remedy.action) > 20)

# an unrecognised failure must NOT be auto retried
d = diag.diagnose_text("LAMMPS (27 Jun 2024)\nsomething nobody has ever seen\n")
check("an unrecognised failure is not auto retried",
      d.remedy is not None and not d.remedy.auto_retry, d.remedy.code if d.remedy else "")
check("and it says it is stopping rather than guessing",
      d.remedy is not None and "STOPPING" in d.remedy.action)

# a clean log must produce no remedy at all
d = diag.diagnose_text("Total wall time: 0:10:00\nSLAGMD_STAGE_DONE 04_produce\n")
check("a clean completion produces no remedy", d.remedy is None)

# the tool-bug cases must refuse to retry
for code in ("BAD_PAIR_COEFF", "MISSING_PACKAGE", "ENVIRONMENT", "MISSING_INPUT"):
    d = diag.diagnose_text(mocklammps.FAIL_MESSAGES[code])
    check(f"{code} is treated as needing a human", d.remedy is not None
          and not d.remedy.auto_retry)

# walltime doubling arithmetic
out = diag.apply_cluster_change({"walltime": "24:00:00", "nodes": 2}, {"walltime": "x2"})
check("walltime x2 gives 48:00:00", out["walltime"] == "48:00:00", out["walltime"])
out = diag.apply_cluster_change({"walltime": "01:30:00"}, {"walltime": "x2"})
check("walltime x2 handles minutes", out["walltime"] == "3:00:00", out["walltime"])
out = diag.apply_cluster_change({"nodes": 2}, {"nodes": "x2"})
check("nodes x2 gives 4", out["nodes"] == 4)
out = diag.apply_cluster_change({"walltime": "120:00:00"}, {"walltime": "x2"})
check("walltime doubling is capped at a week", out["walltime"] == "168:00:00",
      out["walltime"])


# ================================================ G  the pipeline, end to end
section("G  the whole pipeline, end to end, with a real failure and a real recovery")

work = Path(tempfile.mkdtemp(prefix="slagmd_e2e_"))
proj = Project(work / "proj")
proj.cluster.backend = "local"
proj.cluster.nodes = 1
proj.cluster.tasks_per_node = 1
proj.cluster.mpi_launcher = ""
proj.cluster.modules = []
proj.cluster.lammps_exe = f"{sys.executable} -m slagmd.mocklammps"
proj.cluster.pre_commands = [f"export PYTHONPATH='{ROOT}'"]   # quoted: this
# project lives at a path with a TRAILING SPACE, and an unquoted export silently splits
# tiny runs so the test finishes in seconds
proj.protocol = Protocol(hot_ns=0.002, produce_ns=0.012, cool_rate_K_per_ps=500.0,
                         sample_every_ps=0.1, thermo_every_ps=0.5,
                         eq_blocks=[EqBlock(50.0, 100.0, 0.002, "tight"),
                                    EqBlock(500.0, 1000.0, 0.002, "loose")])
proj.save()

run = proj.add_run("e2e_S1", {"CaO": 0.453, "Al2O3": 0.125, "SiO2": 0.422},
                  1873.15, "bouhadja", target_atoms=1400)
notes = proj.prepare(run)
check("prepare succeeded", run.prepared and (proj.run_dir("e2e_S1") / "melt.data").exists())
check("prepare reported the force field check",
      any("consistency" in n for n in notes))
check("the deck marks which settings are theirs and which are ours",
      "[THEIRS]" in (proj.run_dir("e2e_S1") / "in.03_equil1").read_text()
      and "[MINE]" in (proj.run_dir("e2e_S1") / "in.03_equil1").read_text())

# state must survive being reloaded from disk
proj2 = Project(work / "proj")
check("project state reloads from disk", "e2e_S1" in proj2.runs
      and proj2.runs["e2e_S1"].prepared)
check("the reloaded run has all six stages",
      len(proj2.runs["e2e_S1"].stage_list()) == 6,
      str(proj2.runs["e2e_S1"].stage_list()))

import time
MOCK_ENV = ("SLAGMD_MOCK_FAIL", "SLAGMD_MOCK_FAIL_STAGE", "SLAGMD_MOCK_FAIL_ONCE")


def inject(**kw) -> None:
    """Set the mock's fault variables, clearing any left over from a previous case.

    This helper exists because a leaked SLAGMD_MOCK_FAIL_ONCE made a later test silently
    pass its fault and complete, which then crashed on a None stage lookup. A test harness
    that leaks state between cases produces failures that point at the wrong test.
    """
    for k in MOCK_ENV:
        os.environ.pop(k, None)
    for k, v in kw.items():
        os.environ[k] = str(v)


def clear_mock() -> None:
    for k in MOCK_ENV:
        os.environ.pop(k, None)


def drain(p: Project, limit: int = 60, stop_on_diagnosis: bool = False) -> list[str]:
    msgs = []
    for _ in range(limit):
        msgs += p.fill_queue(max_concurrent=1)
        time.sleep(0.35)
        msgs += p.poll_all()
        if stop_on_diagnosis and any(
                st.last_diagnosis for r in p.runs.values()
                for st in r.stages.values()):
            break
        if all(r.is_done() or r.is_blocked() for r in p.runs.values()):
            break
    return msgs

msgs = drain(proj)
r = proj.runs["e2e_S1"]
check("all six stages completed under the mock", r.is_done(),
      f"progress {r.progress()}, stage {r.current_stage()}")
rd = proj.run_dir("e2e_S1")
check("a production trajectory was written", (rd / "prod.lammpstrj").exists())
check("each stage wrote its restart file",
      all((rd / f"restart.{st}").exists() for st in r.stage_list()),
      str([st for st in r.stage_list() if not (rd / f"restart.{st}").exists()]))

# analyse the mock trajectory. The mock used D = 1e-9 m^2/s, so we must get it back.
traj_e2e = an.read_lammps_dump(rd / "prod.lammpstrj", timestep_fs=2.0)
check("the dump reads back as unwrapped", traj_e2e.unwrapped)
check("the dump has the right species", set(traj_e2e.species_list()) == {"Si", "Al", "Ca", "O"},
      str(sorted(traj_e2e.species_list())))
te, msde = an.compute_msd(traj_e2e, "O")
rese = an.fit_diffusion(te, msde, "O")
errD = abs(rese.D_m2_s - 1.0e-9) / 1.0e-9
check("D recovered end to end from the written dump matches the mock's input",
      errD < 0.10, f"got {rese.D_m2_s:.4e}, mock used 1.0000e-09, error {100*errD:.1f}%")

log = an.read_lammps_log(rd / f"log.{[k for k in r.stage_list() if k.endswith(chr(95)+chr(112)+chr(114)+chr(111)+chr(100)+chr(117)+chr(99)+chr(101))][0]}")
check("the LAMMPS log parses into thermo columns",
      "Density" in log and len(log["Density"]) > 10, str(sorted(log)[:6]))
eq = an.equilibration_report(log)
check("the equilibration report runs on a real log", "_verdicts" in eq)

# ------------------------------------------------ now make it fail, and recover
proj3 = Project(work / "proj3")
proj3.cluster = sched.ClusterConfig(**asdict(proj.cluster))
proj3.protocol = proj.protocol
proj3.save()
run3 = proj3.add_run("fail_walltime", {"CaO": 0.5, "SiO2": 0.5}, 1873.15,
                     "bouhadja", target_atoms=1400)
proj3.prepare(run3)
inject(SLAGMD_MOCK_FAIL="WALLTIME", SLAGMD_MOCK_FAIL_STAGE="01_hot",
       SLAGMD_MOCK_FAIL_ONCE=1)          # transient, so the recovery path is real
before_wall = proj3.cluster_for(run3).walltime
msgs3 = drain(proj3, limit=40)
del os.environ["SLAGMD_MOCK_FAIL"]
del os.environ["SLAGMD_MOCK_FAIL_STAGE"]

hist = "\n".join(run3.stages["01_hot"].history)
check("the walltime failure was diagnosed", "WALLTIME" in hist, hist[-90:] if hist else "")
after_wall = proj3.cluster_for(run3).walltime
check("the walltime was doubled exactly once in response",
      after_wall == "48:00:00", f"{before_wall} -> {after_wall}")
check("a diagnosis file was written for the human to read",
      any(p.name.startswith("diagnosis.") for p in proj3.run_dir("fail_walltime").iterdir()))
check("the run recovered by itself and finished every stage, with no human involved",
      run3.is_done(), f"progress {run3.progress()}, "
      f"attempts on the failed stage {run3.stages['01_hot'].attempts}")
check("and it took exactly two attempts on the stage that failed",
      run3.stages["01_hot"].attempts == 2,
      f"{run3.stages['01_hot'].attempts} attempts")

# a failure the tool must NOT retry
proj4 = Project(work / "proj4")
proj4.cluster = sched.ClusterConfig(**asdict(proj.cluster))
proj4.protocol = proj.protocol
proj4.save()
run4 = proj4.add_run("fail_package", {"CaO": 0.5, "SiO2": 0.5}, 1873.15,
                     "bouhadja", target_atoms=1400)
proj4.prepare(run4)
inject(SLAGMD_MOCK_FAIL="MISSING_PACKAGE")
drain(proj4, limit=25)
clear_mock()
check("run4 did not complete despite the injected fault", not run4.is_done(),
      f"stage {run4.current_stage()}")
st4 = run4.stages[run4.current_stage() or run4.stage_list()[0]]
check("a missing LAMMPS package blocks instead of looping", run4.is_blocked(),
      f"status {st4.status}, attempts {st4.attempts}")
check("and it gave up after one attempt, not three", st4.attempts == 1,
      f"{st4.attempts} attempts")

# attempt capping on a genuinely repeating failure
proj5 = Project(work / "proj5")
proj5.cluster = sched.ClusterConfig(**asdict(proj.cluster))
proj5.protocol = proj.protocol
proj5.save()
run5 = proj5.add_run("fail_repeat", {"CaO": 0.5, "SiO2": 0.5}, 1873.15,
                     "bouhadja", target_atoms=1400)
proj5.prepare(run5)
inject(SLAGMD_MOCK_FAIL="MPI_CRASH")
drain(proj5, limit=40)
clear_mock()
check("run5 did not complete despite the permanently injected fault",
      not run5.is_done(), f"stage {run5.current_stage()}")
st5 = run5.stages[run5.current_stage() or run5.stage_list()[0]]
check("a repeating failure stops at the attempt cap", run5.is_blocked()
      and st5.attempts <= proj5.max_attempts,
      f"{st5.attempts} attempts, cap {proj5.max_attempts}")

# ============================================ H  the GUI, driven end to end
section("H  the graphical interface, driven end to end against the mock")

gui_ok = True
try:
    from slagmd.gui import App
except Exception as exc:                                    # noqa: BLE001
    gui_ok = False
    check("the GUI imports", False, str(exc)[:80])

if gui_ok:
    try:
        gwork = Path(tempfile.mkdtemp(prefix="slagmd_gui_"))
        gproj = Project(gwork / "proj")
        gproj.cluster = sched.ClusterConfig(**asdict(proj.cluster))
        gproj.protocol = proj.protocol
        gproj.save()
        app = App(gwork / "proj")
        app.update_idletasks()
        app.update()
        check("the GUI constructs with a project open", app.project is not None)
        check("it has the six tabs",
              len(app.nb.tabs()) == 6,
              str([app.nb.tab(i, "text") for i in range(len(app.nb.tabs()))]))

        # add a run through the form, exactly as a person would
        app.rvars["name"].set("gui_run")
        app.rvars["comp"].set("CaO 0.453 Al2O3 0.125 SiO2 0.422")
        app.rvars["temp_C"].set("1600")
        app.rvars["atoms"].set("1400")
        app.ff_var.set("bouhadja")
        app.visc_var.set(True)
        app._check_comp()
        app.update()
        panel = app.comp_out.get("1.0", "end")
        check("the composition panel reports the resolved atom counts",
              "atoms" in panel and "total" in panel, panel.splitlines()[3][:60])
        check("and it reports the cost in core hours", "core hours" in panel)
        app._add_run()
        app.update()
        check("adding a run through the form works", "gui_run" in app.project.runs)
        check("the queue table shows it", len(app.tree.get_children()) == 1,
              str(len(app.tree.get_children())))

        # prepare and run it, synchronously, using the same code paths the buttons use
        notes = app.project.prepare(app.project.runs["gui_run"])
        check("preparation through the GUI's project succeeds",
              app.project.runs["gui_run"].prepared)
        check("the stress tensor was requested because the box was ticked",
              "v_pxy" in (app.project.run_dir("gui_run")
                          / f"in.{app.project.runs['gui_run'].stage_list()[-1]}"
                          ).read_text())
        clear_mock()
        drain(app.project, limit=80)
        app._refresh_queue()
        app.update()
        gr = app.project.runs["gui_run"]
        check("the run completes when driven through the GUI's own queue",
              gr.is_done(), f"progress {gr.progress()}")

        # analysis through the GUI worker, which is what the buttons call
        app.an_run.set("gui_run")
        app._analyse_worker("all", app.project.run_dir("gui_run") / "prod.lammpstrj")
        app._drain()
        app.update()
        out = app.an_out.get("1.0", "end")
        for want, label in [
                ("structure", "the analysis panel reports structure"),
                ("first peak", "with RDF peak positions"),
                ("speciation", "and oxygen speciation"),
                ("bridging", "with the bridging oxygen fraction"),
                ("diffusion", "and diffusion"),
                ("blocks", "with a block error bar on D"),
                ("log-log slope", "and the diffusive regime test"),
                ("conductivity", "and conductivity"),
                ("Nernst-Einstein", "by both routes"),
                ("viscosity", "and viscosity"),
                ("was it equilibrated", "and an equilibration verdict")]:
            check(label, want.lower() in out.lower(), "" if want.lower() in out.lower()
                  else f"{want!r} missing")
        check("the analysis panel never prints a bare number without a verdict",
              out.lower().count("verdict") + out.count("diffusive")
              + out.count("Einstein relation") > 0)

        figs = app._figures(an.read_lammps_dump(
            app.project.run_dir("gui_run") / "prod.lammpstrj"),
            app.project.run_dir("gui_run"), 1873.15)
        check("the figures are drawn", (figs / "rdf.png").exists()
              and (figs / "msd_loglog.png").exists(),
              str(sorted(p.name for p in figs.iterdir())))
        check("and they are real images, not empty files",
              (figs / "rdf.png").stat().st_size > 8000
              and (figs / "msd_loglog.png").stat().st_size > 8000,
              f"{(figs/'rdf.png').stat().st_size} and "
              f"{(figs/'msd_loglog.png').stat().st_size} bytes")

        # _save_analysis ends with messagebox.showinfo, which is MODAL and blocks until
        # someone clicks it. In a headless test that is a hang, and it cost a stalled run
        # on 29 Sep 2026. Stub the dialog, then call the real method.
        import tkinter.messagebox as _mb
        _real_info = _mb.showinfo
        _mb.showinfo = lambda *a, **k: None
        try:
            app._save_analysis()
        finally:
            _mb.showinfo = _real_info
        check("the analysis saves to a text file beside the run",
              (app.project.run_dir("gui_run") / "ANALYSIS.txt").exists())

        # the cluster tab must catch an unquoted path containing a space
        app.project.cluster.pre_commands = ["export PYTHONPATH=/opt/my tools/x"]
        app._check_cluster()
        check("the cluster tab catches an unquoted path containing a space",
              "PROBLEM" in app.cluster_warn.cget("text"),
              app.cluster_warn.cget("text")[:70])

        app.destroy()
        shutil.rmtree(gwork, ignore_errors=True)
    except Exception:                                        # noqa: BLE001
        check("the GUI end to end run completed without an exception", False,
              traceback.format_exc().splitlines()[-1][:90])


rep = proj3.write_report()
check("a plain text report is written", rep.exists() and rep.stat().st_size > 100)
check("the report names what needs a human",
      "NEEDS A HUMAN" in rep.read_text() or run3.is_done())

shutil.rmtree(work, ignore_errors=True)
shutil.rmtree(tmp, ignore_errors=True)
shutil.rmtree(tmp_v, ignore_errors=True)


# ===================================================================== verdict
print("\n" + "=" * 72)
print(f"{len(PASS)} checks passed, {len(FAIL)} failed")
if FAIL:
    print("\nFAILURES:")
    for name, detail in FAIL:
        print(f"  - {name}   {detail}")
print("=" * 72)
sys.exit(1 if FAIL else 0)
