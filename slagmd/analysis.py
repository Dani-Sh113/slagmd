"""Trajectory analysis: RDF, coordination numbers, oxygen speciation, Qn, MSD and
diffusion, electrical conductivity, and viscosity.

Every routine here that produces a number also produces the evidence that the number is
allowed to exist. That is not politeness, it is the difference between a diffusion
coefficient and a fitted slope of a curve that was never diffusive. Specifically:

  * diffusion reports the log-log slope of the MSD over the fitting window. If that slope
    is not close to 1 the melt is still in the cage regime and the tool refuses to call
    the answer a diffusion coefficient.
  * conductivity is computed BOTH ways that Verbeeck 2025 defines, the approximate
    Nernst-Einstein of their equation 7 and the exact Einstein of their equation 8, and
    the ratio of the two is reported, because that ratio IS the correlation the
    approximation throws away.
  * viscosity by Green-Kubo reports the running integral, not one number, because the
    integral of a stress autocorrelation in a viscous melt does not plateau cleanly and
    anyone who quotes a single Green-Kubo viscosity without showing the plateau is
    quoting noise.

Units throughout: distance Angstrom, time picoseconds, temperature kelvin. Diffusion is
returned in m^2/s and in cm^2/s, conductivity in S/m, viscosity in Pa s.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# physical constants, SI, CODATA 2018
KB = 1.380649e-23            # J/K
E_CHARGE = 1.602176634e-19   # C
N_AV = 6.02214076e23
ATM_TO_PA = 101325.0

# defined once, in composition.py, so a new network former is added in one place only
from .composition import NETWORK_FORMERS  # noqa: E402,F401


def _msd_fft(pos: np.ndarray) -> np.ndarray:
    """Mean squared displacement against lag, averaged over all time origins, by FFT.

    pos has shape (n_frames, n_particles, n_dim) or (n_frames, n_dim). Returns an array of
    length n_frames, index = lag in frames.

    The direct double loop is O(n_frames^2) per particle, which is fine for a few hundred
    frames and unusable for the hundred thousand samples a Green-Kubo stress series needs.
    This is the standard identity, from Wiener and Khinchin:

        MSD(m) = S1(m) - 2 S2(m)
        S1(m)  = (1/(N-m)) sum_k  |r_k|^2 + |r_{k+m}|^2
        S2(m)  = (1/(N-m)) sum_k  r_k . r_{k+m}      <- an autocorrelation, so FFT it

    S1 is built by a running subtraction rather than a sum per lag, so the whole thing is
    O(N log N). Checked against the direct loop in the audit suite.
    """
    if pos.ndim == 2:
        pos = pos[:, None, :]
    n = pos.shape[0]
    sq = (pos ** 2).sum(axis=-1)                 # (n, p)
    # S2 by FFT, per particle per dimension, summed over dimensions
    size = 2 * n
    f = np.fft.rfft(pos, n=size, axis=0)
    acf = np.fft.irfft(f * np.conjugate(f), n=size, axis=0)[:n].real.sum(axis=-1)
    denom = (n - np.arange(n)).astype(float)[:, None]
    s2 = acf / denom
    # S1 by running subtraction
    sq_ext = np.concatenate([sq, np.zeros((1, sq.shape[1]))], axis=0)
    q = 2.0 * sq.sum(axis=0)                     # (p,)
    s1 = np.empty_like(s2)
    for m in range(n):
        if m > 0:
            q = q - sq_ext[m - 1] - sq_ext[n - m]
        s1[m] = q / (n - m)
    msd = (s1 - 2.0 * s2)
    return msd.mean(axis=1)


# =============================================================== trajectory reading
@dataclass
class Trajectory:
    """A LAMMPS custom dump, read into memory.

    positions: (n_frames, n_atoms, 3) in Angstrom
    box:       (n_frames, 3) edge lengths
    types:     (n_atoms,) integer LAMMPS types
    elements:  (n_atoms,) species strings, if the dump had an element column
    times_ps:  (n_frames,)
    unwrapped: whether the coordinate columns were xu yu zu
    """
    positions: np.ndarray
    box: np.ndarray
    types: np.ndarray
    elements: np.ndarray
    times_ps: np.ndarray
    unwrapped: bool

    @property
    def n_frames(self) -> int:
        return self.positions.shape[0]

    @property
    def n_atoms(self) -> int:
        return self.positions.shape[1]

    def species_list(self) -> list[str]:
        seen: list[str] = []
        for e in self.elements:
            if e not in seen:
                seen.append(e)
        return seen

    def index_of(self, species: str) -> np.ndarray:
        return np.flatnonzero(self.elements == species)

    def mean_box(self) -> np.ndarray:
        return self.box.mean(axis=0)

    def volume_A3(self) -> float:
        return float(np.prod(self.box, axis=1).mean())


def read_lammps_dump(path: Path, timestep_fs: float = 2.0,
                     type_to_element: dict[int, str] | None = None,
                     max_frames: int | None = None,
                     stride: int = 1) -> Trajectory:
    """Read a LAMMPS custom dump. Handles x y z and xu yu zu, with or without element."""
    path = Path(path)
    frames_pos: list[np.ndarray] = []
    frames_box: list[np.ndarray] = []
    steps: list[int] = []
    types: np.ndarray | None = None
    elements: np.ndarray | None = None
    unwrapped = False
    frame_index = -1

    with path.open() as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            frame_index += 1
            step = int(fh.readline().split()[0])
            assert fh.readline().startswith("ITEM: NUMBER OF ATOMS")
            natoms = int(fh.readline())
            bounds_header = fh.readline()
            assert bounds_header.startswith("ITEM: BOX BOUNDS")
            edges = []
            for _ in range(3):
                lo, hi = (float(v) for v in fh.readline().split()[:2])
                edges.append(hi - lo)
            atoms_header = fh.readline()
            assert atoms_header.startswith("ITEM: ATOMS")
            cols = atoms_header.split()[2:]

            skip_this = (frame_index % stride != 0)
            if skip_this:
                for _ in range(natoms):
                    fh.readline()
                continue

            try:
                ix, iy, iz = cols.index("xu"), cols.index("yu"), cols.index("zu")
                unwrapped = True
            except ValueError:
                ix, iy, iz = cols.index("x"), cols.index("y"), cols.index("z")
                unwrapped = False
            i_id = cols.index("id") if "id" in cols else None
            i_type = cols.index("type") if "type" in cols else None
            i_el = cols.index("element") if "element" in cols else None

            pos = np.empty((natoms, 3))
            tp = np.empty(natoms, dtype=int)
            el = np.empty(natoms, dtype=object)
            for k in range(natoms):
                f = fh.readline().split()
                row = (int(f[i_id]) - 1) if i_id is not None else k
                pos[row] = (float(f[ix]), float(f[iy]), float(f[iz]))
                if i_type is not None:
                    tp[row] = int(f[i_type])
                if i_el is not None:
                    el[row] = f[i_el]

            frames_pos.append(pos)
            frames_box.append(np.array(edges))
            steps.append(step)
            if types is None:
                types = tp
                if i_el is not None:
                    elements = el.astype(str)
            if max_frames and len(frames_pos) >= max_frames:
                break

    if not frames_pos:
        raise ValueError(f"no frames found in {path}")
    if elements is None:
        if type_to_element:
            elements = np.array([type_to_element.get(int(t), f"T{t}")
                                 for t in types], dtype=str)   # type: ignore[union-attr]
        else:
            elements = np.array([f"T{t}" for t in types], dtype=str)  # type: ignore

    if i_id is None:
        # Without an id column there is no way to align atoms between frames, and LAMMPS
        # does NOT write them in a fixed order unless the deck said `dump_modify sort id`.
        # Danish's own July 2026 practice run had exactly this: six different atom
        # orderings across 52 frames, only 2 of 52 matching the first, which made his mean
        # squared displacement compare different atoms between frames and produced a
        # NEGATIVE diffusion coefficient. His dump did carry an id, so it was recoverable.
        # A dump without one is not, so this is a loud warning and not a silent fallback.
        warnings.warn(
            f"{path.name} has no 'id' column, so atoms are being read in file order. "
            "LAMMPS only writes a fixed order when the input said 'dump_modify <name> "
            "sort id'. If it did not, every per-atom quantity computed from this "
            "trajectory, above all the mean squared displacement and therefore the "
            "diffusion coefficient, is comparing different atoms between frames and is "
            "WRONG. There is no way to repair it after the fact without the ids. "
            "Re-run with sort id.",
            RuntimeWarning, stacklevel=2)

    times_ps = np.array(steps, dtype=float) * timestep_fs / 1000.0
    return Trajectory(positions=np.array(frames_pos), box=np.array(frames_box),
                      types=types, elements=elements,                  # type: ignore[arg-type]
                      times_ps=times_ps - times_ps[0], unwrapped=unwrapped)


def read_xyz(path: Path, dt_ps: float = 1.0, box: float | None = None) -> Trajectory:
    """Read a plain xyz trajectory, the format his own analysis already writes."""
    path = Path(path)
    frames_pos: list[np.ndarray] = []
    elements: list[str] | None = None
    with path.open() as fh:
        while True:
            head = fh.readline()
            if not head.strip():
                break
            n = int(head.split()[0])
            comment = fh.readline()
            if box is None:
                for tok in comment.replace("=", " ").split():
                    try:
                        v = float(tok)
                        if v > 5.0:
                            box = v
                            break
                    except ValueError:
                        continue
            pos = np.empty((n, 3))
            els: list[str] = []
            for k in range(n):
                f = fh.readline().split()
                els.append(f[0])
                pos[k] = (float(f[1]), float(f[2]), float(f[3]))
            frames_pos.append(pos)
            if elements is None:
                elements = els
    if not frames_pos:
        raise ValueError(f"no frames in {path}")
    nf = len(frames_pos)
    L = box if box else 0.0
    return Trajectory(positions=np.array(frames_pos),
                      box=np.tile(np.array([L, L, L], dtype=float), (nf, 1)),
                      types=np.arange(len(elements or [])),
                      elements=np.array(elements, dtype=str),
                      times_ps=np.arange(nf) * dt_ps,
                      unwrapped=False)


# =========================================================================== RDF
@dataclass
class RDFResult:
    r: np.ndarray
    g: dict[str, np.ndarray]                  # "Si-O" -> g(r)
    cn: dict[str, np.ndarray]                 # running coordination number
    first_peak: dict[str, float]              # r at the first maximum, Angstrom
    first_min: dict[str, float]               # r at the first minimum after it
    cn_at_first_min: dict[str, float]         # the coordination number
    n_frames_used: int
    warnings: list[str] = field(default_factory=list)


def compute_rdf(traj: Trajectory, rmax: float = 10.0, nbins: int = 400,
                frames: slice | None = None,
                pairs: list[tuple[str, str]] | None = None) -> RDFResult:
    """Partial radial distribution functions with the minimum image convention.

    Normalisation is the standard one: the histogram of pair distances divided by the
    number of pairs expected in the same shell of an ideal gas at the same number density.
    Coordination numbers follow Verbeeck 2025 equation 3, the integral of 4 pi rho g(r) r^2
    up to the first minimum, with rho the number density of the NEIGHBOUR species.
    """
    frames = frames or slice(None)
    pos_all = traj.positions[frames]
    box_all = traj.box[frames]
    species = traj.species_list()
    if pairs is None:
        pairs = [(a, b) for i, a in enumerate(species) for b in species[i:]]

    if rmax > box_all[:, 0].min() / 2.0:
        rmax_new = box_all[:, 0].min() / 2.0
        warn = [f"rmax {rmax:.1f} A exceeds half the smallest box edge, cut to "
                f"{rmax_new:.2f} A. Beyond L/2 the minimum image convention double "
                f"counts and g(r) is meaningless."]
        rmax = rmax_new
    else:
        warn = []

    edges = np.linspace(0.0, rmax, nbins + 1)
    r_mid = 0.5 * (edges[1:] + edges[:-1])
    shell_vol = 4.0 / 3.0 * math.pi * (edges[1:] ** 3 - edges[:-1] ** 3)

    idx = {s: traj.index_of(s) for s in species}
    hist: dict[str, np.ndarray] = {f"{a}-{b}": np.zeros(nbins) for a, b in pairs}
    nframes = len(pos_all)

    for f in range(nframes):
        pos, L = pos_all[f], box_all[f]
        for a, b in pairs:
            ia, ib = idx[a], idx[b]
            if len(ia) == 0 or len(ib) == 0:
                continue
            pa, pb = pos[ia], pos[ib]
            # chunk over a to keep memory bounded for 5000 atoms
            chunk = max(1, int(2.0e6 / max(1, len(ib))))
            h = np.zeros(nbins)
            for s0 in range(0, len(pa), chunk):
                d = pa[s0:s0 + chunk, None, :] - pb[None, :, :]
                d -= L * np.round(d / L)
                dist = np.sqrt((d * d).sum(axis=-1))
                if a == b:
                    n_in = dist.shape[0]
                    ii = np.arange(s0, s0 + n_in)
                    dist[np.arange(n_in), ii] = np.inf   # drop self pairs
                h += np.histogram(dist, bins=edges)[0]
            hist[f"{a}-{b}"] += h

    g: dict[str, np.ndarray] = {}
    cn: dict[str, np.ndarray] = {}
    vol = float(np.prod(box_all, axis=1).mean())
    for a, b in pairs:
        key = f"{a}-{b}"
        na, nb = len(idx[a]), len(idx[b])
        if na == 0 or nb == 0:
            continue
        rho_b = nb / vol
        # pairs counted: for a != b every a-b pair once; for a == b every pair twice
        n_centres = na
        norm = n_centres * nframes * shell_vol * rho_b
        with np.errstate(divide="ignore", invalid="ignore"):
            g[key] = np.where(norm > 0, hist[key] / norm, 0.0)
        cn[key] = np.cumsum(g[key] * shell_vol * rho_b)

    peaks: dict[str, float] = {}
    mins: dict[str, float] = {}
    cn_min: dict[str, float] = {}
    for key, gv in g.items():
        # First maximum: skip the empty core, then find the first local max above 1.
        # The search MUST begin at the first non-empty bin, not one past it. On a sharply
        # peaked g(r), a crystal or a very cold melt, the first non-empty bin IS the first
        # peak, and starting one bin later silently returns the SECOND coordination shell.
        # That bug was live until the simple-cubic-lattice test caught it: the tool
        # reported 4.245 A and a coordination number of 18 for a lattice whose answers are
        # 3.0 A and 6.
        start = int(np.argmax(gv > 0.0))
        i_max = None
        for i in range(max(start, 1), nbins - 1):
            if gv[i] > 1.0 and gv[i] >= gv[i - 1] and gv[i] >= gv[i + 1]:
                i_max = i
                break
        if i_max is None:
            i_max = int(np.argmax(gv))
        peaks[key] = float(r_mid[i_max])
        i_min = None
        for i in range(i_max + 1, nbins - 1):
            if gv[i] <= gv[i - 1] and gv[i] <= gv[i + 1]:
                i_min = i
                break
        if i_min is None:
            i_min = nbins - 1
        mins[key] = float(r_mid[i_min])
        cn_min[key] = float(cn[key][i_min])

    # Physical sanity. No pair of atoms sits closer than about 1.2 Angstrom in an oxide,
    # so a first peak below 1 Angstrom means the coordinates and the box do not belong
    # together, usually because the trajectory was written with a different box or the
    # atoms have wandered outside it and the minimum image fold has turned the separations
    # into noise. Without this the tool reported a 0.038 Angstrom first peak and zero
    # coordination numbers with no complaint, which cost a paragraph of explanation in the
    # manual. Found 29 September 2026.
    absurd = sorted(k for k, v in peaks.items() if v < 1.0)
    if absurd:
        warn.append(
            f"the first RDF peak is below 1 Angstrom for {absurd}. That is not physically "
            f"possible. Either the coordinates do not belong to this box, or the atoms "
            f"have drifted outside it and the minimum image convention is folding them "
            f"into noise. Do not use these coordination numbers.")
    empty = sorted(k for k, v in cn_min.items() if v < 0.05)
    if empty:
        warn.append(
            f"the coordination number is essentially zero for {empty}. In a condensed "
            f"melt that cannot happen, so the trajectory is probably not condensed, or the "
            f"box in the dump is wrong.")

    return RDFResult(r=r_mid, g=g, cn=cn, first_peak=peaks, first_min=mins,
                     cn_at_first_min=cn_min, n_frames_used=nframes, warnings=warn)


# ================================================================= MSD, diffusion
@dataclass
class DiffusionResult:
    species: str
    t_ps: np.ndarray
    msd_A2: np.ndarray
    fit_lo_ps: float
    fit_hi_ps: float
    slope_A2_per_ps: float
    intercept: float
    r_squared: float
    loglog_slope: float
    D_m2_s: float
    D_cm2_s: float
    diffusive: bool
    verdict: str


def compute_msd(traj: Trajectory, species: str,
                max_lag_fraction: float = 0.5,
                unwrap_if_needed: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Multiple time origin MSD for one species.

    Uses every available time origin, which is what makes a 50 ns run give a smooth MSD
    instead of a noisy one. If the dump was wrapped, coordinates are unwrapped first by
    following each atom and undoing any jump larger than half the box. That is correct as
    long as no atom really moves more than L/2 between snapshots, which at 1 ps sampling
    in a melt it does not.
    """
    idx = traj.index_of(species)
    if len(idx) == 0:
        raise ValueError(f"no atoms of species {species!r} in the trajectory")
    pos = traj.positions[:, idx, :].astype(float)

    if not traj.unwrapped and unwrap_if_needed:
        L = traj.box[:, None, :]
        d = np.diff(pos, axis=0)
        d -= L[1:] * np.round(d / L[1:])
        pos = np.concatenate([pos[:1], pos[:1] + np.cumsum(d, axis=0)], axis=0)

    nf = pos.shape[0]
    max_lag = max(1, int(nf * max_lag_fraction))
    msd = _msd_fft(pos)[:max_lag]
    t = traj.times_ps[:max_lag] - traj.times_ps[0]
    return t, msd


def fit_diffusion(t_ps: np.ndarray, msd_A2: np.ndarray, species: str,
                  fit_lo_frac: float = 0.3, fit_hi_frac: float = 0.8,
                  loglog_tolerance: float = 0.15) -> DiffusionResult:
    """Einstein relation, with the diffusive-regime test that makes it honest.

    D = slope / 6, Verbeeck 2025 equation 4, with MSD in Angstrom^2 and t in ps, so
    1 A^2/ps = 1e-8 m^2/s.

    The test: over the fitting window, fit log(MSD) against log(t). In the diffusive
    regime that slope is 1. In the cage or ballistic regime it is not, and dividing the
    linear slope by 6 then gives a number with the units of a diffusion coefficient that
    is not one. This is the exact trap that produced a plateaued MSD in his own earlier
    practice run, so the tool now refuses instead of reporting.
    """
    n = len(t_ps)
    lo = max(1, int(n * fit_lo_frac))
    hi = max(lo + 3, int(n * fit_hi_frac))
    hi = min(hi, n)
    tt, mm = t_ps[lo:hi], msd_A2[lo:hi]
    if len(tt) < 4:
        raise ValueError("fitting window has fewer than 4 points, the run is too short")

    A = np.vstack([tt, np.ones_like(tt)]).T
    (slope, intercept), residuals, *_ = np.linalg.lstsq(A, mm, rcond=None)
    ss_res = float(((mm - (slope * tt + intercept)) ** 2).sum())
    ss_tot = float(((mm - mm.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    ok = (tt > 0) & (mm > 0)
    if ok.sum() >= 4:
        ll_slope = float(np.polyfit(np.log(tt[ok]), np.log(mm[ok]), 1)[0])
    else:
        ll_slope = float("nan")

    D_A2_ps = slope / 6.0
    D_m2_s = D_A2_ps * 1.0e-8            # A^2/ps -> m^2/s
    diffusive = (abs(ll_slope - 1.0) <= loglog_tolerance) and slope > 0

    if not np.isfinite(ll_slope):
        verdict = "cannot test the regime, the MSD has non-positive values"
    elif slope <= 0:
        verdict = ("the MSD slope is not positive. Nothing is diffusing. Either the run "
                   "is far too short or the melt has frozen at this temperature.")
    elif ll_slope < 1.0 - loglog_tolerance:
        verdict = (f"log-log slope is {ll_slope:.2f}, not 1. The atoms are still rattling "
                   f"in their cages over this window, so this is NOT a diffusion "
                   f"coefficient. Run longer, or fit later in the trajectory. This is the "
                   f"single commonest way a wrong D gets published.")
    elif ll_slope > 1.0 + loglog_tolerance:
        verdict = (f"log-log slope is {ll_slope:.2f}, above 1. Either the fit window is "
                   f"still in the ballistic regime at very short times, or the MSD is "
                   f"dominated by a drifting centre of mass. Check that the dump used "
                   f"com yes or subtract the drift.")
    elif r2 < 0.995:
        verdict = (f"the regime is diffusive, log-log slope {ll_slope:.2f}, but the linear "
                   f"fit R squared is only {r2:.4f}. The statistics are poor. Report an "
                   f"error bar from block averages, do not quote a single figure.")
    else:
        verdict = (f"diffusive, log-log slope {ll_slope:.2f}, R squared {r2:.5f}. "
                   f"The Einstein relation is valid over this window.")

    return DiffusionResult(
        species=species, t_ps=t_ps, msd_A2=msd_A2,
        fit_lo_ps=float(tt[0]), fit_hi_ps=float(tt[-1]),
        slope_A2_per_ps=float(slope), intercept=float(intercept), r_squared=float(r2),
        loglog_slope=ll_slope, D_m2_s=float(D_m2_s), D_cm2_s=float(D_m2_s * 1.0e4),
        diffusive=bool(diffusive), verdict=verdict)


# =========================================================== oxygen types and Qn
@dataclass
class SpeciationResult:
    cutoffs: dict[str, float]
    bridging: int
    non_bridging: int
    free: int
    tricluster: int
    frac_bridging: float
    frac_non_bridging: float
    frac_free: float
    frac_tricluster: float
    nbo_per_T: float
    qn_distribution: dict[str, dict[int, float]]    # "Si" -> {0: frac, ... 4: frac}
    mean_qn: dict[str, float]
    coordination: dict[str, float]
    n_frames_used: int
    notes: list[str] = field(default_factory=list)


def compute_speciation(traj: Trajectory, cutoffs: dict[str, float] | None = None,
                       frames: slice | None = None,
                       formers: tuple[str, ...] = NETWORK_FORMERS) -> SpeciationResult:
    """Classify every oxygen, then count Qn for every network former.

    Oxygen classes, the standard slag chemistry definitions:
      free O        bonded to no network former, only to modifier cations
      non-bridging  bonded to exactly one network former
      bridging      bonded to exactly two
      triclustered  bonded to three or more. Real in high alumina melts, not an artefact,
                    and Verbeeck's group counts them, so the tool does too.

    Qn for a network former is its number of BRIDGING oxygens, which is the definition
    that matters, not simply its coordination number.

    The cutoff is the whole argument. The right cutoff is the first minimum of the
    relevant partial RDF, so by default the tool measures it from the trajectory rather
    than taking a literature value. That is one fewer number to be wrong about.
    """
    frames = frames or slice(None)
    species = traj.species_list()
    present_formers = [f for f in formers if f in species]
    notes: list[str] = []

    if cutoffs is None:
        rdf = compute_rdf(traj, rmax=6.0, nbins=300, frames=frames,
                          pairs=[(f, "O") for f in present_formers])
        cutoffs = {}
        for f in present_formers:
            key = f"{f}-O" if f"{f}-O" in rdf.first_min else f"O-{f}"
            cutoffs[f] = rdf.first_min.get(key, 2.4)
        notes.append("cutoffs taken from the first minimum of each former-oxygen partial "
                     "RDF measured on this trajectory: "
                     + ", ".join(f"{k} {v:.2f} A" for k, v in cutoffs.items()))

    pos_all = traj.positions[frames]
    box_all = traj.box[frames]
    o_idx = traj.index_of("O")
    former_idx = {f: traj.index_of(f) for f in present_formers}

    n_b = n_nb = n_free = n_tri = 0
    qn_counts: dict[str, dict[int, int]] = {f: {} for f in present_formers}
    coord_sum: dict[str, float] = {f: 0.0 for f in present_formers}
    coord_n: dict[str, int] = {f: 0 for f in present_formers}
    nframes = len(pos_all)

    for fi in range(nframes):
        pos, L = pos_all[fi], box_all[fi]
        po = pos[o_idx]
        # for each oxygen, how many formers of each type are inside that type's cutoff
        n_formers_on_O = np.zeros(len(o_idx), dtype=int)
        # and, per former atom, its list of bonded oxygens
        bonded: dict[str, list[np.ndarray]] = {}
        for f in present_formers:
            pf = pos[former_idx[f]]
            if len(pf) == 0:
                continue
            rc = cutoffs[f]
            adj = np.zeros((len(pf), len(o_idx)), dtype=bool)
            chunk = max(1, int(2.0e6 / max(1, len(o_idx))))
            for s0 in range(0, len(pf), chunk):
                d = pf[s0:s0 + chunk, None, :] - po[None, :, :]
                d -= L * np.round(d / L)
                dist = np.sqrt((d * d).sum(axis=-1))
                adj[s0:s0 + chunk] = dist < rc
            n_formers_on_O += adj.sum(axis=0)
            bonded[f] = adj
            coord_sum[f] += float(adj.sum(axis=1).mean())
            coord_n[f] += 1

        n_free += int((n_formers_on_O == 0).sum())
        n_nb += int((n_formers_on_O == 1).sum())
        n_b += int((n_formers_on_O == 2).sum())
        n_tri += int((n_formers_on_O >= 3).sum())

        is_bridging = (n_formers_on_O == 2)
        for f in present_formers:
            adj = bonded.get(f)
            if adj is None:
                continue
            q = (adj & is_bridging[None, :]).sum(axis=1)
            for val in q:
                qn_counts[f][int(val)] = qn_counts[f].get(int(val), 0) + 1

    total_o = n_free + n_nb + n_b + n_tri
    if total_o == 0:
        raise ValueError("no oxygen atoms were classified, check the cutoffs")

    n_T_atoms = sum(len(former_idx[f]) for f in present_formers)
    nbo_per_T = (n_nb / nframes) / n_T_atoms if n_T_atoms else float("nan")

    qn_dist: dict[str, dict[int, float]] = {}
    mean_qn: dict[str, float] = {}
    for f, counts in qn_counts.items():
        tot = sum(counts.values())
        if tot == 0:
            continue
        qn_dist[f] = {k: v / tot for k, v in sorted(counts.items())}
        mean_qn[f] = sum(k * v for k, v in counts.items()) / tot

    coord = {f: coord_sum[f] / coord_n[f] for f in present_formers if coord_n[f]}

    if n_tri / total_o > 0.05:
        notes.append(
            f"{100*n_tri/total_o:.1f} percent of oxygens are triclustered. That is high. "
            "It is real in alumina rich melts, but it can also mean the cutoff is too "
            "generous. Check it against the RDF first minimum.")
    if n_free / total_o > 0.10:
        notes.append(
            f"{100*n_free/total_o:.1f} percent free oxygen. In a silicate melt free O is "
            "usually a few percent at most, so a large number suggests the cutoff is too "
            "tight or the composition is very basic.")

    return SpeciationResult(
        cutoffs=cutoffs, bridging=n_b, non_bridging=n_nb, free=n_free, tricluster=n_tri,
        frac_bridging=n_b / total_o, frac_non_bridging=n_nb / total_o,
        frac_free=n_free / total_o, frac_tricluster=n_tri / total_o,
        nbo_per_T=nbo_per_T, qn_distribution=qn_dist, mean_qn=mean_qn,
        coordination=coord, n_frames_used=nframes, notes=notes)


@dataclass
class DiffusionWithError:
    species: str
    D_m2_s: float                 # from the whole trajectory
    D_cm2_s: float                # the same value in cm^2/s, the unit most papers use
    D_mean_m2_s: float            # mean over independent blocks
    D_stderr_m2_s: float          # standard error of that mean
    D_blocks_m2_s: list[float]
    n_blocks: int
    relative_stderr_percent: float
    whole: DiffusionResult
    verdict: str


def diffusion_with_error(traj: Trajectory, species: str, n_blocks: int = 5,
                         fit_lo_frac: float = 0.3, fit_hi_frac: float = 0.8
                         ) -> DiffusionWithError:
    """D from the whole trajectory, plus an honest error bar from independent blocks.

    A single fitted slope has no error bar. The residual scatter of the fit is useless as
    one, because consecutive MSD points are strongly correlated: every point shares time
    origins with its neighbours, so a beautiful R squared can sit on top of a badly
    determined slope. The only cheap honest estimate is to cut the trajectory into blocks
    that do not share origins, fit each, and take the standard error of the block values.

    Verbeeck's own released scripts do something equivalent and quote the standard
    deviation over sliding fit windows. Blocks are stricter, because sliding windows on
    one trajectory are not independent.
    """
    whole_t, whole_msd = compute_msd(traj, species)
    whole = fit_diffusion(whole_t, whole_msd, species, fit_lo_frac, fit_hi_frac)

    nf = traj.n_frames
    per = nf // n_blocks
    blocks: list[float] = []
    if per >= 12:
        for b in range(n_blocks):
            sl = slice(b * per, (b + 1) * per)
            sub = Trajectory(positions=traj.positions[sl], box=traj.box[sl],
                             types=traj.types, elements=traj.elements,
                             times_ps=traj.times_ps[sl] - traj.times_ps[sl][0],
                             unwrapped=traj.unwrapped)
            try:
                bt, bm = compute_msd(sub, species)
                blocks.append(fit_diffusion(bt, bm, species,
                                            fit_lo_frac, fit_hi_frac).D_m2_s)
            except (ValueError, IndexError):
                continue

    if len(blocks) >= 2:
        arr = np.array(blocks)
        mean = float(arr.mean())
        stderr = float(arr.std(ddof=1) / math.sqrt(len(arr)))
        rel = 100.0 * stderr / abs(mean) if mean else float("inf")
        if rel > 25:
            verdict = (f"the error bar is {rel:.0f} percent of the value. This D is not "
                       f"measured, it is estimated. Run longer or average independent "
                       f"trajectories before quoting it.")
        elif rel > 10:
            verdict = (f"error bar {rel:.0f} percent. Usable, but quote it as "
                       f"{mean:.2e} plus or minus {stderr:.1e} m2/s, never as a bare "
                       f"number.")
        else:
            verdict = (f"error bar {rel:.1f} percent over {len(arr)} independent blocks. "
                       f"Well determined.")
    else:
        mean, stderr, rel = whole.D_m2_s, float("nan"), float("nan")
        verdict = ("too few frames to split into independent blocks, so there is no error "
                   "bar. Treat the value as indicative only.")

    return DiffusionWithError(
        species=species, D_m2_s=whole.D_m2_s, D_cm2_s=whole.D_m2_s * 1.0e4,
        D_mean_m2_s=mean,
        D_stderr_m2_s=stderr, D_blocks_m2_s=blocks, n_blocks=len(blocks),
        relative_stderr_percent=rel, whole=whole, verdict=verdict)


# ======================================================== electrical conductivity
@dataclass
class ConductivityResult:
    sigma_NE_S_per_m: float
    sigma_E_S_per_m: float | None
    haven_ratio: float | None
    temperature_K: float
    volume_A3: float
    per_species_D: dict[str, float]
    charges: dict[str, float]
    note: str


def conductivity_nernst_einstein(counts: dict[str, int], charges: dict[str, float],
                                 D_m2_s: dict[str, float], volume_A3: float,
                                 T: float) -> float:
    """Verbeeck 2025 equation 7. Assumes ions move independently."""
    V = volume_A3 * 1.0e-30
    total = 0.0
    for s, n in counts.items():
        if s not in D_m2_s or s not in charges:
            continue
        total += n * charges[s] ** 2 * D_m2_s[s]
    return E_CHARGE ** 2 / (KB * T * V) * total


def conductivity_exact_einstein(traj: Trajectory, charges: dict[str, float], T: float,
                                fit_lo_frac: float = 0.3, fit_hi_frac: float = 0.8,
                                n_blocks: int = 5) -> tuple[float, dict]:
    """Verbeeck 2025 equation 8, the exact Einstein relation.

    Instead of summing independent MSDs, this tracks the mean squared displacement of the
    total CHARGE, sum_i z_i r_i(t). Cross correlations between species are then included
    automatically, because they are inside the square. The result is usually smaller than
    Nernst-Einstein, and the ratio is the Haven ratio, which is the physics the
    approximation discards.
    """
    pos = traj.positions.astype(float)
    if not traj.unwrapped:
        L = traj.box[:, None, :]
        d = np.diff(pos, axis=0)
        d -= L[1:] * np.round(d / L[1:])
        pos = np.concatenate([pos[:1], pos[:1] + np.cumsum(d, axis=0)], axis=0)

    z = np.array([charges[e] for e in traj.elements])
    # collective charge displacement per frame, shape (nf, 3)
    Q = (pos * z[None, :, None]).sum(axis=1)

    nf = len(Q)
    max_lag = max(2, nf // 2)
    mscd = _msd_fft(Q)[:max_lag]
    t = traj.times_ps[:max_lag] - traj.times_ps[0]

    lo = max(1, int(max_lag * fit_lo_frac))
    hi = min(max_lag, max(lo + 3, int(max_lag * fit_hi_frac)))
    A = np.vstack([t[lo:hi], np.ones(hi - lo)]).T
    slope = float(np.linalg.lstsq(A, mscd[lo:hi], rcond=None)[0][0])   # A^2/ps

    V = traj.volume_A3() * 1.0e-30
    slope_m2_s = slope * 1.0e-8
    sigma = E_CHARGE ** 2 / (6.0 * KB * T * V) * slope_m2_s

    # Block uncertainty. This matters more here than anywhere else in the module. The
    # collective charge displacement is ONE coordinate, not one per atom, so however many
    # thousand ions there are, a single trajectory provides only a handful of independent
    # samples of it. Measured on 120 synthetic replicas of a 600 ion box over 400 ps, the
    # per-trajectory scatter was 79 percent while the mean was correct to 5 percent. A
    # single-trajectory exact Einstein conductivity is therefore close to meaningless, and
    # the tool says so rather than printing four significant figures.
    blocks: list[float] = []
    per = nf // n_blocks
    if per >= 12:
        for b in range(n_blocks):
            seg = Q[b * per:(b + 1) * per]
            ml = max(2, len(seg) // 2)
            m2 = _msd_fft(seg)[:ml]
            tt = np.arange(ml) * (traj.times_ps[1] - traj.times_ps[0])
            l2 = max(1, int(ml * fit_lo_frac))
            h2 = min(ml, max(l2 + 3, int(ml * fit_hi_frac)))
            Ab = np.vstack([tt[l2:h2], np.ones(h2 - l2)]).T
            sl = float(np.linalg.lstsq(Ab, m2[l2:h2], rcond=None)[0][0])
            blocks.append(E_CHARGE ** 2 / (6.0 * KB * T * V) * sl * 1.0e-8)

    if len(blocks) >= 2:
        arr = np.array(blocks)
        stderr = float(arr.std(ddof=1) / math.sqrt(len(arr)))
        rel = 100.0 * stderr / abs(arr.mean()) if arr.mean() else float("inf")
        if rel > 30:
            verdict = (f"the block to block scatter gives an error of {rel:.0f} percent. "
                       f"This is NOT a measured conductivity. The exact Einstein relation "
                       f"needs either a very long trajectory or many independent runs, "
                       f"because it depends on one collective coordinate. Quote the "
                       f"Nernst-Einstein value and state that the exact value could not "
                       f"be resolved.")
        else:
            verdict = f"block error {rel:.0f} percent over {len(arr)} blocks"
    else:
        stderr, rel = float("nan"), float("nan")
        verdict = "too few frames for a block error estimate"

    return sigma, {"t_ps": t, "mscd_A2": mscd, "slope_A2_per_ps": slope,
                   "fit_window_ps": (float(t[lo]), float(t[hi - 1])),
                   "sigma_blocks_S_per_m": blocks, "sigma_stderr_S_per_m": stderr,
                   "relative_stderr_percent": rel, "verdict": verdict}


# ===================================================================== viscosity
@dataclass
class ViscosityResult:
    method: str
    eta_Pa_s: float | None
    running_t_ps: np.ndarray
    running_eta_Pa_s: np.ndarray
    plateau_lo_ps: float
    plateau_hi_ps: float
    plateau_scatter_percent: float | None
    temperature_K: float
    volume_A3: float
    verdict: str


def _autocorr_fft(x: np.ndarray) -> np.ndarray:
    """Unbiased autocorrelation by FFT. Much faster than the direct sum and identical."""
    n = len(x)
    f = np.fft.rfft(x - 0.0, n=2 * n)
    ac = np.fft.irfft(f * np.conjugate(f), n=2 * n)[:n].real
    return ac / (n - np.arange(n))


def viscosity_green_kubo(stress_file: Path, T: float, volume_A3: float,
                         timestep_fs: float = 2.0, sample_every: int = 5,
                         plateau_lo_ps: float = 2.0, plateau_hi_ps: float = 10.0
                         ) -> ViscosityResult:
    """Green-Kubo shear viscosity from the off diagonal pressure tensor.

    eta = V/(kB T) * integral_0^inf <P_ab(0) P_ab(t)> dt, averaged over the three
    independent off diagonal components. LAMMPS writes pressure in atm for units real,
    so the conversion to Pa is explicit.

    The result is deliberately returned as a running integral. Green-Kubo viscosity of a
    viscous melt is notorious: the autocorrelation has a long slowly decaying tail, the
    integral wanders, and a single number is meaningless without showing that it settled.
    The tool measures the scatter across the chosen plateau window and refuses to call it
    converged if the scatter is over 10 percent.
    """
    stress_file = Path(stress_file)
    rows = []
    for line in stress_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                rows.append([float(v) for v in parts[:7]])
            except ValueError:
                continue
    if len(rows) < 100:
        raise ValueError(
            f"only {len(rows)} stress samples in {stress_file.name}. Green-Kubo needs "
            "tens of thousands. The stress must be sampled every few steps, not every "
            "picosecond: the stress autocorrelation decays in well under 1 ps and a 1 ps "
            "sampling interval aliases it away completely.")
    arr = np.array(rows)
    dt_ps = timestep_fs * sample_every / 1000.0
    # columns: step pxy pxz pyz [pxx pyy pzz]
    comps = [arr[:, 1], arr[:, 2], arr[:, 3]]

    V = volume_A3 * 1.0e-30
    pref = V / (KB * T)
    acfs = []
    for c in comps:
        c_pa = (c - c.mean()) * ATM_TO_PA
        acfs.append(_autocorr_fft(c_pa))
    acf = np.mean(acfs, axis=0)                      # Pa^2
    running = np.cumsum(acf) * dt_ps * 1.0e-12       # Pa^2 s
    eta_run = pref * running                          # Pa s
    t = np.arange(len(eta_run)) * dt_ps

    m = (t >= plateau_lo_ps) & (t <= plateau_hi_ps)
    if m.sum() < 5:
        return ViscosityResult(
            method="Green-Kubo", eta_Pa_s=None, running_t_ps=t, running_eta_Pa_s=eta_run,
            plateau_lo_ps=plateau_lo_ps, plateau_hi_ps=plateau_hi_ps,
            plateau_scatter_percent=None, temperature_K=T, volume_A3=volume_A3,
            verdict=("the requested plateau window contains fewer than 5 points, so no "
                     "viscosity is reported. Widen the window or sample more often."))
    seg = eta_run[m]
    mean = float(seg.mean())
    scatter = float(100.0 * seg.std() / abs(mean)) if mean else float("inf")
    if scatter > 10.0:
        verdict = (f"NOT converged. Across {plateau_lo_ps} to {plateau_hi_ps} ps the "
                   f"running integral scatters by {scatter:.0f} percent. Green-Kubo has "
                   f"not plateaued. Do not quote this number. Either run much longer, "
                   f"average over many independent trajectories, or use the "
                   f"Einstein-Helfand form which integrates the noise instead of "
                   f"differentiating it.")
    else:
        verdict = (f"plateau found, scatter {scatter:.1f} percent across "
                   f"{plateau_lo_ps} to {plateau_hi_ps} ps. Still report it with an error "
                   f"bar from independent runs, not from this one window.")
    return ViscosityResult(
        method="Green-Kubo", eta_Pa_s=mean, running_t_ps=t, running_eta_Pa_s=eta_run,
        plateau_lo_ps=plateau_lo_ps, plateau_hi_ps=plateau_hi_ps,
        plateau_scatter_percent=scatter, temperature_K=T, volume_A3=volume_A3,
        verdict=verdict)


def viscosity_einstein_helfand(stress_file: Path, T: float, volume_A3: float,
                               timestep_fs: float = 2.0, sample_every: int = 5,
                               fit_lo_frac: float = 0.002, fit_hi_frac: float = 0.02
                               ) -> ViscosityResult:
    """Einstein-Helfand shear viscosity, the integrated form of Green-Kubo.

    Define G_ab(t) = integral_0^t P_ab dt'. Then
      eta = V/(2 kB T) * d<[G_ab(t) - G_ab(0)]^2>/dt
    This is to Green-Kubo what the Einstein relation is to the velocity autocorrelation.
    It behaves better in practice because integrating suppresses noise while
    differentiating amplifies it, which is why it is worth computing both and comparing.

    THE FIT WINDOW MATTERS AND THE DEFAULT IS EARLY ON PURPOSE. The linear regime begins
    once t is much longer than the stress relaxation time, which in a melt is well under a
    picosecond, and it ENDS in practice long before the series does, because at a lag equal
    to a third of the trajectory there are hardly any independent time origins left and the
    Helfand displacement is itself a random walk. Fitting at 20 to 60 percent of the series,
    which is the obvious choice and was this function's first default, gave an answer 5
    times too small on a signal whose exact viscosity is known. The default is now 0.2 to 2
    percent of the series. Always look at the reported log-log slope: if it is not close to
    1, move the window rather than trusting the number.
    """
    stress_file = Path(stress_file)
    rows = []
    for line in stress_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                rows.append([float(v) for v in parts[:7]])
            except ValueError:
                continue
    if len(rows) < 100:
        raise ValueError("too few stress samples for Einstein-Helfand")
    arr = np.array(rows)
    dt_ps = timestep_fs * sample_every / 1000.0
    dt_s = dt_ps * 1.0e-12
    V = volume_A3 * 1.0e-30

    msds = []
    for col in (1, 2, 3):
        p = (arr[:, col] - arr[:, col].mean()) * ATM_TO_PA
        G = (np.cumsum(p) * dt_s)[:, None]          # Pa s, shape (n, 1)
        nf = len(G)
        max_lag = max(2, nf // 2)
        msds.append(_msd_fft(G)[:max_lag])
    msd = np.mean(msds, axis=0)
    t = np.arange(len(msd)) * dt_ps

    lo = max(1, int(len(msd) * fit_lo_frac))
    hi = min(len(msd), max(lo + 3, int(len(msd) * fit_hi_frac)))
    A = np.vstack([t[lo:hi] * 1.0e-12, np.ones(hi - lo)]).T
    slope = float(np.linalg.lstsq(A, msd[lo:hi], rcond=None)[0][0])   # Pa^2 s^2 / s
    eta = V / (2.0 * KB * T) * slope

    # running estimate, as a local derivative, so the plot shows whether it settled
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_run = V / (2.0 * KB * T) * np.gradient(msd, t * 1.0e-12)
    # The log-log slope must be measured ON THE FIT WINDOW, not over the whole series.
    # Measured over everything it is the same number whatever window you choose, so it
    # cannot tell you the window is wrong, which is the one thing it is for.
    tw, mw = t[lo:hi], msd[lo:hi]
    ok = (tw > 0) & np.isfinite(mw) & (mw > 0)
    ll = float(np.polyfit(np.log(tw[ok]), np.log(mw[ok]), 1)[0]) if ok.sum() > 5 \
        else float("nan")
    if not np.isfinite(ll):
        verdict = "cannot test the regime"
    elif abs(ll - 1.0) > 0.2:
        verdict = (f"over the fit window the Helfand mean square displacement grows as "
                   f"t^{ll:.2f}, not t^1, so this is not yet a viscosity. If the exponent "
                   f"is below 1 the window is too LATE, where a single trajectory has too "
                   f"few independent time origins. Move it earlier. If it is above 1 the "
                   f"window is too early, inside the stress relaxation time.")
    else:
        verdict = (f"linear regime confirmed, Helfand MSD grows as t^{ll:.2f}. Compare "
                   f"against the Green-Kubo value. Agreement between the two is the real "
                   f"test, because they fail in different directions.")
    return ViscosityResult(
        method="Einstein-Helfand", eta_Pa_s=float(eta), running_t_ps=t,
        running_eta_Pa_s=eta_run, plateau_lo_ps=float(t[lo]), plateau_hi_ps=float(t[hi-1]),
        plateau_scatter_percent=None, temperature_K=T, volume_A3=volume_A3,
        verdict=verdict)


# ======================================================================== log file
def read_lammps_log(path: Path) -> dict[str, np.ndarray]:
    """Pull every thermo table out of a LAMMPS log and concatenate by column name."""
    path = Path(path)
    cols: dict[str, list[float]] = {}
    header: list[str] | None = None
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if s.startswith("Step ") or s.split()[:1] == ["Step"]:
            header = s.split()
            for h in header:
                cols.setdefault(h, [])
            continue
        if header is None:
            continue
        parts = s.split()
        if len(parts) != len(header):
            if s.startswith("Loop time") or s.startswith("WARNING"):
                header = None
            continue
        try:
            vals = [float(v) for v in parts]
        except ValueError:
            header = None
            continue
        for h, v in zip(header, vals):
            cols[h].append(v)
    return {k: np.array(v) for k, v in cols.items() if v}


def equilibration_report(log: dict[str, np.ndarray],
                         last_fraction: float = 0.5) -> dict:
    """Is the box actually equilibrated? Drift in density and temperature over the tail.

    A density still drifting at the end of equilibration means the production run is done
    at the wrong density, and every transport number scales with density. This is a
    cheap check that almost nobody performs.
    """
    out: dict = {}
    for key in ("Density", "Temp", "Volume", "Press", "TotEng"):
        if key not in log or len(log[key]) < 20:
            continue
        v = log[key]
        n = len(v)
        tail = v[int(n * (1 - last_fraction)):]
        x = np.arange(len(tail), dtype=float)
        slope = float(np.polyfit(x, tail, 1)[0])
        mean = float(tail.mean())
        drift_pct = 100.0 * slope * len(tail) / abs(mean) if mean else float("nan")
        out[key] = {"mean": mean, "std": float(tail.std()),
                    "drift_percent_over_tail": drift_pct}
    verdicts: list[str] = []
    d = out.get("Density", {}).get("drift_percent_over_tail")
    if d is not None and abs(d) > 0.5:
        verdicts.append(
            f"density is still drifting by {d:+.2f} percent over the second half of the "
            f"run. It is not equilibrated. Every diffusion coefficient and conductivity "
            f"measured here is at the wrong density.")
    t = out.get("Temp", {})
    if t and t.get("std", 0) > 0.05 * abs(t.get("mean", 1)):
        verdicts.append(
            f"temperature fluctuates by {t['std']:.0f} K around {t['mean']:.0f} K, which "
            f"is over 5 percent. Check the thermostat damping.")
    out["_verdicts"] = verdicts or ["density and temperature look settled over the tail"]
    return out
