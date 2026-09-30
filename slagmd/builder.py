"""Build a random, charge neutral starting configuration and write a LAMMPS data file.

Verbeeck 2025 uses Packmol for this step. Packmol is not always installed, and for a
single-atom-per-molecule ionic melt it is doing nothing Packmol is special at, so this
module does both: it calls Packmol when it is on PATH and you ask for it, and otherwise
uses its own grid accelerated random insertion with a minimum separation criterion.

The starting structure does not need to be good. It needs to be overlap free, because an
overlapping pair in a Buckingham or BMH potential gives a force of order 1e6 kcal/mol/A
and LAMMPS loses the atom on step one. Everything after that is the 4000 K randomisation
stage's job.
"""
from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .composition import ATOMIC_MASS, Composition


class BuildError(RuntimeError):
    pass


# Minimum allowed separation in Angstrom during random insertion. These are deliberately
# a little below the real first-shell distances, because the minimiser will open them up
# and being too strict makes insertion fail at high density.
DEFAULT_MIN_SEP = 1.6
CATION_CATION_MIN_SEP = 2.4      # like charges, keep them further apart


@dataclass
class Structure:
    species: list[str]              # per-atom species symbol
    positions: np.ndarray           # (N, 3) in Angstrom
    box: float                      # cubic edge in Angstrom
    type_order: list[str]           # species -> LAMMPS type is index+1
    warnings: list[str] = field(default_factory=list)

    @property
    def n_atoms(self) -> int:
        return len(self.species)

    def min_distance(self) -> float:
        """Smallest periodic pair distance. O(N^2), only for the self-check."""
        pos, L = self.positions, self.box
        n = len(pos)
        best = float("inf")
        for i in range(n):
            d = pos[i + 1:] - pos[i]
            d -= L * np.round(d / L)
            if len(d):
                best = min(best, float(np.sqrt((d * d).sum(axis=1)).min()))
        return best

    def density_g_cm3(self) -> float:
        mass_g = sum(ATOMIC_MASS[s] for s in self.species) / 6.02214076e23
        return mass_g / (self.box ** 3 * 1.0e-24)


def _cell_list_insert(counts: dict[str, int], box: float, rng: np.random.Generator,
                      min_sep: float, cation_cation_sep: float,
                      max_tries_per_atom: int = 4000) -> tuple[list[str], np.ndarray]:
    """Insert atoms one at a time, rejecting overlaps, using a uniform cell grid.

    Grid cell size is the largest separation we ever test, so a candidate only has to be
    compared against the 27 neighbouring cells. Without this the insertion of 5000 atoms
    is O(N^2) per try and takes minutes instead of a second.
    """
    cutoff = max(min_sep, cation_cation_sep)
    ncell = max(1, int(box / cutoff))
    cell = box / ncell
    grid: dict[tuple[int, int, int], list[int]] = {}

    # place the big, awkward species first: the ones with the largest required spacing
    order: list[str] = sorted(
        [s for s in counts if s != "O"], key=lambda s: -counts[s]) + ["O"]

    species: list[str] = []
    pos_list: list[np.ndarray] = []

    def cell_of(p: np.ndarray) -> tuple[int, int, int]:
        return tuple(int(x) % ncell for x in (p / cell))      # type: ignore[return-value]

    for sp in order:
        for _ in range(counts.get(sp, 0)):
            placed = False
            for _try in range(max_tries_per_atom):
                p = rng.random(3) * box
                ci, cj, ck = cell_of(p)
                ok = True
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        for dk in (-1, 0, 1):
                            key = ((ci + di) % ncell, (cj + dj) % ncell,
                                   (ck + dk) % ncell)
                            for idx in grid.get(key, ()):
                                d = p - pos_list[idx]
                                d -= box * np.round(d / box)
                                r = math.sqrt(float(d @ d))
                                other = species[idx]
                                need = (cation_cation_sep
                                        if (sp != "O" and other != "O")
                                        else min_sep)
                                if r < need:
                                    ok = False
                                    break
                            if not ok:
                                break
                        if not ok:
                            break
                    if not ok:
                        break
                if ok:
                    grid.setdefault(cell_of(p), []).append(len(pos_list))
                    pos_list.append(p)
                    species.append(sp)
                    placed = True
                    break
            if not placed:
                raise BuildError(
                    f"could not place a {sp} atom after {max_tries_per_atom} tries at "
                    f"box {box:.2f} A. The box is too small for the separation criteria. "
                    f"Lower the initial density or reduce min_sep.")
    return species, np.array(pos_list)


def build_random(comp: Composition, density_g_cm3: float = 2.6,
                 seed: int = 12345, min_sep: float = DEFAULT_MIN_SEP,
                 cation_cation_sep: float = CATION_CATION_MIN_SEP,
                 cutoff_A: float = 11.0) -> Structure:
    """Random insertion at a requested initial density.

    The density here is only a starting guess. Verbeeck's protocol runs NPT at 1 atm from
    4000 K down, so the box finds its own volume and the guess is forgotten. It still has
    to be in the right region, because starting far too dense makes insertion impossible
    and starting far too dilute wastes the whole randomisation stage collapsing the box.
    2.6 g/cm3 is a reasonable middle for CaO-Al2O3-SiO2 liquids near 1600 C.
    """
    rng = np.random.default_rng(seed)
    box = comp.box_length_A(density_g_cm3)
    # The hard requirement is L > 2 * cutoff, so that L/2 > cutoff and the minimum image
    # convention is single valued for the real space part of the Ewald sum. Below that the
    # answer is not merely noisy, it is wrong. Between 2 and 2.5 times the cutoff it is
    # valid but the box is small enough that finite size effects on transport are
    # measurable, so that range warns instead of refusing.
    hard = 2.0 * cutoff_A
    soft = 2.5 * cutoff_A
    if box <= hard:
        raise BuildError(
            f"box edge {box:.1f} A is not more than twice the {cutoff_A:g} A Coulomb "
            f"cutoff. With {comp.n_atoms} atoms an atom would interact with its own "
            f"periodic image. Increase target_atoms to at least "
            f"{int(comp.n_atoms * (hard/box)**3) + 1}.")
    warnings: list[str] = []
    if box < soft:
        warnings.append(
            f"box edge {box:.1f} A is valid but under {soft:.0f} A, which is 2.5 times "
            f"the cutoff. Finite size effects on diffusion and viscosity are measurable "
            f"at this size. Verbeeck 2025 used about 5000 atoms for this reason.")
    species, positions = _cell_list_insert(
        comp.counts, box, rng, min_sep, cation_cation_sep)
    return Structure(species=species, positions=positions, box=box,
                     type_order=comp.species, warnings=warnings)


def build_with_packmol(comp: Composition, workdir: Path, density_g_cm3: float = 2.6,
                       seed: int = 12345, tolerance: float = 1.8) -> Structure:
    """Use Packmol if it is installed, which is what Verbeeck 2025 did.

    Falls back with a clear message rather than silently doing something else.
    """
    exe = shutil.which("packmol")
    if exe is None:
        raise BuildError("packmol is not on PATH. Use build_random instead.")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    box = comp.box_length_A(density_g_cm3)

    # one single-atom pdb per species
    for sp in comp.species:
        elem = "".join(ch for ch in sp if ch.isalpha())
        (workdir / f"{sp}.pdb").write_text(
            f"HETATM    1 {elem:<4s}{sp[:3]:>3s} X   1       0.000   0.000   0.000"
            f"  1.00  0.00          {elem:>2s}\nEND\n")

    lines = [f"tolerance {tolerance}", "filetype pdb",
             f"output {workdir / 'packed.pdb'}", f"seed {seed}", ""]
    for sp in comp.species:
        lines += [f"structure {workdir / (sp + '.pdb')}",
                  f"  number {comp.counts[sp]}",
                  f"  inside box 0. 0. 0. {box:.4f} {box:.4f} {box:.4f}",
                  "end structure", ""]
    inp = workdir / "pack.inp"
    inp.write_text("\n".join(lines))

    with inp.open() as fh:
        res = subprocess.run([exe], stdin=fh, capture_output=True, text=True,
                             cwd=workdir, timeout=3600)
    packed = workdir / "packed.pdb"
    if not packed.exists():
        raise BuildError(f"packmol produced no output.\n{res.stdout[-2000:]}")

    species: list[str] = []
    pos: list[list[float]] = []
    for line in packed.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            species.append(resname)
            pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if len(species) != comp.n_atoms:
        raise BuildError(
            f"packmol placed {len(species)} atoms, expected {comp.n_atoms}")
    return Structure(species=species, positions=np.array(pos), box=box,
                     type_order=comp.species)


def write_data_file(struct: Structure, charges: dict[str, float], path: Path,
                    title: str = "slagmd generated melt") -> None:
    """LAMMPS data file, atom_style charge."""
    path = Path(path)
    types = {s: i + 1 for i, s in enumerate(struct.type_order)}
    L = struct.box
    out = [f"{title}", "",
           f"{struct.n_atoms} atoms",
           f"{len(struct.type_order)} atom types", "",
           f"0.0 {L:.6f} xlo xhi",
           f"0.0 {L:.6f} ylo yhi",
           f"0.0 {L:.6f} zlo zhi", "",
           "Masses", ""]
    for s in struct.type_order:
        out.append(f"{types[s]} {ATOMIC_MASS[s]:.4f}   # {s}")
    out += ["", "Atoms # charge", ""]
    for i, (s, p) in enumerate(zip(struct.species, struct.positions), start=1):
        out.append(f"{i} {types[s]} {charges[s]:.6f} "
                   f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    out.append("")
    path.write_text("\n".join(out))


def check_structure(struct: Structure, comp: Composition,
                    charges: dict[str, float], cutoff_A: float = 11.0) -> list[str]:
    """Everything that can be wrong with a starting configuration, before it costs
    a node-hour. Returns a list of problems, empty means good."""
    problems: list[str] = []
    if struct.n_atoms != comp.n_atoms:
        problems.append(
            f"atom count {struct.n_atoms} does not match composition {comp.n_atoms}")
    got: dict[str, int] = {}
    for s in struct.species:
        got[s] = got.get(s, 0) + 1
    for s, n in comp.counts.items():
        if got.get(s, 0) != n:
            problems.append(f"species {s}: structure has {got.get(s,0)}, needs {n}")
    q = sum(charges[s] for s in struct.species)
    if abs(q) > 1e-6:
        problems.append(f"net box charge is {q:+.6f}, not zero. Ewald would be wrong.")
    dmin = struct.min_distance()
    if dmin < 1.0:
        problems.append(
            f"closest pair is {dmin:.3f} A. Below 1 A the repulsive term explodes and "
            f"LAMMPS will lose atoms on the first step.")
    if struct.box <= 2.0 * cutoff_A:
        problems.append(
            f"box edge {struct.box:.1f} A is not more than twice the {cutoff_A:g} A "
            f"cutoff, so the minimum image convention is not single valued")
    return problems
