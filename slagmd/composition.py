"""Turn an oxide composition in mole fractions into an exact, charge neutral atom count.

The whole point of this module is that you cannot just multiply and round. Rounding the
oxide counts independently breaks stoichiometry, and broken stoichiometry breaks charge
neutrality, and a LAMMPS box with net charge and an Ewald sum is silently wrong rather
than loudly wrong. So oxygen is never rounded. It is counted from the cations.

Oxides supported are the ones that matter for a CaO-Al2O3-SiO2 slag plus the four that
extend it, FeO, Fe2O3, MgO and TiO2. Those last two matter because Bellemans' open gaps
are iron redox and titanium.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# oxide -> (cation symbol, cations per formula unit, oxygens per formula unit)
OXIDES: dict[str, tuple[str, int, int]] = {
    "CaO":   ("Ca", 1, 1),
    "Al2O3": ("Al", 2, 3),
    "SiO2":  ("Si", 1, 2),
    "MgO":   ("Mg", 1, 1),
    "FeO":   ("Fe2", 1, 1),    # Fe(II), kept as its own species because its charge differs
    "Fe2O3": ("Fe3", 2, 3),    # Fe(III)
    "TiO2":  ("Ti", 1, 2),
    "Na2O":  ("Na", 2, 1),
    "K2O":   ("K", 2, 1),
}

# g/mol, IUPAC 2021 standard atomic weights
ATOMIC_MASS: dict[str, float] = {
    "O": 15.999, "Si": 28.085, "Al": 26.982, "Ca": 40.078, "Mg": 24.305,
    "Fe2": 55.845, "Fe3": 55.845, "Ti": 47.867, "Na": 22.990, "K": 39.098,
}

# formal oxidation state, used only to check that a partial-charge set is a uniform
# scaling of the formal charges. If it is, any stoichiometric mixture is exactly neutral.
FORMAL_CHARGE: dict[str, float] = {
    "O": -2.0, "Si": 4.0, "Al": 3.0, "Ca": 2.0, "Mg": 2.0,
    "Fe2": 2.0, "Fe3": 3.0, "Ti": 4.0, "Na": 1.0, "K": 1.0,
}

NETWORK_FORMERS = ("Si", "Al", "Fe3", "Ti")   # species that can carry bridging oxygen


class CompositionError(ValueError):
    pass


@dataclass
class Composition:
    """An oxide melt composition resolved to exact integer atom counts."""

    mole_fractions: dict[str, float]
    target_atoms: int = 5000
    counts: dict[str, int] = field(default_factory=dict)
    formula_units: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for ox in self.mole_fractions:
            if ox not in OXIDES:
                raise CompositionError(
                    f"oxide {ox!r} is not in the supported list: {sorted(OXIDES)}")
        total = sum(self.mole_fractions.values())
        if total <= 0:
            raise CompositionError("mole fractions sum to zero")
        if abs(total - 1.0) > 1e-6:
            # normalise rather than refuse, but say so
            self.mole_fractions = {k: v / total for k, v in self.mole_fractions.items()}
            self.normalised_from = total
        else:
            self.normalised_from = None
        if self.target_atoms < 50:
            raise CompositionError("target_atoms below 50 is not a melt, it is a cluster")
        self._resolve()

    # ------------------------------------------------------------------ resolve
    def _resolve(self) -> None:
        """Pick integer formula-unit counts that hit target_atoms as closely as possible."""
        atoms_per_unit = {
            ox: OXIDES[ox][1] + OXIDES[ox][2] for ox in self.mole_fractions}
        # n_units_total such that sum(x_ox * n * atoms_per_unit) == target
        denom = sum(self.mole_fractions[ox] * atoms_per_unit[ox]
                    for ox in self.mole_fractions)
        n_total = self.target_atoms / denom

        # largest-remainder rounding on the formula units, so the mole fractions stay
        # as close as integers allow
        raw = {ox: self.mole_fractions[ox] * n_total for ox in self.mole_fractions}
        units = {ox: int(v) for ox, v in raw.items()}
        shortfall = round(sum(raw.values())) - sum(units.values())
        order = sorted(raw, key=lambda ox: raw[ox] - units[ox], reverse=True)
        for i in range(max(0, shortfall)):
            units[order[i % len(order)]] += 1
        for ox, v in units.items():
            if v < 1 and self.mole_fractions[ox] > 0:
                units[ox] = 1

        self.formula_units = units

        # cations from the formula units, oxygen from the cations. Never rounded.
        counts: dict[str, int] = {}
        n_oxygen = 0
        for ox, n_units in units.items():
            cation, n_cat, n_ox = OXIDES[ox]
            counts[cation] = counts.get(cation, 0) + n_cat * n_units
            n_oxygen += n_ox * n_units
        counts["O"] = n_oxygen
        self.counts = counts

    # ------------------------------------------------------------------ queries
    @property
    def n_atoms(self) -> int:
        return sum(self.counts.values())

    @property
    def species(self) -> list[str]:
        """Cations in a stable order, then oxygen last, which is how the data file is typed."""
        cats = [s for s in self.counts if s != "O"]
        cats.sort(key=lambda s: list(ATOMIC_MASS).index(s))
        return cats + ["O"]

    @property
    def achieved_mole_fractions(self) -> dict[str, float]:
        tot = sum(self.formula_units.values())
        return {ox: n / tot for ox, n in self.formula_units.items()}

    @property
    def max_mole_fraction_error(self) -> float:
        got = self.achieved_mole_fractions
        return max(abs(got[ox] - self.mole_fractions[ox]) for ox in self.mole_fractions)

    def molar_mass_g(self) -> float:
        return sum(ATOMIC_MASS[s] * n for s, n in self.counts.items())

    def net_charge(self, charges: dict[str, float]) -> float:
        """Net box charge for a given partial-charge set. Must be zero."""
        missing = [s for s in self.counts if s not in charges]
        if missing:
            raise CompositionError(
                f"the force field has no charge for {missing}. Do not guess one.")
        return sum(charges[s] * n for s, n in self.counts.items())

    def charges_are_scaled_formal(self, charges: dict[str, float]) -> tuple[bool, float]:
        """True if every charge is the same multiple of its formal charge.

        If it is, neutrality holds for ANY stoichiometric mixture, not just this one,
        which is a stronger statement than checking one box.
        """
        ratios = []
        for s in self.counts:
            if s not in charges or FORMAL_CHARGE[s] == 0:
                return False, 0.0
            ratios.append(charges[s] / FORMAL_CHARGE[s])
        lo, hi = min(ratios), max(ratios)
        return (hi - lo) < 1e-6, sum(ratios) / len(ratios)

    def box_length_A(self, density_g_cm3: float) -> float:
        """Cubic box edge in Angstrom for a requested mass density."""
        if density_g_cm3 <= 0:
            raise CompositionError("density must be positive")
        avogadro = 6.02214076e23
        mass_g = self.molar_mass_g() / avogadro
        volume_cm3 = mass_g / density_g_cm3
        volume_A3 = volume_cm3 * 1.0e24
        return volume_A3 ** (1.0 / 3.0)

    # ------------------------------------------------------- structural bookkeeping
    def nbo_per_tetrahedron_ideal(self) -> float | None:
        """Ideal NBO/T from stoichiometry alone, the textbook count.

        T is every network former. Each mole of a modifier oxide supplies one mole of
        O that cannot bridge. This is the number MD should be compared against, not a
        substitute for measuring it, because aluminium can act as a former or a modifier
        depending on how much charge balancing calcium there is.
        """
        n_T = sum(self.counts.get(s, 0) for s in NETWORK_FORMERS)
        if n_T == 0:
            return None
        n_O = self.counts["O"]
        # bridging-oxygen accounting: every T is 4-coordinate in the ideal network
        return (2.0 * n_O - 4.0 * n_T) / n_T

    def summary(self) -> str:
        lines = []
        req = "  ".join(f"{ox} {v:.4f}" for ox, v in sorted(self.mole_fractions.items()))
        got = self.achieved_mole_fractions
        gots = "  ".join(f"{ox} {got[ox]:.4f}" for ox in sorted(got))
        lines.append(f"requested  {req}")
        lines.append(f"achieved   {gots}   (max error {self.max_mole_fraction_error:.5f})")
        lines.append("formula units  " + "  ".join(
            f"{ox}x{n}" for ox, n in sorted(self.formula_units.items())))
        lines.append("atoms          " + "  ".join(
            f"{s} {self.counts[s]}" for s in self.species) + f"   total {self.n_atoms}")
        nbo = self.nbo_per_tetrahedron_ideal()
        if nbo is not None:
            lines.append(f"ideal NBO/T from stoichiometry  {nbo:.3f}")
        return "\n".join(lines)


def parse_composition_line(line: str) -> dict[str, float]:
    """Parse 'CaO 0.453 Al2O3 0.125 SiO2 0.422' or 'CaO=0.453,Al2O3=0.125,SiO2=0.422'."""
    cleaned = line.replace(",", " ").replace("=", " ").replace(":", " ")
    toks = cleaned.split()
    if len(toks) % 2 != 0:
        raise CompositionError(f"cannot pair up oxide and fraction in {line!r}")
    out: dict[str, float] = {}
    for i in range(0, len(toks), 2):
        ox = toks[i]
        try:
            out[ox] = float(toks[i + 1])
        except ValueError as exc:
            raise CompositionError(f"{toks[i+1]!r} is not a number") from exc
    return out


# The ten benchmark compositions of Verbeeck 2025 Table 5, transcribed verbatim.
# Kept here so the tool can reproduce their study as a self-test.
VERBEECK_2025_TABLE5 = [
    ("S1",  {"CaO": 0.453, "Al2O3": 0.125, "SiO2": 0.422}, (1400, 1600)),
    ("S2",  {"CaO": 0.515, "Al2O3": 0.126, "SiO2": 0.359}, (1500, 1600)),
    ("S3",  {"CaO": 0.438, "Al2O3": 0.125, "SiO2": 0.437}, (1400, 1600)),
    ("S4",  {"CaO": 0.470, "Al2O3": 0.060, "SiO2": 0.470}, (1450, 1600)),
    ("S5",  {"CaO": 0.340, "Al2O3": 0.330, "SiO2": 0.330}, (1550, 1600)),
    ("S6",  {"CaO": 0.670, "Al2O3": 0.230, "SiO2": 0.100}, (1550, 1600)),
    ("S7",  {"CaO": 0.630, "Al2O3": 0.210, "SiO2": 0.160}, (1550, 1600)),
    ("S8",  {"CaO": 0.503, "Al2O3": 0.082, "SiO2": 0.415}, (1400, 1600)),
    ("S9",  {"CaO": 0.350, "Al2O3": 0.100, "SiO2": 0.550}, (1400, 1600)),
    ("S10", {"CaO": 0.250, "Al2O3": 0.100, "SiO2": 0.650}, (1500, 1600)),
]
