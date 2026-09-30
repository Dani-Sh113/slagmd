"""Load a verified force field and emit LAMMPS pair coefficients.

Two rules are enforced here and they are the reason this file exists.

1. A force field whose parameters are not marked verified=true is REFUSED. The tool will
   not write a deck from remembered or estimated numbers.
2. Every pair the composition needs must be present. A missing pair is a hard error, not
   a silent zero, because LAMMPS will happily run with A=0 and give you a liquid held
   together by nothing but Coulomb.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

KJ_PER_KCAL = 4.184          # exact by definition of the thermochemical calorie
HERE = Path(__file__).resolve().parent


class ForceFieldError(ValueError):
    pass


def _key(a: str, b: str) -> str:
    return "-".join(sorted((a, b)))


@dataclass
class ForceField:
    name: str
    label: str
    form: str                       # "bmh" or "buck"
    charges: dict[str, float]
    pairs: dict[str, dict[str, float]]
    source: str
    cation_cation_coulomb_only: bool
    meta: dict

    # ------------------------------------------------------------------- loading
    @classmethod
    def load(cls, name: str, path: Path | None = None) -> "ForceField":
        path = path or (HERE / "forcefields.json")
        data = json.loads(path.read_text())
        key = name.lower()
        if key not in data or key.startswith("_"):
            available = [k for k in data if not k.startswith("_")]
            raise ForceFieldError(f"no force field {name!r}. Available: {available}")
        d = data[key]
        if not d.get("verified", False):
            raise ForceFieldError(
                f"force field {name!r} is marked verified=false. Refusing to build a deck "
                "from unverified parameters. Add a printed source first.")
        pairs = {_key(*k.split("-")): v for k, v in d["pairs"].items()}
        return cls(
            name=key,
            label=d["label"],
            form=d["form"],
            charges=dict(d["charges"]),
            pairs=pairs,
            source=d["source"],
            cation_cation_coulomb_only=("cation_cation" in d),
            meta=d,
        )

    @staticmethod
    def available(path: Path | None = None) -> list[str]:
        path = path or (HERE / "forcefields.json")
        data = json.loads(path.read_text())
        return [k for k, v in data.items()
                if not k.startswith("_") and v.get("verified", False)]

    @staticmethod
    def protocol_reference(path: Path | None = None) -> dict:
        path = path or (HERE / "forcefields.json")
        return json.loads(path.read_text())["_protocol_reference"]

    # ---------------------------------------------------------------- validation
    def check_self_consistency(self) -> list[str]:
        """Re-run the algebraic identities recorded in the JSON. Returns the failures.

        This is not decoration. It is the check that catches a typo in a transcribed
        table, which is the single most likely way a wrong number enters this tool.
        """
        problems: list[str] = []

        # A force field added to the JSON later gets no identity check unless someone also
        # writes one here. Until 29 Sep 2026 it would simply pass, while its JSON appeared
        # to declare a check, which quietly undercut the whole "verified parameters"
        # guarantee. So: say so, loudly, in both directions.
        implemented = {"matsui", "guillot", "bouhadja"}
        declares = "identity_check" in self.meta
        if self.name not in implemented:
            if declares:
                problems.append(
                    f"{self.name} declares an identity_check in forcefields.json but no "
                    f"implementation exists in check_self_consistency, so nothing is "
                    f"actually being verified. Either implement it or remove the block.")
            else:
                problems.append(
                    f"{self.name} has no algebraic identity check. Its numbers have been "
                    f"transcribed but never cross-checked against the functional form its "
                    f"authors published. Treat the parameters as unconfirmed until one is "
                    f"added, or confirm them against a released input file.")

        # Ground truth check: compare against the authors' own released pair_coeff lines
        # where we have them. This is stronger than an algebraic identity, because it
        # catches a transcription error even in a parameter the algebra cannot constrain.
        released = self.meta.get("released_pair_coeff_kcal")
        if released:
            for key, vals in released.items():
                if key.startswith("_"):
                    continue
                kk = _key(*key.split("-"))
                if kk not in self.pairs:
                    problems.append(f"{self.name}: released file has {key}, we do not")
                    continue
                got = self.pairs[kk]
                A_t, rho_t, sig_t, C_t = vals
                for name, got_v, want_v in (("A", got["A"] / KJ_PER_KCAL, A_t),
                                            ("C", got["C"] / KJ_PER_KCAL, C_t),
                                            ("rho", got["rho"], rho_t),
                                            ("sigma", got.get("sigma", sig_t), sig_t)):
                    if want_v and abs(got_v - want_v) / abs(want_v) > 1e-3:
                        problems.append(
                            f"{self.name} {key} {name}: ours is {got_v:.6g}, the authors' "
                            f"released input file says {want_v:.6g}")

        if self.name == "matsui":
            for k, p in self.pairs.items():
                expect = KJ_PER_KCAL * p["rho"]
                if abs(p["A"] - expect) > 1e-3:
                    problems.append(
                        f"matsui {k}: A={p['A']} but 4.184*rho={expect:.4f}, "
                        "so either A or rho is mistyped")
        if self.name == "guillot":
            try:
                mat = ForceField.load("matsui")
            except ForceFieldError:
                return problems
            for k in ("Al-O", "Ca-O", "O-Si"):
                kk = _key(*k.split("-"))
                if kk not in self.pairs or kk not in mat.pairs:
                    continue
                m, g = mat.pairs[kk], self.pairs[kk]
                if abs(g["rho"] - m["rho"]) > 1e-6:
                    continue          # genuinely refitted, identity does not apply
                expect = m["A"] * math.exp(m["sigma"] / m["rho"])
                # 0.5 percent, not 0.1. The identity is exact in algebra but Matsui's A
                # is printed to four decimals and exp(sigma/rho) is about 4e6, so a
                # rounding in the last printed digit moves the answer by a few tenths of
                # a percent. Measured gap on 29 Sep 2026 was 0.21 to 0.22 percent for all
                # three cation-oxygen pairs, systematic and therefore rounding. A real
                # typo would show as a digit transposition, which is tens of percent or
                # more, and still fails this test.
                if expect > 0 and abs(g["A"] - expect) / expect > 5e-3:
                    problems.append(
                        f"guillot {kk}: A={g['A']:.1f} but A_matsui*exp(sigma/rho)="
                        f"{expect:.1f}, a {100*abs(g['A']-expect)/expect:.2f} percent gap")
        return problems

    def missing_pairs(self, species: list[str]) -> list[str]:
        """Pairs the composition needs and the force field does not have."""
        cations = [s for s in species if s != "O"]
        needed: list[str] = []
        for i, a in enumerate(species):
            for b in species[i:]:
                k = _key(a, b)
                if k in self.pairs:
                    continue
                if self.cation_cation_coulomb_only and a in cations and b in cations:
                    continue          # deliberately absent, Coulomb only
                needed.append(k)
        return needed

    def check_ready(self, species: list[str]) -> None:
        missing_charge = [s for s in species if s not in self.charges]
        if missing_charge:
            raise ForceFieldError(
                f"{self.label} has no effective charge for {missing_charge}. "
                "This force field was not parameterised for those species. "
                "Do not invent a charge, choose a different force field or add a "
                "verified parameter set.")
        missing = self.missing_pairs(species)
        if missing:
            raise ForceFieldError(
                f"{self.label} has no parameters for the pairs {missing}. "
                "LAMMPS would run with zero short-range repulsion for them and the "
                "melt would be physically meaningless. Refusing.")
        bad = self.check_self_consistency()
        if bad:
            raise ForceFieldError(
                "the stored parameters failed their own consistency identity:\n  "
                + "\n  ".join(bad))

    # ------------------------------------------------------------ LAMMPS emission
    def pair_style(self, cutoff: float = 11.0) -> str:
        base = "born/coul/long" if self.form == "bmh" else "buck/coul/long"
        return f"pair_style {base} {cutoff:g}"

    def pair_coeff_lines(self, species: list[str]) -> list[str]:
        """One pair_coeff per type pair, energies converted to kcal/mol for units real.

        LAMMPS born:  E = A exp((sigma - r)/rho) - C/r^6 + D/r^8
        LAMMPS buck:  E = A exp(-r/rho) - C/r^6
        """
        types = {s: i + 1 for i, s in enumerate(species)}
        lines: list[str] = []
        for i, a in enumerate(species):
            for b in species[i:]:
                ta, tb = sorted((types[a], types[b]))
                k = _key(a, b)
                if k not in self.pairs:
                    # only reachable for a deliberately Coulomb-only cation pair
                    if self.form == "bmh":
                        lines.append(f"pair_coeff {ta} {tb} 0.0 1.0 0.0 0.0 0.0"
                                     f"    # {a}-{b} Coulomb only")
                    else:
                        lines.append(f"pair_coeff {ta} {tb} 0.0 1.0 0.0"
                                     f"    # {a}-{b} Coulomb only")
                    continue
                p = self.pairs[k]
                A = p["A"] / KJ_PER_KCAL
                C = p["C"] / KJ_PER_KCAL
                if self.form == "bmh":
                    lines.append(
                        f"pair_coeff {ta} {tb} {A:.8g} {p['rho']:.6g} "
                        f"{p['sigma']:.6g} {C:.8g} 0.0    # {a}-{b}")
                else:
                    lines.append(
                        f"pair_coeff {ta} {tb} {A:.8g} {p['rho']:.6g} {C:.8g}"
                        f"    # {a}-{b}")
        return lines

    def energy_kJ(self, a: str, b: str, r: float) -> float:
        """Short-range part only, kJ/mol, for plotting the curves of Verbeeck Fig 3."""
        k = _key(a, b)
        if k not in self.pairs:
            return 0.0
        p = self.pairs[k]
        if self.form == "bmh":
            rep = p["A"] * math.exp((p["sigma"] - r) / p["rho"])
        else:
            rep = p["A"] * math.exp(-r / p["rho"])
        return rep - p["C"] / r ** 6

    def coulomb_kJ(self, a: str, b: str, r: float) -> float:
        """q_i q_j e^2 / (4 pi eps0 r) in kJ/mol with r in Angstrom."""
        const = 1389.35458  # kJ mol^-1 Angstrom, e^2/(4 pi eps0) in these units
        return const * self.charges[a] * self.charges[b] / r

    def total_kJ(self, a: str, b: str, r: float) -> float:
        return self.energy_kJ(a, b, r) + self.coulomb_kJ(a, b, r)

    def well_depth(self, a: str, b: str,
                   rmin: float = 0.8, rmax: float = 8.0, n: int = 20000
                   ) -> tuple[float, float]:
        """(r at the minimum, depth in kJ/mol) of the total pair energy.

        Used to reproduce the claim in Verbeeck 2025 that Bouhadja's cation-oxygen wells
        are about 66 percent deeper for Ca-O and 57 percent for Al-O and Si-O.
        """
        best_r, best_u = rmin, float("inf")
        for i in range(n):
            r = rmin + (rmax - rmin) * i / (n - 1)
            u = self.total_kJ(a, b, r)
            if u < best_u:
                best_u, best_r = u, r
        return best_r, best_u
