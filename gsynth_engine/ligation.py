"""Ligation arithmetic — how much of each thing goes in the tube.

The calculation is short and gets done wrong constantly, because the
quantity that matters is molar and the quantity you pipette is mass. A
5.4 kb vector and a 150 bp insert at equal mass are at a 36:1 molar ratio
in the vector's favour, which is why the plate comes back empty.

    mass_insert = mass_vector × (length_insert / length_vector) × ratio

Everything else here is the surrounding advice: what ratio to use for the
kind of ends you have, how much DNA a reaction can take, and what a
molar amount means in molecules when the numbers stop feeling real.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.sequence import SequenceError

#: Average mass of one base pair of sodium-salt double-stranded DNA, in
#: g/mol. The figure quoted by every supplier's calculator.
DALTONS_PER_BP = 650.0

AVOGADRO = 6.02214076e23

#: What to use when nothing else is known. Sticky ends ligate efficiently,
#: so a modest excess of insert is enough; blunt ends need more insert to
#: compete with the vector closing on itself.
RECOMMENDED_RATIOS: dict[str, float] = {
    "5'": 3.0,
    "3'": 3.0,
    "blunt": 5.0,
}


def femtomoles(nanograms: float, length_bp: int) -> float:
    """How many fmol a mass of double-stranded DNA of this length is."""
    if length_bp <= 0:
        raise SequenceError("Length must be a positive number of base pairs.")
    # ng / (bp × 650 g/mol) → nmol, ×1e6 → fmol.
    return (nanograms / (length_bp * DALTONS_PER_BP)) * 1e6


def nanograms(femtomoles_wanted: float, length_bp: int) -> float:
    """The inverse: what mass gives this many fmol."""
    if length_bp <= 0:
        raise SequenceError("Length must be a positive number of base pairs.")
    return femtomoles_wanted * length_bp * DALTONS_PER_BP / 1e6


@dataclass
class LigationPlan:
    """One ligation reaction, as amounts to pipette."""

    vector_length: int
    insert_length: int
    vector_ng: float
    insert_ng: float
    ratio: float                #: insert : vector, molar
    vector_fmol: float
    insert_fmol: float
    total_volume_uL: float
    ends: str = "5'"
    warnings: list[str] = field(default_factory=list)

    @property
    def insert_molecules(self) -> float:
        """Insert molecules in the reaction — a count, not a concentration.

        Converted from femtomoles, so a typical reaction is of the order of
        10^9. This is the figure the ratio is really about: at equal *mass* a
        5.4 kb vector is outnumbered by a 150 bp insert thirty-six to one,
        which is why the protocol is worked out in fmol rather than in ng.
        """
        return self.insert_fmol * 1e-15 * AVOGADRO

    @property
    def total_ng(self) -> float:
        return round(self.vector_ng + self.insert_ng, 1)

    def as_rows(self) -> list[dict[str, str]]:
        """The reaction as a table, for a protocol or a spreadsheet."""
        return [
            {
                "Component": "vector",
                "Length": f"{self.vector_length:,} bp",
                "Amount": f"{self.vector_ng:.1f} ng",
                "Moles": f"{self.vector_fmol:.1f} fmol",
            },
            {
                "Component": "insert",
                "Length": f"{self.insert_length:,} bp",
                "Amount": f"{self.insert_ng:.1f} ng",
                "Moles": f"{self.insert_fmol:.1f} fmol",
            },
        ]


def plan_ligation(
    *,
    vector_length: int,
    insert_length: int,
    vector_ng: float = 50.0,
    ratio: float | None = None,
    ends: str = "5'",
    total_volume_uL: float = 20.0,
) -> LigationPlan:
    """Work out how much insert to add to a given amount of vector.

    Args:
        ratio: insert:vector, molar. Left out, the recommendation for the
            kind of ends is used — sticky ends ligate efficiently and need
            only a modest excess, blunt ends need more.

    Raises:
        SequenceError: for lengths or masses that cannot describe a reaction.
    """
    if vector_length <= 0 or insert_length <= 0:
        raise SequenceError("Both lengths must be a positive number of base pairs.")
    if vector_ng <= 0:
        raise SequenceError("The amount of vector must be greater than zero.")

    chosen = RECOMMENDED_RATIOS.get(ends, 3.0) if ratio is None else ratio
    if chosen <= 0:
        raise SequenceError("The molar ratio must be greater than zero.")

    vector_fmol = femtomoles(vector_ng, vector_length)
    insert_fmol = vector_fmol * chosen
    insert_ng = nanograms(insert_fmol, insert_length)

    warnings: list[str] = []
    total = vector_ng + insert_ng
    if total > 500:
        warnings.append(
            f"{total:.0f} ng of DNA in {total_volume_uL:.0f} µL is a lot for one "
            f"ligation. Scale the vector down rather than the ratio."
        )
    if insert_ng < 1:
        warnings.append(
            f"{insert_ng:.2f} ng of insert is below what most people can "
            f"pipette accurately. Dilute less, or use more vector."
        )
    if insert_length > vector_length:
        warnings.append(
            "The insert is longer than the vector. Check the two lengths are "
            "the right way round."
        )
    if ends == "blunt":
        warnings.append(
            "Blunt ends: dephosphorylate the vector, or it will mostly "
            "re-close on itself."
        )

    return LigationPlan(
        vector_length=vector_length,
        insert_length=insert_length,
        vector_ng=round(vector_ng, 2),
        insert_ng=round(insert_ng, 2),
        ratio=chosen,
        vector_fmol=round(vector_fmol, 2),
        insert_fmol=round(insert_fmol, 2),
        total_volume_uL=total_volume_uL,
        ends=ends,
        warnings=warnings,
    )


def ratio_of(
    *, vector_length: int, insert_length: int, vector_ng: float, insert_ng: float,
) -> float:
    """The molar ratio two masses actually represent.

    The other direction, for checking a reaction someone has already set up
    — usually after it failed.
    """
    vector_fmol = femtomoles(vector_ng, vector_length)
    if vector_fmol <= 0:
        raise SequenceError("The amount of vector must be greater than zero.")
    return round(femtomoles(insert_ng, insert_length) / vector_fmol, 2)


def ligation_series(
    *,
    vector_length: int,
    insert_length: int,
    vector_ng: float = 50.0,
    ratios: tuple[float, ...] = (1.0, 3.0, 5.0),
    ends: str = "5'",
) -> list[LigationPlan]:
    """The same reaction at several ratios, which is how it is set up.

    Nobody runs one ligation: they run a small series and pick whichever
    plate gives colonies. Handing back the whole series is what the bench
    actually needs.
    """
    return [
        plan_ligation(
            vector_length=vector_length, insert_length=insert_length,
            vector_ng=vector_ng, ratio=ratio, ends=ends,
        )
        for ratio in ratios
    ]
