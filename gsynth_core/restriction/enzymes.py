"""Restriction-enzyme database — 80 curated enzymes; REBASE fallback via Bio.Restriction."""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Literal

from gsynth_core.sequence.ops import find_motif, reverse_complement


class OverhangType(str, Enum):
    BLUNT = "blunt"
    FIVE_PRIME = "5_prime"
    THREE_PRIME = "3_prime"


@dataclass(frozen=True, slots=True)
class Enzyme:
    name: str
    site: str
    cut_fwd: int
    cut_rev: int
    supplier: str = ""

    @property
    def is_palindromic(self) -> bool:
        s = self.site.upper()
        return s == reverse_complement(s)

    @property
    def overhang_type(self) -> OverhangType:
        if self.cut_fwd == self.cut_rev:
            return OverhangType.BLUNT
        return OverhangType.FIVE_PRIME if self.cut_fwd < self.cut_rev else OverhangType.THREE_PRIME

    @property
    def overhang_length(self) -> int:
        return abs(self.cut_fwd - self.cut_rev)

    @property
    def overhang_sequence(self) -> str:
        if self.overhang_type == OverhangType.BLUNT:
            return ""
        lo, hi = min(self.cut_fwd, self.cut_rev), max(self.cut_fwd, self.cut_rev)
        return self.site[lo:hi]


@dataclass(frozen=True, slots=True)
class RestrictionSite:
    enzyme: Enzyme
    position: int
    strand: Literal["+", "-"]
    cut_position_fwd: int
    cut_position_rev: int


_CURATED_ENZYMES: list[Enzyme] = [
    Enzyme("AatII", "GACGTC", 5, 1), Enzyme("AccI", "GTMKAC", 2, 4),
    Enzyme("AclI", "AACGTT", 2, 4), Enzyme("AflII", "CTTAAG", 1, 5),
    Enzyme("AflIII", "ACRYGT", 1, 5), Enzyme("AgeI", "ACCGGT", 1, 5),
    Enzyme("ApaI", "GGGCCC", 5, 1), Enzyme("ApaLI", "GTGCAC", 1, 5),
    Enzyme("ApoI", "RAATTY", 1, 5), Enzyme("AscI", "GGCGCGCC", 2, 6),
    Enzyme("AseI", "ATTAAT", 2, 4), Enzyme("AvaI", "CYCGRG", 1, 5),
    Enzyme("AvrII", "CCTAGG", 1, 5), Enzyme("BamHI", "GGATCC", 1, 5),
    Enzyme("BbsI", "GAAGAC", 8, 12), Enzyme("BclI", "TGATCA", 1, 5),
    Enzyme("BglII", "AGATCT", 1, 5),
    Enzyme("BmtI", "GCTAGC", 5, 1), Enzyme("BsaI", "GGTCTC", 7, 11),
    Enzyme("BsiWI", "CGTACG", 1, 5), Enzyme("BsmBI", "CGTCTC", 7, 11),
    Enzyme("BspEI", "TCCGGA", 1, 5), Enzyme("BspHI", "TCATGA", 1, 5),
    Enzyme("BsrGI", "TGTACA", 1, 5), Enzyme("BssHII", "GCGCGC", 1, 5),
    Enzyme("BstBI", "TTCGAA", 2, 4), Enzyme("BstEII", "GGTNACC", 1, 6),
    Enzyme("BstZ17I", "GTATAC", 3, 3), Enzyme("ClaI", "ATCGAT", 2, 4),
    Enzyme("DpnI", "GATC", 2, 2), Enzyme("DraI", "TTTAAA", 3, 3),
    Enzyme("DraIII", "CACNNNGTG", 6, 3),
    Enzyme("EagI", "CGGCCG", 1, 5), Enzyme("EcoRI", "GAATTC", 1, 5),
    Enzyme("EcoRV", "GATATC", 3, 3), Enzyme("FseI", "GGCCGGCC", 6, 2),
    Enzyme("FspI", "TGCGCA", 3, 3), Enzyme("HaeIII", "GGCC", 2, 2),
    Enzyme("HindIII", "AAGCTT", 1, 5), Enzyme("HpaI", "GTTAAC", 3, 3),
    Enzyme("KasI", "GGCGCC", 1, 5), Enzyme("KpnI", "GGTACC", 5, 1),
    Enzyme("MfeI", "CAATTG", 1, 5), Enzyme("MluI", "ACGCGT", 1, 5),
    Enzyme("MscI", "TGGCCA", 3, 3), Enzyme("MseI", "TTAA", 1, 3),
    Enzyme("NaeI", "GCCGGC", 3, 3), Enzyme("NcoI", "CCATGG", 1, 5),
    Enzyme("NdeI", "CATATG", 2, 4), Enzyme("NgoMIV", "GCCGGC", 1, 5),
    Enzyme("NheI", "GCTAGC", 1, 5), Enzyme("NotI", "GCGGCCGC", 2, 6),
    Enzyme("NruI", "TCGCGA", 3, 3), Enzyme("NsiI", "ATGCAT", 5, 1),
    Enzyme("PacI", "TTAATTAA", 5, 3), Enzyme("PciI", "ACATGT", 1, 5),
    Enzyme("PmeI", "GTTTAAAC", 4, 4), Enzyme("PmlI", "CACGTG", 3, 3),
    Enzyme("PstI", "CTGCAG", 5, 1), Enzyme("PvuI", "CGATCG", 4, 2),
    Enzyme("PvuII", "CAGCTG", 3, 3), Enzyme("SacI", "GAGCTC", 5, 1),
    Enzyme("SacII", "CCGCGG", 4, 2), Enzyme("SalI", "GTCGAC", 1, 5),
    Enzyme("SapI", "GCTCTTC", 8, 11), Enzyme("SbfI", "CCTGCAGG", 6, 2),
    Enzyme("ScaI", "AGTACT", 3, 3), Enzyme("SfiI", "GGCCNNNNNGGCC", 8, 5),
    Enzyme("SmaI", "CCCGGG", 3, 3),
    Enzyme("SnaBI", "TACGTA", 3, 3), Enzyme("SpeI", "ACTAGT", 1, 5),
    Enzyme("SphI", "GCATGC", 5, 1), Enzyme("SrfI", "GCCCGGGC", 4, 4),
    Enzyme("SspI", "AATATT", 3, 3), Enzyme("StuI", "AGGCCT", 3, 3),
    Enzyme("SwaI", "ATTTAAAT", 4, 4), Enzyme("XbaI", "TCTAGA", 1, 5),
    Enzyme("XhoI", "CTCGAG", 1, 5), Enzyme("XmaI", "CCCGGG", 1, 5),
    Enzyme("ZraI", "GACGTC", 3, 3),
]
_CURATED_BY_NAME: dict[str, Enzyme] = {e.name: e for e in _CURATED_ENZYMES}


def _biopython_available() -> bool:
    try:
        import Bio.Restriction  # noqa: F401
        return True
    except ImportError:
        return False


def list_enzymes(source: Literal["curated", "rebase", "auto"] = "auto") -> list[str]:
    if source == "curated":
        return sorted(_CURATED_BY_NAME)
    if source in ("rebase", "auto") and _biopython_available():
        from Bio.Restriction import AllEnzymes
        return sorted(e.__name__ for e in AllEnzymes)
    return sorted(_CURATED_BY_NAME)


def get_enzyme(name: str) -> Enzyme:
    if name in _CURATED_BY_NAME:
        return _CURATED_BY_NAME[name]
    if _biopython_available():
        from Bio import Restriction as BR
        cls = getattr(BR, name, None)
        if cls is None:
            raise KeyError(f"Unknown enzyme: {name}")
        site = str(cls.site)
        try:
            fwd = cls.fst5
            rev = cls.fst5 + cls.ovhgseq_len() if hasattr(cls, "ovhgseq_len") else cls.fst5
        except Exception:
            fwd, rev = 0, 0
        return Enzyme(name=name, site=site, cut_fwd=fwd, cut_rev=rev)
    raise KeyError(f"Unknown enzyme: {name}")


def find_sites(sequence: str, enzymes: list[str | Enzyme] | None = None,
               *, both_strands: bool = True) -> list[RestrictionSite]:
    enz_list: list[Enzyme] = (list(_CURATED_ENZYMES) if enzymes is None
                              else [e if isinstance(e, Enzyme) else get_enzyme(e) for e in enzymes])
    hits: list[RestrictionSite] = []
    seq_up = sequence.upper()
    for enz in enz_list:
        for pos in find_motif(seq_up, enz.site, allow_iupac=True):
            hits.append(RestrictionSite(
                enzyme=enz, position=pos, strand="+",
                cut_position_fwd=pos + enz.cut_fwd,
                cut_position_rev=pos + enz.cut_rev,
            ))
        if both_strands and not enz.is_palindromic:
            rc_site = reverse_complement(enz.site)
            for pos in find_motif(seq_up, rc_site, allow_iupac=True):
                site_end = pos + len(enz.site)
                hits.append(RestrictionSite(
                    enzyme=enz, position=pos, strand="-",
                    cut_position_fwd=site_end - enz.cut_rev,
                    cut_position_rev=site_end - enz.cut_fwd,
                ))
    hits.sort(key=lambda h: (h.position, h.enzyme.name))
    return hits


def suggest_compatible_ends(enzyme_name: str) -> list[str]:
    """Enzymes whose overhangs ligate compatibly (e.g. BamHI↔BglII↔BclI)."""
    target = get_enzyme(enzyme_name)
    if target.overhang_type == OverhangType.BLUNT:
        return sorted(e.name for e in _CURATED_ENZYMES
                      if e.overhang_type == OverhangType.BLUNT and e.name != target.name)
    target_over = target.overhang_sequence.upper()
    return sorted(e.name for e in _CURATED_ENZYMES
                  if e.name != target.name
                  and e.overhang_type == target.overhang_type
                  and e.overhang_sequence.upper() == target_over)
