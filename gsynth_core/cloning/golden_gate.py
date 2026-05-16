"""Golden Gate / MoClo assembly with type-IIS enzymes."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.restriction.enzymes import get_enzyme
from gsynth_core.sequence.ops import ensure_dna, reverse_complement


@dataclass(frozen=True, slots=True)
class GoldenGatePart:
    name: str
    insert: str
    left_fusion: str
    right_fusion: str
    full_sequence: str
    order: int


def _standard_fusion_sites(n_parts: int) -> list[tuple[str, str]]:
    """CIDAR orthogonal 4-nt overhang pool (Hamming distance ≥ 2 pairwise)."""
    pool = [
        "GGAG", "TACT", "AATG", "AGGT", "GCTT", "CGCT", "TGCC", "ACTA",
        "TTAC", "CAGA", "TAGG", "GACG", "TTCG", "ATAG", "CGAA", "AATA",
    ]
    if n_parts > len(pool):
        raise ValueError(f"Need {n_parts} orthogonal overhangs, pool has {len(pool)}")
    overhangs = pool[:n_parts]
    return [(overhangs[i], overhangs[(i + 1) % n_parts]) for i in range(n_parts)]


def design_golden_gate(parts: list[tuple[str, str]], *, enzyme: str = "BsaI",
                       fusion_sites: list[tuple[str, str]] | None = None,
                       flanking_nt: str = "A") -> list[GoldenGatePart]:
    """Design parts with type-IIS enzyme + 4-nt overhangs."""
    if not parts:
        raise ValueError("No parts provided")
    enz = get_enzyme(enzyme)
    if enz.overhang_length != 4 or enz.overhang_type.value != "5_prime":
        raise ValueError(f"{enzyme} is not a Golden Gate enzyme (need 4-nt 5' overhang)")
    if fusion_sites is None:
        fusion_sites = _standard_fusion_sites(len(parts))
    if len(fusion_sites) != len(parts):
        raise ValueError("fusion_sites length must equal parts length")
    enz_site = enz.site; rc_site = reverse_complement(enz_site)
    designed: list[GoldenGatePart] = []
    for i, ((name, insert), (lf, rf)) in enumerate(zip(parts, fusion_sites)):
        insert = ensure_dna(insert)
        full = enz_site + flanking_nt + lf + insert + rf + flanking_nt + rc_site
        designed.append(GoldenGatePart(
            name=name, insert=insert, left_fusion=lf, right_fusion=rf,
            full_sequence=full, order=i,
        ))
    return designed


def simulate_golden_gate_assembly(parts: list[GoldenGatePart], *, circular: bool = True) -> str:
    """Simulate cut-and-ligate. After digestion: lf + insert + rf per part."""
    if not parts:
        return ""
    assembled: list[str] = []
    for i, p in enumerate(parts):
        if i == 0:
            assembled.append(p.left_fusion + p.insert + p.right_fusion)
        else:
            if p.left_fusion != parts[i-1].right_fusion:
                raise ValueError(
                    f"Overhang mismatch between parts {i-1} and {i}: "
                    f"{parts[i-1].right_fusion!r} vs {p.left_fusion!r}"
                )
            assembled.append(p.insert + p.right_fusion)
    joined = "".join(assembled)
    if circular and len(parts) > 1:
        if parts[-1].right_fusion != parts[0].left_fusion:
            raise ValueError(
                f"Circular assembly impossible: {parts[-1].right_fusion!r} vs {parts[0].left_fusion!r}"
            )
        joined = joined[:-len(parts[-1].right_fusion)]
    return joined
