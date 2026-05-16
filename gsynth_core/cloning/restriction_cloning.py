"""Classical restriction + ligase cloning."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.restriction.digest import DigestFragment, simulate_double_digestion
from gsynth_core.restriction.enzymes import OverhangType, get_enzyme


@dataclass(frozen=True, slots=True)
class CloneProduct:
    sequence: str
    is_circular: bool
    directional: bool
    vector_fragment: DigestFragment
    insert_fragment: DigestFragment


def check_ligation_compatibility(left_overhang: str, left_type: OverhangType,
                                 right_overhang: str, right_type: OverhangType) -> bool:
    """Two ends ligate when types match and overhangs are RC of each other."""
    if left_type != right_type:
        return False
    if left_type == OverhangType.BLUNT:
        return True
    from gsynth_core.sequence.ops import reverse_complement
    return left_overhang.upper() == reverse_complement(right_overhang.upper())


def simulate_restriction_cloning(vector: str, insert: str, *,
                                 enzyme_a: str, enzyme_b: str,
                                 vector_circular: bool = True) -> CloneProduct | None:
    """Simulate directional restriction cloning."""
    ea = get_enzyme(enzyme_a); eb = get_enzyme(enzyme_b)
    vec_frags = simulate_double_digestion(vector, ea, eb, circular=vector_circular)
    if not vec_frags:
        return None
    backbone = max(vec_frags, key=lambda f: f.length)
    ins_frags = simulate_double_digestion(insert, ea, eb, circular=False)
    if not ins_frags:
        return None
    for ins in ins_frags:
        if (check_ligation_compatibility(backbone.right_overhang, backbone.right_overhang_type,
                                         ins.left_overhang, ins.left_overhang_type)
            and check_ligation_compatibility(ins.right_overhang, ins.right_overhang_type,
                                             backbone.left_overhang, backbone.left_overhang_type)):
            return CloneProduct(
                sequence=backbone.sequence + ins.sequence, is_circular=True,
                directional=(backbone.left_overhang != backbone.right_overhang),
                vector_fragment=backbone, insert_fragment=ins,
            )
    return None
