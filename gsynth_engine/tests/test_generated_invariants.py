"""Deterministic generated cases for invariants hand-picked examples miss.

The seed is reported by pytest through the parametrization, so every failure
is reproducible. These tests deliberately vary sequence content and length
while asserting molecular round trips, not implementation details.
"""
from __future__ import annotations

import random

import pytest

from gsynth_engine.cloning import clone, find_sites, linearise
from gsynth_engine.constants import RESTRICTION_ENZYMES
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.preflight import assembly_preflight, cloning_preflight, verification_state
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.verify import verify


def site_free_dna(length: int, seed: int) -> str:
    """Grow DNA that cannot accidentally complete a supported enzyme site."""
    rng = random.Random(seed)
    sites = [str(info["recognition"]) for info in RESTRICTION_ENZYMES.values()]
    sites += [reverse_complement(site) for site in sites]
    longest = max(map(len, sites))
    bases: list[str] = []
    while len(bases) < length:
        for base in rng.sample("ACGT", 4):
            tail = "".join(bases[-longest:]) + base
            if not any(tail.endswith(site) for site in sites):
                bases.append(base)
                break
        else:
            bases.pop()
    return "".join(bases)


def unique_site_vector(left: str, right: str, seed: int) -> str:
    vector = (
        site_free_dna(240, seed)
        + str(RESTRICTION_ENZYMES[left]["recognition"])
        + site_free_dna(37, seed + 100)
        + str(RESTRICTION_ENZYMES[right]["recognition"])
        + site_free_dna(180, seed + 200)
    )
    assert len(find_sites(vector, left)) == 1
    assert len(find_sites(vector, right)) == 1
    return vector


@pytest.mark.parametrize("seed", range(20))
def test_generated_assemblies_reconstruct_both_strands(seed: int):
    length = random.Random(seed).randint(45, 520)
    insert = site_free_dna(length, seed + 1000)
    plan = design_merzoug_assembly(
        insert,
        enzyme_pair="NdeI / XhoI",
        target_oligo_length=random.Random(seed + 1).choice([55, 70, 90]),
    )
    assert plan.verify() == []
    assert assembly_preflight(plan).verdict != "blocked"


@pytest.mark.parametrize(
    ("left", "right"),
    [
        ("NdeI", "XhoI"),       # two 5' ends
        ("KpnI", "SacI"),       # two 3' ends
        ("NdeI", "KpnI"),       # mixed polarity
        ("EcoRV", "SmaI"),      # blunt ends
        ("BamHI", "EcoRI"),
        ("ApaI", "PstI"),
    ],
)
@pytest.mark.parametrize("seed", range(5))
def test_generated_clone_digest_round_trip(left: str, right: str, seed: int):
    insert = site_free_dna(60 + seed * 17, seed + 3000)
    designed = design_small_sequence(insert, enzyme_pair=f"{left} / {right}")
    vector = unique_site_vector(left, right, seed + 4000)
    result = clone(
        vector,
        designed.forward,
        insert_reverse=designed.reverse,
        left_enzyme=left,
        right_enzyme=right,
    )
    assert result.is_clonable, result.problems
    assert cloning_preflight(result).verdict != "blocked"
    recut = linearise(result.plasmid, left_enzyme=left, right_enzyme=right)
    assert recut.removed_length == len(designed.forward)
    assert recut.length == result.backbone_length


@pytest.mark.parametrize("seed", range(20))
def test_generated_single_base_changes_never_report_verified(seed: int):
    design = site_free_dna(180, seed + 7000)
    changed = list(design)
    position = random.Random(seed).randrange(len(changed))
    changed[position] = next(base for base in "ACGT" if base != changed[position])
    report = verify(design, {f"read-{seed}": "".join(changed)}, trim=0)
    assert not report.is_verified
    assert verification_state(report) == "differences_detected"
    assert any(difference.position == position for difference in report.differences)
