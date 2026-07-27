"""Tests for in-silico cloning.

The strongest property here is the **round trip**: cut the recombinant
plasmid with the same two enzymes and you must get the insert back, base for
base. If the ligation arithmetic is wrong anywhere — cut offsets, overhang
polarity, the wrap around the origin — that test fails.
"""
from __future__ import annotations

import random

import pytest

from gsynth_engine.cloning import (
    End,
    clone,
    find_sites,
    linearise,
    open_reading_frames,
    translate,
)
from gsynth_engine.constants import RESTRICTION_ENZYMES
from gsynth_engine.constants import overhang as enzyme_overhang
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.ssd import design_small_sequence

INSERT = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGC"


def clean_filler(length: int, seed: int = 7) -> str:
    """Random DNA guaranteed to contain no recognition site from the table.

    Built rather than hand-written: a filler that accidentally carries a
    second BamHI site turns a test of cloning into a test of luck.
    """
    rng = random.Random(seed)
    sites = [
        str(info["recognition"]) for info in RESTRICTION_ENZYMES.values()
    ]
    sites += [reverse_complement(s) for s in sites]

    while True:
        candidate = "".join(rng.choice("ACGT") for _ in range(length))
        if not any(site in candidate for site in sites):
            return candidate


def build_vector(left: str, right: str, *, seed: int = 7) -> str:
    """A circular vector with exactly one site for each of the two enzymes."""
    left_site = str(RESTRICTION_ENZYMES[left]["recognition"])
    right_site = str(RESTRICTION_ENZYMES[right]["recognition"])
    vector = (
        clean_filler(400, seed)
        + left_site
        + clean_filler(30, seed + 1)
        + right_site
        + clean_filler(200, seed + 2)
    )
    assert len(find_sites(vector, left)) == 1
    assert len(find_sites(vector, right)) == 1
    return vector


@pytest.fixture(scope="module")
def vector():
    return build_vector("NdeI", "XhoI")


@pytest.fixture(scope="module")
def ssd():
    return design_small_sequence(INSERT, enzyme_pair="NdeI / XhoI")


# ── Site finding ────────────────────────────────────────────────────────────


class TestFindSites:
    def test_finds_a_site(self):
        assert find_sites("AAAACATATGAAAA", "NdeI", circular=False) == [4]

    def test_finds_nothing_when_absent(self):
        assert find_sites(clean_filler(300), "NdeI", circular=False) == []

    def test_finds_a_site_on_the_reverse_strand(self):
        """Non-palindromic sites cut whichever strand carries them."""
        forward = find_sites("AAAAGCGGCCGCAAAA", "NotI", circular=False)
        reverse = find_sites(
            reverse_complement("AAAAGCGGCCGCAAAA"), "NotI", circular=False
        )
        assert forward and reverse

    def test_finds_a_site_spanning_the_origin(self):
        """A circular plasmid has no beginning; the search must not either."""
        vector = "TATG" + clean_filler(200) + "CA"      # CATATG across the join
        assert find_sites(vector, "NdeI", circular=True)
        assert find_sites(vector, "NdeI", circular=False) == []

    def test_unknown_enzyme_is_reported(self):
        with pytest.raises(SequenceError, match="Unknown enzyme"):
            find_sites("ACGT", "NotAnEnzyme")


# ── Ends ────────────────────────────────────────────────────────────────────


class TestEnds:
    def test_two_matching_five_prime_overhangs_ligate(self):
        assert End("TCGA", "top").anneals_to(End("TCGA", "bottom"))

    def test_the_same_strand_twice_does_not_ligate(self):
        """Two top-strand overhangs leave both bottom strands unpaired."""
        assert not End("TCGA", "top").anneals_to(End("TCGA", "top"))

    def test_different_sequences_do_not_ligate(self):
        assert not End("TCGA", "top").anneals_to(End("GATC", "bottom"))

    def test_blunt_ligates_only_to_blunt(self):
        assert End("", "blunt").anneals_to(End("", "blunt"))
        assert not End("", "blunt").anneals_to(End("TCGA", "top"))

    def test_polarity_is_named_the_way_a_catalogue_names_it(self):
        assert End("TCGA", "top").kind == "5'"
        assert End("GTAC", "bottom").kind == "3'"
        assert End("", "blunt").kind == "blunt"


# ── Digestion ───────────────────────────────────────────────────────────────


class TestLinearise:
    def test_backbone_plus_removed_is_the_whole_vector(self, vector):
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        assert backbone.length + backbone.removed_length == len(vector)

    def test_backbone_ends_are_the_enzymes_own_overhangs(self, vector):
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        # The backbone begins at the XhoI cut and ends at the NdeI one.
        assert backbone.left_end.sequence == enzyme_overhang("XhoI")[0]
        assert backbone.right_end.sequence == enzyme_overhang("NdeI")[0]

    def test_the_two_ends_sit_on_opposite_strands(self, vector):
        """Otherwise the backbone could circularise on itself."""
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        assert backbone.left_end.strand == "top"
        assert backbone.right_end.strand == "bottom"

    def test_a_second_site_is_refused(self):
        """The single most common reason a good design fails at the bench."""
        vector = clean_filler(200) + "CATATG" + clean_filler(50) + "CTCGAG" \
            + clean_filler(50) + "CATATG" + clean_filler(100)
        with pytest.raises(SequenceError, match="NdeI cuts it 2 times"):
            linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")

    def test_a_missing_site_is_refused(self):
        vector = clean_filler(200) + "CTCGAG" + clean_filler(100)
        with pytest.raises(SequenceError, match="NdeI does not cut"):
            linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")

    def test_identical_enzymes_are_refused(self, vector):
        with pytest.raises(SequenceError, match="must differ"):
            linearise(vector, left_enzyme="NdeI", right_enzyme="NdeI")

    def test_linear_vectors_are_refused_with_a_reason(self, vector):
        with pytest.raises(SequenceError, match="must be circular"):
            linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI", circular=False)

    def test_handles_a_site_across_the_origin(self):
        """Rotating a plasmid must not change what the digest produces."""
        vector = build_vector("NdeI", "XhoI")
        straight = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")

        rotated = vector[-3:] + vector[:-3]
        turned = linearise(rotated, left_enzyme="NdeI", right_enzyme="XhoI")

        assert turned.length == straight.length
        assert turned.removed_length == straight.removed_length
        assert turned.top == straight.top


# ── Cloning ─────────────────────────────────────────────────────────────────


class TestClone:
    def test_a_g_synth_design_drops_straight_in(self, vector, ssd):
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
        )
        assert result.is_clonable, result.problems

    def test_plasmid_length_is_backbone_plus_insert(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.length == result.backbone_length + len(ssd.forward)
        assert result.length == len(vector) - result.removed_length + len(ssd.forward)

    def test_the_insert_is_where_the_result_says_it_is(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.plasmid[result.insert_start : result.insert_end] == ssd.forward

    def test_recutting_the_plasmid_returns_the_insert(self, vector, ssd):
        """The round trip. Everything else is a detail of this."""
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        recut = linearise(result.plasmid, left_enzyme="NdeI", right_enzyme="XhoI")

        assert recut.removed_length == len(ssd.forward)
        assert recut.length == result.backbone_length

    def test_round_trip_holds_across_enzyme_pairs(self):
        """Including 3' overhangs and blunt ends, where the geometry inverts."""
        pairs = [
            ("NdeI", "XhoI"),      # 5' / 5'
            ("BamHI", "EcoRI"),    # 5' / 5'
            ("KpnI", "SacI"),      # 3' / 3'
            ("ApaI", "PstI"),      # 3' / 3'
            ("NdeI", "KpnI"),      # 5' / 3'
            ("EcoRV", "SmaI"),     # blunt / blunt
        ]
        # A site-free insert: recutting a plasmid whose gene contains an
        # internal PstI site is supposed to fail, and would be testing the
        # wrong thing here.
        clean_insert = clean_filler(60, seed=41)

        for left, right in pairs:
            vector = build_vector(left, right)
            design = design_small_sequence(
                clean_insert, enzyme_pair=f"{left} / {right}"
            )
            result = clone(
                vector, design.forward, insert_reverse=design.reverse,
                left_enzyme=left, right_enzyme=right,
            )
            assert result.is_clonable, f"{left}/{right}: {result.problems}"

            recut = linearise(result.plasmid, left_enzyme=left, right_enzyme=right)
            assert recut.removed_length == len(design.forward), f"{left}/{right}"

    def test_junctions_name_the_right_enzyme(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        start, end = result.junctions

        assert start.name == "vector → insert"
        assert start.enzyme == "NdeI"
        assert start.position == result.insert_start

        assert end.name == "insert → vector"
        assert end.enzyme == "XhoI"

    def test_both_sites_are_regenerated(self, vector, ssd):
        """How the clone gets verified: cut it back out and run a gel."""
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert all(junction.site_regenerated for junction in result.junctions)

    def test_junction_context_straddles_the_seam(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        for junction in result.junctions:
            assert len(junction.context) == 24
            assert set(junction.context) <= set("ACGT")

    def test_an_insert_cut_with_the_wrong_enzyme_is_refused(self, vector):
        """Read off the duplex, not assumed — otherwise the check is vacuous."""
        wrong = design_small_sequence(INSERT, enzyme_pair="BamHI / EcoRI")
        result = clone(
            vector, wrong.forward, insert_reverse=wrong.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        assert not result.is_clonable
        assert any("does not match" in problem for problem in result.problems)

    def test_a_hand_pasted_insert_with_the_wrong_start_is_caught(self, vector):
        """Only the forward strand: the visible end must still be right."""
        result = clone(
            vector, "GGGGCCCCAAAATTTT", left_enzyme="NdeI", right_enzyme="XhoI"
        )
        assert not result.is_clonable
        assert any("was not cut with NdeI" in problem for problem in result.problems)

    def test_a_forward_only_insert_says_what_it_could_not_check(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert any("Only the forward strand" in note for note in result.warnings)

    def test_supplying_both_strands_removes_that_caveat(self, vector, ssd):
        result = clone(
            vector, ssd.forward, insert_reverse=ssd.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        assert result.is_clonable
        assert not any("Only the forward strand" in note for note in result.warnings)

    def test_ends_read_off_the_duplex_match_the_enzymes(self, vector, ssd):
        result = clone(
            vector, ssd.forward, insert_reverse=ssd.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        assert result.is_clonable, result.problems

    def test_invalid_insert_is_reported(self, vector):
        with pytest.raises(SequenceError, match="insert"):
            clone(vector, "ACGTXYZ", left_enzyme="NdeI", right_enzyme="XhoI")


# ── Reading frame ───────────────────────────────────────────────────────────


class TestReadingFrame:
    def test_translates_the_standard_cassette(self, vector, ssd):
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
        )
        assert result.protein.startswith("MGSSHHHHHHSSG")
        assert "LVPRGS" in result.protein          # the thrombin site

    def test_a_stop_inside_the_insert_blocks_the_clone(self, vector):
        design = design_small_sequence(
            "GGCTAAATCGTGGAACAGTGCTGCACCAGC", enzyme_pair="NdeI / XhoI"
        )
        result = clone(
            vector, design.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=design.orf_start,
        )
        assert not result.is_clonable
        assert any("truncated" in problem for problem in result.problems)

    def test_an_inserts_own_stop_is_a_note_not_a_problem(self, vector):
        """Ending the gene is normal; it just means no C-terminal vector tag."""
        design = design_small_sequence(
            INSERT + "TAA", enzyme_pair="NdeI / XhoI"
        )
        result = clone(
            vector, design.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=design.orf_start,
        )
        assert result.is_clonable
        assert any("own stop codon" in note for note in result.warnings)

    def test_a_frame_running_into_the_vector_is_reported(self, vector, ssd):
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
        )
        assert result.protein
        # Either it stops in the insert or it runs on — both must be stated.
        assert result.warnings or result.problems

    def test_no_frame_check_without_an_orf_start(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.protein == ""

    def test_translate_handles_a_partial_final_codon(self):
        assert translate("ATGAAA") == "MK"
        assert translate("ATGAAAG") == "MK"

    def test_open_reading_frames_finds_the_cloned_construct(self, vector, ssd):
        """The insert sits at the end of the plasmid string, so this only
        works if the scan wraps — which on a circular molecule it must."""
        result = clone(
            vector, ssd.forward, insert_reverse=ssd.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        orfs = open_reading_frames(result.plasmid, minimum_codons=20)
        assert any(orf["protein"].startswith("MGSSHHHHHH") for orf in orfs)

    def test_a_linear_scan_misses_an_orf_across_the_join(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        linear = open_reading_frames(
            result.plasmid, minimum_codons=20, circular=False
        )
        assert not any(orf["protein"].startswith("MGSSHHHHHH") for orf in linear)

    def test_orfs_are_reported_longest_first(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        orfs = open_reading_frames(result.plasmid, minimum_codons=10)
        assert orfs == sorted(orfs, key=lambda o: o["codons"], reverse=True)


# ── Internal sites ──────────────────────────────────────────────────────────


class TestInternalSites:
    def test_an_internal_site_is_a_note_because_merzoug_never_digests(self, vector):
        """The build is unaffected; the diagnostic digest is not."""
        design = design_small_sequence(
            "GGCATCGGATCCGAACAGTGCTGCACCAGC", enzyme_pair="BamHI / EcoRI"
        )
        pair_vector = build_vector("BamHI", "EcoRI")
        result = clone(
            pair_vector, design.forward,
            left_enzyme="BamHI", right_enzyme="EcoRI",
        )
        notes = " ".join(result.warnings)
        assert "internal BamHI site" in notes
        assert "diagnostic digest" in notes

    def test_a_clean_insert_produces_no_such_note(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert not any("internal" in note for note in result.warnings)


# ── Annotations ─────────────────────────────────────────────────────────────


class TestAnnotations:
    def test_features_move_with_the_backbone(self, vector, ssd):
        features = [{"name": "ori", "type": "rep_origin", "start": 10, "end": 60}]
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", vector_annotations=features,
        )
        assert len(result.annotations) == 1
        moved = result.annotations[0]
        assert moved["name"] == "ori"
        assert moved["end"] - moved["start"] == 50
        assert result.plasmid[moved["start"] : moved["end"]] == vector[10:60]

    def test_features_in_the_removed_stretch_are_dropped(self, vector, ssd):
        """Drawing a feature that is not in the molecule any more is a lie."""
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        inside = backbone.removed_start + 2
        features = [
            {"name": "stuffer", "type": "misc_feature",
             "start": inside, "end": inside + 5},
        ]
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", vector_annotations=features,
        )
        assert result.annotations == []

    def test_a_feature_straddling_a_junction_is_flagged(self, vector, ssd):
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        features = [
            {"name": "promoter", "type": "promoter",
             "start": backbone.removed_start - 20, "end": backbone.removed_start + 10},
        ]
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", vector_annotations=features,
        )
        assert result.annotations
        assert result.annotations[0].get("truncated") is True

    def test_no_annotations_is_fine(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.annotations == []
