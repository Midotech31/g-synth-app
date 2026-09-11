from __future__ import annotations

import random

import pytest

from gsynth_engine import vectors
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

    rng = random.Random(seed)
    sites = [str(info["recognition"]) for info in RESTRICTION_ENZYMES.values()]
    sites += [reverse_complement(site) for site in sites]
    longest = max(len(site) for site in sites)

    out: list[str] = []
    while len(out) < length:
        for base in rng.sample("ACGT", 4):
            tail = "".join(out[-(longest - 1):]) + base
            if not any(tail.endswith(site) for site in sites):
                out.append(base)
                break
        else:
            out.pop()
    return "".join(out)


def build_vector(left: str, right: str, *, seed: int = 7) -> str:

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


class TestFindSites:
    def test_finds_a_site(self):
        assert find_sites("AAAACATATGAAAA", "NdeI", circular=False) == [4]

    def test_finds_nothing_when_absent(self):
        assert find_sites(clean_filler(300), "NdeI", circular=False) == []

    def test_finds_a_site_on_the_reverse_strand(self):

        forward = find_sites("AAAAGCGGCCGCAAAA", "NotI", circular=False)
        reverse = find_sites(
            reverse_complement("AAAAGCGGCCGCAAAA"), "NotI", circular=False
        )
        assert forward and reverse

    def test_finds_a_site_spanning_the_origin(self):

        vector = "TATG" + clean_filler(200) + "CA"
        assert find_sites(vector, "NdeI", circular=True)
        assert find_sites(vector, "NdeI", circular=False) == []

    def test_unknown_enzyme_is_reported(self):
        with pytest.raises(SequenceError, match="Unknown enzyme"):
            find_sites("ACGT", "NotAnEnzyme")


class TestEnds:
    def test_two_matching_five_prime_overhangs_ligate(self):
        assert End("TCGA", "top").anneals_to(End("TCGA", "bottom"))

    def test_the_same_strand_twice_does_not_ligate(self):

        assert not End("TCGA", "top").anneals_to(End("TCGA", "top"))

    def test_different_sequences_do_not_ligate(self):
        assert not End("TCGA", "top").anneals_to(End("GATC", "bottom"))

    def test_blunt_ligates_only_to_blunt(self):
        assert End("", "blunt").anneals_to(End("", "blunt"))
        assert not End("", "blunt").anneals_to(End("TCGA", "top"))

    def test_polarity_depends_on_the_end_as_well_as_the_strand(self):

        assert End("TCGA", "top", "left").kind == "5'"
        assert End("TCGA", "top", "right").kind == "3'"
        assert End("GTAC", "bottom", "left").kind == "3'"
        assert End("GTAC", "bottom", "right").kind == "5'"
        assert End("", "blunt").kind == "blunt"


class TestLinearise:
    def test_backbone_plus_removed_is_the_whole_vector(self, vector):
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        assert backbone.length + backbone.removed_length == len(vector)

    def test_backbone_ends_are_the_enzymes_own_overhangs(self, vector):
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")

        assert backbone.left_end.sequence == enzyme_overhang("XhoI")[0]
        assert backbone.right_end.sequence == enzyme_overhang("NdeI")[0]

    def test_the_two_ends_sit_on_opposite_strands(self, vector):

        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        assert backbone.left_end.strand == "top"
        assert backbone.right_end.strand == "bottom"

    def test_a_second_site_is_refused(self):

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

        vector = build_vector("NdeI", "XhoI")
        straight = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")

        rotated = vector[-3:] + vector[:-3]
        turned = linearise(rotated, left_enzyme="NdeI", right_enzyme="XhoI")

        assert turned.length == straight.length
        assert turned.removed_length == straight.removed_length
        assert turned.top == straight.top


class TestClone:
    def test_a_g_synth_design_drops_straight_in(self, vector, ssd):
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
        )
        assert result.is_clonable, result.problems

    def test_an_internal_mismatch_between_insert_strands_blocks_cloning(self, vector, ssd):
        reverse = list(ssd.reverse)
        at = len(reverse) // 2
        reverse[at] = next(base for base in "ACGT" if base != reverse[at])
        result = clone(
            vector,
            ssd.forward,
            insert_reverse="".join(reverse),
            left_enzyme="NdeI",
            right_enzyme="XhoI",
        )
        assert not result.is_clonable
        assert any("do not pair" in problem for problem in result.problems)

    def test_plasmid_length_is_backbone_plus_insert(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.length == result.backbone_length + len(ssd.forward)
        assert result.length == len(vector) - result.removed_length + len(ssd.forward)

    def test_the_insert_is_where_the_result_says_it_is(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.plasmid[result.insert_start : result.insert_end] == ssd.forward

    def test_recutting_the_plasmid_returns_the_insert(self, vector, ssd):

        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        recut = linearise(result.plasmid, left_enzyme="NdeI", right_enzyme="XhoI")

        assert recut.removed_length == len(ssd.forward)
        assert recut.length == result.backbone_length

    def test_round_trip_holds_across_enzyme_pairs(self):

        pairs = [
            ("NdeI", "XhoI"),
            ("BamHI", "EcoRI"),
            ("KpnI", "SacI"),
            ("ApaI", "PstI"),
            ("NdeI", "KpnI"),
            ("EcoRV", "SmaI"),
        ]


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

    def test_junction_polarity_matches_the_enzyme_catalogue(self):

        for left, right in [
            ("NdeI", "XhoI"),
            ("KpnI", "SacI"),
            ("NdeI", "KpnI"),
            ("EcoRV", "SmaI"),
        ]:
            design = design_small_sequence(
                clean_filler(60, seed=11), enzyme_pair=f"{left} / {right}"
            )
            result = clone(
                build_vector(left, right), design.forward,
                insert_reverse=design.reverse,
                left_enzyme=left, right_enzyme=right,
            )
            reported = {j.enzyme: (j.kind, j.overhang) for j in result.junctions}
            for enzyme in (left, right):
                expected_sequence, expected_kind = enzyme_overhang(enzyme)
                assert reported[enzyme] == (expected_kind, expected_sequence), enzyme

    def test_both_sites_are_regenerated(self, vector, ssd):

        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert all(junction.site_regenerated for junction in result.junctions)

    def test_junction_context_straddles_the_seam(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        for junction in result.junctions:
            assert len(junction.context) == 24
            assert set(junction.context) <= set("ACGT")

    def test_an_insert_cut_with_the_wrong_enzyme_is_refused(self, vector):

        wrong = design_small_sequence(INSERT, enzyme_pair="BamHI / EcoRI")
        result = clone(
            vector, wrong.forward, insert_reverse=wrong.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        assert not result.is_clonable
        assert any("does not match" in problem for problem in result.problems)

    def test_a_hand_pasted_insert_with_the_wrong_start_is_caught(self, vector):

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


class TestReadingFrame:
    def test_translates_the_standard_cassette(self, vector, ssd):
        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
        )
        assert result.protein.startswith("MGSSHHHHHHSSG")
        assert "LVPRGS" in result.protein
        assert result.translation_start == result.insert_start + ssd.orf_start

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

        assert result.warnings or result.problems

    def test_no_frame_check_without_an_orf_start(self, vector, ssd):
        result = clone(vector, ssd.forward, left_enzyme="NdeI", right_enzyme="XhoI")
        assert result.protein == ""
        assert result.translation_start is None
        assert result.reading_frame.status == "review"

    def test_pre_digested_insert_can_expose_one_start_candidate(self):
        vector_record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(INSERT, enzyme_pair="NdeI / XhoI")
        result = clone(
            vector_record["sequence"],
            design.forward,
            insert_reverse=design.reverse,
            left_enzyme="NdeI",
            right_enzyme="XhoI",
            vector_annotations=vector_record["annotations"],
            vector_spec=vectors.get("pET-21a"),
            auto_detect_frame=True,
        )

        assert result.protein
        assert result.translation_start == result.insert_start + design.orf_start
        assert result.reading_frame.start_source == "sequence_candidate"
        assert result.reading_frame.status == "review"

    def test_translation_origin_is_retained_on_a_non_expression_backbone(self, vector, ssd):

        from dataclasses import replace

        result = clone(
            vector, ssd.forward,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=ssd.orf_start,
            vector_spec=replace(vectors.get("pET-21a"), expression_capable=False),
        )
        assert result.reading_frame.status == "not_applicable"
        assert result.translation_start == result.insert_start + ssd.orf_start
        assert result.protein.startswith("MGSSHHHHHHSSG")

    def test_real_pet21a_expression_context_is_confirmed(self):
        vector_record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(INSERT, enzyme_pair="NdeI / XhoI")
        result = clone(
            vector_record["sequence"],
            design.forward,
            insert_reverse=design.reverse,
            left_enzyme="NdeI",
            right_enzyme="XhoI",
            orf_start=design.orf_start,
            vector_annotations=vector_record["annotations"],
            vector_spec=vectors.get("pET-21a"),
        )

        assert result.reading_frame.confirmed
        assert result.reading_frame.start_codon == "ATG"
        assert result.reading_frame.rbs_name == "RBS"
        assert result.reading_frame.rbs_spacing_nt == 8
        assert result.reading_frame.promoter_name == "T7 promoter"
        assert result.reading_frame.stop_context == "vector"
        assert all(check.status == "pass" for check in result.reading_frame.checks)

    def test_wrong_declared_start_blocks_expression_frame(self):
        vector_record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(INSERT, enzyme_pair="NdeI / XhoI")
        result = clone(
            vector_record["sequence"],
            design.forward,
            insert_reverse=design.reverse,
            left_enzyme="NdeI",
            right_enzyme="XhoI",
            orf_start=design.orf_start + 1,
            vector_annotations=vector_record["annotations"],
            vector_spec=vectors.get("pET-21a"),
        )

        assert result.reading_frame.status == "block"
        start = next(check for check in result.reading_frame.checks if check.code == "FRAME_START_CODON")
        assert start.status == "block"

    def test_storage_vector_does_not_claim_expression_validation(self):
        design = design_small_sequence(INSERT, enzyme_pair="EcoRI / HindIII")
        result = clone(
            build_vector("EcoRI", "HindIII"),
            design.forward,
            insert_reverse=design.reverse,
            left_enzyme="EcoRI",
            right_enzyme="HindIII",
            orf_start=design.orf_start,
            vector_spec=vectors.get("pUC19"),
        )

        assert result.reading_frame.status == "not_applicable"
        assert not result.reading_frame.confirmed

    def test_expression_vector_without_rbs_is_blocked(self):
        vector_record = vectors.sequence_of("pET-21")
        design = design_small_sequence(INSERT, enzyme_pair="BamHI / XhoI")
        result = clone(
            vector_record["sequence"],
            design.forward,
            insert_reverse=design.reverse,
            left_enzyme="BamHI",
            right_enzyme="XhoI",
            orf_start=design.orf_start,
            vector_annotations=vector_record["annotations"],
            vector_spec=vectors.get("pET-21"),
        )

        assert result.reading_frame.status == "block"
        rbs = next(check for check in result.reading_frame.checks if check.code == "FRAME_RBS_CONTEXT")
        assert rbs.status == "block"

    def test_translate_handles_a_partial_final_codon(self):
        assert translate("ATGAAA") == "MK"
        assert translate("ATGAAAG") == "MK"

    def test_open_reading_frames_finds_the_cloned_construct(self, vector, ssd):

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


class TestInternalSites:
    def test_an_internal_site_is_a_note_because_esd_never_digests(self, vector):

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


class TestOrientation:


    def test_a_real_pet21a_keeps_its_backbone(self):

        from gsynth_engine import vectors

        record = vectors.sequence_of("pET-21a")
        backbone = linearise(
            record["sequence"], left_enzyme="NdeI", right_enzyme="XhoI"
        )

        assert backbone.reversed_insert, "NdeI cuts after XhoI in pET-21a"
        assert backbone.length > 5000, "the backbone is nearly the whole vector"
        assert backbone.removed_length < 200, "only the stuffer comes out"
        assert backbone.length + backbone.removed_length == record["length"]

    def test_the_insert_lands_against_the_c_terminal_tag(self):

        from gsynth_engine import vectors

        backbone = linearise(
            vectors.sequence_of("pET-21a")["sequence"],
            left_enzyme="NdeI", right_enzyme="XhoI",
        )

        assert backbone.top.startswith("TCGAGCACCACCACCACCACCACTGA")

    def test_round_trip_holds_on_a_real_vector(self):
        from gsynth_engine import vectors

        record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(
            clean_filler(90, 31), enzyme_pair="NdeI / XhoI"
        )
        result = clone(
            record["sequence"], design.forward, insert_reverse=design.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
        )
        assert result.is_clonable, result.problems
        assert result.length == record["length"] - result.removed_length + len(
            design.forward
        )

        recut = linearise(result.plasmid, left_enzyme="NdeI", right_enzyme="XhoI")
        assert recut.removed_length == len(design.forward)

    def test_a_forward_oriented_vector_is_not_flipped(self):

        vector = build_vector("NdeI", "XhoI")
        backbone = linearise(vector, left_enzyme="NdeI", right_enzyme="XhoI")
        assert not backbone.reversed_insert

    def test_vector_features_follow_the_flip(self):

        from gsynth_engine import vectors

        record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(clean_filler(60, 21), enzyme_pair="NdeI / XhoI")
        result = clone(
            record["sequence"], design.forward, insert_reverse=design.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
            vector_annotations=record["annotations"],
        )
        by_name = {a["name"]: a for a in result.annotations}
        assert "AmpR" in by_name and "ori" in by_name

        amp = by_name["AmpR"]
        original = [a for a in record["annotations"] if a["name"] == "AmpR"][0]
        assert amp["end"] - amp["start"] == original["end"] - original["start"]

        assert result.plasmid[amp["start"]:amp["end"]] == reverse_complement(
            record["sequence"][original["start"]:original["end"]]
        )
