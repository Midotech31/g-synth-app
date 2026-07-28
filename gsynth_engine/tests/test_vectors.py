"""Tests for the vector catalogue.

The catalogue carries no sequences, so what is worth testing is what it
*decides*: that a supplied sequence really is the vector that was selected,
and that the vector's tags are reported by looking at the protein rather
than by a rule that could be wrong for a lab's own copy of a backbone.
"""
from __future__ import annotations

import pytest

from gsynth_engine import vectors
from gsynth_engine.cloning import clone
from gsynth_engine.constants import RESTRICTION_ENZYMES
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.tests.test_cloning import clean_filler


def fake_vector(spec: vectors.VectorSpec, *, seed: int = 3) -> str:
    """A sequence that passes as `spec`: right length, right motifs, one site each.

    Standing in for the real backbone, which the user supplies. It exists so
    the checks can be tested without shipping 5 443 bases nobody verified.
    """
    sites = "".join(
        str(RESTRICTION_ENZYMES[e]["recognition"])
        for e in spec.unique_sites
        if e in RESTRICTION_ENZYMES
    )
    motifs = "".join(motif for _, motif in spec.motifs)
    fixed = sites + motifs
    padding = spec.length - len(fixed)
    assert padding > 0, f"{spec.name} is too short for its own features"

    # Split the filler so no site sits at the very edge of the sequence.
    head = clean_filler(padding // 2, seed)
    tail = clean_filler(padding - len(head), seed + 1)
    return head + motifs + sites + tail


# ── The catalogue itself ────────────────────────────────────────────────────


class TestCatalogue:
    def test_pet21a_is_the_default(self):
        """What the app offers before the user chooses anything.

        pET-21a(+) leads because it is the vector the G-Synth cassette is
        built around: it supplies the ribosome binding site and the ATG
        inside an NdeI site, so the default NdeI / XhoI pair works on it.
        """
        assert vectors.DEFAULT_VECTOR.key == "pET-21a"
        assert vectors.CATALOGUE[0] is vectors.DEFAULT_VECTOR
        assert vectors.DEFAULT_VECTOR.has_sequence
        assert vectors.DEFAULT_VECTOR.recommended_pairs[0] == "NdeI / XhoI"

    def test_the_bundled_sequence_validates_against_its_own_entry(self):
        """A bundled sequence that failed its own checks would be worthless."""
        for spec in vectors.CATALOGUE:
            record = vectors.sequence_of(spec.key)
            if record is None:
                continue
            check = vectors.validate(record["sequence"], spec)
            assert check.matches, f"{spec.name}: {check.problems}"
            assert record["length"] == spec.length

    def test_pet21_has_no_ndei_site(self):
        """The distinction from pET-21a(+), and the reason the default pair
        for this vector is not the G-Synth default pair."""
        record = vectors.sequence_of("pET-21")
        from gsynth_engine.cloning import find_sites

        assert find_sites(record["sequence"], "NdeI") == []
        assert not any(
            "NdeI" in pair for pair in vectors.get("pET-21").recommended_pairs
        )

    def test_pet21_does_not_claim_to_supply_a_start(self):
        assert vectors.get("pET-21").supplies_translation_start is False
        assert vectors.get("pET-28a").supplies_translation_start is True

    def test_every_entry_is_complete_enough_to_act_on(self):
        for spec in vectors.CATALOGUE:
            assert spec.length > 0
            assert spec.resistance and spec.promoter and spec.summary
            assert spec.recommended_pairs, f"{spec.name} suggests no enzyme pair"

    def test_recommended_pairs_use_known_enzymes_that_the_vector_has(self):
        """A suggested pair that the vector cannot be cut with is a trap."""
        for spec in vectors.CATALOGUE:
            for pair in spec.recommended_pairs:
                left, right = (part.strip() for part in pair.split("/"))
                assert left in RESTRICTION_ENZYMES, f"{spec.name}: {left}"
                assert right in RESTRICTION_ENZYMES, f"{spec.name}: {right}"
                assert left != right
                assert left in spec.unique_sites, f"{spec.name} has no {left}"
                assert right in spec.unique_sites, f"{spec.name} has no {right}"

    def test_keys_and_names_are_unique(self):
        keys = [spec.key for spec in vectors.CATALOGUE]
        assert len(keys) == len(set(keys))

    def test_lookup_tolerates_how_people_write_vector_names(self):
        for written in ("pET-21a", "pET21a", "PET-21A(+)", "pet 21 a", "pET-21a(+)"):
            assert vectors.get(written) is vectors.DEFAULT_VECTOR

    def test_pet21_and_pet21a_are_separate_entries(self):
        """They differ by 74 bp, and only one of them has NdeI."""
        plain, letter = vectors.get("pET-21"), vectors.get("pET-21a")
        assert plain is not letter
        assert letter.length - plain.length == 74
        assert "NdeI" in letter.unique_sites
        assert "NdeI" not in plain.unique_sites

    def test_unknown_vector_returns_none(self):
        assert vectors.get("pNotARealVector") is None

    def test_tags_are_declared_as_peptides(self):
        """Codon-independent, so a lab's own copy still matches."""
        for spec in vectors.CATALOGUE:
            for tag in spec.tags:
                assert tag.end in {"N", "C"}
                assert set(tag.motif) <= set("ACDEFGHIKLMNPQRSTVWY"), tag.name

    def test_pet21a_supplies_a_c_terminal_his_tag(self):
        """The fact that decides whether the user should add one themselves."""
        spec = vectors.get("pET-21a")
        his = [tag for tag in spec.tags if tag.name == "His-tag"]
        assert len(his) == 1
        assert his[0].end == "C"

    def test_pet28a_supplies_the_same_cassette_g_synth_builds(self):
        """Which is exactly why it needs a warning rather than silence."""
        spec = vectors.get("pET-28a")
        n_terminal = {tag.name for tag in spec.tags if tag.end == "N"}
        assert {"His-tag", "thrombin site"} <= n_terminal


# ── Validation ──────────────────────────────────────────────────────────────


class TestValidate:
    def test_a_matching_sequence_passes(self):
        spec = vectors.DEFAULT_VECTOR
        check = vectors.validate(vectors.sequence_of(spec.key)["sequence"], spec)
        assert check.matches, check.problems
        assert check.length == spec.length

    def test_the_wrong_length_is_caught(self):
        """pET-28a is 74 bp shorter than pET-21a — a silent substitution."""
        spec = vectors.get("pET-21a")
        supplied = fake_vector(vectors.get("pET-28a"))
        check = vectors.validate(supplied, spec)

        assert not check.matches
        assert any("5,369 bp" in problem for problem in check.problems)
        assert any("pET-21a(+)" in problem for problem in check.problems)

    def test_two_vectors_of_the_same_length_are_told_apart(self):
        """pET-21(+) and pET-28a(+) are both 5,369 bp with the same promoter.

        Length cannot separate them. What can: pET-28a has NdeI and NcoI,
        pET-21 has neither; pET-21 is AmpR, pET-28a is KanR.
        """
        real = vectors.sequence_of("pET-21")["sequence"]
        assert vectors.validate(real, vectors.get("pET-21")).matches

        wrong = vectors.validate(real, vectors.get("pET-28a"))
        assert not wrong.matches
        assert any("NdeI does not cut" in p for p in wrong.problems)

    def test_a_missing_promoter_is_caught(self):
        spec = vectors.DEFAULT_VECTOR
        supplied = clean_filler(spec.length, 9)
        check = vectors.validate(supplied, spec)

        assert not check.matches
        assert "T7 promoter" in check.missing_motifs

    def test_motifs_are_found_on_either_strand(self):
        """A plasmid has no canonical orientation."""
        spec = vectors.DEFAULT_VECTOR
        flipped = reverse_complement(vectors.sequence_of(spec.key)["sequence"])
        check = vectors.validate(flipped, spec)
        assert "T7 promoter" in check.found_motifs

    def test_an_empty_sequence_is_refused(self):
        check = vectors.validate("", vectors.DEFAULT_VECTOR)
        assert not check.matches

    def test_a_duplicated_site_is_a_note_not_a_refusal(self):
        """A lab's modified copy is still their vector — but they must know."""
        spec = vectors.DEFAULT_VECTOR
        supplied = vectors.sequence_of(spec.key)["sequence"]
        supplied = supplied[:100] + "GGATCC" + supplied[100:-6]
        check = vectors.validate(supplied, spec)

        assert any("BamHI cuts 2 times" in note for note in check.notes)
        assert not any("BamHI" in problem for problem in check.problems)

    def test_tolerance_allows_a_known_variant(self):
        spec = vectors.DEFAULT_VECTOR
        supplied = vectors.sequence_of(spec.key)["sequence"] + "ACGT"
        assert not vectors.validate(supplied, spec).matches
        assert vectors.validate(supplied, spec, length_tolerance=10).matches


class TestIdentify:
    def test_recognises_a_bundled_vector_exactly(self):
        spec = vectors.DEFAULT_VECTOR
        assert vectors.identify(vectors.sequence_of(spec.key)["sequence"]) is spec

    def test_recognises_a_catalogue_vector_by_length_and_motifs(self):
        spec = vectors.get("pET-22b")
        assert vectors.identify(fake_vector(spec)) is spec

    def test_refuses_to_guess_between_two_that_fit_equally(self):
        """Better to ask than to pick the first of two 5,369 bp entries."""
        ambiguous = fake_vector(vectors.get("pET-28a"))
        assert vectors.identify(ambiguous) in (vectors.get("pET-28a"), None)

    def test_returns_none_for_something_unknown(self):
        assert vectors.identify(clean_filler(1234, 5)) is None

    def test_returns_none_for_an_empty_sequence(self):
        assert vectors.identify("") is None


# ── Tags on the finished protein ────────────────────────────────────────────


class TestTagOutcomes:
    """Read off the translated product, not asserted from the enzyme."""

    def build(self, vector_key: str, *, insert: str, his: bool = True, stop: bool = False):
        """Clone into the real backbone when one is bundled."""
        spec = vectors.get(vector_key)
        left, right = (p.strip() for p in spec.recommended_pairs[0].split("/"))
        design = design_small_sequence(
            insert + ("TAA" if stop else ""),
            enzyme_pair=f"{left} / {right}",
            include_his_tag=his,
        )
        record = vectors.sequence_of(spec.key)
        backbone = record["sequence"] if record else fake_vector(spec)
        return spec, clone(
            backbone, design.forward, insert_reverse=design.reverse,
            left_enzyme=left, right_enzyme=right,
            vector_spec=spec, orf_start=design.orf_start,
        )

    def test_reports_one_outcome_per_declared_tag(self):
        spec, result = self.build("pET-21a", insert=clean_filler(60, 21))
        assert len(result.tags) == len(spec.tags)
        assert {t.name for t in result.tags} == {t.name for t in spec.tags}

    def test_a_vector_tag_is_only_credited_in_the_vectors_own_stretch(self):
        """Searching the whole protein finds the *insert's* 6xHis and calls it
        the vector's C-terminal tag — a construct that will not bind the
        column, described as one that will."""
        _, result = self.build("pET-21a", insert=clean_filler(60, 21), his=True)
        his = [t for t in result.tags if t.name == "His-tag"][0]

        assert result.protein.startswith("MGSSHHHHHH")   # the cassette's
        assert his.present
        assert his.position is not None and his.position > 30   # the vector's
        assert result.protein.endswith("HHHHHH")

    def test_the_vector_his_tag_follows_the_xhoi_leucine_glutamate(self):
        """CTCGAG reads LE in frame, then the six histidines."""
        _, result = self.build("pET-21a", insert=clean_filler(60, 21), his=False)
        assert result.protein.endswith("LEHHHHHH")

    def test_ndei_cloning_removes_the_vectors_t7_tag(self):
        """It sits just downstream of the ATG, so the insert replaces it."""
        _, result = self.build("pET-21a", insert=clean_filler(60, 21))
        t7 = [t for t in result.tags if t.name == "T7·Tag"][0]
        assert not t7.present
        assert any("T7·Tag is not on this protein" in n for n in result.warnings)

    def test_an_insert_that_stops_loses_the_vectors_c_terminal_tag(self):
        """The classic silent failure: the protein does not bind the column."""
        _, result = self.build("pET-21a", insert=clean_filler(60, 21),
                               his=False, stop=True)
        his = [t for t in result.tags if t.name == "His-tag"][0]
        assert not his.present
        note = " ".join(result.warnings)
        assert "His-tag is not on this protein" in note
        assert "in frame through XhoI" in note

    def test_two_his_tags_are_flagged(self):
        """The cassette adds one N-terminally and pET-21a adds one after XhoI."""
        _, result = self.build("pET-21a", insert=clean_filler(60, 21), his=True)
        assert result.protein.count("HHHHHH") == 2
        assert any("appears more than once" in n for n in result.warnings)

    def test_one_tag_when_the_cassette_leaves_it_to_the_vector(self):
        _, result = self.build("pET-21a", insert=clean_filler(60, 21), his=False)
        assert result.protein.count("HHHHHH") == 1
        assert not any("appears more than once" in n for n in result.warnings)

    def test_no_spec_means_no_tag_reporting(self):
        """Cloning into an unknown backbone must not invent tag claims."""
        design = design_small_sequence(clean_filler(60, 21), enzyme_pair="NdeI / XhoI")
        result = clone(
            vectors.sequence_of("pET-21a")["sequence"], design.forward,
            insert_reverse=design.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=design.orf_start,
        )
        assert result.tags == []
        assert not any("not on this protein" in n for n in result.warnings)
