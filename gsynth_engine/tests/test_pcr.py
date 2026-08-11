"""PCR: what the two oligos would actually produce, and what cutting it gives.

These tests are claims about molecules. The product is checked against the
template it was copied from rather than against the design that built it; the
insert's ends are read off the cut duplex rather than looked up from the
enzyme table; and the whole path — amplify, digest, ligate — is run into a
real vector, because each of those steps can be individually right and still
not compose.
"""
import pytest

from gsynth_engine import cloning, vectors
from gsynth_engine.cloning import _observed_insert_ends, translate
from gsynth_engine.pcr import (
    DEFAULT_CLAMP,
    MAX_ANNEAL,
    MIN_ANNEAL,
    _cross_dimer,
    _primer_warnings,
    design_pcr,
)
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.thermo import ANNEALING, PCR, melting_temperature

#: A gene with no NdeI, XhoI, BamHI or EcoRI site of its own, and a length
#: that is a whole number of codons — so a frame warning here means the code
#: put it there, not the fixture.
GENE = (
    "ATGAAAGGTGAAGAATTGTTCACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTT"
    "AACGGTCACAAATTCTCTGTTTCTGGTGAAGGTGAAGGTGATGCTACCTACGGTAAACTGACC"
    "CTGAAA"
)


class TestConventionalPcr:
    """No enzymes named: copy a region and give back exactly that."""

    def test_the_product_is_the_region_and_nothing_else(self):
        r = design_pcr(GENE)
        assert r.product == GENE
        assert r.amplified_region == GENE
        assert r.digest is None

    def test_a_sub_region_is_amplified_to_its_own_bounds(self):
        r = design_pcr(GENE, target_start=30, target_end=120)
        assert r.product == GENE[30:120]
        assert r.product_length == 90

    def test_both_primers_bind_the_template_they_were_designed_from(self):
        """The forward primer reads off the top strand; the reverse primer is
        written 5'→3' in its own direction, so it is the reverse complement of
        the stretch it covers. Getting this backwards is invisible on a
        palindromic test fixture, which is why the check is against the gene."""
        r = design_pcr(GENE)
        assert r.forward.anneals in GENE
        assert reverse_complement(r.reverse.anneals) in GENE

    def test_the_primers_sit_at_the_ends_of_the_requested_region(self):
        r = design_pcr(GENE, target_start=12, target_end=126)
        assert r.forward.start == 12
        assert r.reverse.end == 126

    def test_primer_lengths_stay_inside_the_stated_bounds(self):
        r = design_pcr(GENE)
        for primer in (r.forward, r.reverse):
            assert MIN_ANNEAL <= primer.anneal_length <= MAX_ANNEAL

    def test_a_region_too_short_for_two_primers_is_refused(self):
        with pytest.raises(SequenceError, match="too short"):
            design_pcr(GENE, target_start=0, target_end=20)

    def test_an_empty_template_is_refused(self):
        with pytest.raises(SequenceError, match="template"):
            design_pcr("")

    def test_a_region_outside_the_template_is_refused(self):
        with pytest.raises(SequenceError, match="inside the template"):
            design_pcr(GENE, target_start=0, target_end=len(GENE) + 50)


class TestAnnealingTemperature:
    """Ta comes from the part that binds in cycle one, not from the whole
    oligo. This is the difference between a reaction that works and one that
    produces nothing, and on a tailed primer the two figures are ten degrees
    or more apart."""

    def test_tm_is_reported_for_the_annealing_region_not_the_whole_primer(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.forward.tm == pytest.approx(
            melting_temperature(r.forward.anneals, conditions=PCR), abs=0.1
        )
        assert r.forward.tm_full == pytest.approx(
            melting_temperature(r.forward.sequence, conditions=PCR), abs=0.1
        )

    def test_a_tail_raises_the_full_tm_well_above_the_annealing_tm(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.forward.tm_full > r.forward.tm + 5

    def test_ta_follows_the_lower_annealing_tm(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.annealing_temperature == pytest.approx(
            min(r.forward.tm, r.reverse.tm) - 5.0, abs=0.1
        )

    def test_ta_ignores_the_tail(self):
        """The tailed and untailed designs share their annealing regions, so
        they must share Ta. If the tail leaked into the calculation this is
        where it shows."""
        plain = design_pcr(GENE)
        tailed = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert plain.annealing_temperature == tailed.annealing_temperature

    def test_pcr_conditions_are_not_the_annealing_reaction(self):
        """A Tm quoted under the wrong reaction is not actionable. PCR carries
        Mg²⁺, which the G-Synth annealing buffer does not, so the same oligo
        melts at a different temperature in each."""
        oligo = GENE[:24]
        assert melting_temperature(oligo, conditions=PCR) != pytest.approx(
            melting_temperature(oligo, conditions=ANNEALING), abs=0.5
        )


class TestCloningTails:
    def test_each_primer_carries_clamp_then_site(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.forward.tail.endswith("CATATG")
        assert len(r.forward.tail) == DEFAULT_CLAMP + len("CATATG")
        assert r.forward.sequence == r.forward.tail + r.forward.anneals

    def test_the_reverse_tail_carries_the_site_reverse_complemented(self):
        """The reverse primer is written in its own direction, so the site it
        adds must appear reverse-complemented in the oligo and the right way
        round in the product. EcoRI is not palindromic-equivalent under this
        mistake the way XhoI is, so it is the one that catches it."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="EcoRI")
        assert r.reverse.tail.endswith(reverse_complement("GAATTC"))
        assert "GAATTC" in r.product

    def test_the_product_carries_both_sites_and_the_region_verbatim(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert "CATATG" in r.product
        assert "CTCGAG" in r.product
        assert r.amplified_region in r.product
        assert r.amplified_region == GENE

    def test_the_clamp_can_be_shortened_and_says_what_that_costs(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI", clamp=2)
        assert len(r.forward.tail) == 2 + len("CATATG")
        assert any("clamp" in w for w in r.warnings)

    def test_naming_one_enzyme_alone_is_refused(self):
        with pytest.raises(SequenceError, match="both ends or for neither"):
            design_pcr(GENE, left_enzyme="NdeI")

    def test_an_unknown_enzyme_is_refused_by_name(self):
        with pytest.raises(SequenceError, match="not an enzyme"):
            design_pcr(GENE, left_enzyme="NotAnEnzyme", right_enzyme="XhoI")


class TestTheDigest:
    """What cutting the product actually leaves."""

    def test_the_insert_presents_the_ends_the_enzymes_leave(self):
        """Read off the cut molecule, not looked up from the table — a value
        copied from the table agrees with the table whatever the bases spell."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.digest.left_end.sequence == "TA"
        assert r.digest.left_end.kind == "5'"
        assert r.digest.right_end.sequence == "TCGA"
        assert r.digest.right_end.kind == "5'"

    def test_polarity_is_read_per_end_not_per_strand(self):
        """Both of these overhangs sit on the top strand, and both are 5' —
        but at opposite ends of the fragment, which is the case that catches
        polarity derived from the strand alone."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.digest.left_end.strand == "top"
        assert r.digest.right_end.strand == "bottom"

    def test_the_digest_keeps_the_gene_and_discards_only_the_stubs(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert GENE in r.digest.top
        assert r.digest.length < r.product_length
        assert r.digest.trimmed_left + r.digest.trimmed_right == (
            r.product_length - r.digest.length
        )

    def test_the_two_strands_are_staggered_not_a_plain_reverse_complement(self):
        """The stagger between the strands *is* the overhang. Returning the
        bottom strand as the plain reverse complement of the top describes a
        blunt fragment whatever the enzymes leave — and it stays hidden for
        exactly as long as the ends are passed around instead of measured."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.digest.bottom != reverse_complement(r.digest.top)
        # NdeI leaves 2 and XhoI leaves 4, both on the strand that protrudes.
        assert len(r.digest.bottom) - len(r.digest.top) == 4 - 2

    def test_the_ends_can_be_measured_back_off_the_two_strands(self):
        """`clone()` re-reads the ends from the strands it is handed rather
        than trusting the ones it is told. This asserts the strands alone
        carry the same answer, because that is the path a caller who passes
        only `top` and `bottom` will take."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        left, right = _observed_insert_ends(r.digest.top, r.digest.bottom, "NdeI")
        assert (left.sequence, left.strand) == (
            r.digest.left_end.sequence, r.digest.left_end.strand
        )
        assert (right.sequence, right.strand) == (
            r.digest.right_end.sequence, r.digest.right_end.strand
        )

    def test_a_product_the_enzyme_does_not_cut_is_refused(self):
        with pytest.raises(SequenceError, match="does not cut"):
            cloning.digest_linear(GENE, left_enzyme="NdeI", right_enzyme="XhoI")

    def test_enzymes_the_wrong_way_round_are_refused(self):
        """Swapping them would return the stub that was meant to be thrown
        away, at a length that still looks plausible on a gel."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        with pytest.raises(SequenceError, match="wrong way round"):
            cloning.digest_linear(r.product, left_enzyme="XhoI", right_enzyme="NdeI")


class TestSitesInsideTheGene:
    """The check that decides whether a strategy works at all."""

    def test_an_internal_site_is_a_problem_not_a_warning(self):
        """The digest that opens the ends cuts the middle too, and what goes
        into the ligation is a piece of the gene. No reaction condition
        rescues that, so it blocks rather than advises."""
        gene = "ATG" + "GCTAGC" + GENE[3:]     # NheI site planted at position 3
        r = design_pcr(gene, left_enzyme="NheI", right_enzyme="XhoI")
        assert not r.is_clean
        assert any("cuts inside" in p for p in r.problems)

    def test_no_insert_is_offered_when_the_strategy_is_broken(self):
        gene = "ATG" + "GCTAGC" + GENE[3:]
        r = design_pcr(gene, left_enzyme="NheI", right_enzyme="XhoI")
        assert r.digest is None

    def test_a_clean_gene_raises_no_problem(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        assert r.is_clean
        assert r.problems == []


class TestReadingFrame:
    def test_an_atg_supplying_enzyme_anchors_the_frame_on_its_own_atg(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI", keep_frame=True)
        start = r.insert_orf_start
        assert r.digest.top[start:start + 3] == "ATG"

    def test_a_duplicated_start_codon_is_reported(self):
        """NdeI's CATATG carries an ATG. A gene that also begins with one
        yields Met-Met, which is a real extra residue on the product."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI", keep_frame=True)
        start = r.insert_orf_start
        assert translate(r.digest.top[start:start + 6]) == "MM"
        assert any("Met-Met" in w for w in r.warnings)

    def test_a_non_atg_enzyme_puts_the_region_on_a_codon_boundary(self):
        """Here the vector supplies the start codon upstream, so the insert
        must begin a whole number of codons from its own 5' end."""
        r = design_pcr(GENE, left_enzyme="BamHI", right_enzyme="XhoI", keep_frame=True)
        assert r.insert_orf_start % 3 == 0
        start = r.insert_orf_start
        assert translate(r.digest.top[start:start + 6]) == "MK"

    def test_a_region_that_is_not_whole_codons_is_flagged(self):
        r = design_pcr(GENE[:-2], left_enzyme="BamHI", right_enzyme="XhoI",
                       keep_frame=True)
        assert any("whole number of codons" in w for w in r.warnings)

    def test_a_well_formed_design_is_not_flagged_for_frame(self):
        """A warning that fires on every design is one the reader learns to
        skip past. The insert carries an overhang at each end that belongs to
        no codon, so its own length is never a multiple of three — testing
        that instead of the region would flag everything."""
        r = design_pcr(GENE, left_enzyme="BamHI", right_enzyme="XhoI",
                       keep_frame=True)
        assert not any("whole number of codons" in w for w in r.warnings)

    def test_an_internal_stop_is_reported(self):
        gene = GENE[:30] + "TAA" + GENE[33:]
        r = design_pcr(gene, left_enzyme="BamHI", right_enzyme="XhoI",
                       keep_frame=True)
        assert any("stop codon" in w for w in r.warnings)


class TestPrimerQuality:
    """Each check is asserted on an input that trips it *and* on one that does
    not. A quality check only ever tested against the case it fires on cannot
    be distinguished from one that fires on everything."""

    def test_a_weak_three_prime_end_is_reported(self):
        notes = _primer_warnings("GGCACCGGTGTTGTTCCGATTCTAA", label="forward")
        assert any("ends in A or T" in n for n in notes)

    def test_a_strong_three_prime_end_is_not_reported(self):
        notes = _primer_warnings("GAAGAATTGTTCACCGGTGTTGTTC", label="forward")
        assert not any("ends in A or T" in n for n in notes)

    def test_a_gc_clamp_that_is_too_strong_is_reported(self):
        notes = _primer_warnings("ATGAAAGGTGAAGAATTGGCGCC", label="forward")
        assert any("last five bases" in n for n in notes)

    def test_an_at_rich_region_is_reported(self):
        notes = _primer_warnings("ATATATATATATATATATAT", label="forward")
        assert any("Below 40%" in n for n in notes)

    def test_a_gc_rich_region_is_reported(self):
        notes = _primer_warnings("GCGCGGCCGCGGCCGCGGCC", label="forward")
        assert any("Above 65%" in n for n in notes)

    def test_a_long_homopolymer_is_reported(self):
        notes = _primer_warnings("ATGAAAAAAAAGGTCACCGGTGTT", label="forward")
        assert any("identical bases" in n for n in notes)

    def test_a_balanced_primer_draws_no_complaint(self):
        assert _primer_warnings("GAAGAATTGTTCACCGGTGTTGTTC", label="forward") == []

    def test_primer_dimer_between_the_pair_is_detected(self):
        """Two primers complementary at their 3' ends extend one another into
        a short product that then amplifies far better than the template."""
        forward = "ATGCGTACGTTGACCGGGCCCAAAGGGCC"
        reverse = "GGCCCTTTGGGCCCGGTCAACGTACGCAT"
        assert _cross_dimer(forward, reverse) >= 5

    def test_an_unrelated_pair_is_not_called_a_dimer(self):
        assert _cross_dimer("ATGAAAGGTGAAGAATTGTT", "TTTCAGGGTCAGTTTACCGT") < 5

    def test_the_gc_clamp_property_reads_the_three_prime_end(self):
        r = design_pcr(GENE)
        assert r.forward.has_gc_clamp == any(b in "GC" for b in r.forward.anneals[-2:])


class TestTheWholeJourney:
    """Amplify, cut, ligate. Each step can be right on its own and still not
    compose — this is the only test that says the three of them do."""

    def test_a_pcr_product_clones_into_pet21a(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI", name="gene")
        assert r.is_clean

        record = vectors.sequence_of("pET-21a")
        result = cloning.clone(
            record["sequence"], r.digest.top,
            left_enzyme="NdeI", right_enzyme="XhoI",
            vector_spec=vectors.get("pET-21a"),
            insert_reverse=r.digest.bottom,
            insert_left_end=r.digest.left_end,
            insert_right_end=r.digest.right_end,
            name="pET21a-gene",
        )
        assert result.is_clonable
        assert result.problems == []

    def test_the_cloned_insert_is_the_gene_that_was_amplified(self):
        """Round trip: what comes back out of the plasmid is what went into
        the reaction."""
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        record = vectors.sequence_of("pET-21a")
        result = cloning.clone(
            record["sequence"], r.digest.top,
            left_enzyme="NdeI", right_enzyme="XhoI",
            insert_reverse=r.digest.bottom,
            insert_left_end=r.digest.left_end,
            insert_right_end=r.digest.right_end,
        )
        plasmid = result.plasmid
        # The cassette reads on the minus strand in pET-21a, so look on both.
        assert GENE in plasmid + plasmid or GENE in reverse_complement(plasmid) * 2

    def test_recutting_the_plasmid_returns_an_insert_of_the_right_length(self):
        r = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        record = vectors.sequence_of("pET-21a")
        result = cloning.clone(
            record["sequence"], r.digest.top,
            left_enzyme="NdeI", right_enzyme="XhoI",
            insert_reverse=r.digest.bottom,
            insert_left_end=r.digest.left_end,
            insert_right_end=r.digest.right_end,
        )
        assert result.insert_length == r.digest.length
