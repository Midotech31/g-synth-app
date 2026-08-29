"""Six-frame translation and ORF coordinates."""
from gsynth_engine.cloning import translate
from gsynth_engine.constants import CLEAVAGE_SITES
from gsynth_engine.orf import analyse_sequence
from gsynth_engine.sequence import reverse_complement


def test_all_six_frames_are_reported_in_conventional_order():
    result = analyse_sequence("ATGAAATAA" + "CCC" * 4, minimum_codons=1)
    assert [frame.frame for frame in result.frames] == [1, 2, 3, -1, -2, -3]


def test_forward_orf_has_half_open_top_strand_coordinates():
    result = analyse_sequence("CCCATGAAATAAGGG", minimum_codons=2)
    orf = next(record for record in result.orfs if record.strand == "forward")
    assert (orf.start, orf.end, orf.frame) == (3, 12, 1)
    assert orf.dna == "ATGAAATAA"
    assert orf.protein == "MK"
    assert orf.stop_codon == "TAA"


def test_reverse_orf_maps_back_to_original_coordinates():
    reverse_orf = "ATGCCCTAG"
    sequence = "AAAA" + reverse_complement(reverse_orf) + "GG"
    result = analyse_sequence(sequence, minimum_codons=2)
    orf = next(record for record in result.orfs if record.strand == "reverse")
    assert (orf.start, orf.end) == (4, 13)
    assert orf.dna == reverse_orf
    assert orf.protein == "MP"


def test_nested_starts_are_kept_as_distinct_candidates():
    result = analyse_sequence("ATGAAAATGCCCTAAGGG", minimum_codons=1)
    proteins = [record.protein for record in result.orfs if record.strand == "forward"]
    assert proteins == ["MKMP", "MP"]


def test_first_atg_translation_stays_in_the_selected_frame():
    result = analyse_sequence("CATGAAATAA", minimum_codons=1)
    plus_two = next(frame for frame in result.frames if frame.frame == 2)
    assert plus_two.first_atg == 1
    assert plus_two.protein_from_first_atg == "MK*"


def test_factor_xa_legacy_and_modern_dna_encode_the_same_site():
    assert translate(CLEAVAGE_SITES["Factor Xa"]) == "IEGR"
    assert translate(CLEAVAGE_SITES["Factor Xa (legacy DNA)"]) == "IEGR"
    assert CLEAVAGE_SITES["Factor Xa"] != CLEAVAGE_SITES["Factor Xa (legacy DNA)"]
