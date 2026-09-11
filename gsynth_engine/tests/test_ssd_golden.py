import pytest

from gsynth_engine.cloning import translate
from gsynth_engine.constants import CLEAVAGE_SITES
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence

EX1_INPUT = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA"
EX1_FORWARD = (
    "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAAC"
)
EX1_REVERSE = (
    "TCGAGTTAGCCGCAGTAGTTTTCCAGCTGGTACAGGCTGCAGATGCTGGTGCAGCACTGTTCCACGATGCC"
    "AGAACCACGCGGCACCAGACCAGAAGAGTGGTGGTGGTGGTGGTGAGAAGAACCCA"
)


EX2_INPUT = (
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCAAAAACCCGCCGCTAA"
)
EX2_FORWARD = (
    "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCAAAAACCCGCCGCTAAC"
)
EX2_REVERSE = (
    "TCGAGTTAGCGGCGGGTTTTTGGGGTGTAAAAAAAGCCGCGTTCGCCGCACACCAGGTACAGCGCTTCC"
    "ACCAGATGGCTGCCGCACAGATGCTGGTTCACAAAAGAACCACGCGGCACCAGACCAGAAGAGTGGTGG"
    "TGGTGGTGGTGAGAAGAACCCA"
)

GOLDEN = [
    pytest.param(EX1_INPUT, EX1_FORWARD, EX1_REVERSE, id="example-1"),
    pytest.param(EX2_INPUT, EX2_FORWARD, EX2_REVERSE, id="example-2"),
]


@pytest.mark.parametrize(("insert", "expected_forward", "expected_reverse"), GOLDEN)
def test_forward_oligo_matches_specification(insert, expected_forward, expected_reverse):
    result = design_small_sequence(
        insert, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    assert result.forward == expected_forward


@pytest.mark.parametrize(("insert", "expected_forward", "expected_reverse"), GOLDEN)
def test_reverse_oligo_matches_specification(insert, expected_forward, expected_reverse):
    result = design_small_sequence(
        insert, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    assert result.reverse == expected_reverse


@pytest.mark.parametrize(("insert", "expected_forward", "expected_reverse"), GOLDEN)
def test_duplex_presents_the_expected_sticky_ends(insert, expected_forward, expected_reverse):

    result = design_small_sequence(
        insert, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    forward = result.forward
    bottom_in_top_sense = reverse_complement(result.reverse)


    offset = 2
    assert forward[:offset] == "TA", "NdeI 5' overhang"
    assert result.left_overhang == "TA"


    core_top = forward[offset:]
    core_bottom = bottom_in_top_sense[: len(core_top)]
    assert core_top == core_bottom, "forward and reverse must anneal without mismatches"


    tail = bottom_in_top_sense[len(core_top):]
    assert reverse_complement(tail) == "TCGA", "XhoI 5' overhang"
    assert result.right_overhang == "TCGA"


def test_ndei_start_codon_overlaps_its_own_overhang():

    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    assert result.forward.startswith("TATG")
    assert result.left_overhang == "TA"
    assert result.orf_start == 1
    assert result.coding_region.startswith("ATG")

    assert "CA" + result.forward[:4] == "CATATG"


def test_his_tag_and_thrombin_site_are_present_in_frame():

    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    coding = result.coding_region
    assert coding.startswith("ATG")
    assert "CACCACCACCACCACCAC" in coding, "6×His"
    assert "CTGGTGCCGCGTGGTTCT" in coding, "Thrombin site"

    assert coding.index("CACCACCACCACCACCAC") % 3 == 0
    assert coding.index("CTGGTGCCGCGTGGTTCT") % 3 == 0
    assert coding.index(EX1_INPUT) % 3 == 0, "the insert itself must be in frame"


def test_factor_xa_uses_the_standard_iegr_coding_sequence():

    factor_xa = CLEAVAGE_SITES["Factor Xa"]
    assert factor_xa == "ATCGAAGGTCGT"
    assert translate(factor_xa) == "IEGR"

    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False,
        cleavage_site="Factor Xa",
    )
    assert factor_xa in result.coding_region
    assert result.coding_region.index(factor_xa) % 3 == 0


def test_cassette_translates_to_the_expected_protein():

    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
        "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
        "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
        "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    coding = result.coding_region
    protein = "".join(
        codon_table.get(coding[i : i + 3], "X") for i in range(0, len(coding) - 2, 3)
    )
    assert protein.startswith("MGSSHHHHHHSSG"), protein[:20]
    assert "LVPRGS" in protein, "Thrombin recognition LVPR/GS"
    assert protein.endswith("*") or "*" in protein, "the insert's stop must survive"


def test_segments_describe_the_whole_forward_oligo():

    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    assert result.segments, "segments should be populated for display and reports"
    rebuilt = "".join(segment.sequence for segment in result.segments)
    assert rebuilt == result.forward
    cursor = 0
    for segment in result.segments:
        assert segment.start == cursor
        assert segment.end == cursor + len(segment.sequence)
        cursor = segment.end
