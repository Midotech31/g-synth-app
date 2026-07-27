"""Golden tests — the worked examples from the G-Synth specification.

These are the most important tests in the repository. They pin the exact
oligo sequences the lab orders. If one of them fails, the engine would send
someone to the bench with the wrong DNA: fix the code, never the expectation
— unless the protocol itself has genuinely changed, in which case the new
expectation needs a signed-off example to replace the old one.

Source: "Detailed and Rigorous Breakdown of the Logic for Each Tab in the
G-Synth Application", Sequence Design tab, Result Validation, Examples 1
and 2 — Non-Coding, NdeI/XhoI, Thrombin.
"""
import pytest

from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence

# ── Example 1 ────────────────────────────────────────────────────────────────
EX1_INPUT = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA"
EX1_FORWARD = (
    "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAAC"
)
EX1_REVERSE = (
    "TCGAGTTAGCCGCAGTAGTTTTCCAGCTGGTACAGGCTGCAGATGCTGGTGCAGCACTGTTCCACGATGCC"
    "AGAACCACGCGGCACCAGACCAGAAGAGTGGTGGTGGTGGTGGTGAGAAGAACCCA"
)

# ── Example 2 ────────────────────────────────────────────────────────────────
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
    """After annealing, the duplex must show 5'-TA (NdeI) and 5'-TCGA (XhoI).

    Checked structurally rather than by string comparison: the reverse
    oligo's own reverse complement is laid against the forward oligo, and the
    single-stranded ends that remain are the overhangs.
    """
    result = design_small_sequence(
        insert, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    forward = result.forward
    bottom_in_top_sense = reverse_complement(result.reverse)

    # NdeI leaves TATG on the forward and CA on the reverse: the duplex is
    # offset by len("TATG") - len("CA") = 2, leaving 5'-TA unpaired.
    offset = 2
    assert forward[:offset] == "TA", "NdeI 5' overhang"
    assert result.left_overhang == "TA"

    # The double-stranded core must pair perfectly.
    core_top = forward[offset:]
    core_bottom = bottom_in_top_sense[: len(core_top)]
    assert core_top == core_bottom, "forward and reverse must anneal without mismatches"

    # XhoI leaves TCGA hanging off the bottom strand at the right end.
    tail = bottom_in_top_sense[len(core_top):]
    assert reverse_complement(tail) == "TCGA", "XhoI 5' overhang"
    assert result.right_overhang == "TCGA"


def test_ndei_start_codon_overlaps_its_own_overhang():
    """NdeI's subtlety, pinned so nobody 'fixes' it later.

    The forward oligo starts TATG. The 5' overhang is TA (indices 0–1) and
    the start codon is ATG (indices 1–3) — the A is shared. Ligation restores
    CATATG because the cut vector supplies the CA.
    """
    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    assert result.forward.startswith("TATG")
    assert result.left_overhang == "TA"
    assert result.orf_start == 1
    assert result.coding_region.startswith("ATG")
    # The vector's CA plus the oligo's TATG must rebuild the NdeI site.
    assert "CA" + result.forward[:4] == "CATATG"


def test_his_tag_and_thrombin_site_are_present_in_frame():
    """The cassette must be translatable: ATG…linker…His…linker…LVPRGS…insert."""
    result = design_small_sequence(
        EX1_INPUT, enzyme_pair="NdeI / XhoI", is_coding=False, cleavage_site="Thrombin",
    )
    coding = result.coding_region
    assert coding.startswith("ATG")
    assert "CACCACCACCACCACCAC" in coding, "6×His"
    assert "CTGGTGCCGCGTGGTTCT" in coding, "Thrombin site"
    # Everything must sit in frame — offsets from the ATG divisible by 3.
    assert coding.index("CACCACCACCACCACCAC") % 3 == 0
    assert coding.index("CTGGTGCCGCGTGGTTCT") % 3 == 0
    assert coding.index(EX1_INPUT) % 3 == 0, "the insert itself must be in frame"


def test_cassette_translates_to_the_expected_protein():
    """Translate the coding region and check the tag actually reads MGSSHHHHHH…"""
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
    """The labelled segments must tile the forward oligo exactly, in order."""
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
