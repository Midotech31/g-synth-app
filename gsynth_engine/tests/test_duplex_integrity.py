"""Every enzyme must produce a duplex that actually anneals.

G-Synth 2.x stored one pair of "what's left after the cut" values per enzyme,
but what the forward and reverse oligos must carry depends on whether the
enzyme sits at the left or the right end of the insert. The stored values
happened to be right for NdeI-as-left and XhoI-as-right — the validated
default pair — and wrong for everything else, which produced duplexes with a
mismatch at the junction. The oligos would not have annealed cleanly.

These tests check the property directly, for every enzyme, in both roles.
"""
import pytest

from gsynth_engine.constants import RESTRICTION_ENZYMES, left_remainders, overhang
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence

INSERT = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGTAA"
ENZYMES = sorted(RESTRICTION_ENZYMES)


def duplex_layout(result):
    """Lay the two oligos against each other and return the paired region.

    Returns (offset, top, bottom_in_top_sense).
    """
    forward_remainder, reverse_remainder = left_remainders(result.left_enzyme)
    offset = len(forward_remainder) - len(reverse_remainder)
    return offset, result.forward, reverse_complement(result.reverse)


@pytest.mark.parametrize("enzyme", ENZYMES)
def test_left_position_duplex_anneals_without_mismatch(enzyme):
    result = design_small_sequence(
        INSERT, enzyme_pair=f"{enzyme} / XhoI", cleavage_site=None,
        include_his_tag=False, include_linkers=False,
    )
    offset, top, bottom = duplex_layout(result)
    start, end = max(0, offset), min(len(top), offset + len(bottom))
    assert end > start, f"{enzyme}: the two oligos do not overlap at all"
    assert top[start:end] == bottom[start - offset : end - offset], (
        f"{enzyme} at the left end produces a mismatched duplex — the oligos "
        f"would not anneal"
    )


@pytest.mark.parametrize("enzyme", ENZYMES)
def test_right_position_duplex_anneals_without_mismatch(enzyme):
    result = design_small_sequence(
        INSERT, enzyme_pair=f"NdeI / {enzyme}", cleavage_site=None,
        include_his_tag=False, include_linkers=False,
    )
    offset, top, bottom = duplex_layout(result)
    start, end = max(0, offset), min(len(top), offset + len(bottom))
    assert end > start, f"{enzyme}: the two oligos do not overlap at all"
    assert top[start:end] == bottom[start - offset : end - offset], (
        f"{enzyme} at the right end produces a mismatched duplex"
    )


@pytest.mark.parametrize("enzyme", ENZYMES)
def test_ligation_restores_the_recognition_site(enzyme):
    """The point of the remainders: the site must reappear after ligation.

    The vector supplies whatever the insert does not, so remainder + counter-
    remainder must rebuild the recognition sequence exactly.
    """
    site = str(RESTRICTION_ENZYMES[enzyme]["recognition"])
    top_cut = int(RESTRICTION_ENZYMES[enzyme]["cut_top"])  # type: ignore[call-overload]
    forward_remainder, _ = left_remainders(enzyme)
    vector_side = site[:top_cut]
    assert vector_side + forward_remainder == site


@pytest.mark.parametrize("enzyme", ENZYMES)
def test_overhang_matches_the_cut_geometry(enzyme):
    """A 5' overhang means the top is cut before the bottom, and vice versa."""
    top_cut = int(RESTRICTION_ENZYMES[enzyme]["cut_top"])      # type: ignore[call-overload]
    bottom_cut = int(RESTRICTION_ENZYMES[enzyme]["cut_bottom"])  # type: ignore[call-overload]
    sequence, kind = overhang(enzyme)
    if top_cut == bottom_cut:
        assert kind == "blunt" and sequence == ""
    elif top_cut < bottom_cut:
        assert kind == "5'" and len(sequence) == bottom_cut - top_cut
    else:
        assert kind == "3'" and len(sequence) == top_cut - bottom_cut


def test_ndei_xhoi_overhangs_are_the_documented_ones():
    """The validated pair, stated explicitly so a table edit cannot drift."""
    assert overhang("NdeI") == ("TA", "5'")
    assert overhang("XhoI") == ("TCGA", "5'")
    assert left_remainders("NdeI") == ("TATG", "CA")
