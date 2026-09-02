import pytest

from gsynth_engine.gel import GEL_LADDERS, recommended_ladder, restriction_digest_sizes
from gsynth_engine.sequence import SequenceError


def test_generic_ladders_publish_all_marker_sizes():
    assert GEL_LADDERS["100-bp"]["bands"] == (
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500,
    )
    assert GEL_LADDERS["1-kb"]["bands"][-1] == 10000


def test_recommends_a_ladder_that_brackets_the_fragments():
    assert recommended_ladder([248, 900]) == "100-bp"
    assert recommended_ladder([1200, 5443]) == "1-kb"
    assert recommended_ladder([85, 5443]) == "broad-range"


def test_complete_circular_double_digest_returns_all_fragments():
    sequence = "A" * 10 + "AAGCTT" + "C" * 20 + "GAATTC" + "G" * 30
    sizes = restriction_digest_sizes(sequence, ["HindIII", "EcoRI"])

    assert len(sizes) == 2
    assert sum(sizes) == len(sequence)


def test_refuses_an_enzyme_that_does_not_cut():
    with pytest.raises(SequenceError, match="HindIII does not cut"):
        restriction_digest_sizes("ACGT" * 20, ["HindIII"])
