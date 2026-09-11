import random

import pytest
from Bio.Align import PairwiseAligner

from gsynth_engine import codon, esd, pcr, primers
from gsynth_engine.align import Scoring, align
from gsynth_engine.sequence import SequenceError
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.tests.test_cloning import clean_filler
from gsynth_engine.tests.test_preflight import GENE


def test_semiglobal_penalizes_missing_query_prefix():
    result = align('CCCCAAAA', 'AAAA', mode='semi-global', try_reverse=False,
        scoring=Scoring(gap_open=2, gap_extend=1))
    oracle = PairwiseAligner(mode='global', match_score=5, mismatch_score=-4, open_gap_score=-3, extend_gap_score=-1)
    oracle.left_insertion_score = oracle.right_insertion_score = 0
    assert result.score == oracle.score('CCCCAAAA', 'AAAA') == 14
    assert result.top == 'CCCCAAAA'
    assert result.bottom == '----AAAA'


@pytest.mark.parametrize('mode', ['global', 'local', 'semi-global'])
@pytest.mark.parametrize('seed', range(20))
def test_alignment_score_matches_independent_implementation(mode, seed):
    rng = random.Random(seed)
    first = ''.join(rng.choices('ACGT', k=20))
    second = ''.join(rng.choices('ACGT', k=25))
    oracle = PairwiseAligner(mode='local' if mode == 'local' else 'global', match_score=5,
        mismatch_score=-4, open_gap_score=-11, extend_gap_score=-1)
    if mode == 'semi-global':
        oracle.left_insertion_score = oracle.right_insertion_score = 0
    result = align(first, second, mode=mode, try_reverse=False)
    assert result.score == oracle.score(first, second)


@pytest.mark.parametrize('codons,protein,index,expected', [
    (['GCT'], 'A', -1, False), (['ATG'], 'M', 0, False),
    (['TTC'], 'F', 0, False), (['GCT'], 'A', 0, True),
])
def test_synonymous_swap_preserves_translation(codons, protein, index, expected):
    assert codon._swap(codons, index, protein, codon.ECOLI, codon.Constraints(avoid_rare=False), random.Random(0)) is expected
    assert codon.translate(''.join(codons)) == protein


def test_synonymous_swap_removes_forbidden_motif():
    codons = ['GCT']
    assert codon._swap(codons, 0, 'A', codon.ECOLI, codon.Constraints(avoid_motifs=('GCT',)), random.Random(0))
    assert codons != ['GCT'] and codon.translate(''.join(codons)) == 'A'


def test_empty_insert_after_start_and_stop_removal_is_rejected():
    with pytest.raises(SequenceError, match='left nothing'):
        design_small_sequence('ATGTAA', is_coding=True, remove_stop=True)


def test_primer_candidates_with_ambiguous_bases_are_excluded():
    assert primers._pick('N' * 100, anchor=50, direction=1, circular=False,
        tm_min=0, tm_max=100, tm_target=50, length_range=(18, 20), search=2) is None


def test_pcr_tm_difference_and_dimer_are_reported():
    result = pcr.design_pcr('A' * 40 + clean_filler(100) + 'G' * 40)
    assert any('differ by' in warning for warning in result.warnings)
    complementary = pcr.design_pcr('ACGTTGCAAGGCTTAGCCAT' + clean_filler(100) + 'ACGTTGCAAGGCTTAGCCAT')
    assert any('complementary over' in warning for warning in complementary.warnings)


def test_custom_primer_frame_shift_is_reported():
    result = pcr.design_pcr(GENE, left_enzyme='BamHI', right_enzyme='XhoI', keep_frame=True)
    edited = pcr.design_pcr(GENE, left_enzyme='BamHI', right_enzyme='XhoI', keep_frame=True,
        forward_primer=result.forward.tail + 'A' + result.forward.anneals,
        reverse_primer=result.reverse.sequence)
    assert any('out of the requested' in warning for warning in edited.warnings)


def test_junction_search_rejects_boundary_and_overlapping_candidates():
    with pytest.raises(SequenceError, match='Could not place junction'):
        esd._choose_junctions('ACGT' * 8, count=3, overhang_length=8,
            ds_start=0, ds_end=12, forbidden=set(), search_window=12)
    with pytest.raises(SequenceError, match='Could not place junction'):
        esd._choose_junctions('AAGCTGAC', count=3, overhang_length=4,
            ds_start=0, ds_end=30, forbidden=set(), search_window=2)


def test_junction_search_exhaustion_is_distinct_from_low_complexity(monkeypatch):
    def unavailable(*args, **kwargs):
        raise SequenceError('No compatible junction')
    monkeypatch.setattr(esd, '_choose_junctions', unavailable)
    with pytest.raises(SequenceError, match='Tried'):
        esd._place_junctions(clean_filler(500), count=1, overhang_length=4,
            ds_start=0, ds_end=500, forbidden=set(), search_window=2)
