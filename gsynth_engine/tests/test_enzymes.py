"""Tests for the restriction tables.

Two sets, for two different questions.

`RESTRICTION_ENZYMES` is what this lab keeps in the freezer — the list a
cloning pair is *picked* from, and the one whose golden designs are pinned.
It is typed by hand, which is why the first test here checks it against
REBASE: a transposed cut position is invisible, and would silently produce
oligos with the wrong overhang for one enzyme while every other test passes.

`ALL_ENZYMES` answers a different question — "what else cuts here" — and
breadth is the point. pET-21a has sixty single-cutters; a plasmid map drawn
from nineteen shows fourteen of them, and a diagnostic digest planned
against that map can fail on a site nobody was ever shown.
"""
import pytest

from gsynth_engine.constants import (
    ALL_ENZYMES,
    RESTRICTION_ENZYMES,
    left_remainders,
    overhang,
    right_remainders,
)
from gsynth_engine.sequence import reverse_complement

biopython = pytest.importorskip("Bio.Restriction", reason="regeneration source")


class TestTheCuratedSetIsCorrect:
    """Hand-typed data, checked against the authority it came from."""

    @pytest.mark.parametrize("name", sorted(RESTRICTION_ENZYMES))
    def test_geometry_matches_rebase(self, name):
        """Recognition sequence and both cut positions, against Biopython.

        This is the whole reason the curated set is allowed to be hand-typed:
        it is verified against REBASE on every run. Without this, a single
        transposed digit gives one enzyme the wrong overhang for ever.
        """
        enzyme = getattr(biopython, name)
        site = str(enzyme.site)
        ours = RESTRICTION_ENZYMES[name]

        assert ours["recognition"] == site
        assert ours["cut_top"] == enzyme.fst5
        assert ours["cut_bottom"] == len(site) + enzyme.fst3


class TestTheWideTable:
    def test_contains_every_curated_enzyme(self):
        """The two must not diverge; the curated names are the canonical ones."""
        missing = sorted(set(RESTRICTION_ENZYMES) - set(ALL_ENZYMES))
        assert missing == []
        for name, spec in RESTRICTION_ENZYMES.items():
            for key in ("recognition", "cut_top", "cut_bottom"):
                assert ALL_ENZYMES[name][key] == spec[key], name

    def test_is_substantially_wider_than_the_freezer(self):
        """If this collapses back to nineteen the generator has silently
        stopped running and every plasmid map is under-annotated again."""
        assert len(ALL_ENZYMES) >= 100

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_every_entry_is_a_real_type_IIP_specification(self, name):
        """No ambiguity codes, and the cut falls inside the site.

        The scanner matches literally, so an N in a site reports positions
        the enzyme does not cut. And the remainders are slices of the site,
        so a Type IIS enzyme cutting past its own end yields empty strings —
        wrong ends rather than an error, which is worse.
        """
        spec = ALL_ENZYMES[name]
        site = str(spec["recognition"])
        top, bottom = int(spec["cut_top"]), int(spec["cut_bottom"])

        assert site and not set(site) - set("ACGT"), site
        assert 0 <= top <= len(site)
        assert 0 <= bottom <= len(site)

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_the_derived_ends_reassemble_the_site(self, name):
        """Cutting and re-joining must give the recognition sequence back.

        This is the property the whole cut model rests on: what each oligo
        carries is *derived* from the cut positions per role. If the
        derivation is wrong for an enzyme, its two remainders will not
        reconstitute the site it came from.
        """
        site = str(ALL_ENZYMES[name]["recognition"])
        left_fwd, left_rev = left_remainders(name)
        right_fwd, right_rev = right_remainders(name)

        assert right_fwd + left_fwd == site
        assert left_rev + right_rev == reverse_complement(site)

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_overhang_agrees_with_the_cut_positions(self, name):
        """`overhang()` is what the compatibility checks compare against, so
        it must be derived from the same two numbers and not stored."""
        spec = ALL_ENZYMES[name]
        site = str(spec["recognition"])
        top, bottom = int(spec["cut_top"]), int(spec["cut_bottom"])
        sequence, kind = overhang(name)

        if top == bottom:
            assert (sequence, kind) == ("", "blunt")
        else:
            low, high = min(top, bottom), max(top, bottom)
            assert sequence == site[low:high]
            assert kind == ("5'" if top < bottom else "3'")


class TestIsoschizomersAreCollapsed:
    def test_one_name_per_cut_specification(self):
        """Five names for one site turns a plasmid map into noise."""
        seen: dict[tuple, str] = {}
        for name, spec in ALL_ENZYMES.items():
            key = (spec["recognition"], spec["cut_top"], spec["cut_bottom"])
            assert key not in seen, f"{name} duplicates {seen.get(key)}"
            seen[key] = name

    def test_the_familiar_name_survives(self):
        """NheI, not AsuNHI — the canonical choice must favour the name a
        biologist will recognise, or the map is unreadable."""
        assert ALL_ENZYMES["NheI"]["recognition"] == "GCTAGC"
        assert "AsuNHI" not in ALL_ENZYMES, "a true isoschizomer of NheI"

    def test_a_neoschizomer_is_not_collapsed_into_its_twin(self):
        """Same site, opposite cut, and both must survive.

        NheI and BmtI both recognise GCTAGC. NheI leaves G^CTAG_C — a 5'
        overhang; BmtI leaves G_CTAG^C — a 3' one. Grouping by recognition
        sequence would merge them and hand back the wrong sticky end for
        whichever name lost, which is the same class of error as reading
        polarity off the strand instead of the cut.
        """
        nhe, bmt = ALL_ENZYMES["NheI"], ALL_ENZYMES["BmtI"]
        assert nhe["recognition"] == bmt["recognition"] == "GCTAGC"
        assert (nhe["cut_top"], nhe["cut_bottom"]) == (1, 5)
        assert (bmt["cut_top"], bmt["cut_bottom"]) == (5, 1)
        assert overhang("NheI") == ("CTAG", "5'")
        assert overhang("BmtI") == ("CTAG", "3'")
