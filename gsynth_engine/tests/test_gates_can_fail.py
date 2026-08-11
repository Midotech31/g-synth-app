"""The checks that are supposed to stop bad work — deliberately broken.

A guard nothing has ever tripped is indistinguishable from a guard that does
not work. This project has already shipped one: a compatibility check in the
cloning module that passed for every input because it compared a label with
itself.

`AssemblyPlan.verify()` is the gate the whole application rests on — nothing
downloads until it returns empty. Every existing test asserts it returns
empty. None had ever made it speak.
"""
from dataclasses import replace

import pytest

from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.verify import Difference, verify_read

INSERT = ("GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
          "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCTAA")


@pytest.fixture
def plan():
    return design_merzoug_assembly(INSERT, target_oligo_length=60)


class TestVerifyCatchesWhatItPromises:
    """Break one thing at a time and check the right sentence comes back."""

    def test_a_clean_plan_says_nothing(self, plan):
        assert plan.verify() == []

    def test_a_forward_oligo_that_lost_a_base_is_caught(self, plan):
        """The commonest way a fragment builder goes wrong is an off-by-one
        slice, which leaves oligos that look entirely plausible."""
        first = plan.fragments[0]
        plan.fragments[0] = replace(first, forward=first.forward[:-1])

        problems = plan.verify()
        assert any("forward oligos do not reassemble" in p for p in problems)
        # The message must carry both lengths, or nobody can tell by how much.
        assert any("nt rebuilt vs" in p for p in problems)

    def test_a_reverse_oligo_that_does_not_match_its_strand_is_caught(self, plan):
        """The bottom strand is the half a design is most likely to get wrong,
        because it is the half nobody reads."""
        first = plan.fragments[0]
        plan.fragments[0] = replace(
            first, reverse=reverse_complement(reverse_complement(first.reverse)[:-1]),
        )
        assert any("reverse oligos do not reassemble the bottom strand" in p
                   for p in plan.verify())

    def test_two_fragments_that_do_not_share_a_junction_are_caught(self, plan):
        """This is the failure that produces a scrambled construct which runs
        correctly on a gel — the one that is only found at sequencing."""
        first = plan.fragments[0]
        wrong = "AAAA" if first.right_overhang != "AAAA" else "TTTT"
        plan.fragments[0] = replace(first, right_overhang=wrong)

        problems = plan.verify()
        assert any("do not share a junction" in p for p in problems)
        # Naming both fragments is the difference between a usable report and
        # a shrug: at 27 fragments "somewhere" is not an answer.
        assert any(f"Fragments {first.index} and" in p for p in problems)

    def test_the_gate_reports_every_break_at_once(self, plan):
        """Fixing one at a time and re-running is how a person loses an
        afternoon; the report has to be complete on the first pass."""
        first = plan.fragments[0]
        plan.fragments[0] = replace(first, forward=first.forward[:-1],
                                    right_overhang="AAAA")
        assert len(plan.verify()) >= 2


class TestDifferencesAreDescribedForAPerson:
    """`description` is the sentence the user actually reads. Only the
    substitution branch had a test; a read with an indel is routine."""

    DESIGN = ("ATGGCTAGCAAAGAACTGGTTACCGCTCTGTATCTGGTGTGCGGCGAACGCGGCTTTTTCTAC"
              "ACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGCGGC"
              "CCGGGCGCGGGCAGCCTGCAGCCGCTGGCGCTGGAAGGCAGCCTGCAGAAACGCGGCATCGTG")

    def test_an_inserted_base_is_named_as_inserted(self):
        read = list(self.DESIGN[20:180])
        read.insert(80, "G")
        result = verify_read(self.DESIGN, "".join(read), trim=0)

        inserted = [d for d in result.differences if d.kind == "insertion"]
        assert inserted, [d.description for d in result.differences]
        assert "inserted at position" in inserted[0].description

    def test_a_missing_base_is_named_as_missing(self):
        read = list(self.DESIGN[20:180])
        del read[80]
        result = verify_read(self.DESIGN, "".join(read), trim=0)

        deleted = [d for d in result.differences if d.kind == "deletion"]
        assert deleted, [d.description for d in result.differences]
        assert "missing at position" in deleted[0].description

    def test_positions_are_counted_from_one(self):
        """Biologists count from 1. Reporting a 0-based index sends someone
        to the wrong base, and the two look equally plausible."""
        d = Difference(kind="substitution", position=0, expected="A", found="G")
        assert "position 1" in d.description

    def test_a_silent_change_says_the_residue_is_unaffected(self):
        """Severity follows consequence: a silent change is real but harmless,
        and saying so is what stops it being chased."""
        d = Difference(kind="substitution", position=11, expected="T", found="C",
                       residue=4, from_residue="Leu", to_residue="Leu", silent=True)
        assert "silent" in d.description
        assert "Leu" in d.description


class TestTheOrderSheetAdvisesTheRightScale:
    """Past ~60 nt an oligo accumulates truncated products, and a truncated
    oligo in a ligation is a construct short of a base with nothing to show
    for it. Only the short branch had a test."""

    @pytest.mark.parametrize("target,expect_page", [(50, False), (90, True), (150, True)])
    def test_longer_oligos_are_sent_for_purification(self, target, expect_page):
        from gsynth_engine.protocol import order_sheet

        plan = design_merzoug_assembly(INSERT, target_oligo_length=target)
        rows = order_sheet(plan)
        purifications = {row.purification for row in rows}
        assert ("PAGE" in purifications) is expect_page, purifications
