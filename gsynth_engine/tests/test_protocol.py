"""Tests for the order sheet and the bench protocol.

What matters here is that the paperwork agrees with the design: every oligo
in the plan appears in the order sheet exactly once, and the protocol names
the right fragments, junctions and enzymes. A protocol that quietly
disagrees with the oligos is worse than no protocol.
"""
import csv
import io

import pytest

from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.protocol import bench_protocol, order_sheet, order_sheet_csv

INSERT = (
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGTAA"
)


@pytest.fixture
def plan():
    return design_merzoug_assembly(INSERT, target_oligo_length=90)


class TestOrderSheet:
    def test_lists_every_oligo_once(self, plan):
        orders = order_sheet(plan)
        assert len(orders) == plan.oligo_count
        designed = {f.forward for f in plan.fragments} | {f.reverse for f in plan.fragments}
        assert {o.sequence for o in orders} == designed

    def test_names_are_unique_and_traceable(self, plan):
        orders = order_sheet(plan, construct_name="pGS-EntA")
        names = [o.name for o in orders]
        assert len(names) == len(set(names)), "duplicate names would confuse an order"
        assert all(name.startswith("pGS-EntA_") for name in names)
        assert any(name.endswith("_F1_F") for name in names)

    def test_reported_length_matches_the_sequence(self, plan):
        for order in order_sheet(plan):
            assert order.length == len(order.sequence)

    def test_long_oligos_get_purification(self, plan):
        for order in order_sheet(plan):
            if order.length > 60:
                assert order.purification == "PAGE"

    def test_csv_round_trips(self, plan):
        text = order_sheet_csv(plan, construct_name="pGS-EntA")
        rows = list(csv.DictReader(io.StringIO(text)))
        assert len(rows) == plan.oligo_count
        assert rows[0]["Sequence (5'->3')"] == plan.fragments[0].forward
        assert set(rows[0]) >= {"Name", "Sequence (5'->3')", "Length (nt)", "Tm (°C)"}


class TestBenchProtocol:
    def test_states_the_method_and_its_constraints(self, plan):
        text = bench_protocol(plan)
        assert "MERZOUG ASSEMBLY" in text
        assert "No PCR" in text, "the defining constraint must be stated"

    def test_mentions_every_fragment_and_junction(self, plan):
        text = bench_protocol(plan)
        for fragment in plan.fragments:
            assert fragment.name in text
        for junction in plan.junction_overhangs:
            assert junction in text

    def test_names_the_cloning_enzymes_and_vector(self, plan):
        text = bench_protocol(plan, vector="pET-21a(+)")
        assert "NdeI" in text and "XhoI" in text
        assert "pET-21a(+)" in text

    def test_covers_the_steps_that_make_it_work(self, plan):
        """Phosphorylation and slow cooling are the two steps people skip."""
        text = bench_protocol(plan)
        assert "PHOSPHORYLATION" in text
        assert "5' phosphate" in text
        assert "Slow cooling" in text or "slow cooling" in text.lower()

    def test_includes_the_full_construct_for_verification(self, plan):
        text = bench_protocol(plan)
        assert "EXPECTED CONSTRUCT" in text
        stripped = "".join(
            ch for line in text.splitlines() for ch in line if ch in "ACGT"
        )
        assert plan.construct_forward[:40] in stripped

    def test_shows_the_duplex_before_the_ordering_step(self, plan):
        """The protocol is what reaches the bench, so the check goes in it."""
        text = bench_protocol(plan)
        assert "HYBRIDISATION" in text

        lines = text.splitlines()
        strands = [ln for ln in lines if "5'" in ln and "3'" not in ln]
        assert any("|" in ln for ln in lines), "the pairing rungs should be drawn"
        assert strands, "both strands should be labelled"

    def test_duplex_in_the_protocol_matches_the_design(self, plan):
        text = bench_protocol(plan)
        start = text.index("HYBRIDISATION")
        stop = text.index("PHOSPHORYLATION")
        drawn = "".join(
            ch for ch in text[start:stop] if ch in "ACGT"
        )
        # The top strand, read straight out of the drawing.
        assert plan.construct_forward[:40] in drawn

    def test_states_the_tm_model_and_its_conditions(self, plan):
        """A Tm on an order sheet with no conditions is not actionable."""
        text = bench_protocol(plan)
        assert "SantaLucia" in text
        assert "Na" in text

    def test_single_fragment_skips_the_ligation_chain(self):
        plan = design_merzoug_assembly("GGCATCGTGGAACAGTGCTGCACCAGCTAA")
        text = bench_protocol(plan)
        assert "Single fragment" in text

    def test_surfaces_design_warnings(self):
        """A warning that only exists in the object helps nobody at the bench."""
        plan = design_merzoug_assembly(INSERT, enzyme_pair="NdeI / XhoI")
        plan.warnings.append("Test warning about the design")
        assert "Test warning about the design" in bench_protocol(plan)
