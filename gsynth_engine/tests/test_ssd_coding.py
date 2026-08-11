"""The `is_coding=True` path — an insert that is already a complete gene.

Every other test in this suite passes `is_coding=False`, which builds the
expression cassette. But the interface offers a checkbox — "Insert already
has its own ATG" — and ticking it takes a different branch entirely: no
cassette, no tag, no linkers, just the gene between the two sticky ends.

That branch removes a start codon and can remove a stop codon. Both change
the protein. It was reachable from the interface with no test behind it.
"""
import pytest

from gsynth_engine.constants import overhang
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.ssd import design_small_sequence

GENE = ("ATGGCTAGCAAAGAACTGGTTACCGCTCTGTATCTGGTGTGCGGCGAACGCGGCTTTTTCTAC"
        "ACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGTAA")


class TestNoCassetteIsAdded:
    """The point of the flag: the gene is already complete."""

    def test_the_construct_is_only_the_ends_and_the_insert(self):
        result = design_small_sequence(GENE, enzyme_pair="BamHI / EcoRI",
                                       is_coding=True)
        assert [s.name for s in result.segments] == [
            "BamHI overhang", "insert", "EcoRI overhang",
        ]

    def test_no_tag_linker_or_protease_site_is_introduced(self):
        """Adding a 6×His to a gene that already carries its own gives the
        protein two, and the second is not where anyone expects it."""
        result = design_small_sequence(GENE, enzyme_pair="BamHI / EcoRI",
                                       is_coding=True, include_his_tag=True,
                                       include_linkers=True,
                                       cleavage_site="Thrombin")
        assert "CACCACCACCACCACCAC" not in result.forward, "a His tag was added"
        assert result.forward.count(GENE) == 1

    def test_the_insert_survives_verbatim(self):
        result = design_small_sequence(GENE, enzyme_pair="BamHI / EcoRI",
                                       is_coding=True)
        assert GENE in result.forward


class TestTheStartCodon:
    """NdeI's site is CA^TATG — it supplies the ATG itself."""

    def test_ndeI_removes_the_inserts_own_atg(self):
        """Keeping both gives ATG-ATG: one extra methionine, in frame, and
        the protein is a residue longer than the one that was designed."""
        result = design_small_sequence(GENE, enzyme_pair="NdeI / XhoI",
                                       is_coding=True)
        assert GENE[3:] in result.forward
        assert not result.forward.startswith("TATGATG")
        assert any("ATG was removed" in w for w in result.warnings)

    def test_any_other_enzyme_keeps_it(self):
        """Only NdeI supplies a start codon. Stripping the ATG for BamHI
        would leave the gene with no way to begin."""
        result = design_small_sequence(GENE, enzyme_pair="BamHI / EcoRI",
                                       is_coding=True)
        assert GENE in result.forward
        assert not any("ATG was removed" in w for w in result.warnings)

    def test_a_gene_not_starting_with_atg_is_left_alone(self):
        gene = GENE[3:]
        result = design_small_sequence(gene, enzyme_pair="NdeI / XhoI",
                                       is_coding=True)
        assert gene in result.forward
        assert not any("ATG was removed" in w for w in result.warnings)


class TestTheStopCodon:
    """Removing it is how a C-terminal vector tag gets translated."""

    def test_removing_the_stop_takes_the_last_in_frame_one(self):
        result = design_small_sequence(GENE, enzyme_pair="BamHI / EcoRI",
                                       is_coding=True, remove_stop=True)
        assert GENE[:-3] in result.forward
        assert GENE not in result.forward

    def test_a_gene_without_a_stop_says_so_rather_than_failing(self):
        """The user may have trimmed it already; that is not an error."""
        result = design_small_sequence(GENE[:-3], enzyme_pair="BamHI / EcoRI",
                                       is_coding=True, remove_stop=True)
        assert any("No in-frame stop" in w for w in result.warnings)

    def test_a_gene_that_is_only_a_stop_codon_is_refused(self):
        """Nothing left to clone, and the message has to say why rather than
        hand back a construct made of two sticky ends."""
        with pytest.raises(SequenceError, match="reading frame"):
            design_small_sequence("TAA", enzyme_pair="BamHI / EcoRI",
                                  is_coding=True, remove_stop=True)


class TestTheInvariantsStillHold:
    """This branch must not be exempt from the properties everything else is
    held to — it is the branch that reaches the vector unmodified."""

    def test_both_strands_pair(self):
        result = design_small_sequence(GENE, enzyme_pair="NdeI / XhoI",
                                       is_coding=True)
        offset = len("TATG") - len("CA")          # NdeI's left remainders
        bottom = reverse_complement(result.reverse)
        overlap = result.forward[offset:offset + len(bottom)]
        assert overlap == bottom[:len(overlap)]

    @pytest.mark.parametrize("pair", ["NdeI / XhoI", "BamHI / EcoRI", "KpnI / SacI"])
    def test_the_terminal_ends_are_what_the_enzymes_leave(self, pair):
        left, right = (e.strip() for e in pair.split("/"))
        plan = design_merzoug_assembly(GENE, enzyme_pair=pair, is_coding=True,
                                       target_oligo_length=90)
        assert plan.terminal_ends == (overhang(left), overhang(right))

    @pytest.mark.parametrize("pair", ["NdeI / XhoI", "BamHI / EcoRI"])
    def test_the_fragments_still_re_ligate_into_the_construct(self, pair):
        plan = design_merzoug_assembly(GENE, enzyme_pair=pair, is_coding=True,
                                       target_oligo_length=60)
        assert plan.verify() == []
        assert "".join(f.forward for f in plan.fragments) == plan.construct_forward
