"""Tests for the design endpoints.

The engine's own suite proves the biology. These tests prove the HTTP layer
does not corrupt it: the API must return exactly what the engine computed,
guard access, save designs to the right owner, and turn engine errors into
messages a user can act on.
"""
import csv
import io

import pytest
from django.urls import reverse
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.ssd import design_small_sequence

from apps.projects.models import Project

# The specification's Example 1 — the API must reproduce it end to end.
GOLDEN_INSERT = "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA"
GOLDEN_FORWARD = (
    "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAAC"
)
GOLDEN_REVERSE = (
    "TCGAGTTAGCCGCAGTAGTTTTCCAGCTGGTACAGGCTGCAGATGCTGGTGCAGCACTGTTCCACGATGCC"
    "AGAACCACGCGGCACCAGACCAGAAGAGTGGTGGTGGTGGTGGTGAGAAGAACCCA"
)

LONG_INSERT = (
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGTAA"
)


@pytest.mark.django_db
class TestEnzymeCatalogue:
    def test_is_public(self, api_client):
        """A reference table — the UI may need it before sign-in."""
        response = api_client.get(reverse("design-enzymes"))
        assert response.status_code == 200

    def test_lists_enzymes_with_their_overhangs(self, api_client):
        data = api_client.get(reverse("design-enzymes")).data
        by_name = {e["name"]: e for e in data["enzymes"]}
        assert by_name["NdeI"]["overhang"] == "TA"
        assert by_name["NdeI"]["overhang_type"] == "5'"
        assert by_name["NdeI"]["supplies_start_codon"] is True
        assert by_name["XhoI"]["overhang"] == "TCGA"
        assert by_name["KpnI"]["overhang_type"] == "3'"
        assert by_name["SmaI"]["overhang_type"] == "blunt"

    def test_offers_cleavage_sites_and_common_pairs(self, api_client):
        data = api_client.get(reverse("design-enzymes")).data
        assert "NdeI / XhoI" in data["common_pairs"]
        names = {c["name"] for c in data["cleavage_sites"]}
        assert {"Thrombin", "TEV", "Factor Xa"} <= names


@pytest.mark.django_db
class TestSSDEndpoint:
    url_name = "design-ssd"

    def test_requires_authentication(self, api_client):
        response = api_client.post(reverse(self.url_name), {"sequence": GOLDEN_INSERT})
        assert response.status_code == 401

    def test_reproduces_the_specification_example(self, auth_client):
        """The golden example, through HTTP."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": GOLDEN_INSERT,
            "left_enzyme": "NdeI", "right_enzyme": "XhoI",
            "cleavage_site": "Thrombin", "is_coding": False,
        })
        assert response.status_code == 200, response.data
        assert response.data["forward"] == GOLDEN_FORWARD
        assert response.data["reverse"] == GOLDEN_REVERSE
        assert response.data["left_overhang"] == "TA"
        assert response.data["right_overhang"] == "TCGA"

    def test_matches_the_engine_exactly(self, auth_client):
        """The API must not reshape, round or truncate the design."""
        expected = design_small_sequence(GOLDEN_INSERT)
        response = auth_client.post(reverse(self.url_name), {"sequence": GOLDEN_INSERT})
        assert response.data["forward"] == expected.forward
        assert response.data["reverse"] == expected.reverse
        assert response.data["orf_start"] == expected.orf_start

    def test_returns_labelled_segments_for_display(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": GOLDEN_INSERT})
        names = [s["name"] for s in response.data["segments"]]
        assert "6×His tag" in names
        assert any("Thrombin" in name for name in names)
        rebuilt = "".join(s["sequence"] for s in response.data["segments"])
        assert rebuilt == response.data["forward"]

    def test_saves_as_a_project_when_asked(self, auth_client, user):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": GOLDEN_INSERT, "name": "EntA construct", "save_as_project": True,
        })
        assert "project_id" in response.data
        project = Project.objects.get(id=response.data["project_id"])
        assert project.user == user
        assert project.module == "ssd"
        assert project.sequence == GOLDEN_FORWARD

    def test_does_not_save_by_default(self, auth_client):
        auth_client.post(reverse(self.url_name), {"sequence": GOLDEN_INSERT})
        assert Project.objects.count() == 0

    def test_rejects_identical_enzymes(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": GOLDEN_INSERT, "left_enzyme": "NdeI", "right_enzyme": "NdeI",
        })
        assert response.status_code == 400
        assert "differ" in str(response.data)

    def test_invalid_bases_produce_a_usable_message(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": "ATGXYZ"})
        assert response.status_code == 400
        assert "not A, C, G or T" in response.data["detail"]


@pytest.mark.django_db
class TestAssemblyEndpoint:
    url_name = "design-assembly"

    def test_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {"sequence": LONG_INSERT}).status_code == 401

    def test_returns_a_verified_plan(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "target_oligo_length": 90, "overhang_length": 4,
        })
        assert response.status_code == 200, response.data
        # Empty verification means the oligos re-ligate to the design.
        assert response.data["verification"] == []
        assert response.data["fragment_count"] >= 2
        assert response.data["oligo_count"] == 2 * response.data["fragment_count"]

    def test_fragments_reassemble_into_the_construct(self, auth_client):
        """The property the whole method rests on, checked over HTTP."""
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        rebuilt = "".join(f["forward"] for f in response.data["fragments"])
        assert rebuilt == response.data["construct_forward"]

    def test_junctions_are_the_requested_length_and_unique(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "overhang_length": 6, "target_oligo_length": 60,
        })
        junctions = response.data["junction_overhangs"]
        assert junctions
        assert all(len(j) == 6 for j in junctions)
        assert len(junctions) == len(set(junctions))

    def test_terminal_ends_carry_the_enzyme_overhangs(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        fragments = response.data["fragments"]
        assert fragments[0]["left_overhang"] == "TA"
        assert fragments[-1]["right_overhang"] == "TCGA"

    def test_includes_the_order_sheet(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "name": "pGS-EntA",
        })
        oligos = response.data["oligos"]
        assert len(oligos) == response.data["oligo_count"]
        assert all(row["Name"].startswith("pGS-EntA_") for row in oligos)

    def test_matches_the_engine_exactly(self, auth_client):
        expected = design_merzoug_assembly(LONG_INSERT, target_oligo_length=90)
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "target_oligo_length": 90,
        })
        assert response.data["construct_forward"] == expected.construct_forward
        assert [f["forward"] for f in response.data["fragments"]] == [
            f.forward for f in expected.fragments
        ]

    def test_saves_as_a_project_when_asked(self, auth_client, user):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "name": "EntA assembly", "save_as_project": True,
        })
        project = Project.objects.get(id=response.data["project_id"])
        assert project.user == user
        assert project.module == "merzoug_assembly"
        assert project.data["fragment_count"] == response.data["fragment_count"]

    def test_rejects_an_overhang_outside_the_method(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "overhang_length": 2,
        })
        assert response.status_code == 400


@pytest.mark.django_db
class TestDownloads:
    def test_order_sheet_is_csv(self, auth_client):
        response = auth_client.post(reverse("design-order-sheet"), {
            "sequence": LONG_INSERT, "name": "pGS-EntA",
        })
        assert response.status_code == 200
        assert response["Content-Type"].startswith("text/csv")
        assert "pGS-EntA_oligos.csv" in response["Content-Disposition"]
        rows = list(csv.DictReader(io.StringIO(response.content.decode())))
        assert rows and rows[0]["Name"].startswith("pGS-EntA_F1")
        assert "Sequence (5'->3')" in rows[0]

    def test_protocol_is_text_and_states_the_method(self, auth_client):
        response = auth_client.post(reverse("design-protocol"), {
            "sequence": LONG_INSERT, "name": "pGS-EntA",
        })
        assert response.status_code == 200
        text = response.content.decode()
        assert "MERZOUG ASSEMBLY" in text
        assert "No PCR" in text
        assert "PHOSPHORYLATION" in text
        assert "pGS-EntA_protocol.txt" in response["Content-Disposition"]

    def test_downloads_require_authentication(self, api_client):
        for name in ("design-order-sheet", "design-protocol"):
            assert api_client.post(reverse(name), {"sequence": LONG_INSERT}).status_code == 401
