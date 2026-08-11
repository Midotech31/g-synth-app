"""The PCR endpoint.

What these assert is that the response carries what the engine computed —
the same primers, the same product, the same ends — rather than something
that merely looks like a PCR. The biology itself is pinned in
`gsynth_engine/tests/test_pcr.py`; duplicating those claims here would mean
two places to update and one of them silently going stale.
"""
import pytest
from django.urls import reverse

from gsynth_engine.pcr import design_pcr

GENE = (
    "ATGAAAGGTGAAGAATTGTTCACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTT"
    "AACGGTCACAAATTCTCTGTTTCTGGTGAAGGTGAAGGTGATGCTACCTACGGTAAACTGACC"
    "CTGAAA"
)


@pytest.mark.django_db
class TestAccess:
    def test_signed_out_is_refused(self, api_client):
        r = api_client.post(reverse("design-pcr"), {"template": GENE}, format="json")
        assert r.status_code == 401


@pytest.mark.django_db
class TestRequestIsBounded:
    def test_a_missing_template_is_rejected(self, auth_client):
        r = auth_client.post(reverse("design-pcr"), {}, format="json")
        assert r.status_code == 400

    def test_an_oversized_template_is_rejected(self, auth_client):
        r = auth_client.post(
            reverse("design-pcr"), {"template": "A" * 200_001}, format="json"
        )
        assert r.status_code == 400

    def test_one_enzyme_alone_is_rejected(self, auth_client):
        """One site cannot open both ends of a product."""
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI"},
            format="json",
        )
        assert r.status_code == 400
        assert "right_enzyme" in r.data

    def test_an_unknown_enzyme_is_rejected(self, auth_client):
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NotAnEnzyme", "right_enzyme": "XhoI"},
            format="json",
        )
        assert r.status_code == 400

    def test_a_backwards_region_is_rejected(self, auth_client):
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "target_start": 90, "target_end": 30},
            format="json",
        )
        assert r.status_code == 400

    def test_an_unorderable_clamp_is_rejected(self, auth_client):
        """The clamp is synthesised into every primer, so it is bounded."""
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI",
             "clamp": 500},
            format="json",
        )
        assert r.status_code == 400

    def test_a_region_too_short_becomes_a_readable_400(self, auth_client):
        """`SequenceError` messages are shown to the user verbatim, so the
        endpoint must pass one through rather than raise a 500."""
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "target_start": 0, "target_end": 20},
            format="json",
        )
        assert r.status_code == 400
        assert "too short" in r.data["detail"]


@pytest.mark.django_db
class TestTheResponseMatchesTheEngine:
    def test_conventional_pcr_returns_the_engine_primers(self, auth_client):
        expected = design_pcr(GENE)
        r = auth_client.post(reverse("design-pcr"), {"template": GENE}, format="json")

        assert r.status_code == 200
        assert r.data["forward"]["sequence"] == expected.forward.sequence
        assert r.data["reverse"]["sequence"] == expected.reverse.sequence
        assert r.data["product"] == expected.product
        assert r.data["annealing_temperature"] == expected.annealing_temperature
        assert r.data["digest"] is None

    def test_both_tm_figures_survive_serialisation(self, auth_client):
        """Collapsing the two into one is what makes a tailed primer look as
        though it should anneal ten degrees hotter than it does."""
        expected = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI"},
            format="json",
        )
        forward = r.data["forward"]
        assert forward["tm"] == expected.forward.tm
        assert forward["tm_full"] == expected.forward.tm_full
        assert forward["tm_full"] > forward["tm"]
        assert r.data["annealing_temperature"] == expected.annealing_temperature

    def test_the_tail_is_returned_separately_from_the_annealing_part(self, auth_client):
        """The interface draws them differently, so it needs them apart —
        and the two must still concatenate to the oligo that gets ordered."""
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI"},
            format="json",
        )
        forward = r.data["forward"]
        assert forward["tail"] + forward["anneals"] == forward["sequence"]
        assert forward["tail"].endswith("CATATG")

    def test_the_digest_carries_the_ends_the_engine_measured(self, auth_client):
        expected = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI")
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI"},
            format="json",
        )
        digest = r.data["digest"]
        assert digest["top"] == expected.digest.top
        assert digest["bottom"] == expected.digest.bottom
        assert digest["left_end"]["sequence"] == expected.digest.left_end.sequence
        assert digest["left_end"]["kind"] == expected.digest.left_end.kind
        assert digest["right_end"]["sequence"] == expected.digest.right_end.sequence
        assert digest["right_end"]["kind"] == expected.digest.right_end.kind

    def test_the_insert_can_be_posted_straight_to_the_clone_endpoint(self, auth_client):
        """The whole point of the feature: what comes out of the PCR page
        goes into the Clone page without the caller reshaping it."""
        pcr = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI"},
            format="json",
        )
        assert pcr.data["is_clean"]

        clone = auth_client.post(
            reverse("design-clone"),
            {"sequence": pcr.data["digest"]["top"],
             "insert_reverse": pcr.data["digest"]["bottom"],
             "pre_digested": True,
             "vector_key": "pET-21a",
             "left_enzyme": "NdeI", "right_enzyme": "XhoI",
             "is_coding": True, "remove_stop": False, "cleavage_site": None,
             "include_his_tag": False, "include_linkers": False},
            format="json",
        )
        assert clone.status_code == 200, clone.data
        assert clone.data["is_clonable"]


@pytest.mark.django_db
class TestSitesInsideTheGene:
    def test_an_internal_site_comes_back_as_a_problem_with_no_insert(self, auth_client):
        """A problem, not a warning: the digest would cut the gene in two as
        well as opening its ends, and no reaction condition rescues that."""
        gene = "ATG" + "GCTAGC" + GENE[3:]
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": gene, "left_enzyme": "NheI", "right_enzyme": "XhoI"},
            format="json",
        )
        assert r.status_code == 200
        assert not r.data["is_clean"]
        assert any("cuts inside" in p for p in r.data["problems"])
        assert r.data["digest"] is None

    def test_a_clean_design_reports_no_problems(self, auth_client):
        r = auth_client.post(
            reverse("design-pcr"),
            {"template": GENE, "left_enzyme": "NdeI", "right_enzyme": "XhoI"},
            format="json",
        )
        assert r.data["is_clean"]
        assert r.data["problems"] == []
        assert r.data["digest"] is not None
