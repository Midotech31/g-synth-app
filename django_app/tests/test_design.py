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

    def test_includes_the_hybridisation_view(self, auth_client):
        """The client draws the duplex from coordinates, not from prose."""
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        duplex = response.data["duplex"]

        assert duplex["mismatches"] == []
        assert len(duplex["top"]) == len(duplex["bottom"]) == duplex["width"]
        assert len(duplex["pairs"]) == duplex["width"]
        assert duplex["top"].strip() == response.data["construct_forward"]

    def test_duplex_spans_cover_both_strands(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        duplex = response.data["duplex"]

        assert len(duplex["top_fragments"]) == response.data["fragment_count"]
        assert len(duplex["bottom_fragments"]) == response.data["fragment_count"]
        assert duplex["segments"], "the cassette should be labelled for colouring"

        for span, fragment in zip(duplex["top_fragments"], response.data["fragments"]):
            assert duplex["top"][span["start"]:span["end"]] == fragment["forward"]

    def test_reports_which_strand_carries_each_overhang(self, auth_client):
        """NdeI and XhoI both leave 5' overhangs, so both sit on the outer top/bottom."""
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        fragments = response.data["fragments"]

        assert fragments[0]["left_overhang_strand"] == "top"
        assert fragments[-1]["right_overhang_strand"] == "bottom"

    def test_states_the_conditions_every_tm_refers_to(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": LONG_INSERT})
        conditions = response.data["tm_conditions"]

        assert "SantaLucia" in conditions["model"]
        assert "µM" in conditions["summary"] and "Na" in conditions["summary"]

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


# ── Cloning ─────────────────────────────────────────────────────────────────


def build_vector(left: str = "NdeI", right: str = "XhoI") -> str:
    """A circular vector with exactly one site for each enzyme.

    Built from the enzyme table rather than hand-written, so a filler that
    happens to carry a second site cannot turn a cloning test into a test of
    luck.
    """
    from gsynth_engine.cloning import find_sites
    from gsynth_engine.tests.test_cloning import build_vector as engine_build

    vector = engine_build(left, right)
    assert len(find_sites(vector, left)) == 1
    return vector


@pytest.mark.django_db
class TestCloneEndpoint:
    url_name = "design-clone"

    def test_requires_authentication(self, api_client):
        response = api_client.post(reverse(self.url_name), {
            "sequence": GOLDEN_INSERT, "vector": build_vector(),
        })
        assert response.status_code == 401

    def test_returns_the_recombinant_plasmid(self, auth_client):
        vector = build_vector()
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": vector, "name": "pGS-EntA",
        })
        assert response.status_code == 200, response.data
        assert response.data["is_clonable"], response.data["problems"]

        plasmid = response.data["plasmid"]
        assert len(plasmid) == response.data["length"]
        assert response.data["topology"] == "circular"
        assert response.data["backbone_length"] < len(vector)

    def test_the_insert_sits_where_the_response_says(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(),
        })
        data = response.data
        placed = data["plasmid"][data["insert_start"]:data["insert_end"]]
        assert placed == data["assembly"]["construct_forward"]

    def test_reports_both_junctions_with_their_enzymes(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(),
        })
        junctions = response.data["junctions"]
        assert [j["enzyme"] for j in junctions] == ["NdeI", "XhoI"]
        assert all(j["site_regenerated"] for j in junctions)

    def test_translates_the_protein_that_will_be_expressed(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(),
        })
        assert response.data["protein"].startswith("MGSSHHHHHHSSG")
        assert response.data["protein_length"] == len(response.data["protein"])

    def test_the_insert_comes_back_as_a_drawable_annotation(self, auth_client):
        """So the client draws the whole map from one list."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(), "name": "EntA",
        })
        data = response.data
        insert = [a for a in data["annotations"] if a["name"] == "EntA"]
        assert len(insert) == 1
        assert insert[0]["start"] == data["insert_start"]
        assert insert[0]["end"] == data["insert_end"]

    def test_vector_annotations_are_carried_over(self, auth_client):
        vector = build_vector()
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT,
            "vector": vector,
            "vector_annotations": [
                {"name": "ori", "type": "rep_origin", "start": 5, "end": 80,
                 "direction": 1, "color": "#6A4C93"},
            ],
        }, format="json")
        assert response.status_code == 200, response.data
        names = [a["name"] for a in response.data["annotations"]]
        assert "ori" in names

    def test_a_vector_cut_twice_is_refused_with_a_reason(self, auth_client):
        vector = build_vector()
        doubled = vector[:100] + "CATATG" + vector[100:]
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": doubled,
        })
        assert response.status_code == 400
        assert "NdeI cuts it 2 times" in response.data["detail"]

    def test_an_unclonable_design_returns_200_with_problems(self, auth_client):
        """The user needs to see what does not fit, not an error page."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": "GGCTAAATCGTGGAACAGTGCTGCACCAGCTGCAGCCTGTACCAGCTGGAA",
            "vector": build_vector(),
        })
        assert response.status_code == 200, response.data
        assert not response.data["is_clonable"]
        assert response.data["problems"]

    def test_cloning_without_fragmenting_skips_the_assembly(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": GOLDEN_INSERT, "vector": build_vector(), "fragment": False,
        })
        assert response.status_code == 200, response.data
        assert response.data["assembly"] is None
        assert response.data["insert"]["forward"] == GOLDEN_FORWARD

    def test_saves_as_a_project_when_asked(self, auth_client, user):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(),
            "name": "pGS-EntA", "save_as_project": True,
        })
        assert response.data["project_id"]
        project = Project.objects.get(pk=response.data["project_id"])
        assert project.user == user
        assert project.module == "cloning"
        assert project.sequence == response.data["plasmid"]

    def test_does_not_save_by_default(self, auth_client):
        auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": build_vector(),
        })
        assert not Project.objects.filter(module="cloning").exists()

    def test_matches_the_engine_exactly(self, auth_client):
        """The HTTP layer must not reinterpret the biology."""
        from gsynth_engine.cloning import clone

        vector = build_vector()
        plan = design_merzoug_assembly(LONG_INSERT)
        expected = clone(
            vector, plan.construct_forward, insert_reverse=plan.construct_reverse,
            left_enzyme="NdeI", right_enzyme="XhoI", orf_start=plan.ssd.orf_start,
        )
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector": vector,
        })
        assert response.data["plasmid"] == expected.plasmid
        assert response.data["protein"] == expected.protein


@pytest.mark.django_db
class TestVectorCatalogue:
    def test_is_public(self, api_client):
        """The cloning page builds its dropdown before anything is designed."""
        assert api_client.get(reverse("design-vectors")).status_code == 200

    def test_pet21a_is_the_default(self, api_client):
        data = api_client.get(reverse("design-vectors")).data
        assert data["default"] == "pET-21a"
        assert data["vectors"][0]["key"] == "pET-21a"

    def test_entries_carry_what_the_ui_shows(self, api_client):
        data = api_client.get(reverse("design-vectors")).data
        entry = {v["key"]: v for v in data["vectors"]}["pET-21a"]
        assert entry["length"] == 5443
        assert entry["resistance"] == "Ampicillin"
        assert entry["recommended_pairs"][0] == "NdeI / XhoI"
        assert entry["has_sequence"] is True
        assert any(tag["name"] == "His-tag" and tag["end"] == "C"
                   for tag in entry["tags"])

    def test_more_than_one_vector_is_offered(self, api_client):
        data = api_client.get(reverse("design-vectors")).data
        keys = {v["key"] for v in data["vectors"]}
        assert {"pET-21a", "pET-21", "pET-28a"} <= keys

    def test_a_bundled_sequence_can_be_fetched(self, api_client):
        response = api_client.get(reverse("design-vector", args=["pET-21a"]))
        assert response.status_code == 200
        assert response.data["length"] == 5443
        assert len(response.data["sequence"]) == 5443
        assert response.data["annotations"]

    def test_a_vector_without_a_sequence_says_so(self, api_client):
        response = api_client.get(reverse("design-vector", args=["pGEX-4T-1"]))
        assert response.status_code == 404
        assert "Import your own copy" in response.data["detail"]

    def test_an_unknown_vector_is_a_404(self, api_client):
        assert api_client.get(
            reverse("design-vector", args=["pNope"])
        ).status_code == 404


@pytest.mark.django_db
class TestCloningIntoACatalogueVector:
    url_name = "design-clone"

    def test_a_vector_key_alone_is_enough(self, auth_client):
        """No sequence needed: pET-21a's ships with G-Synth."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA",
        })
        assert response.status_code == 200, response.data
        assert response.data["is_clonable"], response.data["problems"]
        assert response.data["vector_name"] == "pET-21a(+)"
        assert response.data["length"] > 5000

    def test_the_backbone_survives(self, auth_client):
        """pET-21a's cassette reads on the minus strand; getting that wrong
        keeps the 78 bp stuffer and discards the origin and the marker."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a",
        })
        data = response.data
        assert data["reversed_insert"] is True
        assert data["backbone_length"] > 5000
        assert data["removed_length"] < 200

    def test_the_vectors_features_come_across(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a",
        })
        names = {a["name"] for a in response.data["annotations"]}
        assert {"AmpR", "ori", "lacI", "T7 promoter"} <= names

    def test_it_reports_which_vector_tags_land_on_the_protein(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a",
        })
        tags = {t["name"]: t for t in response.data["tags"]}
        assert "His-tag" in tags and "T7·Tag" in tags
        # NdeI cloning replaces the T7·Tag with the insert.
        assert tags["T7·Tag"]["present"] is False

    def test_the_supplied_sequence_is_checked_against_the_entry(self, auth_client):
        """Pasting pET-21(+) while pET-21a(+) is selected must be caught."""
        from gsynth_engine import vectors

        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT,
            "vector_key": "pET-21a",
            "vector": vectors.sequence_of("pET-21")["sequence"],
            "left_enzyme": "BamHI", "right_enzyme": "XhoI",
        })
        assert response.status_code == 200, response.data
        check = response.data["vector"]["check"]
        assert check["matches"] is False
        assert any("74 bp shorter" in p for p in check["problems"])

    def test_a_matching_sequence_passes_its_check(self, auth_client):
        from gsynth_engine import vectors

        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT,
            "vector_key": "pET-21a",
            "vector": vectors.sequence_of("pET-21a")["sequence"],
        })
        assert response.data["vector"]["check"]["matches"] is True

    def test_an_unrecognised_sequence_is_still_cloned_into(self, auth_client):
        """A lab's own backbone is not in any catalogue, and still works."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT,
            "vector_key": "",
            "vector": build_vector(),
        })
        assert response.status_code == 200, response.data
        assert response.data["is_clonable"]
        assert response.data["vector"]["recognised"] is False

    def test_pet21_has_no_ndei_so_the_default_pair_is_refused(self, auth_client):
        """The distinction that matters between pET-21(+) and pET-21a(+)."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21",
        })
        assert response.status_code == 400
        assert "NdeI does not cut" in response.data["detail"]

    def test_pet21_works_with_its_own_pair(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pET-21",
            "left_enzyme": "BamHI", "right_enzyme": "XhoI",
        })
        assert response.status_code == 200, response.data
        assert response.data["is_clonable"], response.data["problems"]

    def test_a_vector_without_a_bundled_sequence_asks_for_one(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "vector_key": "pGEX-4T-1",
            "left_enzyme": "BamHI", "right_enzyme": "EcoRI",
        })
        assert response.status_code == 400
        assert "paste or import your own" in str(response.data)


@pytest.mark.django_db
class TestOptimiseEndpoint:
    url_name = "design-optimise"

    # An enterocin-like gene: AT-rich, slow codons, internal NdeI site.
    DONOR = (
        "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA"
        "TTAAATTTAGCATTATTAGGATTAGCAAGTTTATTAGGTAAAGGTATTAGTAAATTAGGA"
    )
    DONOR_PROTEIN = "MTTSKLGKGLGYIGNNGAHMGLNLALLGLASLLGKGISKLG"

    def test_requires_authentication(self, api_client):
        assert api_client.post(
            reverse(self.url_name), {"sequence": self.DONOR}
        ).status_code == 401

    def test_the_protein_is_unchanged(self, auth_client):
        """The invariant, checked over HTTP as well as in the engine."""
        from gsynth_engine.cloning import translate

        response = auth_client.post(reverse(self.url_name), {"sequence": self.DONOR})
        assert response.status_code == 200, response.data
        assert response.data["protein"] == self.DONOR_PROTEIN
        assert translate(response.data["sequence"]).rstrip("*") == self.DONOR_PROTEIN

    def test_it_adapts_the_gene_to_the_host(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": self.DONOR})
        data = response.data
        assert data["cai_after"] > data["cai_before"]
        assert data["rare_codons_after"] < data["rare_codons_before"]
        assert 40 <= data["gc_after"] <= 60

    def test_the_cloning_sites_are_removed(self, auth_client):
        """A gene with an internal NdeI site cannot be cloned NdeI/XhoI."""
        from gsynth_engine.cloning import find_sites

        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR, "avoid_enzymes": ["NdeI", "XhoI"],
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["is_clean"], response.data["problems"]
        assert find_sites(response.data["sequence"], "NdeI", circular=False) == []
        assert "CATATG" in response.data["sites_removed"]

    def test_a_protein_can_be_reverse_translated(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR_PROTEIN, "is_protein": True,
        })
        assert response.status_code == 200, response.data
        assert response.data["protein"] == self.DONOR_PROTEIN
        assert response.data["cai_before"] is None
        assert response.data["length"] == 3 * len(self.DONOR_PROTEIN) + 3

    def test_the_stop_codon_can_be_left_off(self, auth_client):
        """For an insert going into a C-terminal vector tag."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR_PROTEIN, "is_protein": True, "keep_stop": False,
        })
        assert response.data["length"] == 3 * len(self.DONOR_PROTEIN)

    def test_a_reference_set_replaces_the_shipped_table(self, auth_client):
        """The honest way to get a CAI: measure it against genes you chose."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR,
            "reference_genes": ["ATGTTATTATTAAAAAAA" * 3],
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["table"] == "your reference set"
        assert "1 genes you supplied" in response.data["table_source"]

    def test_a_partial_codon_is_a_usable_message(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {"sequence": "ATGAAAA"})
        assert response.status_code == 400
        assert "multiple of three" in response.data["detail"]

    def test_an_inverted_gc_window_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR, "gc_min": 70, "gc_max": 40,
        })
        assert response.status_code == 400
        assert "above the lower one" in str(response.data)

    def test_it_matches_the_engine_exactly(self, auth_client):
        from gsynth_engine.codon import Constraints, optimise

        expected = optimise(
            self.DONOR, constraints=Constraints(avoid_enzymes=("NdeI",)),
        )
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR, "avoid_enzymes": ["NdeI"],
        }, format="json")
        assert response.data["sequence"] == expected.sequence
