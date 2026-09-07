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

from apps.projects.models import Project
from gsynth_engine.constants import ALL_ENZYMES
from gsynth_engine.esd import design_extended_sequence
from gsynth_engine.ssd import design_small_sequence

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
        assert by_name["CviAII"]["supplies_start_codon"] is True
        assert by_name["FatI"]["supplies_start_codon"] is True
        assert by_name["NcoI"]["supplies_start_codon"] is False
        assert by_name["NsiI"]["supplies_start_codon"] is False
        assert by_name["XhoI"]["overhang"] == "TCGA"
        assert by_name["KpnI"]["overhang_type"] == "3'"
        assert by_name["SmaI"]["overhang_type"] == "blunt"
        assert by_name["HindIII"]["recognition"] == "AAGCTT"
        assert by_name["HindIII"]["overhang"] == "AGCT"
        assert by_name["HindIII"]["overhang_type"] == "5'"
        assert "AsuNHI" in by_name["NheI"]["aliases"]
        assert "AvrII" in by_name["BlnI"]["aliases"]
        assert set(by_name) == set(ALL_ENZYMES)
        assert data["canonical_geometries"] == 109
        assert data["selectable_names"] == 289

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


def _random_gene(length: int, seed: int) -> str:
    """A coding insert long enough to need more junctions than 4 nt can supply."""
    import random

    codons = ("GCG", "TGC", "GAT", "GAA", "TTT", "GGC", "CAT", "ATT", "AAA", "CTG",
              "ATG", "AAC", "CCG", "CAG", "CGC", "AGC", "ACC", "GTG", "TGG", "TAT")
    rng = random.Random(seed)
    return "".join(rng.choice(codons) for _ in range(length // 3)) + "TAA"

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

    def test_reports_the_ends_measured_off_the_assembly(self, auth_client):
        """The response must say what the fragments present, not what was asked.

        Every terminal value in the plan is copied from the SSD, so a payload
        built from those labels would agree with itself whatever the oligos
        actually spell. These come from the assembled duplex, and carry the
        polarity — which follows the side, not the strand.
        """
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "left_enzyme": "NdeI", "right_enzyme": "XhoI",
        })
        assert response.data["terminal_ends"] == [
            {"side": "left", "enzyme": "NdeI", "overhang": "TA", "kind": "5'"},
            {"side": "right", "enzyme": "XhoI", "overhang": "TCGA", "kind": "5'"},
        ]

    def test_a_three_prime_pair_is_reported_as_three_prime(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "left_enzyme": "KpnI", "right_enzyme": "SacI",
        })
        assert [e["kind"] for e in response.data["terminal_ends"]] == ["3'", "3'"]
        assert [e["overhang"] for e in response.data["terminal_ends"]] == ["GTAC", "AGCT"]

    def test_the_ends_match_what_the_engine_measured(self, auth_client):
        expected = design_extended_sequence(LONG_INSERT, target_oligo_length=90)
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "target_oligo_length": 90,
        })
        left, right = expected.terminal_ends
        assert [
            (e["overhang"], e["kind"]) for e in response.data["terminal_ends"]
        ] == [left, right]

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

        for span, fragment in zip(
            duplex["top_fragments"], response.data["fragments"], strict=True
        ):
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
        expected = design_extended_sequence(LONG_INSERT, target_oligo_length=90)
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "target_oligo_length": 90,
        })
        assert response.data["construct_forward"] == expected.construct_forward
        assert [f["forward"] for f in response.data["fragments"]] == [
            f.forward for f in expected.fragments
        ]

    def test_saves_as_a_project_when_asked(self, auth_client, user):
        expected = design_extended_sequence(LONG_INSERT)
        insert = next(segment for segment in expected.ssd.segments if segment.name == "insert")
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "name": "EntA assembly", "save_as_project": True,
        })
        project = Project.objects.get(id=response.data["project_id"])
        assert project.user == user
        assert project.module == "extended_sequence_design"
        assert project.data["fragment_count"] == response.data["fragment_count"]
        assert project.data["insert_start"] == response.data["insert_start"] == insert.start
        assert project.data["insert_end"] == response.data["insert_end"] == insert.end
        assert project.data["topology"] == response.data["topology"] == "linear"
        assert "backbone_length" not in project.data

    def test_rejects_an_overhang_outside_the_method(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT, "overhang_length": 2,
        })
        assert response.status_code == 400

    def test_a_whole_gene_designs_over_http(self, auth_client):
        """The endpoint accepts 200 kb, but was only ever exercised on peptides.

        A 2.4 kb gene needs 26 junctions and there are only 22 sets of 4 nt
        overhangs that ligate in one order, so this used to come back as a 400
        naming a junction the user had never heard of.
        """
        gene = _random_gene(2_400, seed=7)
        response = auth_client.post(reverse(self.url_name), {
            "sequence": gene, "target_oligo_length": 90, "overhang_length": 4,
        })
        assert response.status_code == 200, response.data
        assert response.data["verification"] == []
        assert response.data["fragment_count"] > 25

    def test_reports_the_overhang_it_actually_used(self, auth_client):
        """The client shows what was built, not what was asked for.

        The design widens the overhang when 4 nt cannot supply the junctions.
        A response that echoed the request would have the interface — and the
        bench protocol printed from it — describing oligos nobody ordered.
        """
        gene = _random_gene(2_400, seed=7)
        expected = design_extended_sequence(
            gene, target_oligo_length=90, overhang_length=4,
        )
        response = auth_client.post(reverse(self.url_name), {
            "sequence": gene, "target_oligo_length": 90, "overhang_length": 4,
        })
        assert response.data["overhang_length"] == expected.overhang_length > 4
        assert all(
            len(j) == expected.overhang_length
            for j in response.data["junction_overhangs"]
        )
        assert any("widened" in note for note in response.data["warnings"])


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
        assert "EXTENDED SEQUENCE DESIGN" in text
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

    def test_cassette_parts_are_coordinate_level_annotations(self, auth_client):
        """The viewer can show the tag, cleavage site and target without motif guesses."""
        response = auth_client.post(reverse(self.url_name), {
            "sequence": LONG_INSERT,
            "vector": build_vector(),
            "name": "Insulin cassette",
            "cleavage_site": "Thrombin",
        })
        assert response.status_code == 200, response.data
        data = response.data
        by_name = {}
        for annotation in data["annotations"]:
            by_name.setdefault(annotation["name"], []).append(annotation)

        assert "6×His tag" in by_name
        assert "Thrombin site" in by_name
        assert "Insulin cassette target" in by_name

        cassette = next(
            annotation for annotation in by_name["Insulin cassette"]
            if annotation["type"] == "CDS"
        )
        assert cassette["translation_start"] == data["insert_start"] + data["insert"]["orf_start"]
        assert cassette["translation_end"] == data["insert_end"]

        tag = by_name["6×His tag"][0]
        assert data["plasmid"][tag["start"]:tag["end"]] == "CACCACCACCACCACCAC"

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
        plan = design_extended_sequence(LONG_INSERT)
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
        assert response.data["preflight"]["can_export"] is True
        assert response.data["provenance"]["workflow"] == "cloning"
        assert len(response.data["provenance"]["output_sha256"]) == 64

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
        assert data["rare_codons_after"] <= data["rare_codons_before"]
        assert data["changed_codons"] > 0
        assert 40 <= data["gc_after"] <= 60

    def test_host_catalogue_is_public_and_identifies_sources(self, api_client):
        response = api_client.get(reverse("design-codon-hosts"))
        assert response.status_code == 200
        assert response.data["default"] == "ecoli"
        hosts = {host["key"]: host for host in response.data["hosts"]}
        assert len(hosts) == 15
        assert {
            "ecoli", "b_subtilis", "p_putida", "l_lactis",
            "c_glutamicum", "s_coelicolor", "s_cerevisiae", "k_phaffii",
            "k_lactis", "y_lipolytica", "h_sapiens", "c_griseus",
            "s_frugiperda", "d_melanogaster", "n_benthamiana",
        } == set(hosts)
        assert response.data["dataset"]["name"] == "FDA HIVE-CUTs / CoCoPUTs"
        assert response.data["dataset"]["release"] == "September 2021"
        assert len(response.data["dataset"]["sha256"]) == 64
        yeast = hosts["s_cerevisiae"]
        assert yeast["dataset"] == "RefSeq"
        assert yeast["taxon_id"] == 4932
        assert yeast["coding_sequences"] == 5983
        assert yeast["codon_count"] == 2929341
        assert yeast["metric_label"] == "Profile-relative CAI"
        assert "HIVE-CUTs/CoCoPUTs" in yeast["source"]

    def test_selected_host_controls_the_usage_table(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": "MLRLKFY",
            "is_protein": True,
            "keep_stop": False,
            "host": "s_cerevisiae",
            "gc_min": 10,
            "gc_max": 90,
            "avoid_rare": False,
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["host"] == "s_cerevisiae"
        assert response.data["table"] == "Saccharomyces cerevisiae"
        assert response.data["protein"] == "MLRLKFY"
        assert response.data["metric_label"] == "Profile-relative CAI"
        assert response.data["expression_yield_predicted"] is False

    def test_custom_reference_genes_are_identified_as_context_specific(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": "MLRLKFY",
            "is_protein": True,
            "keep_stop": False,
            "reference_genes": ["ATGCTGCGTCTGAAA", "ATGCTGCGTCTGAAG"],
            "gc_min": 10,
            "gc_max": 90,
            "avoid_rare": False,
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["host"] == "custom"
        assert response.data["metric_label"] == "CAI (custom reference set)"
        assert response.data["expression_yield_predicted"] is False

    def test_unknown_host_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": self.DONOR,
            "host": "unknown-host",
        })
        assert response.status_code == 400

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
        assert response.data["protein_context"] == "complete_orf"
        assert response.data["initiator_methionine_added"] is False
        assert response.data["recommended_design_is_coding"] is True

    def test_a_mature_peptide_is_preserved_without_an_artificial_methionine(
        self, auth_client,
    ):
        peptide = "GIVEQCCTSICSLYQLENYCG"
        response = auth_client.post(reverse(self.url_name), {
            "sequence": peptide,
            "is_protein": True,
            "protein_context": "auto",
            "keep_stop": False,
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["input_protein"] == peptide
        assert response.data["protein"] == peptide
        assert response.data["protein_context"] == "mature_peptide"
        assert response.data["initiator_methionine_added"] is False
        assert response.data["recommended_design_is_coding"] is False

    def test_a_complete_orf_can_receive_one_initiator_methionine(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "sequence": "GIVEQ",
            "is_protein": True,
            "protein_context": "complete_orf",
            "keep_stop": False,
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["protein"] == "MGIVEQ"
        assert response.data["initiator_methionine_added"] is True
        assert response.data["recommended_design_is_coding"] is True

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


@pytest.mark.django_db
class TestExport:
    """A design that only exists inside G-Synth is not finished."""

    def parse_genbank(self, response):
        import io

        from Bio import SeqIO

        return SeqIO.read(io.StringIO(response.content.decode()), "genbank")

    def test_the_recombinant_plasmid_exports_as_genbank(self, auth_client):
        response = auth_client.post(reverse("design-clone-export"), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "pGS-EntA",
        })
        assert response.status_code == 200, response.content[:300]
        assert 'filename="pGS-EntA.gb"' in response["Content-Disposition"]

        record = self.parse_genbank(response)
        assert record.annotations["topology"] == "circular"
        assert len(record.seq) > 5000
        labels = {
            f.qualifiers.get("label", [""])[0]
            for f in record.features if f.type != "source"
        }
        assert {"AmpR", "ori", "pGS-EntA"} <= labels

    def test_the_exported_plasmid_is_the_one_the_api_returned(self, auth_client):
        """Two endpoints, one molecule — they must not drift."""
        payload = {"sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA"}
        designed = auth_client.post(reverse("design-clone"), payload).data
        exported = self.parse_genbank(
            auth_client.post(reverse("design-clone-export"), payload)
        )
        assert str(exported.seq).upper() == designed["plasmid"]

    def test_reviewed_product_annotations_are_saved_and_exported(self, auth_client):
        payload = {
            "sequence": LONG_INSERT,
            "vector_key": "pET-21a",
            "name": "EntA",
            "product_annotations": [{
                "name": "reviewed therapeutic insert",
                "type": "mat_peptide",
                "start": 100,
                "end": 160,
                "direction": 1,
                "color": "#3F7A52",
                "translation_start": 100,
                "translation_end": 160,
                "truncated": False,
            }],
            "save_as_project": True,
        }
        designed = auth_client.post(reverse("design-clone"), payload)
        assert designed.status_code == 200, designed.data
        assert designed.data["annotations"] == payload["product_annotations"]

        project = Project.objects.get(pk=designed.data["project_id"])
        assert project.data["annotations"] == payload["product_annotations"]

        exported = self.parse_genbank(
            auth_client.post(reverse("design-clone-export"), payload)
        )
        labels = {
            feature.qualifiers.get("label", [""])[0]
            for feature in exported.features if feature.type != "source"
        }
        assert "reviewed therapeutic insert" in labels

    def test_a_reviewed_annotation_cannot_escape_the_product_coordinates(self, auth_client):
        response = auth_client.post(reverse("design-clone-export"), {
            "sequence": LONG_INSERT,
            "vector_key": "pET-21a",
            "name": "EntA",
            "product_annotations": [{
                "name": "invalid feature",
                "type": "misc_feature",
                "start": 999_999,
                "end": 1_000_020,
                "direction": 1,
                "color": "#3F7A52",
            }],
        })
        assert response.status_code == 400
        assert "outside" in response.data["detail"]

    def test_reviewed_cds_translation_must_stay_inside_its_feature(self, auth_client):
        response = auth_client.post(reverse("design-clone-export"), {
            "sequence": LONG_INSERT,
            "vector_key": "pET-21a",
            "name": "EntA",
            "product_annotations": [{
                "name": "invalid CDS",
                "type": "CDS",
                "start": 100,
                "end": 160,
                "direction": 1,
                "color": "#3F7A52",
                "translation_start": 99,
                "translation_end": 160,
            }],
        })
        assert response.status_code == 400
        assert "Translation coordinates" in str(response.data)

    def test_the_plasmid_exports_as_fasta(self, auth_client):
        response = auth_client.post(
            reverse("design-clone-export") + "?filetype=fasta",
            {"sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA"},
        )
        assert response.status_code == 200
        assert response.content.decode().startswith(">EntA")

    def test_the_plasmid_exports_as_valid_sbol3(self, auth_client):
        from apps.sequences.parsing import parse_sequence_file

        response = auth_client.post(
            reverse("design-clone-export") + "?filetype=sbol3",
            {"sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA"},
        )
        record = parse_sequence_file(response.content, "EntA.sbol.json")

        assert response.status_code == 200
        assert record.topology == "circular"
        assert record.length > 5000
        assert {annotation.name for annotation in record.annotations} >= {"AmpR", "EntA"}

    def test_cloning_worksheet_links_preflight_bands_ligation_and_primers(self, auth_client):
        response = auth_client.post(
            reverse("design-clone-worksheet"),
            {"sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA"},
        )
        assert response.status_code == 200, response.content[:300]
        text = response.content.decode()
        assert "PREFLIGHT RELEASE" in text
        assert "Expected gel bands" in text
        assert "3:1 insert:vector" in text
        assert "SEQUENCING PRIMERS" in text
        assert "Output SHA-256" in text

    def test_the_construct_exports_with_its_cassette_labelled(self, auth_client):
        response = auth_client.post(reverse("design-construct-export"), {
            "sequence": LONG_INSERT, "name": "EntA",
        })
        assert response.status_code == 200, response.content[:300]

        record = self.parse_genbank(response)
        assert record.annotations["topology"] == "linear"
        labels = {
            f.qualifiers.get("label", [""])[0]
            for f in record.features if f.type != "source"
        }
        assert "6×His tag" in labels
        assert "F1" in labels          # the fragments are drawn too

    def test_the_linear_construct_exports_as_valid_sbol3(self, auth_client):
        from apps.sequences.parsing import parse_sequence_file

        response = auth_client.post(
            reverse("design-construct-export") + "?filetype=sbol3",
            {"sequence": LONG_INSERT, "name": "EntA"},
        )
        record = parse_sequence_file(response.content, "EntA.sbol.json")

        assert response.status_code == 200
        assert record.topology == "linear"
        assert any(annotation.name == "6×His tag" for annotation in record.annotations)

    def test_the_oligos_export_as_one_fasta_per_oligo(self, auth_client):
        """Suppliers take a FASTA upload; retyping thirty names is where
        transcription errors come from."""
        import io

        from Bio import SeqIO

        designed = auth_client.post(reverse("design-assembly"), {
            "sequence": LONG_INSERT, "name": "EntA",
        }).data
        response = auth_client.post(
            reverse("design-construct-export") + "?filetype=oligos",
            {"sequence": LONG_INSERT, "name": "EntA"},
        )
        records = list(SeqIO.parse(io.StringIO(response.content.decode()), "fasta"))

        assert len(records) == designed["oligo_count"]
        assert records[0].id.startswith("EntA_")
        assert str(records[0].seq) == designed["oligos"][0]["Sequence (5'->3')"]

    def test_all_sequences_export_contains_the_assembled_duplex_and_every_oligo(
        self, auth_client,
    ):
        import io

        from Bio import SeqIO

        designed = auth_client.post(reverse("design-assembly"), {
            "sequence": LONG_INSERT, "name": "EntA",
        }).data
        response = auth_client.post(
            reverse("design-construct-export") + "?filetype=all-sequences",
            {"sequence": LONG_INSERT, "name": "EntA"},
        )
        records = list(SeqIO.parse(io.StringIO(response.content.decode()), "fasta"))

        assert response.status_code == 200
        assert "EntA_all_sequences.fasta" in response["Content-Disposition"]
        assert len(records) == designed["oligo_count"] + 2
        assert records[0].id == "EntA_assembled_forward"
        assert records[1].id == "EntA_assembled_reverse"
        assert str(records[0].seq) == designed["construct_forward"]
        assert str(records[1].seq) == designed["construct_reverse"]
        assert [record.id for record in records[2:]] == [
            row["Name"] for row in designed["oligos"]
        ]

    def test_export_requires_authentication(self, api_client):
        for name in ("design-clone-export", "design-construct-export"):
            assert api_client.post(
                reverse(name), {"sequence": LONG_INSERT}
            ).status_code == 401

    def test_a_bad_design_still_errors_usefully(self, auth_client):
        response = auth_client.post(reverse("design-construct-export"), {
            "sequence": "ATGXYZ",
        })
        assert response.status_code == 400
        assert "not A, C, G or T" in response.data["detail"]


@pytest.mark.django_db
class TestLigationEndpoint:
    url_name = "design-ligation"

    def test_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {
            "vector_length": 5443, "insert_length": 150,
        }).status_code == 401

    def test_equal_mass_is_nowhere_near_equal_moles(self, auth_client):
        """The mistake the module exists to prevent, over HTTP."""
        response = auth_client.post(reverse(self.url_name), {
            "vector_length": 5443, "insert_length": 150, "vector_ng": 50, "ratio": 3,
        })
        assert response.status_code == 200, response.data
        reaction = response.data["reactions"][0]
        assert reaction["insert_ng"] < 5      # not 50
        assert reaction["insert_fmol"] == pytest.approx(
            reaction["vector_fmol"] * 3, rel=1e-2
        )

    def test_a_series_comes_back_when_ratios_are_given(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "vector_length": 5443, "insert_length": 150, "ratios": [1, 3, 5],
        }, format="json")
        assert [r["ratio"] for r in response.data["reactions"]] == [1.0, 3.0, 5.0]

    def test_blunt_ends_carry_their_advice(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "vector_length": 5000, "insert_length": 500, "ends": "blunt",
        })
        notes = " ".join(response.data["reactions"][0]["warnings"])
        assert "dephosphorylate" in notes

    def test_a_zero_length_is_refused(self, auth_client):
        assert auth_client.post(reverse(self.url_name), {
            "vector_length": 0, "insert_length": 150,
        }).status_code == 400


@pytest.mark.django_db
class TestPrimerEndpoint:
    url_name = "design-primers"

    def clone_something(self, auth_client, insert=None):
        return auth_client.post(reverse("design-clone"), {
            "sequence": insert or LONG_INSERT, "vector_key": "pET-21a", "name": "EntA",
        }).data

    def test_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {
            "template": "ACGT" * 100, "target_start": 10, "target_end": 50,
        }).status_code == 401

    def test_it_designs_primers_that_read_the_insert(self, auth_client):
        cloned = self.clone_something(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
            "name": "EntA",
        })
        assert response.status_code == 200, response.data
        assert response.data["covers_target"], response.data["gaps"]
        assert {p["direction"] for p in response.data["primers"]} == {1, -1}

    def test_every_primer_is_unique_in_the_plasmid(self, auth_client):
        """One that binds twice gives a superimposed trace and no data."""
        from gsynth_engine.sequence import reverse_complement

        cloned = self.clone_something(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
        })
        plasmid = cloned["plasmid"]
        for primer in response.data["primers"]:
            hits = plasmid.count(primer["sequence"]) + plasmid.count(
                reverse_complement(primer["sequence"])
            )
            assert hits == 1, primer["name"]

    def test_the_rows_carry_what_a_supplier_needs(self, auth_client):
        cloned = self.clone_something(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
        })
        row = response.data["rows"][0]
        assert set(row) >= {"Name", "Sequence (5'->3')", "Tm (°C)"}

    def test_a_margin_inside_the_dead_zone_is_refused(self, auth_client):
        cloned = self.clone_something(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
            "margin": 50,
        })
        # 50 is the floor the serializer allows and the engine's dead zone,
        # so this is the boundary rather than an error.
        assert response.status_code == 200

    def test_an_inverted_region_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "template": "ACGT" * 100, "target_start": 90, "target_end": 20,
        })
        assert response.status_code == 400
        assert "end after it starts" in str(response.data)


@pytest.mark.django_db
class TestVerifyEndpoint:
    url_name = "design-verify"

    def build(self, auth_client):
        return auth_client.post(reverse("design-clone"), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA",
        }).data

    def read_of(self, plasmid, start, end, noise=30):
        length = len(plasmid)
        body = "".join(plasmid[i % length] for i in range(start, end))
        return "A" * noise + body + "T" * noise

    def test_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {
            "design": "ACGT" * 100, "reads": {"a": "ACGT" * 20},
        }, format="json").status_code == 401

    def test_a_matching_read_verifies_the_construct(self, auth_client):
        cloned = self.build(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "design": cloned["plasmid"],
            "reads": {
                "T7-F": self.read_of(
                    cloned["plasmid"],
                    cloned["insert_start"] - 200, cloned["insert_end"] + 200,
                ),
            },
            "region_start": cloned["insert_start"],
            "region_end": cloned["insert_end"],
        }, format="json")
        assert response.status_code == 200, response.data
        assert response.data["is_verified"]
        assert response.data["coverage"] == 100.0
        assert response.data["differences"] == []
        assert response.data["verification_state"] == "fully_verified"
        assert response.data["region_start"] == cloned["insert_start"]
        assert response.data["region_end"] == cloned["insert_end"]
        assert response.data["preflight"]["verdict"] == "ready"
        assert response.data["provenance"]["workflow"] == "sequence_verification"

    def test_a_reversed_read_is_handled(self, auth_client):
        """Half of all Sanger reads come back on the other strand."""
        from gsynth_engine.sequence import reverse_complement

        cloned = self.build(auth_client)
        forward = self.read_of(
            cloned["plasmid"], cloned["insert_start"] - 200, cloned["insert_end"] + 200,
        )
        response = auth_client.post(reverse(self.url_name), {
            "design": cloned["plasmid"],
            "reads": {"T7-R": reverse_complement(forward)},
            "region_start": cloned["insert_start"],
            "region_end": cloned["insert_end"],
        }, format="json")
        assert response.data["reads"][0]["reverse_complemented"] is True
        assert response.data["is_verified"]

    def test_a_point_mutation_is_reported_with_its_effect(self, auth_client):
        cloned = self.build(auth_client)
        plasmid = cloned["plasmid"]
        at = cloned["insert_start"] + 60
        replacement = "G" if plasmid[at] != "G" else "C"
        mutated = plasmid[:at] + replacement + plasmid[at + 1:]

        response = auth_client.post(reverse(self.url_name), {
            "design": plasmid,
            "reads": {"T7-F": self.read_of(
                mutated, cloned["insert_start"] - 200, cloned["insert_end"],
            )},
            "coding_start": cloned["insert_start"] + cloned["insert"]["orf_start"],
            "coding_end": cloned["insert_end"],
        }, format="json")
        assert not response.data["is_verified"]

        difference = response.data["differences"][0]
        assert difference["kind"] == "substitution"
        assert difference["position"] == at
        assert difference["found"] == replacement
        assert difference["silent"] in (True, False)
        assert "position" in difference["description"]

    def test_a_gap_in_coverage_is_reported(self, auth_client):
        """Half the insert read is not the insert verified."""
        cloned = self.build(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "design": cloned["plasmid"],
            "reads": {"short": self.read_of(
                cloned["plasmid"],
                cloned["insert_start"] - 100, cloned["insert_start"] + 60,
            )},
            "region_start": cloned["insert_start"],
            "region_end": cloned["insert_end"],
        }, format="json")
        assert not response.data["fully_covered"]
        assert not response.data["is_verified"]
        assert response.data["coverage"] < 100

    def test_one_unplaceable_read_does_not_sink_the_rest(self, auth_client):
        cloned = self.build(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "design": cloned["plasmid"],
            "reads": {
                "good": self.read_of(
                    cloned["plasmid"],
                    cloned["insert_start"] - 200, cloned["insert_end"],
                ),
                "junk": "ACGTACGTGGCCTTAA" * 25,
            },
        }, format="json")
        assert len(response.data["reads"]) == 1
        assert response.data["warnings"]

    def test_no_reads_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "design": "ACGT" * 100, "reads": {},
        }, format="json")
        assert response.status_code == 400
        assert "at least one read" in str(response.data)


@pytest.mark.django_db
class TestAlignEndpoint:
    url_name = "design-align"
    GENE = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA"

    def test_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": self.GENE,
        }).status_code == 401

    def test_identical_sequences_align_perfectly(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": self.GENE,
        })
        assert response.status_code == 200, response.data
        assert response.data["identity"] == 100.0
        assert response.data["gaps"] == 0

    def test_the_alignment_does_not_invent_or_lose_bases(self, auth_client):
        other = self.GENE[:20] + "GGG" + self.GENE[25:]
        response = auth_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": other,
        })
        assert response.data["top"].replace("-", "") == self.GENE
        assert response.data["bottom"].replace("-", "") == other

    def test_a_deletion_aligns_as_one_gap(self, auth_client):
        """Affine penalties exist to stop one event becoming four."""
        deleted = self.GENE[:20] + self.GENE[32:]
        response = auth_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": deleted,
        })
        assert response.data["gaps"] == 12
        runs = [run for run in response.data["bottom"].split("-") if run]
        assert len(runs) == 2

    def test_a_reversed_sequence_is_recognised(self, auth_client):
        from gsynth_engine.sequence import reverse_complement

        response = auth_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": reverse_complement(self.GENE),
        })
        assert response.data["reverse_complemented"] is True
        assert response.data["identity"] == 100.0
        assert any("other way round" in w for w in response.data["warnings"])

    def test_local_finds_only_the_shared_stretch(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "T" * 12 + self.GENE[:40] + "G" * 12,
            "second": "C" * 12 + self.GENE[:40] + "A" * 12,
            "mode": "local",
        })
        assert response.data["length"] == 40
        assert response.data["identity"] == 100.0

    def test_protein_alignment_uses_the_published_matrix(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "MTTSKLGKGLGYIGNN", "second": "MTTSRLGKGLGYVGNN",
            "is_protein": True,
        })
        assert response.data["identity"] < 100
        assert response.data["similarity"] == 100.0
        assert response.data["marks"].count(":") == 2

    def test_the_rows_are_ready_to_draw(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": self.GENE, "second": self.GENE,
        })
        row = response.data["rows"][0]
        assert set(row) >= {"top", "marks", "bottom", "top_start", "bottom_start"}
        assert "".join(r["top"] for r in response.data["rows"]) == response.data["top"]

    def test_an_empty_sequence_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "", "second": self.GENE,
        })
        assert response.status_code == 400

    def test_a_pair_too_large_is_refused_with_a_pointer(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "A" * 2500, "second": "A" * 2500,
        })
        assert response.status_code == 400
        assert "verification tool" in response.data["detail"]


@pytest.mark.django_db
class TestHybridizationEndpoint:
    url_name = "design-hybridize"

    def test_requires_authentication(self, api_client):
        response = api_client.post(reverse(self.url_name), {
            "first": "AATTATGC", "second": "GGCCGCAT",
        })
        assert response.status_code == 401

    def test_draws_antiparallel_strands_and_both_five_prime_overhangs(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "AATTATGC", "second": "GGCCGCAT",
        }, format="json")

        assert response.status_code == 200, response.data
        assert response.data["top"] == "AATTATGC    "
        assert response.data["bottom"] == "    TACGCCGG"
        assert response.data["marks"] == "    ||||    "
        assert response.data["complementarity"] == "exact"
        assert response.data["left_end"]["polarity"] == "5′"
        assert response.data["left_end"]["sequence"] == "AATT"
        assert response.data["right_end"]["polarity"] == "5′"
        assert response.data["right_end"]["sequence"] == "GGCC"

    def test_a_mismatch_is_explicit_and_not_given_a_perfect_match_tm(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "ATGCCGTA", "second": "TACGACAT",
        }, format="json")

        assert response.status_code == 200, response.data
        assert response.data["mismatches"] == 1
        assert "×" in response.data["marks"]
        assert response.data["tm_c"] is None
        assert response.data["predicted_state"] == "mismatches_not_thermodynamically_scored"

    def test_conditions_are_echoed_with_the_thermodynamic_result(self, auth_client):
        sequence = "ATGCCGTAGCTAGCTA"
        response = auth_client.post(reverse(self.url_name), {
            "first": sequence,
            "second": "TAGCTAGCTACGGCAT",
            "analysis_temperature_c": 37,
            "oligo_nM": 20_000,
            "na_mM": 75,
            "mg_mM": 1.5,
        }, format="json")

        assert response.status_code == 200, response.data
        assert response.data["tm_c"] is not None
        assert response.data["conditions"]["summary"] == (
            "20 µM total strand, 75 mM Na⁺, 1.5 mM Mg²⁺"
        )

    def test_invalid_dna_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "first": "ATUG", "second": "CAT",
        }, format="json")
        assert response.status_code == 400
        assert "not A, C, G or T" in response.data["detail"]


@pytest.mark.django_db
class TestPrimerExport:
    url_name = "design-primers-export"

    def build(self, auth_client):
        return auth_client.post(reverse("design-clone"), {
            "sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA",
        }).data

    def test_primers_export_as_csv(self, auth_client):
        import csv
        import io

        cloned = self.build(auth_client)
        response = auth_client.post(reverse(self.url_name), {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
            "name": "EntA",
        })
        assert response.status_code == 200
        assert 'filename="EntA_primers.csv"' in response["Content-Disposition"]

        rows = list(csv.DictReader(io.StringIO(response.content.decode())))
        assert rows
        assert set(rows[0]) >= {"Name", "Sequence (5'->3')", "Tm (°C)"}

    def test_primers_export_as_fasta(self, auth_client):
        import io

        from Bio import SeqIO

        cloned = self.build(auth_client)
        response = auth_client.post(
            reverse(self.url_name) + "?filetype=fasta",
            {
                "template": cloned["plasmid"],
                "target_start": cloned["insert_start"],
                "target_end": cloned["insert_end"],
                "name": "EntA",
            },
        )
        records = list(SeqIO.parse(io.StringIO(response.content.decode()), "fasta"))
        assert records
        assert records[0].id.startswith("EntA_")

    def test_the_export_matches_what_the_design_endpoint_returned(self, auth_client):
        """Two endpoints, one primer set — they must not drift."""
        import csv
        import io

        cloned = self.build(auth_client)
        payload = {
            "template": cloned["plasmid"],
            "target_start": cloned["insert_start"],
            "target_end": cloned["insert_end"],
            "name": "EntA",
        }
        designed = auth_client.post(reverse("design-primers"), payload).data
        exported = list(csv.DictReader(io.StringIO(
            auth_client.post(reverse(self.url_name), payload).content.decode()
        )))
        assert [row["Name"] for row in exported] == [
            p["name"] for p in designed["primers"]
        ]

    def test_export_requires_authentication(self, api_client):
        assert api_client.post(reverse(self.url_name), {
            "template": "ACGT" * 100, "target_start": 10, "target_end": 50,
        }).status_code == 401


@pytest.mark.django_db
class TestValidationAndJunctionViews:
    """The checks the user asked to be able to *see* rather than infer."""

    def cloned(self, auth_client, **extra):
        payload = {
            "sequence": LONG_INSERT, "vector_key": "pET-21a", "name": "EntA",
        }
        payload.update(extra)
        return auth_client.post(reverse("design-clone"), payload).data

    def test_each_seam_comes_back_as_a_drawable_duplex(self, auth_client):
        data = self.cloned(auth_client)
        views = data["junction_views"]
        assert len(views) == 2

        for view in views:
            assert view["compatible"], view["reason"]
            assert len(view["joined_top"]) == len(view["joined_bottom"])
            assert view["joined_pairs"].count("|") == len(view["joined_top"])

    def test_the_overhang_is_located_in_the_drawing(self, auth_client):
        """So "the overhangs match" is checkable rather than asserted."""
        data = self.cloned(auth_client)
        for view in data["junction_views"]:
            low, high = view["overhang_span"]
            assert view["joined_top"][low:high] == view["overhang"]
            assert high - low == len(view["overhang"])

    def test_both_ends_carry_the_overhang_before_ligation(self, auth_client):
        """On opposite strands — that is what lets them anneal."""
        view = self.cloned(auth_client)["junction_views"][0]
        assert view["left_top"].rstrip() != view["left_bottom"].rstrip()
        assert view["right_top"].lstrip() != view["right_bottom"].lstrip()

    def test_the_validation_list_states_each_check_separately(self, auth_client):
        """One banner collapses a dozen questions into a colour."""
        data = self.cloned(auth_client)
        checks = {row["check"]: row for row in data["validation"]}

        assert {"Overhangs are compatible", "Both strands pair everywhere",
                "Each enzyme cuts the vector once", "Orientation is forced",
                "Sites are regenerated at both seams",
                "Expression reading frame"} <= set(checks)
        assert all(row["passed"] for row in data["validation"]), checks
        assert all(row["status"] == "pass" for row in data["validation"]), checks
        assert all(row["detail"] for row in data["validation"])

    def test_a_failing_check_says_which_one(self, auth_client):
        """A design that fails on the frame must not look like one that
        fails on the ends."""
        data = self.cloned(
            auth_client,
            sequence="GGCTAAATCGTGGAACAGTGCTGCACCAGCTGCAGCCTGTACCAGCTGGAA",
        )
        failed = [row["check"] for row in data["validation"] if not row["passed"]]
        assert failed == ["Expression reading frame"]

    def test_a_missing_reading_frame_is_review_not_pass(self, auth_client):
        designed = self.cloned(auth_client)
        assembly = designed["assembly"]
        data = self.cloned(
            auth_client,
            sequence=assembly["construct_forward"],
            insert_reverse=assembly["construct_reverse"],
            pre_digested=True,
        )
        frame = next(
            row for row in data["validation"]
            if row["check"] == "Expression reading frame"
        )

        assert frame["status"] == "review"
        assert frame["passed"] is False
        assert data["reading_frame"]["start_source"] == "sequence_candidate"
        assert "not confirmed" in frame["detail"].lower()

    def test_design_handoff_preserves_and_confirms_the_frame(self, auth_client):
        designed = self.cloned(auth_client)
        assembly = designed["assembly"]
        data = self.cloned(
            auth_client,
            sequence=assembly["construct_forward"],
            insert_reverse=assembly["construct_reverse"],
            pre_digested=True,
            orf_start=assembly["ssd"]["orf_start"],
        )

        assert data["reading_frame"]["status"] == "pass"
        assert data["reading_frame"]["start_source"] == "declared"

    @pytest.mark.parametrize("declare_start", [True, False])
    def test_handoff_annotation_translates_the_same_n_terminus(self, auth_client, declare_start):
        """A duplex handoff must display Met and six His in the engine's frame."""
        from Bio.Seq import Seq

        designed = self.cloned(auth_client)
        assembly = designed["assembly"]
        extra = {"orf_start": assembly["ssd"]["orf_start"]} if declare_start else {}
        data = self.cloned(
            auth_client,
            sequence=assembly["construct_forward"],
            insert_reverse=assembly["construct_reverse"],
            pre_digested=True,
            **extra,
        )
        cassette = next(a for a in data["annotations"] if a["name"] == "EntA")
        start = cassette.get("translation_start", cassette["start"])
        end = cassette.get("translation_end", cassette["end"])
        coding = data["plasmid"][start:end]
        displayed = str(Seq(coding[:len(coding) // 3 * 3]).translate()).split("*")[0]

        assert displayed.startswith("MGSSHHHHHHSSG")
        assert data["protein"].startswith(displayed)
        assert start == data["reading_frame"]["translation_start"]
        assert data["reading_frame"]["start_source"] == (
            "declared" if declare_start else "sequence_candidate"
        )

    def test_uploaded_annotated_vector_is_validated_without_catalogue_identity(self, auth_client):
        from gsynth_engine import vectors

        record = vectors.sequence_of("pET-21a")
        custom_vector = "AAA" + record["sequence"]
        shifted = [
            {**feature, "start": feature["start"] + 3, "end": feature["end"] + 3}
            for feature in record["annotations"]
        ]
        data = self.cloned(
            auth_client,
            vector_key="",
            vector=custom_vector,
            vector_name="Imported expression vector",
            vector_annotations=shifted,
        )

        assert data["vector"]["recognised"] is False
        assert data["reading_frame"]["status"] == "pass"
        assert data["reading_frame"]["rbs_source"] == "annotation"
        assert data["reading_frame"]["promoter_source"] == "annotation"

    def test_unannotated_uploaded_vector_gets_reviewable_sequence_annotations(self, auth_client):
        from gsynth_engine import vectors

        record = vectors.sequence_of("pET-21a")
        data = self.cloned(
            auth_client,
            vector_key="",
            vector="AAA" + record["sequence"],
            vector_name="Unannotated expression vector",
            vector_annotations=[],
        )

        assert data["vector"]["recognised"] is False
        assert data["reading_frame"]["status"] == "review"
        assert data["reading_frame"]["rbs_source"] == "sequence_motif"
        assert data["reading_frame"]["promoter_source"] == "sequence_motif"
        assert any(feature.get("inferred") for feature in data["annotations"])

    def test_expression_frame_exposes_vector_context_evidence(self, auth_client):
        frame = self.cloned(auth_client)["reading_frame"]

        assert frame["confirmed"] is True
        assert frame["status"] == "pass"
        assert frame["start_codon"] == "ATG"
        assert frame["rbs_name"] == "RBS"
        assert frame["rbs_spacing_nt"] == 8
        assert frame["promoter_name"] == "T7 promoter"
        assert frame["stop_context"] in {"insert", "vector"}
        assert {check["code"] for check in frame["checks"]} == {
            "FRAME_START_CODON",
            "FRAME_RBS_CONTEXT",
            "FRAME_PROMOTER_CONTEXT",
            "FRAME_JUNCTION_PHASE",
            "FRAME_TRANSLATED_TERMINUS",
        }

    def test_restriction_sites_are_annotated_on_the_map(self, auth_client):
        data = self.cloned(auth_client)
        sites = {s["name"]: s for s in data["restriction_sites"]}

        assert "NdeI" in sites and "XhoI" in sites
        assert sites["NdeI"]["used"] is True
        assert any(s["cuts"] == 1 and not s["used"] for s in data["restriction_sites"])
        assert any(s["cuts"] > 1 and not s["used"] for s in data["restriction_sites"])

    def test_hindiii_is_returned_when_it_is_used_for_cloning(self, auth_client):
        data = self.cloned(
            auth_client,
            left_enzyme="NdeI",
            right_enzyme="HindIII",
        )
        sites = [site for site in data["restriction_sites"] if site["name"] == "HindIII"]

        assert data["is_clonable"], data["problems"]
        assert len(sites) == 1
        assert sites[0]["used"] is True
        assert sites[0]["recognition"] == "AAGCTT"

    def test_cloning_response_includes_the_complete_diagnostic_digest_gel(self, auth_client):
        data = self.cloned(auth_client, right_enzyme="HindIII")
        gel = data["gel"]
        lane = gel["lanes"][0]
        sizes = [band["size_bp"] for band in lane["bands"]]

        assert gel["prediction_only"] is True
        assert lane["name"] == "NdeI + HindIII"
        assert len(sizes) == 2
        assert sum(sizes) == data["length"]

    def test_every_catalogued_restriction_site_has_drawable_coordinates(self, auth_client):
        data = self.cloned(auth_client)
        for site in data["restriction_sites"]:
            assert site["start"] >= 0
            assert site["end"] > site["start"]
            assert len(site["recognition"]) == site["end"] - site["start"]

    def test_an_annotated_site_really_is_where_it_says(self, auth_client):
        data = self.cloned(auth_client)
        from gsynth_engine.sequence import reverse_complement

        plasmid = data["plasmid"]
        length = len(plasmid)
        for site in data["restriction_sites"]:
            # A circular molecule has no beginning: a site can straddle it.
            found = "".join(
                plasmid[i % length] for i in range(site["start"], site["end"])
            )
            assert found in (site["recognition"],
                             reverse_complement(site["recognition"])), site["name"]
            assert site["wraps"] == (site["end"] > length)

    def test_multi_cutters_are_returned_for_the_explicit_all_sites_filter(self, auth_client):
        """The client may hide these by default, but the API must not erase them."""
        data = self.cloned(auth_client)
        assert any(site["cuts"] > 1 and not site["used"] for site in data["restriction_sites"])


@pytest.mark.django_db
class TestTraceVerifyEndpoint:
    """Uploading ABIF/SCF traces rather than pasting the letters off them.

    The trace is what separates a mutation from a bad call, so the response
    has to carry the confidence — not just the difference.
    """

    url_name = "design-verify-traces"

    DESIGN = (
        "ATGGCTAGCAAAGAACTGGTTACCGCTCTGTATCTGGTGTGCGGCGAACGCGGCTTTTTCTACACCCCG"
        "AAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGCGGCCCGGGCGCGGGC"
        "AGCCTGCAGCCGCTGGCGCTGGAAGGCAGCCTGCAGAAACGCGGCATCGTGGAACAGTGCTGCACCAGC"
        "ATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGCTTTGTGAACCAGCATCTGTGCGGCAGC"
    )

    def _upload(self, sequence, quality=None, name="fwd.ab1"):
        from django.core.files.uploadedfile import SimpleUploadedFile

        from gsynth_engine.tests.test_chromatogram import build_ab1

        blob = build_ab1(sequence, quality if quality is not None else [45] * len(sequence))
        return SimpleUploadedFile(name, blob, content_type="application/octet-stream")

    def _upload_scf(self, sequence, quality=None, name="facility-export.ab1"):
        from django.core.files.uploadedfile import SimpleUploadedFile

        from gsynth_engine.tests.test_chromatogram import build_scf

        blob = build_scf(sequence, quality if quality is not None else [45] * len(sequence))
        return SimpleUploadedFile(name, blob, content_type="application/octet-stream")

    def _changed_read(self, at=100):
        piece = list(self.DESIGN[30:230])
        piece[at] = {"A": "G", "G": "A", "C": "T", "T": "C"}[piece[at]]
        return "".join(piece)

    def test_requires_authentication(self, api_client):
        response = api_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "traces": [self._upload(self.DESIGN[30:230])],
        }, format="multipart")
        assert response.status_code == 401

    def test_a_clean_trace_verifies_the_design(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(self.DESIGN[30:230])],
            "region_start": 30, "region_end": 230,
        }, format="multipart")
        assert response.status_code == 200, response.data
        assert response.data["is_verified"] is True
        assert response.data["differences"] == []
        assert response.data["traces"][0]["mean_quality"] == 45.0

    def test_an_scf_trace_with_an_ab1_name_is_detected_by_content(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload_scf(self.DESIGN[30:230])],
            "region_start": 30, "region_end": 230,
        }, format="multipart")
        assert response.status_code == 200, response.data
        assert response.data["is_verified"] is True
        assert response.data["differences"] == []
        assert response.data["traces"][0]["mean_quality"] == 45.0

    def test_a_poor_peak_is_returned_as_unconfident(self, auth_client):
        """Same letters as a real mutation; the response must not conflate them."""
        read = self._changed_read()
        quality = [45] * len(read)
        quality[100] = 7

        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read, quality)],
        }, format="multipart")
        assert response.status_code == 200, response.data
        difference = response.data["differences"][0]
        assert difference["quality"] == 7
        assert difference["confident"] is False
        assert "Q7" in difference["description"]

    def test_a_clean_peak_is_returned_as_confident(self, auth_client):
        read = self._changed_read()
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read)],
        }, format="multipart")
        difference = response.data["differences"][0]
        assert difference["quality"] == 45
        assert difference["confident"] is True

    def test_the_peaks_around_each_difference_come_back(self, auth_client):
        """So it can be looked at, not taken on trust. Only the window
        travels — whole traces would be megabytes per read."""
        read = self._changed_read()
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read)],
        }, format="multipart")

        windows = response.data["trace_windows"]
        assert len(windows) == 1
        window = windows[0]
        assert window["read"] == "fwd.ab1"
        assert set(window["traces"]) == set("ACGT")
        assert window["bases"], "the window must name the bases it spans"
        assert any(b["index"] == window["centre"] for b in window["bases"])

    def test_reference_aligned_trace_track_uses_only_verified_bases(self, auth_client):
        """The overview must not present quality-discarded ends as evidence."""
        read = self.DESIGN[30:230]
        quality = [2] * 12 + [45] * (len(read) - 24) + [2] * 12
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read, quality)],
        }, format="multipart")

        assert response.status_code == 200, response.data
        track = response.data["trace_tracks"][0]
        aligned = response.data["reads"][0]
        assert track["read"] == "fwd.ab1"
        assert track["reference_start"] == aligned["start"]
        assert track["sequence"] == read[12:-12]
        assert len(track["qualities"]) == len(track["sequence"])
        assert len(track["peaks"]) == len(track["sequence"])
        assert set(track["traces"]) == set("ACGT")

    def test_requested_quality_cutoff_controls_the_verified_alignment(self, auth_client):
        """The API must not display one cutoff while silently using Q13."""
        read = self.DESIGN[30:230]
        quality = [2] * 12 + [45] * (len(read) - 24) + [2] * 12

        strict = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read, quality)], "trim_quality": 13,
        }, format="multipart")
        permissive = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read, quality)], "trim_quality": 0,
        }, format="multipart")

        assert strict.status_code == permissive.status_code == 200
        assert strict.data["reads"][0]["trimmed_start"] == 12
        assert permissive.data["reads"][0]["trimmed_start"] == 0
        assert permissive.data["coverage"] > strict.data["coverage"]

    def test_reverse_track_is_complemented_into_reference_orientation(self, auth_client):
        from gsynth_engine.sequence import reverse_complement

        forward = self.DESIGN[30:230]
        reverse_read = reverse_complement(forward)
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(reverse_read, name="rev.ab1")],
        }, format="multipart")

        assert response.status_code == 200, response.data
        track = response.data["trace_tracks"][0]
        assert track["reverse_complemented"] is True
        assert track["sequence"] == forward
        assert track["peaks"] == sorted(track["peaks"])

    def test_forward_reverse_calls_are_assembled_before_coverage(self, auth_client):
        from gsynth_engine.sequence import reverse_complement

        forward = self.DESIGN[:190]
        reverse_read = reverse_complement(self.DESIGN[120:])
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [
                self._upload(forward, name="forward.ab1"),
                self._upload(reverse_read, name="reverse.ab1"),
            ],
            "region_start": 0, "region_end": len(self.DESIGN),
        }, format="multipart")

        assert response.status_code == 200, response.data
        consensus = response.data["raw_consensus"]
        assert consensus["coverage"] == 100.0
        assert consensus["identity"] == 100.0
        assert consensus["fully_covered"] is True
        assert consensus["bidirectional_overlap"] > 0
        assert consensus["bidirectional_agreement"] == 100.0

    def test_matches_what_the_engine_measured(self, auth_client):
        """The response reports the engine's numbers, not its own."""
        from gsynth_engine.chromatogram import read_ab1
        from gsynth_engine.tests.test_chromatogram import build_ab1
        from gsynth_engine.verify import verify

        read = self._changed_read()
        quality = [45] * len(read)
        quality[100] = 9
        blob = build_ab1(read, quality)

        trace = read_ab1(blob, name="fwd.ab1")
        expected = verify(self.DESIGN, {"fwd.ab1": read}, circular=False,
                          traces={"fwd.ab1": trace})

        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN, "circular": False,
            "traces": [self._upload(read, quality)],
        }, format="multipart")

        assert [d["position"] for d in response.data["differences"]] == [
            d.position for d in expected.differences
        ]
        assert [d["quality"] for d in response.data["differences"]] == [
            d.quality for d in expected.differences
        ]

    def test_a_file_that_is_not_a_trace_is_refused_in_plain_words(self, auth_client):
        from django.core.files.uploadedfile import SimpleUploadedFile

        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN,
            "traces": [SimpleUploadedFile("scan.pdf", b"%PDF-1.4" + b"\x00" * 500)],
        }, format="multipart")
        assert response.status_code == 400
        assert "not a supported Sanger trace" in response.data["detail"]

    def test_no_trace_at_all_is_refused(self, auth_client):
        response = auth_client.post(reverse(self.url_name), {
            "design": self.DESIGN,
        }, format="multipart")
        assert response.status_code == 400
