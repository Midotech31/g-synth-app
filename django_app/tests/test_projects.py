"""Project CRUD, and — crucially — user isolation."""
import pytest
from django.urls import reverse

from apps.projects.models import Project


@pytest.mark.django_db
class TestProjectCrud:
    def test_requires_auth(self, api_client):
        assert api_client.get("/api/projects/").status_code == 401

    def test_list_empty(self, auth_client):
        r = auth_client.get("/api/projects/")
        assert r.status_code == 200
        assert r.data["count"] == 0

    def test_create_and_retrieve(self, auth_client, user):
        r = auth_client.post("/api/projects/", {
            "name": "Insulin v1", "module": "extended_sequence_design",
            "sequence": "ATGAAACGT",
            "notes": "first test",
            "data": {"guides": []},
        }, format="json")
        assert r.status_code == 201, r.data
        pid = r.data["id"]

        got = auth_client.get(f"/api/projects/{pid}/")
        assert got.status_code == 200
        assert got.data["name"] == "Insulin v1"
        assert got.data["sequence"] == "ATGAAACGT"
        assert got.data["provenance"]["schema"] == "gsynth.provenance/v1"
        assert len(got.data["provenance"]["output_sha256"]) == 64
        assert Project.objects.filter(id=pid, user=user).exists()

    def test_provenance_is_server_generated_and_read_only(self, auth_client):
        response = auth_client.post("/api/projects/", {
            "name": "traceable",
            "module": "ssd",
            "sequence": "ATGAAATAA",
            "provenance": {"output_sha256": "forged"},
        }, format="json")
        assert response.status_code == 201
        assert response.data["provenance"]["output_sha256"] != "forged"

        patched = auth_client.patch(
            f"/api/projects/{response.data['id']}/",
            {"provenance": {"output_sha256": "changed"}},
            format="json",
        )
        assert patched.data["provenance"]["output_sha256"] == response.data["provenance"]["output_sha256"]

    def test_update(self, auth_client, user):
        p = Project.objects.create(user=user, name="A", notes="")
        r = auth_client.patch(f"/api/projects/{p.id}/", {"notes": "updated"}, format="json")
        assert r.status_code == 200
        p.refresh_from_db()
        assert p.notes == "updated"

    def test_delete(self, auth_client, user):
        p = Project.objects.create(user=user, name="A")
        r = auth_client.delete(f"/api/projects/{p.id}/")
        assert r.status_code == 204
        assert not Project.objects.filter(id=p.id).exists()


@pytest.mark.django_db
class TestEditableAnnotations:
    def test_replaces_annotations_without_overwriting_other_project_data(self, auth_client, user):
        project = Project.objects.create(
            user=user,
            name="new insert",
            sequence="ATGAAACCCGGGTAA",
            data={"topology": "linear", "preflight": {"verdict": "ready"}},
        )
        annotation = {
            "name": "My insulin insert",
            "type": "CDS",
            "start": 0,
            "end": 15,
            "direction": 1,
            "color": "#0E6E77",
            "translation_start": 0,
            "translation_end": 15,
        }

        response = auth_client.patch(
            f"/api/projects/{project.id}/annotations/",
            {"annotations": [annotation]},
            format="json",
        )

        assert response.status_code == 200, response.data
        assert response.data["data"]["annotations"] == [annotation]
        assert response.data["data"]["preflight"] == {"verdict": "ready"}

    def test_rejects_invalid_linear_and_wrapped_coordinates(self, auth_client, user):
        project = Project.objects.create(
            user=user,
            name="linear",
            sequence="A" * 20,
            data={"topology": "linear"},
        )
        invalid = {
            "name": "past end", "type": "misc_feature", "start": 18, "end": 25,
            "direction": 1, "color": "#0E6E77",
        }

        response = auth_client.patch(
            f"/api/projects/{project.id}/annotations/",
            {"annotations": [invalid]},
            format="json",
        )

        assert response.status_code == 400
        assert "linear feature" in str(response.data).lower()

    def test_accepts_a_feature_crossing_the_circular_origin(self, auth_client, user):
        project = Project.objects.create(
            user=user,
            name="plasmid",
            sequence="A" * 20,
            data={"topology": "circular"},
        )
        wrapped = {
            "name": "origin insert", "type": "CDS", "start": 17, "end": 24,
            "direction": -1, "color": "#3F7A52",
        }

        response = auth_client.patch(
            f"/api/projects/{project.id}/annotations/",
            {"annotations": [wrapped]},
            format="json",
        )

        assert response.status_code == 200, response.data
        assert response.data["data"]["annotations"][0]["end"] == 24

    def test_annotation_updates_are_scoped_to_the_owner(self, auth_client, other_user):
        project = Project.objects.create(
            user=other_user, name="private", sequence="ACGT", data={"topology": "linear"}
        )
        response = auth_client.patch(
            f"/api/projects/{project.id}/annotations/", {"annotations": []}, format="json",
        )
        assert response.status_code == 404

    def test_detects_common_features_without_saving_them(self, auth_client, user):
        project = Project.objects.create(
            user=user,
            name="unannotated insert",
            sequence="AAATAATACGACTCACTATAGGGCC",
            data={"topology": "linear", "annotations": []},
        )

        response = auth_client.get(f"/api/projects/{project.id}/detect-common-features/")

        assert response.status_code == 200, response.data
        promoter = next(
            item for item in response.data["matches"]
            if item["annotation"]["name"] == "T7 promoter"
        )
        assert promoter["annotation"]["start"] == 3
        assert promoter["annotation"]["direction"] == 1
        project.refresh_from_db()
        assert project.data["annotations"] == []

    def test_does_not_propose_a_known_feature_already_annotated(self, auth_client, user):
        motif = "TAATACGACTCACTATAGGG"
        project = Project.objects.create(
            user=user,
            name="annotated",
            sequence=motif,
            data={
                "topology": "linear",
                "annotations": [{
                    "name": "T7 promoter", "type": "promoter", "start": 0,
                    "end": len(motif), "direction": 1, "color": "#C97634",
                }],
            },
        )

        response = auth_client.get(f"/api/projects/{project.id}/detect-common-features/")

        assert not [
            item for item in response.data["matches"]
            if item["annotation"]["name"] == "T7 promoter"
        ]


@pytest.mark.django_db
class TestUserIsolation:
    """A user must never see or touch another user's projects."""

    def test_list_excludes_other_users_projects(self, auth_client, user, other_user):
        Project.objects.create(user=other_user, name="hidden")
        Project.objects.create(user=user, name="mine")
        r = auth_client.get("/api/projects/")
        names = [p["name"] for p in r.data["results"]]
        assert names == ["mine"]

    def test_cannot_retrieve_other_users_project(self, auth_client, other_user):
        theirs = Project.objects.create(user=other_user, name="hidden")
        assert auth_client.get(f"/api/projects/{theirs.id}/").status_code == 404

    def test_cannot_delete_other_users_project(self, auth_client, other_user):
        theirs = Project.objects.create(user=other_user, name="hidden")
        r = auth_client.delete(f"/api/projects/{theirs.id}/")
        assert r.status_code == 404
        assert Project.objects.filter(id=theirs.id).exists()


@pytest.mark.django_db
class TestProjectExport:
    """Saved work that cannot be taken out is not really saved."""

    def test_a_project_exports_as_genbank(self, auth_client, user):
        project = Project.objects.create(
            user=user, name="pGS EntA", module="cloning",
            sequence="ATGGGTTCTTCTCACCACCACCACCACCACTAA",
            notes="EntA in pET-21a(+)",
            data={
                "topology": "circular",
                "annotations": [
                    {"name": "6xHis", "type": "CDS", "start": 12, "end": 30,
                     "direction": 1},
                ],
            },
        )
        response = auth_client.get(reverse("project-export", args=[project.id]))
        assert response.status_code == 200
        assert 'filename="pGS_EntA.gb"' in response["Content-Disposition"]

        import io

        from Bio import SeqIO

        record = SeqIO.read(io.StringIO(response.content.decode()), "genbank")
        assert str(record.seq).upper() == project.sequence
        assert record.annotations["topology"] == "circular"
        labels = [
            f.qualifiers.get("label", [""])[0]
            for f in record.features if f.type != "source"
        ]
        assert labels == ["6xHis"]

    def test_a_manually_named_insert_survives_into_genbank(self, auth_client, user):
        project = Project.objects.create(
            user=user,
            name="editable",
            sequence="ATGAAACCCGGGTAA",
            data={"topology": "linear", "annotations": []},
        )
        annotation = {
            "name": "Insulin glargine insert",
            "type": "CDS",
            "start": 0,
            "end": 15,
            "direction": 1,
            "color": "#3F7A52",
            "translation_start": 0,
            "translation_end": 15,
        }
        updated = auth_client.patch(
            f"/api/projects/{project.id}/annotations/",
            {"annotations": [annotation]},
            format="json",
        )
        assert updated.status_code == 200

        response = auth_client.get(reverse("project-export", args=[project.id]))
        assert '/label="Insulin glargine insert"' in response.content.decode()

    def test_it_can_export_as_fasta(self, auth_client, user):
        project = Project.objects.create(
            user=user, name="insert", module="ssd", sequence="ATGAAATAA",
        )
        response = auth_client.get(
            reverse("project-export", args=[project.id]) + "?filetype=fasta"
        )
        assert response.status_code == 200
        assert response.content.decode().startswith(">insert")

    def test_export_is_scoped_to_the_owner(self, auth_client, other_user):
        """The same rule as every other project endpoint."""
        theirs = Project.objects.create(
            user=other_user, name="theirs", module="ssd", sequence="ATG",
        )
        response = auth_client.get(reverse("project-export", args=[theirs.id]))
        assert response.status_code == 404

    def test_export_requires_authentication(self, api_client, user):
        mine = Project.objects.create(
            user=user, name="mine", module="ssd", sequence="ATG",
        )
        assert api_client.get(
            reverse("project-export", args=[mine.id])
        ).status_code == 401
