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
            "name": "Insulin v1", "module": "merzoug_assembly",
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
        assert Project.objects.filter(id=pid, user=user).exists()

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
