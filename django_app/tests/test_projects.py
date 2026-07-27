"""Project CRUD, and — crucially — user isolation."""
import pytest

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
            "name": "Insulin v1", "module": "crispr_designer",
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
