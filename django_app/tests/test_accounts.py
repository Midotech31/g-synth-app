"""Auth flow: register, login, refresh, me, change-password."""
import pytest
from django.contrib.auth import get_user_model
from django.urls import reverse

User = get_user_model()


@pytest.mark.django_db
class TestRegister:
    def test_creates_user(self, api_client):
        r = api_client.post(reverse("auth-register"), {
            "email": "new@example.com", "name": "New User",
            "password": "GoodPass123", "password2": "GoodPass123",
        })
        assert r.status_code == 201
        assert r.data["email"] == "new@example.com"
        assert User.objects.filter(email="new@example.com").exists()

    def test_rejects_short_password(self, api_client):
        r = api_client.post(reverse("auth-register"), {
            "email": "a@b.com", "name": "X",
            "password": "short", "password2": "short",
        })
        assert r.status_code == 400

    def test_rejects_password_mismatch(self, api_client):
        r = api_client.post(reverse("auth-register"), {
            "email": "a@b.com", "name": "X",
            "password": "GoodPass123", "password2": "DifferentPass123",
        })
        assert r.status_code == 400

    def test_rejects_duplicate_email(self, api_client, user):
        r = api_client.post(reverse("auth-register"), {
            "email": user.email, "name": "X",
            "password": "GoodPass123", "password2": "GoodPass123",
        })
        assert r.status_code == 400


@pytest.mark.django_db
class TestLogin:
    def test_returns_tokens(self, api_client, user):
        r = api_client.post(reverse("auth-login"), {
            "email": user.email, "password": "StrongPass123",
        })
        assert r.status_code == 200
        assert "access" in r.data and "refresh" in r.data

    def test_wrong_password_rejected(self, api_client, user):
        r = api_client.post(reverse("auth-login"), {
            "email": user.email, "password": "WrongPassword",
        })
        assert r.status_code == 401


@pytest.mark.django_db
class TestMe:
    def test_requires_auth(self, api_client):
        assert api_client.get(reverse("auth-me")).status_code == 401

    def test_returns_profile(self, auth_client, user):
        r = auth_client.get(reverse("auth-me"))
        assert r.status_code == 200
        assert r.data["email"] == user.email
        assert r.data["name"] == user.name


@pytest.mark.django_db
class TestChangePassword:
    def test_flow(self, auth_client, user):
        r = auth_client.post(reverse("auth-change-password"), {
            "current_password": "StrongPass123",
            "new_password": "AnotherGood456",
            "new_password2": "AnotherGood456",
        })
        assert r.status_code == 200
        user.refresh_from_db()
        assert user.check_password("AnotherGood456")

    def test_wrong_current_rejected(self, auth_client):
        r = auth_client.post(reverse("auth-change-password"), {
            "current_password": "Wrong",
            "new_password": "AnotherGood456",
            "new_password2": "AnotherGood456",
        })
        assert r.status_code == 400
