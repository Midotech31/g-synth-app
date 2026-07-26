"""Shared pytest fixtures."""
import pytest
from django.contrib.auth import get_user_model
from rest_framework.test import APIClient

User = get_user_model()


@pytest.fixture
def api_client():
    return APIClient()


@pytest.fixture
def user(db):
    return User.objects.create_user(
        email="alice@example.com",
        password="StrongPass123",
        name="Alice Test",
    )


@pytest.fixture
def other_user(db):
    return User.objects.create_user(
        email="bob@example.com",
        password="StrongPass123",
        name="Bob Other",
    )


@pytest.fixture
def auth_client(api_client, user):
    """APIClient pre-authenticated as `user`.

    Mints via VersionedRefreshToken, the same path the login endpoint uses —
    a bare `RefreshToken.for_user()` would omit the `ver` claim and the
    request would be (correctly) rejected.
    """
    from apps.accounts.tokens import VersionedRefreshToken
    token = str(VersionedRefreshToken.for_user(user).access_token)
    api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {token}")
    return api_client
