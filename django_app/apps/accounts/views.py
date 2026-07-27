"""Auth endpoints: register, login, logout, me, change-password.

Refresh is provided by simplejwt directly (wired in urls.py).
"""
from rest_framework import generics, permissions, status
from rest_framework.response import Response
from rest_framework.throttling import ScopedRateThrottle
from rest_framework.views import APIView
from rest_framework_simplejwt.views import TokenObtainPairView

from apps.accounts.serializers import (
    ChangePasswordSerializer,
    LogoutSerializer,
    RegisterSerializer,
    UserSerializer,
)


class RegisterView(generics.CreateAPIView):
    """POST /api/auth/register — create a new account. Public, rate-limited."""

    permission_classes = (permissions.AllowAny,)
    serializer_class = RegisterSerializer
    throttle_classes = (ScopedRateThrottle,)
    throttle_scope = "register"

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        user = serializer.save()
        return Response(UserSerializer(user).data, status=status.HTTP_201_CREATED)


class LoginView(TokenObtainPairView):
    """POST /api/auth/login — email + password → access & refresh tokens.

    Rate-limited on its own scope: this is the endpoint an attacker hammers
    to guess passwords.
    """

    throttle_classes = (ScopedRateThrottle,)
    throttle_scope = "login"


class LogoutView(APIView):
    """POST /api/auth/logout — blacklist the supplied refresh token."""

    def post(self, request):
        serializer = LogoutSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response(status=status.HTTP_205_RESET_CONTENT)


class LogoutAllView(APIView):
    """POST /api/auth/logout-all — revoke every session for this user.

    The "I think someone has my password" button: bumps the token version
    and blacklists all outstanding refresh tokens.
    """

    def post(self, request):
        revoked = request.user.revoke_all_tokens()
        return Response({"detail": "All sessions revoked.", "sessions_revoked": revoked})


class MeView(APIView):
    """GET /api/auth/me — profile of the currently signed-in user."""

    def get(self, request):
        return Response(UserSerializer(request.user).data)


class ChangePasswordView(APIView):
    """POST /api/auth/change-password — requires current password.

    Succeeding here revokes every existing token for the account, so the
    client must sign in again with the new password.
    """

    def post(self, request):
        serializer = ChangePasswordSerializer(
            data=request.data, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response({
            "detail": "Password updated. All existing sessions were signed out.",
        })
