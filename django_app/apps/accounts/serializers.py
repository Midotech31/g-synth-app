"""DRF serializers for user registration, profile, and token issuance."""
from django.contrib.auth import get_user_model
from django.contrib.auth.password_validation import validate_password
from rest_framework import serializers
from rest_framework_simplejwt.serializers import TokenObtainPairSerializer

from apps.accounts.tokens import VersionedRefreshToken

User = get_user_model()


class VersionedTokenObtainPairSerializer(TokenObtainPairSerializer):
    """Issues tokens carrying the user's current `token_version`.

    simplejwt copies custom claims from a refresh token onto the access
    tokens minted from it, so a refresh token issued before a password
    change produces access tokens carrying the stale version — which
    `VersionedJWTAuthentication` then rejects.
    """

    token_class = VersionedRefreshToken

    @classmethod
    def get_token(cls, user):
        return VersionedRefreshToken.for_user(user)


class UserSerializer(serializers.ModelSerializer):
    """Read-only projection returned by /me and /register."""

    class Meta:
        model = User
        fields = ("id", "email", "name", "date_joined")
        read_only_fields = ("id", "date_joined")


class RegisterSerializer(serializers.ModelSerializer):
    """Sign-up payload: email + name + password (validated) + confirm."""

    password = serializers.CharField(
        write_only=True, required=True, validators=[validate_password],
    )
    password2 = serializers.CharField(write_only=True, required=True)

    class Meta:
        model = User
        fields = ("email", "name", "password", "password2")

    def validate(self, attrs):
        if attrs["password"] != attrs["password2"]:
            raise serializers.ValidationError(
                {"password2": "The two passwords do not match."}
            )
        return attrs

    def create(self, validated_data):
        validated_data.pop("password2")
        password = validated_data.pop("password")
        return User.objects.create_user(password=password, **validated_data)


class ChangePasswordSerializer(serializers.Serializer):
    """Change password: needs the current password + new + confirm."""

    current_password = serializers.CharField(write_only=True, required=True)
    new_password = serializers.CharField(
        write_only=True, required=True, validators=[validate_password],
    )
    new_password2 = serializers.CharField(write_only=True, required=True)

    def validate(self, attrs):
        if attrs["new_password"] != attrs["new_password2"]:
            raise serializers.ValidationError(
                {"new_password2": "The two new passwords do not match."}
            )
        user = self.context["request"].user
        if not user.check_password(attrs["current_password"]):
            raise serializers.ValidationError(
                {"current_password": "Current password is incorrect."}
            )
        if attrs["new_password"] == attrs["current_password"]:
            raise serializers.ValidationError(
                {"new_password": "The new password must differ from the current one."}
            )
        return attrs

    def save(self, **kwargs):
        user = self.context["request"].user
        user.set_password(self.validated_data["new_password"])
        user.save(update_fields=["password"])
        # Cut off every token issued under the old password.
        user.revoke_all_tokens()
        return user


class LogoutSerializer(serializers.Serializer):
    """Blacklist one refresh token (single-device sign-out)."""

    refresh = serializers.CharField(required=True)

    def save(self, **kwargs):
        from rest_framework_simplejwt.exceptions import TokenError
        from rest_framework_simplejwt.tokens import RefreshToken

        try:
            RefreshToken(self.validated_data["refresh"]).blacklist()
        except TokenError as exc:
            raise serializers.ValidationError(
                {"refresh": "Token is invalid or already revoked."}
            ) from exc
