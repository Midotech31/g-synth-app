from rest_framework_simplejwt.authentication import JWTAuthentication
from rest_framework_simplejwt.exceptions import AuthenticationFailed

TOKEN_VERSION_CLAIM = "ver"


class VersionedJWTAuthentication(JWTAuthentication):


    def get_user(self, validated_token):
        user = super().get_user(validated_token)


        token_version = validated_token.get(TOKEN_VERSION_CLAIM)
        if token_version is None or int(token_version) != user.token_version:
            raise AuthenticationFailed(
                "This token was revoked. Please sign in again.",
                code="token_revoked",
            )
        return user
