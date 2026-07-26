"""JWT authentication that honours credential revocation.

A JWT is self-contained: once signed, it stays valid until it expires, no
matter what happens to the account. Plain `JWTAuthentication` therefore lets
a stolen access token keep working *after* the victim changes their password
— exactly when it must not.

Every token we issue carries a `ver` claim holding the user's
`token_version` at issue time (see `VersionedTokenObtainPairSerializer`).
Changing credentials bumps that counter, so previously issued tokens no
longer match and are refused here.
"""
from rest_framework_simplejwt.authentication import JWTAuthentication
from rest_framework_simplejwt.exceptions import AuthenticationFailed

TOKEN_VERSION_CLAIM = "ver"


class VersionedJWTAuthentication(JWTAuthentication):
    """`JWTAuthentication` + a token-version check against the user record."""

    def get_user(self, validated_token):
        user = super().get_user(validated_token)

        # Tokens minted before this claim existed have no version and are
        # treated as stale — fail closed rather than trusting them.
        token_version = validated_token.get(TOKEN_VERSION_CLAIM)
        if token_version is None or int(token_version) != user.token_version:
            raise AuthenticationFailed(
                "This token was revoked. Please sign in again.",
                code="token_revoked",
            )
        return user
