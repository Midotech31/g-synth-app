"""Token minting with the revocation claim attached.

Every token must carry `ver` (the user's `token_version` at issue time) or
`VersionedJWTAuthentication` will refuse it. Putting that in the login
serializer alone is a footgun: any other call site — a management command,
an invite flow, a test fixture — that reaches for
`RefreshToken.for_user(user)` mints a token that is silently dead on
arrival. Minting goes through this class instead.
"""
from rest_framework_simplejwt.tokens import RefreshToken

from apps.accounts.authentication import TOKEN_VERSION_CLAIM


class VersionedRefreshToken(RefreshToken):
    """`RefreshToken` that stamps the user's current token version."""

    @classmethod
    def for_user(cls, user):
        token = super().for_user(user)
        token[TOKEN_VERSION_CLAIM] = user.token_version
        return token
