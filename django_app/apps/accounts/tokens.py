from rest_framework_simplejwt.tokens import RefreshToken

from apps.accounts.authentication import TOKEN_VERSION_CLAIM


class VersionedRefreshToken(RefreshToken):


    @classmethod
    def for_user(cls, user):
        token = super().for_user(user)
        token[TOKEN_VERSION_CLAIM] = user.token_version
        return token
