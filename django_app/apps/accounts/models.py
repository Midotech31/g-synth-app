"""Custom user model — email as the login identifier, no separate username."""
from django.contrib.auth.models import AbstractBaseUser, BaseUserManager, PermissionsMixin
from django.db import models
from django.utils import timezone


class UserManager(BaseUserManager):
    """Manager that normalises the email + hashes the password."""

    use_in_migrations = True

    def _create_user(self, email: str, password: str | None, **extra_fields):
        if not email:
            raise ValueError("Users must have an email address.")
        email = self.normalize_email(email).lower()
        user = self.model(email=email, **extra_fields)
        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_user(self, email: str, password: str | None = None, **extra_fields):
        extra_fields.setdefault("is_staff", False)
        extra_fields.setdefault("is_superuser", False)
        return self._create_user(email, password, **extra_fields)

    def create_superuser(self, email: str, password: str | None = None, **extra_fields):
        extra_fields.setdefault("is_staff", True)
        extra_fields.setdefault("is_superuser", True)
        if not extra_fields["is_staff"] or not extra_fields["is_superuser"]:
            raise ValueError("Superuser must have is_staff=True and is_superuser=True.")
        return self._create_user(email, password, **extra_fields)


class User(AbstractBaseUser, PermissionsMixin):
    """G-Synth user account."""

    email = models.EmailField(unique=True, db_index=True)
    name = models.CharField(max_length=255, blank=True, default="")
    is_active = models.BooleanField(default=True)
    is_staff = models.BooleanField(default=False)
    date_joined = models.DateTimeField(default=timezone.now)

    # Bumped whenever credentials change. Every issued JWT carries the value
    # it was minted with; VersionedJWTAuthentication rejects any token whose
    # version is stale. This is what makes "change my password" actually cut
    # off an attacker holding a previously stolen access token — JWTs are
    # stateless and cannot otherwise be revoked before they expire.
    token_version = models.PositiveIntegerField(default=0)

    objects = UserManager()

    USERNAME_FIELD = "email"
    REQUIRED_FIELDS: list[str] = []

    class Meta:
        db_table = "gsynth_user"
        ordering = ("-date_joined",)

    def __str__(self) -> str:
        return self.email

    @property
    def display_name(self) -> str:
        return self.name or self.email.split("@")[0]

    def revoke_all_tokens(self) -> int:
        """Invalidate every JWT previously issued to this user.

        Two mechanisms, because access and refresh tokens fail differently:

        1. `token_version` is bumped, so any *access* token still in flight
           carries a stale version and is rejected on its next request.
        2. Every outstanding *refresh* token is blacklisted, so it can no
           longer be exchanged for a fresh access token.

        Returns the number of refresh tokens blacklisted.
        """
        from rest_framework_simplejwt.token_blacklist.models import (
            BlacklistedToken,
            OutstandingToken,
        )

        self.token_version += 1
        self.save(update_fields=["token_version"])

        blacklisted = 0
        for token in OutstandingToken.objects.filter(user=self):
            _, created = BlacklistedToken.objects.get_or_create(token=token)
            blacklisted += int(created)
        return blacklisted
