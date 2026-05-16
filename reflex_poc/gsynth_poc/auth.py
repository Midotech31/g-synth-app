"""Authentication state.

Local SQLite + PBKDF2-SHA256 password hashing. When you wire up Supabase,
swap _hash_password/_verify and the DB calls for the Supabase SDK — the
rest of the app doesn't need to change.
"""
from __future__ import annotations
import hashlib
import hmac
import os
import re

import reflex as rx
from sqlmodel import select

from gsynth_poc.models import User


_ITER = 200_000


def _hash_password(password: str, *, salt: bytes | None = None) -> str:
    """PBKDF2-SHA256 → 'pbkdf2_sha256$<iter>$<salt_hex>$<hash_hex>'."""
    if salt is None:
        salt = os.urandom(16)
    dk = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, _ITER)
    return f"pbkdf2_sha256${_ITER}${salt.hex()}${dk.hex()}"


def _verify_password(password: str, stored: str) -> bool:
    try:
        scheme, iterations, salt_hex, hash_hex = stored.split("$")
        if scheme != "pbkdf2_sha256":
            return False
        salt = bytes.fromhex(salt_hex)
        expected = bytes.fromhex(hash_hex)
        dk = hashlib.pbkdf2_hmac(
            "sha256", password.encode("utf-8"), salt, int(iterations)
        )
        return hmac.compare_digest(dk, expected)
    except (ValueError, AttributeError):
        return False


_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


class AuthState(rx.State):
    """Session-scoped authentication state."""

    email_input: str = ""
    password_input: str = ""
    password_confirm: str = ""
    error: str = ""
    info: str = ""

    # The signed-in user, persisted across reruns via Reflex client-state
    user_id: int = 0
    user_email: str = ""

    @rx.var
    def is_authenticated(self) -> bool:
        return self.user_id > 0

    def _reset_form(self) -> None:
        self.email_input = ""
        self.password_input = ""
        self.password_confirm = ""

    @rx.event
    def set_email_input(self, v: str):
        self.email_input = v

    @rx.event
    def set_password_input(self, v: str):
        self.password_input = v

    @rx.event
    def set_password_confirm(self, v: str):
        self.password_confirm = v

    @rx.event
    def sign_up(self):
        self.error = ""
        self.info = ""
        email = self.email_input.strip().lower()
        password = self.password_input
        if not _EMAIL_RE.match(email):
            self.error = "Please enter a valid email address."
            return
        if len(password) < 8:
            self.error = "Password must be at least 8 characters."
            return
        if password != self.password_confirm:
            self.error = "Passwords do not match."
            return
        with rx.session() as session:
            existing = session.exec(select(User).where(User.email == email)).first()
            if existing:
                self.error = "An account with that email already exists."
                return
            user = User(email=email, password_hash=_hash_password(password))
            session.add(user)
            session.commit()
            session.refresh(user)
            self.user_id = user.id
            self.user_email = user.email
        self._reset_form()
        return rx.redirect("/dashboard")

    @rx.event
    def sign_in(self):
        self.error = ""
        self.info = ""
        email = self.email_input.strip().lower()
        password = self.password_input
        if not email or not password:
            self.error = "Email and password are required."
            return
        with rx.session() as session:
            user = session.exec(select(User).where(User.email == email)).first()
            if user is None or not _verify_password(password, user.password_hash):
                self.error = "Invalid email or password."
                return
            self.user_id = user.id
            self.user_email = user.email
        self._reset_form()
        return rx.redirect("/dashboard")

    @rx.event
    def sign_out(self):
        self.user_id = 0
        self.user_email = ""
        self._reset_form()
        return rx.redirect("/login")

    @rx.event
    def require_auth(self):
        """Page-load guard for protected routes."""
        if not self.is_authenticated:
            return rx.redirect("/login")
