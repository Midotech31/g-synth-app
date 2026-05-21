# -*- coding: utf-8 -*-
"""
utils/auth.py

Authentication for G-Synth: PBKDF2-SHA256 password hashing (Python
standard library — no extra dependency) plus Streamlit-session-based
login state.

The signed-in user is tracked in `st.session_state`:
    auth_user_id, auth_user_email, auth_user_name
"""

import hashlib
import hmac
import logging
import os
import re

import streamlit as st

from utils import database as db

logger = logging.getLogger("G-Synth:AUTH")

_PBKDF2_ITERATIONS = 240_000
_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
_MIN_PASSWORD_LEN = 8


# ───────────────────────────────────────────────────────────────────────────────
# PASSWORD HASHING
# ───────────────────────────────────────────────────────────────────────────────
def hash_password(password: str, *, salt: bytes | None = None) -> str:
    """Return 'pbkdf2_sha256$<iterations>$<salt_hex>$<hash_hex>'."""
    if salt is None:
        salt = os.urandom(16)
    dk = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt,
                             _PBKDF2_ITERATIONS)
    return f"pbkdf2_sha256${_PBKDF2_ITERATIONS}${salt.hex()}${dk.hex()}"


def verify_password(password: str, stored: str) -> bool:
    """Constant-time verification against a stored hash."""
    try:
        scheme, iterations, salt_hex, hash_hex = stored.split("$")
        if scheme != "pbkdf2_sha256":
            return False
        dk = hashlib.pbkdf2_hmac(
            "sha256", password.encode("utf-8"),
            bytes.fromhex(salt_hex), int(iterations),
        )
        return hmac.compare_digest(dk, bytes.fromhex(hash_hex))
    except (ValueError, AttributeError):
        return False


# ───────────────────────────────────────────────────────────────────────────────
# VALIDATION
# ───────────────────────────────────────────────────────────────────────────────
def validate_email(email: str) -> bool:
    return bool(_EMAIL_RE.match(email.strip()))


def validate_password(password: str) -> tuple[bool, str]:
    if len(password) < _MIN_PASSWORD_LEN:
        return False, f"Password must be at least {_MIN_PASSWORD_LEN} characters."
    if password.isdigit() or password.isalpha():
        return False, "Password must mix letters and numbers."
    return True, ""


# ───────────────────────────────────────────────────────────────────────────────
# SIGN-UP / SIGN-IN / SIGN-OUT
# ───────────────────────────────────────────────────────────────────────────────
def sign_up(email: str, name: str, password: str, password_confirm: str) -> tuple[bool, str]:
    """Create a new account. Returns (ok, message)."""
    email = email.strip().lower()
    if not validate_email(email):
        return False, "Please enter a valid email address."
    if not name.strip():
        return False, "Please enter your name."
    ok, msg = validate_password(password)
    if not ok:
        return False, msg
    if password != password_confirm:
        return False, "The two passwords do not match."
    try:
        if db.get_user_by_email(email) is not None:
            return False, "An account with this email already exists."
        user = db.create_user(email, name, hash_password(password))
    except Exception as e:  # noqa: BLE001 — surface DB errors to the UI
        logger.error("Sign-up failed: %s", e)
        return False, f"Could not create the account: {e}"
    _set_session(user.id, user.email, user.name)
    return True, "Account created."


def sign_in(email: str, password: str) -> tuple[bool, str]:
    """Authenticate an existing account. Returns (ok, message)."""
    email = email.strip().lower()
    if not email or not password:
        return False, "Email and password are required."
    try:
        user = db.get_user_by_email(email)
    except Exception as e:  # noqa: BLE001
        logger.error("Sign-in DB error: %s", e)
        return False, f"Database error: {e}"
    if user is None or not verify_password(password, user.password_hash):
        return False, "Invalid email or password."
    _set_session(user.id, user.email, user.name)
    return True, "Signed in."


def sign_out() -> None:
    for key in ("auth_user_id", "auth_user_email", "auth_user_name"):
        st.session_state.pop(key, None)


# ───────────────────────────────────────────────────────────────────────────────
# SESSION HELPERS
# ───────────────────────────────────────────────────────────────────────────────
def _set_session(user_id: int, email: str, name: str) -> None:
    st.session_state["auth_user_id"] = user_id
    st.session_state["auth_user_email"] = email
    st.session_state["auth_user_name"] = name


def is_authenticated() -> bool:
    return bool(st.session_state.get("auth_user_id"))


def current_user() -> dict | None:
    if not is_authenticated():
        return None
    return {
        "id": st.session_state.get("auth_user_id"),
        "email": st.session_state.get("auth_user_email", ""),
        "name": st.session_state.get("auth_user_name", ""),
    }


def current_user_id() -> int | None:
    return st.session_state.get("auth_user_id")
