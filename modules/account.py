# -*- coding: utf-8 -*-
"""
modules/account.py

Login / sign-up gate and the account-management page for G-Synth.

`render_auth_gate()` is shown by app.py when no user is signed in.
`main()` is the in-app "Account" page (profile + change password).
"""

import logging

import streamlit as st

from utils import auth
from utils import database as db

logger = logging.getLogger("G-Synth:ACCOUNT")


# ───────────────────────────────────────────────────────────────────────────────
# AUTH GATE (shown before login)
# ───────────────────────────────────────────────────────────────────────────────
def render_auth_gate():
    """Full-page login / sign-up screen. Called by app.py when logged out."""
    col_left, col_mid, col_right = st.columns([1, 1.4, 1])

    with col_mid:
        try:
            st.image("assets/logo.png", use_container_width=True)
        except Exception:
            st.markdown("## 🧬 G-Synth")

        st.markdown(
            "<p style='text-align:center; color:#475569; font-size:1.05rem; "
            "margin-top:-0.5rem;'>Advanced Gene Synthesis Platform</p>",
            unsafe_allow_html=True,
        )

        if db.using_ephemeral_sqlite():
            st.info(
                "⚙️ **Demo mode** — running on a temporary local database. "
                "Accounts persist only until the app restarts. To make data "
                "permanent, set `DATABASE_URL` in the app secrets "
                "(see SETUP_MULTIUSER.md).",
                icon="ℹ️",
            )

        tab_login, tab_signup = st.tabs(["🔑  Sign in", "✨  Create account"])

        with tab_login:
            _login_form()

        with tab_signup:
            _signup_form()

        st.markdown(
            "<p style='text-align:center; color:#94a3b8; font-size:0.8rem; "
            "margin-top:2rem;'>© 2025 G-Synth — Dr. Mohamed Merzoug, ESSBO</p>",
            unsafe_allow_html=True,
        )


def _login_form():
    with st.form("login_form", clear_on_submit=False):
        email = st.text_input("Email", key="login_email",
                              placeholder="you@example.com")
        password = st.text_input("Password", type="password", key="login_pw",
                                 placeholder="Your password")
        submitted = st.form_submit_button("Sign in", use_container_width=True,
                                          type="primary")
    if submitted:
        ok, msg = auth.sign_in(email, password)
        if ok:
            st.success(msg)
            st.rerun()
        else:
            st.error(msg)


def _signup_form():
    with st.form("signup_form", clear_on_submit=False):
        name = st.text_input("Full name", key="su_name",
                             placeholder="Dr. Jane Doe")
        email = st.text_input("Email", key="su_email",
                              placeholder="you@example.com")
        password = st.text_input("Password", type="password", key="su_pw",
                                 placeholder="At least 8 characters, letters + numbers")
        password2 = st.text_input("Confirm password", type="password",
                                  key="su_pw2", placeholder="Repeat your password")
        submitted = st.form_submit_button("Create account",
                                          use_container_width=True, type="primary")
    if submitted:
        ok, msg = auth.sign_up(email, name, password, password2)
        if ok:
            st.success(msg)
            st.rerun()
        else:
            st.error(msg)


# ───────────────────────────────────────────────────────────────────────────────
# ACCOUNT PAGE (shown inside the app, after login)
# ───────────────────────────────────────────────────────────────────────────────
def main():
    """In-app account page — module entry point."""
    user = auth.current_user()
    if not user:
        st.warning("You are not signed in.")
        return

    st.title("⚡ My Account")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("#### Profile")
        st.markdown(f"**Name:** {user['name'] or '—'}")
        st.markdown(f"**Email:** {user['email']}")
    with c2:
        st.markdown("#### Storage")
        if db.using_ephemeral_sqlite():
            st.warning("Temporary database — data resets on app restart.")
        else:
            st.success("Persistent database connected.")
        try:
            n = len(db.list_projects(user["id"]))
            st.markdown(f"**Saved projects:** {n}")
        except Exception as e:  # noqa: BLE001
            st.error(f"Could not read projects: {e}")

    st.divider()
    st.markdown("#### Change password")
    with st.form("change_pw_form", clear_on_submit=True):
        current = st.text_input("Current password", type="password")
        new1 = st.text_input("New password", type="password")
        new2 = st.text_input("Confirm new password", type="password")
        submitted = st.form_submit_button("Update password", type="primary")
    if submitted:
        record = db.get_user_by_email(user["email"])
        if record is None or not auth.verify_password(current, record.password_hash):
            st.error("Current password is incorrect.")
        else:
            ok, msg = auth.validate_password(new1)
            if not ok:
                st.error(msg)
            elif new1 != new2:
                st.error("The two new passwords do not match.")
            else:
                try:
                    db.update_user_password(user["id"], auth.hash_password(new1))
                    st.success("Password updated.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Could not update password: {e}")


def app():
    """Alternate entry point for modular integration."""
    main()


if __name__ == "__main__":
    main()
