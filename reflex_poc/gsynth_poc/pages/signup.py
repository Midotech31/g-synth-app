"""Sign-up page."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page
from gsynth_poc.pages.login import _auth_card


def signup_form() -> rx.Component:
    return _auth_card(
        rx.vstack(
            rx.heading("Create your account", size="6"),
            rx.text("Free, no credit card. You stay in control of your data.", color="gray"),
            rx.input(
                placeholder="email@example.com",
                value=AuthState.email_input,
                on_change=AuthState.set_email_input,
                type="email",
                width="100%",
            ),
            rx.input(
                placeholder="Password (8+ characters)",
                value=AuthState.password_input,
                on_change=AuthState.set_password_input,
                type="password",
                width="100%",
            ),
            rx.input(
                placeholder="Confirm password",
                value=AuthState.password_confirm,
                on_change=AuthState.set_password_confirm,
                type="password",
                width="100%",
            ),
            rx.cond(
                AuthState.error != "",
                rx.callout(AuthState.error, icon="triangle_alert", color_scheme="red", width="100%"),
            ),
            rx.button(
                "Create account",
                on_click=AuthState.sign_up,
                color_scheme="cyan",
                size="3",
                width="100%",
            ),
            rx.hstack(
                rx.text("Already have an account?", color="gray", size="2"),
                rx.link("Sign in", href="/login", color="#1e3a8a"),
                spacing="2",
                justify="center",
                width="100%",
            ),
            spacing="3",
            width="100%",
        ),
    )


def signup_page() -> rx.Component:
    return page("Create account", signup_form())
