"""Sign-in page."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page


def _auth_card(children: rx.Component) -> rx.Component:
    return rx.box(
        children,
        max_width="420px",
        margin="0 auto",
        padding="2em",
        background="white",
        border_radius="12px",
        box_shadow="0 4px 12px rgba(0,0,0,0.08)",
        border="1px solid #e2e8f0",
    )


def login_form() -> rx.Component:
    return _auth_card(
        rx.vstack(
            rx.heading("Welcome back", size="6"),
            rx.text("Sign in to access your saved projects.", color="gray"),
            rx.input(
                placeholder="email@example.com",
                value=AuthState.email_input,
                on_change=AuthState.set_email_input,
                type="email",
                width="100%",
            ),
            rx.input(
                placeholder="Password",
                value=AuthState.password_input,
                on_change=AuthState.set_password_input,
                type="password",
                width="100%",
            ),
            rx.cond(
                AuthState.error != "",
                rx.callout(AuthState.error, icon="triangle_alert", color_scheme="red", width="100%"),
            ),
            rx.button(
                "Sign in",
                on_click=AuthState.sign_in,
                color_scheme="cyan",
                size="3",
                width="100%",
            ),
            rx.hstack(
                rx.text("No account?", color="gray", size="2"),
                rx.link("Sign up", href="/signup", color="#1e3a8a"),
                spacing="2",
                justify="center",
                width="100%",
            ),
            spacing="3",
            width="100%",
        ),
    )


def login_page() -> rx.Component:
    return page("Sign in", login_form())
