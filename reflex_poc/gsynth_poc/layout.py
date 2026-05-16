"""Shared layout (navbar + page wrapper)."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState


_BRAND_COLOR = "#1e3a8a"


def navbar() -> rx.Component:
    return rx.hstack(
        rx.hstack(
            rx.icon("dna", size=22, color="white"),
            rx.heading("G-Synth", size="5", color="white"),
            rx.badge("Reflex POC", color_scheme="cyan", variant="solid"),
            spacing="3",
            align="center",
        ),
        rx.spacer(),
        rx.hstack(
            rx.cond(
                AuthState.is_authenticated,
                rx.hstack(
                    rx.link(
                        rx.button("Dashboard", variant="ghost", color="white"),
                        href="/dashboard",
                    ),
                    rx.link(
                        rx.button("CRISPR Designer", variant="ghost", color="white"),
                        href="/crispr",
                    ),
                    rx.badge(AuthState.user_email, color_scheme="gray", variant="soft"),
                    rx.button(
                        "Sign out",
                        on_click=AuthState.sign_out,
                        variant="outline",
                        color="white",
                        border_color="white",
                    ),
                    spacing="3",
                    align="center",
                ),
                rx.hstack(
                    rx.link(rx.button("Sign in", variant="outline", color="white", border_color="white"), href="/login"),
                    rx.link(rx.button("Sign up", variant="solid", color_scheme="cyan"), href="/signup"),
                    spacing="2",
                ),
            ),
        ),
        width="100%",
        padding="0.75em 2em",
        background=_BRAND_COLOR,
        align="center",
    )


def page(title: str, *children: rx.Component) -> rx.Component:
    return rx.vstack(
        navbar(),
        rx.box(
            rx.vstack(
                rx.heading(title, size="7", color=_BRAND_COLOR),
                *children,
                spacing="4",
                align_items="stretch",
                width="100%",
            ),
            max_width="1100px",
            margin="2em auto",
            padding="0 2em",
            width="100%",
        ),
        spacing="0",
        width="100%",
        min_height="100vh",
        background="#f8fafc",
    )
