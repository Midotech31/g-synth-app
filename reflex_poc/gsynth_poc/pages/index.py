"""Landing page (/) — hero + CTA."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page


def _hero() -> rx.Component:
    return rx.center(
        rx.vstack(
            rx.icon("dna", size=64, color="#1e3a8a"),
            rx.heading("G-Synth — multi-user proof of concept", size="8", text_align="center"),
            rx.text(
                "Reflex + SQLite, demonstrating the architecture that replaces "
                "Streamlit: real accounts, persistent projects, and an "
                "interactive UI that compiles to React under the hood.",
                size="4",
                color="gray",
                text_align="center",
                max_width="700px",
            ),
            rx.hstack(
                rx.cond(
                    AuthState.is_authenticated,
                    rx.link(rx.button("Go to dashboard", color_scheme="cyan", size="3"), href="/dashboard"),
                    rx.hstack(
                        rx.link(rx.button("Get started — sign up", color_scheme="cyan", size="3"), href="/signup"),
                        rx.link(rx.button("Sign in", size="3", variant="outline"), href="/login"),
                        spacing="3",
                    ),
                ),
                margin_top="1em",
            ),
            spacing="4",
            align="center",
        ),
        padding="4em 2em",
        width="100%",
    )


def _feature(icon: str, title: str, desc: str) -> rx.Component:
    return rx.box(
        rx.vstack(
            rx.icon(icon, size=32, color="#06b6d4"),
            rx.heading(title, size="4"),
            rx.text(desc, color="gray", size="2"),
            spacing="2",
            align_items="start",
        ),
        padding="1.5em",
        background="white",
        border_radius="12px",
        border="1px solid #e2e8f0",
        flex="1",
    )


def _features() -> rx.Component:
    return rx.hstack(
        _feature("user-check", "Per-user accounts",
                 "Email/password sign-in. Each user sees only their own data."),
        _feature("database", "Persistent projects",
                 "Save your CRISPR runs to the database and reopen them later."),
        _feature("zap", "Interactive UI",
                 "Reflex compiles to React + FastAPI — no page reloads."),
        spacing="4",
        width="100%",
    )


def index_page() -> rx.Component:
    return page("", _hero(), _features())
