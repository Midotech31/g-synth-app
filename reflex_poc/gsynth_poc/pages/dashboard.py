"""User dashboard — saved-projects list."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page
from gsynth_poc.state import ProjectState


def _empty_state() -> rx.Component:
    return rx.center(
        rx.vstack(
            rx.icon("inbox", size=48, color="gray"),
            rx.text("No saved projects yet.", color="gray"),
            rx.link(
                rx.button("Design your first guides", color_scheme="cyan", size="3"),
                href="/crispr",
            ),
            spacing="3",
            align="center",
        ),
        padding="3em",
        width="100%",
        background="white",
        border_radius="12px",
        border="1px solid #e2e8f0",
    )


def _project_row(project: dict) -> rx.Component:
    return rx.hstack(
        rx.vstack(
            rx.heading(project["name"], size="4"),
            rx.hstack(
                rx.badge(project["module"], color_scheme="cyan", variant="soft"),
                rx.text(f"{project['guide_count']} guides", size="2", color="gray"),
                rx.text("·", size="2", color="gray"),
                rx.text(project["created_at"], size="2", color="gray"),
                spacing="2",
            ),
            align_items="start",
            spacing="1",
        ),
        rx.spacer(),
        rx.hstack(
            rx.button(
                "Open",
                on_click=ProjectState.open_project(project["id"]),
                color_scheme="cyan",
                variant="solid",
            ),
            rx.button(
                rx.icon("trash-2", size=16),
                on_click=ProjectState.delete_project(project["id"]),
                color_scheme="red",
                variant="ghost",
            ),
            spacing="2",
        ),
        padding="1em 1.25em",
        background="white",
        border_radius="10px",
        border="1px solid #e2e8f0",
        align="center",
        width="100%",
    )


def dashboard_body() -> rx.Component:
    return rx.vstack(
        rx.hstack(
            rx.text("Your saved CRISPR designs.", color="gray", size="3"),
            rx.spacer(),
            rx.link(
                rx.button("New design", color_scheme="cyan", size="2"),
                href="/crispr",
            ),
            width="100%",
            align="center",
        ),
        rx.cond(
            ProjectState.saved_projects.length() > 0,
            rx.vstack(
                rx.foreach(ProjectState.saved_projects, _project_row),
                spacing="3",
                width="100%",
            ),
            _empty_state(),
        ),
        spacing="4",
        width="100%",
    )


def dashboard_page() -> rx.Component:
    return page(f"Dashboard", dashboard_body())


dashboard_page.on_load = [AuthState.require_auth, ProjectState.load_projects]
