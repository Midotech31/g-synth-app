"""CRISPR sgRNA designer page."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page
from gsynth_poc.state import ProjectState


def _input_card() -> rx.Component:
    return rx.box(
        rx.vstack(
            rx.heading("Input", size="5"),
            rx.input(
                placeholder="Project name",
                value=ProjectState.project_name,
                on_change=ProjectState.set_project_name,
                width="100%",
            ),
            rx.text_area(
                placeholder="Paste DNA sequence (A/C/G/T only — IUPAC will be stripped)…",
                value=ProjectState.input_sequence,
                on_change=ProjectState.set_input_sequence,
                rows="8",
                width="100%",
                font_family="ui-monospace, monospace",
            ),
            rx.hstack(
                rx.vstack(
                    rx.text("Min GC %", size="2", color="gray"),
                    rx.input(
                        type="number",
                        value=ProjectState.min_gc.to_string(),
                        on_change=ProjectState.set_min_gc,
                        width="100%",
                    ),
                    align_items="start",
                    width="100%",
                ),
                rx.vstack(
                    rx.text("Max GC %", size="2", color="gray"),
                    rx.input(
                        type="number",
                        value=ProjectState.max_gc.to_string(),
                        on_change=ProjectState.set_max_gc,
                        width="100%",
                    ),
                    align_items="start",
                    width="100%",
                ),
                rx.vstack(
                    rx.text("Min on-target", size="2", color="gray"),
                    rx.input(
                        type="number",
                        value=ProjectState.min_on_target.to_string(),
                        on_change=ProjectState.set_min_on_target,
                        step="0.05",
                        width="100%",
                    ),
                    align_items="start",
                    width="100%",
                ),
                spacing="3",
                width="100%",
            ),
            rx.hstack(
                rx.button(
                    rx.icon("play", size=18),
                    "Find guides",
                    on_click=ProjectState.run_design,
                    color_scheme="cyan",
                    size="3",
                ),
                rx.button(
                    rx.icon("save", size=18),
                    "Save to dashboard",
                    on_click=ProjectState.save_project,
                    color_scheme="green",
                    size="3",
                    variant="soft",
                ),
                spacing="3",
            ),
            rx.cond(
                ProjectState.last_error != "",
                rx.callout(ProjectState.last_error, icon="triangle_alert", color_scheme="red", width="100%"),
            ),
            rx.cond(
                ProjectState.last_message != "",
                rx.callout(ProjectState.last_message, icon="info", color_scheme="blue", width="100%"),
            ),
            spacing="3",
            width="100%",
        ),
        padding="1.5em",
        background="white",
        border_radius="12px",
        border="1px solid #e2e8f0",
        width="100%",
    )


def _guide_row(guide: dict) -> rx.Component:
    # Badge colors are precomputed server-side (see crispr_logic.Guide.to_dict).
    return rx.table.row(
        rx.table.cell(rx.code(guide["protospacer"])),
        rx.table.cell(rx.code(guide["pam"])),
        rx.table.cell(guide["start"].to_string()),
        rx.table.cell(rx.badge(guide["strand"], color_scheme=guide["strand_color"])),
        rx.table.cell(guide["gc_percent"].to_string()),
        rx.table.cell(rx.badge(guide["on_target"].to_string(), color_scheme=guide["score_color"])),
    )


def _results_card() -> rx.Component:
    return rx.box(
        rx.vstack(
            rx.heading("Results", size="5"),
            rx.cond(
                ProjectState.last_guides.length() > 0,
                rx.table.root(
                    rx.table.header(
                        rx.table.row(
                            rx.table.column_header_cell("Protospacer (20 nt)"),
                            rx.table.column_header_cell("PAM"),
                            rx.table.column_header_cell("Start"),
                            rx.table.column_header_cell("Strand"),
                            rx.table.column_header_cell("GC %"),
                            rx.table.column_header_cell("On-target"),
                        )
                    ),
                    rx.table.body(
                        rx.foreach(ProjectState.last_guides, _guide_row),
                    ),
                    variant="surface",
                    width="100%",
                ),
                rx.text("Run a design to see candidate guides here.", color="gray"),
            ),
            spacing="3",
            width="100%",
        ),
        padding="1.5em",
        background="white",
        border_radius="12px",
        border="1px solid #e2e8f0",
        width="100%",
    )


def crispr_body() -> rx.Component:
    return rx.vstack(
        _input_card(),
        _results_card(),
        spacing="4",
        width="100%",
    )


def crispr_page() -> rx.Component:
    return page("CRISPR sgRNA Designer", crispr_body())


crispr_page.on_load = [AuthState.require_auth]
