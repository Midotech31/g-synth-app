"""Read-only project detail page (loaded via /projects/<id>)."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.auth import AuthState
from gsynth_poc.layout import page
from gsynth_poc.state import ProjectState


def _guide_row(guide: dict) -> rx.Component:
    return rx.table.row(
        rx.table.cell(rx.code(guide["protospacer"])),
        rx.table.cell(rx.code(guide["pam"])),
        rx.table.cell(guide["start"].to_string()),
        rx.table.cell(rx.badge(guide["strand"], color_scheme=guide["strand_color"])),
        rx.table.cell(guide["gc_percent"].to_string()),
        rx.table.cell(rx.badge(guide["on_target"].to_string(), color_scheme=guide["score_color"])),
    )


def project_body() -> rx.Component:
    return rx.vstack(
        rx.hstack(
            rx.vstack(
                rx.heading(ProjectState.opened_name, size="6"),
                rx.text(f"Saved on {ProjectState.opened_created_at}", color="gray", size="2"),
                align_items="start",
                spacing="1",
            ),
            rx.spacer(),
            rx.link(rx.button("← Back to dashboard", variant="ghost"), href="/dashboard"),
            width="100%",
            align="center",
        ),
        rx.box(
            rx.vstack(
                rx.heading("Input sequence", size="4"),
                rx.code(ProjectState.opened_sequence, width="100%", style={"white_space": "pre_wrap"}),
                spacing="2",
                width="100%",
            ),
            padding="1.5em",
            background="white",
            border_radius="12px",
            border="1px solid #e2e8f0",
            width="100%",
        ),
        rx.box(
            rx.vstack(
                rx.heading(f"Guides", size="4"),
                rx.cond(
                    ProjectState.opened_guides.length() > 0,
                    rx.table.root(
                        rx.table.header(
                            rx.table.row(
                                rx.table.column_header_cell("Protospacer"),
                                rx.table.column_header_cell("PAM"),
                                rx.table.column_header_cell("Start"),
                                rx.table.column_header_cell("Strand"),
                                rx.table.column_header_cell("GC %"),
                                rx.table.column_header_cell("On-target"),
                            )
                        ),
                        rx.table.body(rx.foreach(ProjectState.opened_guides, _guide_row)),
                        variant="surface",
                        width="100%",
                    ),
                    rx.text("This project has no saved guides.", color="gray"),
                ),
                spacing="2",
                width="100%",
            ),
            padding="1.5em",
            background="white",
            border_radius="12px",
            border="1px solid #e2e8f0",
            width="100%",
        ),
        spacing="4",
        width="100%",
    )


def project_detail_page() -> rx.Component:
    return page("Project", project_body())


project_detail_page.on_load = [AuthState.require_auth]
