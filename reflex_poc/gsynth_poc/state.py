"""Application state: CRISPR designer + project list."""
from __future__ import annotations
import json
from datetime import datetime

import reflex as rx
from sqlmodel import select

from gsynth_poc.auth import AuthState
from gsynth_poc.crispr_logic import Guide, clean_dna, design_guides
from gsynth_poc.models import Project


SAMPLE_SEQ = (
    "ATGAAACGCCTGGCTGTTTTTGTGCTGCTGTTCCTGGCTTGGCTGGGCGCTGCTCAGCGCAGCGCA"
    "GAATTCGCGCGCCGCTGACTAGGGCCAGGCCAAGGTCACCTCCAATGACTAGGGTGGTTATGCGAA"
    "TCCATGGAGTACACCACGGAGTGCTGAGCATGGAGCTGTAGTCCAGGATGCGTAGGAGTAGATGCG"
)


class ProjectState(rx.State):
    """CRISPR analysis + saved projects."""

    # Input form
    project_name: str = "Untitled design"
    input_sequence: str = SAMPLE_SEQ
    min_gc: int = 30
    max_gc: int = 70
    min_on_target: float = 0.4

    # Computation results (transient — not persisted until user clicks Save)
    last_guides: list[dict] = []
    last_message: str = ""
    last_error: str = ""

    # Saved projects loaded from DB for the dashboard
    saved_projects: list[dict] = []

    # Project currently opened
    opened_id: int = 0
    opened_name: str = ""
    opened_sequence: str = ""
    opened_guides: list[dict] = []
    opened_created_at: str = ""

    # Explicit setters (Reflex 0.9 no longer auto-generates them)
    @rx.event
    def set_project_name(self, v: str):
        self.project_name = v

    @rx.event
    def set_input_sequence(self, v: str):
        self.input_sequence = v

    @rx.event
    def set_min_gc(self, v: str):
        try:
            self.min_gc = max(0, min(100, int(float(v))))
        except (ValueError, TypeError):
            pass

    @rx.event
    def set_max_gc(self, v: str):
        try:
            self.max_gc = max(0, min(100, int(float(v))))
        except (ValueError, TypeError):
            pass

    @rx.event
    def set_min_on_target(self, v: str):
        try:
            self.min_on_target = max(0.0, min(1.0, float(v)))
        except (ValueError, TypeError):
            pass

    @rx.event
    def run_design(self):
        self.last_error = ""
        self.last_message = ""
        cleaned = clean_dna(self.input_sequence)
        if len(cleaned) < 23:
            self.last_error = "Sequence must contain at least 23 valid DNA bases."
            self.last_guides = []
            return
        try:
            guides = design_guides(
                cleaned,
                min_gc=float(self.min_gc),
                max_gc=float(self.max_gc),
                min_on_target=self.min_on_target,
            )
        except Exception as e:
            self.last_error = f"Run failed: {e}"
            self.last_guides = []
            return
        self.last_guides = [g.to_dict() for g in guides]
        self.last_message = (
            f"Found {len(guides)} guide{'s' if len(guides) != 1 else ''}."
            if guides else
            "No guides matched these filters — try relaxing GC bounds or min on-target."
        )

    @rx.event
    async def save_project(self):
        auth = await self.get_state(AuthState)
        if not auth.is_authenticated:
            return rx.redirect("/login")
        if not self.last_guides:
            self.last_error = "Run a design before saving — there's nothing to save yet."
            return
        with rx.session() as session:
            project = Project(
                user_id=auth.user_id,
                name=self.project_name or "Untitled design",
                module="crispr",
                input_sequence=clean_dna(self.input_sequence),
                result_json=json.dumps({"guides": self.last_guides}),
            )
            session.add(project)
            session.commit()
            session.refresh(project)
        self.last_message = f"Saved as project #{project.id}."
        return rx.redirect("/dashboard")

    @rx.event
    async def load_projects(self):
        auth = await self.get_state(AuthState)
        if not auth.is_authenticated:
            self.saved_projects = []
            return
        with rx.session() as session:
            rows = session.exec(
                select(Project).where(Project.user_id == auth.user_id).order_by(Project.created_at.desc())
            ).all()
        self.saved_projects = [
            {
                "id": p.id,
                "name": p.name,
                "module": p.module,
                "created_at": _fmt(p.created_at),
                "guide_count": _count_guides(p.result_json),
            }
            for p in rows
        ]

    @rx.event
    async def open_project(self, project_id: int):
        auth = await self.get_state(AuthState)
        if not auth.is_authenticated:
            return rx.redirect("/login")
        with rx.session() as session:
            project = session.exec(
                select(Project).where(Project.id == project_id, Project.user_id == auth.user_id)
            ).first()
        if project is None:
            return rx.redirect("/dashboard")
        payload = json.loads(project.result_json or "{}")
        self.opened_id = project.id
        self.opened_name = project.name
        self.opened_sequence = project.input_sequence
        self.opened_guides = payload.get("guides", [])
        self.opened_created_at = _fmt(project.created_at)
        return rx.redirect(f"/projects/{project.id}")

    @rx.event
    async def delete_project(self, project_id: int):
        auth = await self.get_state(AuthState)
        if not auth.is_authenticated:
            return rx.redirect("/login")
        with rx.session() as session:
            project = session.exec(
                select(Project).where(Project.id == project_id, Project.user_id == auth.user_id)
            ).first()
            if project is not None:
                session.delete(project)
                session.commit()
        return ProjectState.load_projects


def _fmt(dt: datetime | None) -> str:
    if not dt:
        return ""
    return dt.strftime("%Y-%m-%d %H:%M")


def _count_guides(result_json: str) -> int:
    try:
        return len(json.loads(result_json or "{}").get("guides", []))
    except Exception:
        return 0
