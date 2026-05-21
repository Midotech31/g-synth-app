# -*- coding: utf-8 -*-
"""
modules/my_projects.py

Per-user saved-projects workspace.

A project is a lightweight, persistent container — a name, an optional
DNA/protein sequence, free-text notes, and a JSON `data` blob. Any
G-Synth module can write into it via `utils.database.create_project`
/ `update_project`, so this page doubles as the hub where a user finds
everything they have saved.
"""

import logging

import streamlit as st

from utils import auth
from utils import database as db

logger = logging.getLogger("G-Synth:PROJECTS")


def main():
    user = auth.current_user()
    if not user:
        st.warning("Please sign in to use projects.")
        return

    st.title("📁 My Projects")
    st.caption(f"Signed in as **{user['email']}** — your projects are private to you.")

    if db.using_ephemeral_sqlite():
        st.info(
            "Demo mode: projects are stored in a temporary database and will "
            "reset when the app restarts. Connect a persistent database "
            "(`DATABASE_URL` secret) to keep them.",
            icon="ℹ️",
        )

    _new_project_form(user["id"])
    st.divider()
    _project_list(user["id"])


def _new_project_form(user_id: int):
    with st.expander("➕  New project", expanded=False):
        with st.form("new_project_form", clear_on_submit=True):
            name = st.text_input("Project name", placeholder="e.g. Insulin construct v1")
            sequence = st.text_area(
                "Sequence (optional)", height=110,
                placeholder="Paste a DNA or protein sequence to keep with this project…",
            )
            notes = st.text_area("Notes (optional)", height=80,
                                 placeholder="Anything you want to remember…")
            submitted = st.form_submit_button("Create project", type="primary")
        if submitted:
            if not name.strip():
                st.error("Please give the project a name.")
            else:
                try:
                    pid = db.create_project(
                        user_id, name,
                        sequence="".join(sequence.split()),
                        notes=notes.strip(),
                    )
                    st.success(f"Project created (#{pid}).")
                    st.rerun()
                except Exception as e:  # noqa: BLE001
                    st.error(f"Could not create project: {e}")


def _project_list(user_id: int):
    try:
        projects = db.list_projects(user_id)
    except Exception as e:  # noqa: BLE001
        st.error(f"Could not load projects: {e}")
        return

    if not projects:
        st.markdown(
            "<div style='text-align:center; padding:2.5rem; color:#64748b;'>"
            "<div style='font-size:2.5rem;'>📭</div>"
            "<p>No projects yet. Create your first one above.</p></div>",
            unsafe_allow_html=True,
        )
        return

    st.markdown(f"**{len(projects)} project{'s' if len(projects) != 1 else ''}**")

    for p in projects:
        with st.container(border=True):
            head, actions = st.columns([4, 1])
            with head:
                st.markdown(f"### {p['name']}")
                meta = f"module: `{p['module']}`"
                if p["updated_at"]:
                    meta += f" · updated {p['updated_at']:%Y-%m-%d %H:%M}"
                st.caption(meta)
            with actions:
                if st.button("🗑 Delete", key=f"del_{p['id']}",
                             use_container_width=True):
                    try:
                        db.delete_project(user_id, p["id"])
                        st.rerun()
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Delete failed: {e}")

            with st.expander("View / edit"):
                _edit_project(user_id, p)


def _edit_project(user_id: int, project: dict):
    pid = project["id"]
    with st.form(f"edit_{pid}", clear_on_submit=False):
        name = st.text_input("Name", value=project["name"], key=f"name_{pid}")
        sequence = st.text_area("Sequence", value=project["sequence"],
                                height=110, key=f"seq_{pid}")
        notes = st.text_area("Notes", value=project["notes"],
                             height=80, key=f"notes_{pid}")
        saved = st.form_submit_button("💾 Save changes", type="primary")
    if saved:
        try:
            db.update_project(
                user_id, pid,
                name=name, sequence="".join(sequence.split()), notes=notes,
            )
            st.success("Saved.")
            st.rerun()
        except Exception as e:  # noqa: BLE001
            st.error(f"Save failed: {e}")

    if project["data"]:
        st.markdown("**Stored data**")
        st.json(project["data"], expanded=False)


def app():
    """Alternate entry point for modular integration."""
    main()


if __name__ == "__main__":
    main()
