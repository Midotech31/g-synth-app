# -*- coding: utf-8 -*-
"""
utils/database.py

Database layer for G-Synth multi-user support.

Uses SQLAlchemy. The connection string comes from Streamlit secrets
(`DATABASE_URL`) when available, otherwise falls back to a local SQLite
file.

IMPORTANT — persistence on Streamlit Community Cloud:
    The container filesystem is EPHEMERAL. A SQLite file is wiped every
    time the app reboots (redeploy, inactivity recycle, ...). For real
    multi-user persistence you MUST set `DATABASE_URL` in the app's
    Streamlit secrets to an external Postgres database (e.g. a free
    Supabase or Neon instance). See README / SETUP_MULTIUSER.md.
"""

import json
import logging
from contextlib import contextmanager
from datetime import datetime, timezone

import streamlit as st
from sqlalchemy import (
    Column, DateTime, ForeignKey, Integer, String, Text, create_engine,
)
from sqlalchemy.orm import declarative_base, relationship, sessionmaker

logger = logging.getLogger("G-Synth:DB")

Base = declarative_base()


# ───────────────────────────────────────────────────────────────────────────────
# MODELS
# ───────────────────────────────────────────────────────────────────────────────
class User(Base):
    """A registered user. Passwords are stored only as PBKDF2-SHA256 hashes."""
    __tablename__ = "users"

    id = Column(Integer, primary_key=True)
    email = Column(String(255), unique=True, index=True, nullable=False)
    name = Column(String(255), default="")
    password_hash = Column(String(255), nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    projects = relationship("Project", back_populates="user",
                            cascade="all, delete-orphan")


class Project(Base):
    """A saved item belonging to one user.

    `data` is a JSON string — keeps the schema flexible so any G-Synth
    module can persist whatever shape of result it needs without a
    migration.
    """
    __tablename__ = "projects"

    id = Column(Integer, primary_key=True)
    user_id = Column(Integer, ForeignKey("users.id"), index=True, nullable=False)
    name = Column(String(255), nullable=False)
    module = Column(String(100), default="general")
    sequence = Column(Text, default="")
    notes = Column(Text, default="")
    data = Column(Text, default="{}")
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    user = relationship("User", back_populates="projects")


# ───────────────────────────────────────────────────────────────────────────────
# ENGINE / SESSION
# ───────────────────────────────────────────────────────────────────────────────
def _database_url() -> str:
    """Resolve the DB URL: Streamlit secret first, SQLite fallback."""
    try:
        url = st.secrets["DATABASE_URL"]
        if url:
            # Supabase / Heroku style "postgres://" → SQLAlchemy needs "postgresql://"
            if url.startswith("postgres://"):
                url = url.replace("postgres://", "postgresql://", 1)
            return url
    except Exception:
        pass
    return "sqlite:///gsynth.db"


@st.cache_resource(show_spinner=False)
def _get_engine():
    """Create (once) the SQLAlchemy engine and ensure tables exist."""
    url = _database_url()
    is_sqlite = url.startswith("sqlite")
    kwargs = {"pool_pre_ping": True}
    if is_sqlite:
        kwargs["connect_args"] = {"check_same_thread": False}
    engine = create_engine(url, **kwargs)
    Base.metadata.create_all(engine)
    logger.info("Database initialised (%s)", "sqlite" if is_sqlite else "postgres")
    return engine


def using_ephemeral_sqlite() -> bool:
    """True when running on the fallback SQLite — data will not persist."""
    return _database_url().startswith("sqlite")


@contextmanager
def get_session():
    """Context-managed SQLAlchemy session: commits on success, rolls back on error."""
    SessionLocal = sessionmaker(bind=_get_engine(), expire_on_commit=False)
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


# ───────────────────────────────────────────────────────────────────────────────
# USER CRUD
# ───────────────────────────────────────────────────────────────────────────────
def get_user_by_email(email: str):
    with get_session() as s:
        return s.query(User).filter(User.email == email.strip().lower()).first()


def create_user(email: str, name: str, password_hash: str) -> User:
    with get_session() as s:
        user = User(email=email.strip().lower(), name=name.strip(),
                    password_hash=password_hash)
        s.add(user)
        s.flush()
        s.refresh(user)
        return user


def update_user_password(user_id: int, password_hash: str) -> None:
    with get_session() as s:
        user = s.query(User).filter(User.id == user_id).first()
        if user:
            user.password_hash = password_hash


# ───────────────────────────────────────────────────────────────────────────────
# PROJECT CRUD
# ───────────────────────────────────────────────────────────────────────────────
def create_project(user_id: int, name: str, *, module: str = "general",
                    sequence: str = "", notes: str = "", data: dict | None = None) -> int:
    with get_session() as s:
        project = Project(
            user_id=user_id, name=name.strip() or "Untitled",
            module=module, sequence=sequence, notes=notes,
            data=json.dumps(data or {}),
        )
        s.add(project)
        s.flush()
        return project.id


def list_projects(user_id: int) -> list[dict]:
    """Return the user's projects as plain dicts (detached from the session)."""
    with get_session() as s:
        rows = (s.query(Project)
                .filter(Project.user_id == user_id)
                .order_by(Project.updated_at.desc())
                .all())
        return [_project_to_dict(p) for p in rows]


def get_project(user_id: int, project_id: int) -> dict | None:
    with get_session() as s:
        p = (s.query(Project)
             .filter(Project.id == project_id, Project.user_id == user_id)
             .first())
        return _project_to_dict(p) if p else None


def update_project(user_id: int, project_id: int, **fields) -> bool:
    allowed = {"name", "module", "sequence", "notes", "data"}
    with get_session() as s:
        p = (s.query(Project)
             .filter(Project.id == project_id, Project.user_id == user_id)
             .first())
        if not p:
            return False
        for key, value in fields.items():
            if key not in allowed:
                continue
            if key == "data" and isinstance(value, dict):
                value = json.dumps(value)
            setattr(p, key, value)
        p.updated_at = datetime.now(timezone.utc)
        return True


def delete_project(user_id: int, project_id: int) -> bool:
    with get_session() as s:
        p = (s.query(Project)
             .filter(Project.id == project_id, Project.user_id == user_id)
             .first())
        if not p:
            return False
        s.delete(p)
        return True


def _project_to_dict(p: Project) -> dict:
    try:
        data = json.loads(p.data or "{}")
    except (ValueError, TypeError):
        data = {}
    return {
        "id": p.id,
        "name": p.name,
        "module": p.module,
        "sequence": p.sequence or "",
        "notes": p.notes or "",
        "data": data,
        "created_at": p.created_at,
        "updated_at": p.updated_at,
    }
