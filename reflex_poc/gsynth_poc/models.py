"""SQLModel tables — persistent storage layer.

`reflex>=0.9.2` deprecated `rx.Model` and recommends using SQLModel
directly with an explicit engine; we follow that recommendation. The
DB URL still comes from `rxconfig.py` so `rx.session()` keeps working.
"""
from __future__ import annotations
from datetime import datetime, timezone
from typing import Optional

import sqlmodel


class User(sqlmodel.SQLModel, table=True):
    """Registered user. Password stored as a salted PBKDF2-SHA256 hash."""
    id: Optional[int] = sqlmodel.Field(default=None, primary_key=True)
    email: str = sqlmodel.Field(unique=True, index=True)
    password_hash: str
    created_at: datetime = sqlmodel.Field(default_factory=lambda: datetime.now(timezone.utc))


class Project(sqlmodel.SQLModel, table=True):
    """A saved analysis. Currently scoped to CRISPR guide-design runs."""
    id: Optional[int] = sqlmodel.Field(default=None, primary_key=True)
    user_id: int = sqlmodel.Field(foreign_key="user.id", index=True)
    name: str
    module: str                    # "crispr", "primers", … (extensible)
    input_sequence: str
    result_json: str               # serialized result payload
    notes: Optional[str] = None
    created_at: datetime = sqlmodel.Field(default_factory=lambda: datetime.now(timezone.utc))
