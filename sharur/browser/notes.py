"""Human notes and flags on browser entities, kept beside (never inside) the dataset.

Flags record curator judgments (verified, suspicious, interesting, follow up);
notes are free text. Both are attributed and timestamped, and deletions are
soft so the history stays auditable. These are human annotations, kept apart
from the evidence-backed labels.
"""

from __future__ import annotations

import os
import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

FLAGS = {
    "verified": ("Verified", "✓"),
    "suspicious": ("Suspicious", "!"),
    "interesting": ("Interesting", "★"),
    "follow_up": ("Follow up", "↻"),
}
KINDS = ("protein", "genome", "system", "domain", "vog", "function", "contig", "clade", "module", "crispr")

SCHEMA = """
CREATE TABLE IF NOT EXISTS notes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    kind TEXT NOT NULL,
    entity TEXT NOT NULL,
    flag TEXT,
    text TEXT,
    author TEXT NOT NULL,
    created_at TEXT NOT NULL,
    deleted_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_notes_entity ON notes(kind, entity);
"""


def default_notes_path(db_path: Path) -> Path:
    """``browser_notes.sqlite`` beside the dataset, or under ~/.sharur/notes when that is read-only."""
    beside = db_path.resolve().parent / "browser_notes.sqlite"
    if os.access(beside.parent, os.W_OK):
        return beside
    fallback = Path.home() / ".sharur" / "notes" / f"{db_path.resolve().parent.name}.sqlite"
    fallback.parent.mkdir(parents=True, exist_ok=True)
    return fallback


class NotesStore:
    def __init__(self, path: Path):
        self.path = Path(path)
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(str(self.path), check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._conn.executescript(SCHEMA)
        self._conn.commit()

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat(timespec="seconds")

    def state(self, kind: str, entity: str) -> dict[str, Any]:
        with self._lock:
            rows = self._conn.execute(
                "SELECT * FROM notes WHERE kind = ? AND entity = ? AND deleted_at IS NULL ORDER BY id",
                (kind, entity)).fetchall()
        flags: dict[str, list[str]] = {}
        notes = []
        for r in rows:
            if r["flag"]:
                flags.setdefault(r["flag"], []).append(r["author"])
            if r["text"]:
                notes.append({"id": r["id"], "text": r["text"], "author": r["author"], "created_at": r["created_at"]})
        return {"kind": kind, "entity": entity, "flags": flags, "notes": notes}

    def toggle_flag(self, kind: str, entity: str, flag: str, author: str) -> dict[str, Any]:
        if flag not in FLAGS or kind not in KINDS:
            raise ValueError("unknown flag or kind")
        with self._lock:
            existing = self._conn.execute(
                "SELECT id FROM notes WHERE kind = ? AND entity = ? AND flag = ? AND author = ? AND deleted_at IS NULL",
                (kind, entity, flag, author)).fetchone()
            if existing:
                self._conn.execute("UPDATE notes SET deleted_at = ? WHERE id = ?", (self._now(), existing["id"]))
            else:
                self._conn.execute(
                    "INSERT INTO notes (kind, entity, flag, author, created_at) VALUES (?, ?, ?, ?, ?)",
                    (kind, entity, flag, author, self._now()))
            self._conn.commit()
        return self.state(kind, entity)

    def add_note(self, kind: str, entity: str, text: str, author: str) -> dict[str, Any]:
        text = text.strip()
        if kind not in KINDS or not text:
            raise ValueError("unknown kind or empty note")
        with self._lock:
            self._conn.execute("INSERT INTO notes (kind, entity, text, author, created_at) VALUES (?, ?, ?, ?, ?)",
                               (kind, entity, text[:5000], author, self._now()))
            self._conn.commit()
        return self.state(kind, entity)

    def delete(self, note_id: int, author: str) -> tuple[str, str] | None:
        with self._lock:
            row = self._conn.execute("SELECT kind, entity, author FROM notes WHERE id = ? AND deleted_at IS NULL",
                                     (note_id,)).fetchone()
            if row is None or row["author"] != author:
                return None
            self._conn.execute("UPDATE notes SET deleted_at = ? WHERE id = ?", (self._now(), note_id))
            self._conn.commit()
        return row["kind"], row["entity"]

    def listing(self, *, flag: str = "", kind: str = "", author: str = "") -> list[dict[str, Any]]:
        where, params = ["deleted_at IS NULL"], []
        if flag:
            where.append("flag = ?")
            params.append(flag)
        if kind:
            where.append("kind = ?")
            params.append(kind)
        if author:
            where.append("author = ?")
            params.append(author)
        with self._lock:
            rows = self._conn.execute(f"SELECT * FROM notes WHERE {' AND '.join(where)} ORDER BY id DESC",
                                      params).fetchall()
        return [dict(r) for r in rows]

    def flagged_entities(self, flag: str, kind: str = "") -> list[tuple[str, str]]:
        """(kind, entity) pairs currently carrying a flag, oldest flag first."""
        where, params = ["deleted_at IS NULL", "flag = ?"], [flag]
        if kind:
            where.append("kind = ?")
            params.append(kind)
        with self._lock:
            rows = self._conn.execute(f"SELECT kind, entity, MIN(id) FROM notes WHERE {' AND '.join(where)} "
                                      "GROUP BY kind, entity ORDER BY MIN(id)", params).fetchall()
        return [(r[0], r[1]) for r in rows]

    def flag_index(self, kind: str, entities: list[str]) -> dict[str, list[str]]:
        """entity -> active flags, for marking rows in lists."""
        if not entities:
            return {}
        out: dict[str, list[str]] = {}
        with self._lock:
            for start in range(0, len(entities), 500):
                chunk = entities[start:start + 500]
                for entity, flag in self._conn.execute(
                        f"SELECT DISTINCT entity, flag FROM notes WHERE kind = ? AND flag IS NOT NULL AND "
                        f"deleted_at IS NULL AND entity IN ({','.join('?' * len(chunk))})", (kind, *chunk)):
                    out.setdefault(entity, []).append(flag)
        return out
