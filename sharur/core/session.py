"""
Exploration session state.

Tracks working sets, focus stack, hypotheses, and provenance, and exposes
handles to the DuckDB store and optional vector store.
"""

from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


from sharur.core.types import (
    Evidence,
    FocusEntity,
    Hypothesis,
    HypothesisStatus,
    ProvenanceEntry,
    WorkingSet,
)
from sharur.storage.duckdb_store import DuckDBStore
from sharur.storage.vector_store import FAISSStore

MAX_FOCUS_DEPTH = 10
FOCUS_EXPIRY_TURNS = 10


class ExplorationSession:
    """
    Manages state for an exploratory analysis session.

    State includes:
    - Working sets: Named protein collections
    - Focus stack: Recently referenced entities (for pronoun resolution)
    - Hypotheses: Tracked claims with evidence
    - Provenance: Complete audit trail
    """

    def __init__(
        self,
        session_id: Optional[uuid.UUID] = None,
        db_path: Optional[Path] = None,
        *,
        read_only: bool = False,
    ):
        self.session_id = session_id or uuid.uuid4()
        self.created_at = datetime.now(timezone.utc)
        self.read_only = read_only

        # State
        self._working_sets: dict[str, WorkingSet] = {}
        self._focus_stack: list[FocusEntity] = []  # Most recent last
        self._hypotheses: dict[uuid.UUID, Hypothesis] = {}
        self._provenance: list[ProvenanceEntry] = []

        # Database connection
        self._db_path = Path(db_path) if db_path else None
        self._db = DuckDBStore(db_path, read_only=read_only) if db_path else None
        self._vector_store = None
        self._vector_store_path: Optional[Path] = None
        self._vector_store_initialized = False
        self._vector_store_error: Optional[str] = None
        if db_path:
            self._discover_vector_store(Path(db_path))

        # Turn counter for focus expiry
        self._turn_counter = 0

    # ------------------------------------------------------------------ #
    # Properties
    # ------------------------------------------------------------------ #
    @property
    def db(self) -> DuckDBStore:
        if self._db is None:
            self._db = DuckDBStore(read_only=self.read_only)
        return self._db

    @property
    def vector_store(self):
        """Load the embedding index on first similarity operation only."""
        if not self._vector_store_initialized:
            self._load_vector_store()
        return self._vector_store

    @vector_store.setter
    def vector_store(self, store) -> None:
        self._vector_store = store
        self._vector_store_initialized = True
        self._vector_store_error = None

    @property
    def vector_store_available(self) -> bool:
        """Return whether an embeddings file was discovered without loading it."""
        return self._vector_store_path is not None

    @property
    def vector_store_path(self) -> Optional[Path]:
        """Return the discovered H5 path without constructing a FAISS index."""
        return self._vector_store_path

    @property
    def vector_store_error(self) -> Optional[str]:
        """Return the most recent vector-store initialization error, if any."""
        return self._vector_store_error

    # ------------------------------------------------------------------ #
    # Working sets
    # ------------------------------------------------------------------ #
    def create_set(self, name: str, protein_ids: list[str], description: Optional[str] = None) -> WorkingSet:
        ws = WorkingSet(name=name, protein_ids=protein_ids, description=description)
        self._working_sets[name.lower()] = ws
        self.push_focus("set", ws.set_id.hex)
        return ws

    def get_set(self, name: str) -> Optional[WorkingSet]:
        return self._working_sets.get(name.lower())

    def add_to_set(self, name: str, protein_ids: list[str]) -> None:
        ws = self.get_set(name)
        if not ws:
            raise KeyError(f"Working set '{name}' not found")
        ws.protein_ids = list(dict.fromkeys(ws.protein_ids + protein_ids))  # dedupe

    def remove_from_set(self, name: str, protein_ids: list[str]) -> None:
        ws = self.get_set(name)
        if not ws:
            raise KeyError(f"Working set '{name}' not found")
        ids_to_remove = set(protein_ids)
        ws.protein_ids = [pid for pid in ws.protein_ids if pid not in ids_to_remove]

    def delete_set(self, name: str) -> None:
        self._working_sets.pop(name.lower(), None)

    def list_sets(self) -> list[WorkingSet]:
        return list(self._working_sets.values())

    # ------------------------------------------------------------------ #
    # Focus stack
    # ------------------------------------------------------------------ #
    def push_focus(self, entity_type: str, entity_id: str) -> None:
        self._focus_stack.append(
            FocusEntity(entity_type=entity_type, entity_id=entity_id, mentioned_at=datetime.now(timezone.utc))
        )
        if len(self._focus_stack) > MAX_FOCUS_DEPTH:
            self._focus_stack = self._focus_stack[-MAX_FOCUS_DEPTH:]

    def get_focus(self, entity_type: Optional[str] = None) -> Optional[FocusEntity]:
        self._prune_focus()
        if not self._focus_stack:
            return None
        if entity_type is None:
            return self._focus_stack[-1]
        for ent in reversed(self._focus_stack):
            if ent.entity_type == entity_type:
                return ent
        return None

    def resolve_reference(self, reference: str) -> Optional[str]:
        """Resolve pronouns and simple phrases to entity_id."""
        ref = reference.lower().strip()

        # Direct set name lookup
        if ref in self._working_sets:
            return self._working_sets[ref].set_id.hex

        pronoun_single = {"it", "this", "that", "the protein"}
        pronoun_plural = {"those", "them", "these", "the proteins", "the set", "the results"}

        if ref in pronoun_single:
            ent = self.get_focus()
            return ent.entity_id if ent else None

        if ref in pronoun_plural:
            ent = self.get_focus("set") or self.get_focus()
            return ent.entity_id if ent else None

        # Type-qualified
        for type_word in ("protein", "locus", "bin", "contig", "set"):
            if type_word in ref:
                ent = self.get_focus(type_word)
                if ent:
                    return ent.entity_id

        # Partial match on working set names
        for name, ws in self._working_sets.items():
            if name in ref or ref in name:
                return ws.set_id.hex

        return None

    def clear_focus(self) -> None:
        self._focus_stack.clear()

    def _prune_focus(self) -> None:
        """Remove focus items that are too old in turn count."""
        if not self._focus_stack:
            return
        if len(self._focus_stack) <= MAX_FOCUS_DEPTH:
            return
        self._focus_stack = self._focus_stack[-MAX_FOCUS_DEPTH :]

    # ------------------------------------------------------------------ #
    # Vector store discovery
    # ------------------------------------------------------------------ #
    def _discover_vector_store(self, db_path: Path) -> None:
        """
        Discover H5 embeddings alongside the DuckDB without loading them.

        Checks both new standard path (embeddings/) and legacy path (stage06_embeddings/).
        """
        for directory_name in ("embeddings", "stage06_embeddings"):
            embeddings_dir = db_path.parent / directory_name
            h5_path = embeddings_dir / "protein_embeddings.h5"
            if h5_path.is_file():
                self._vector_store_path = h5_path
                return

    def _load_vector_store(self) -> None:
        """Construct the FAISS store once, when similarity is first requested."""
        self._vector_store_initialized = True
        if self._vector_store_path is None:
            return

        try:
            self._vector_store = FAISSStore(str(self._vector_store_path))
        except Exception as exc:
            self._vector_store = None
            self._vector_store_error = f"{type(exc).__name__}: {exc}"

    # ------------------------------------------------------------------ #
    # Hypotheses
    # ------------------------------------------------------------------ #
    def propose_hypothesis(self, statement: str) -> Hypothesis:
        hypo = Hypothesis(statement=statement)
        self._hypotheses[hypo.hypothesis_id] = hypo
        return hypo

    def add_evidence(
        self, hypothesis_id: uuid.UUID, query: str, result_summary: str, supports: bool, confidence: float
    ) -> None:
        if hypothesis_id not in self._hypotheses:
            raise KeyError("Hypothesis not found")
        ev = Evidence(query=query, result_summary=result_summary, supports=supports, confidence=confidence)
        self._hypotheses[hypothesis_id].evidence.append(ev)

    def update_hypothesis_status(self, hypothesis_id: uuid.UUID, status: HypothesisStatus) -> None:
        if hypothesis_id not in self._hypotheses:
            raise KeyError("Hypothesis not found")
        self._hypotheses[hypothesis_id].status = status

    def add_evidence_from_provenance(
        self,
        hypothesis_id: uuid.UUID,
        entry_id: uuid.UUID,
        supports: bool,
        confidence: float,
    ) -> None:
        """Add evidence to a hypothesis by referencing a provenance entry.

        Automatically populates query and result_summary from the referenced
        provenance entry, creating an auditable link between the hypothesis
        evidence chain and the actual analytical step that produced it.
        """
        entry = next(
            (p for p in self._provenance if p.entry_id == entry_id), None
        )
        if entry is None:
            raise KeyError(f"Provenance entry {entry_id} not found")
        if hypothesis_id not in self._hypotheses:
            raise KeyError(f"Hypothesis {hypothesis_id} not found")
        ev = Evidence(
            provenance_id=entry_id,
            query=entry.query,
            result_summary=entry.results_summary,
            supports=supports,
            confidence=confidence,
        )
        self._hypotheses[hypothesis_id].evidence.append(ev)

    def list_hypotheses(self, status: Optional[HypothesisStatus] = None) -> list[Hypothesis]:
        if status is None:
            return list(self._hypotheses.values())
        return [h for h in self._hypotheses.values() if h.status == status]

    # ------------------------------------------------------------------ #
    # Provenance
    # ------------------------------------------------------------------ #
    def log_query(
        self,
        query: str,
        tool_calls: list[dict],
        results_summary: str,
        duration_ms: int,
        error: Optional[str] = None,
        parent_ids: Optional[list[uuid.UUID]] = None,
    ) -> ProvenanceEntry:
        entry = ProvenanceEntry(
            query=query,
            tool_calls=tool_calls,
            results_summary=results_summary,
            duration_ms=duration_ms,
            error=error,
            parent_ids=parent_ids or [],
        )
        self._provenance.append(entry)
        self._turn_counter += 1
        self._prune_focus()
        return entry

    def get_provenance(self) -> list[ProvenanceEntry]:
        return self._provenance

    def export_provenance(self, format: str = "json") -> str:
        if format != "json":
            raise ValueError("Only json export is supported")
        return json.dumps([p.model_dump() for p in self._provenance], default=str, indent=2)

    # ------------------------------------------------------------------ #
    # Persistence
    # ------------------------------------------------------------------ #
    def save(self, path: Path) -> None:
        state = {
            "session_id": str(self.session_id),
            "created_at": self.created_at.isoformat(),
            "db_path": str(self._db_path) if self._db_path else None,
            "read_only": self.read_only,
            "turn_counter": self._turn_counter,
            "working_sets": [ws.model_dump() for ws in self._working_sets.values()],
            "focus_stack": [fe.model_dump() for fe in self._focus_stack],
            "hypotheses": [h.model_dump() for h in self._hypotheses.values()],
            "provenance": [p.model_dump() for p in self._provenance],
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        temporary_path.write_text(json.dumps(state, default=str, indent=2))
        temporary_path.replace(path)

    @classmethod
    def load(
        cls,
        path: Path,
        *,
        db_path: Optional[Path] = None,
        read_only: Optional[bool] = None,
    ) -> "ExplorationSession":
        data = json.loads(path.read_text())
        saved_db_path = data.get("db_path")
        resolved_db_path = db_path or (Path(saved_db_path) if saved_db_path else None)
        resolved_read_only = (
            bool(data.get("read_only", False))
            if read_only is None
            else read_only
        )
        session = cls(
            session_id=uuid.UUID(data["session_id"]),
            db_path=resolved_db_path,
            read_only=resolved_read_only,
        )
        session.created_at = datetime.fromisoformat(data["created_at"])
        session._working_sets = {ws["name"].lower(): WorkingSet(**ws) for ws in data["working_sets"]}
        session._focus_stack = [FocusEntity(**fe) for fe in data["focus_stack"]]
        session._hypotheses = {uuid.UUID(h["hypothesis_id"]): Hypothesis(**h) for h in data["hypotheses"]}
        session._provenance = [ProvenanceEntry(**p) for p in data["provenance"]]
        session._turn_counter = int(data.get("turn_counter", len(session._provenance)))
        return session

    def to_notebook(self) -> str:
        """Generate Jupyter notebook from session provenance."""
        cells = [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# Sharur Exploration Session\n",
                    f"**Session ID**: {self.session_id}\n",
                    f"**Created**: {self.created_at.isoformat()}\n",
                ],
            }
        ]
        for entry in self._provenance:
            cells.append(
                {
                    "cell_type": "markdown",
                    "metadata": {},
                    "source": [f"## Query: {entry.query}\n\nResults: {entry.results_summary}\n"],
                }
            )
        nb = {
            "nbformat": 4,
            "nbformat_minor": 5,
            "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}},
            "cells": cells,
        }
        return json.dumps(nb, indent=2)


__all__ = ["ExplorationSession"]

# Backwards compatibility alias
SessionState = ExplorationSession
