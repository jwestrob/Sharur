"""
Persistent hypothesis registry for cross-session tracking.

Hypotheses registered here survive across exploration sessions and
subagent runs, enabling multi-agent hypothesis development where
different subagents contribute evidence to shared hypotheses.

Storage: JSON file at {dataset_dir}/exploration/hypotheses.json
"""

from __future__ import annotations

import json
from pathlib import Path
from uuid import UUID

from sharur.core.types import Evidence, Hypothesis, HypothesisStatus


class HypothesisRegistry:
    """Persistent hypothesis store across exploration sessions."""

    def __init__(self, path: Path):
        self.path = path
        self._hypotheses: dict[UUID, Hypothesis] = {}
        if path.exists():
            self._load()

    def _load(self) -> None:
        text = self.path.read_text().strip()
        if not text:
            return
        data = json.loads(text)
        for h in data:
            hypo = Hypothesis(**h)
            self._hypotheses[hypo.hypothesis_id] = hypo

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps(
                [h.model_dump() for h in self._hypotheses.values()],
                default=str,
                indent=2,
            )
        )

    def register(self, hypothesis: Hypothesis) -> None:
        """Register a new hypothesis. Saves immediately."""
        self._hypotheses[hypothesis.hypothesis_id] = hypothesis
        self.save()

    def get(self, hypothesis_id: UUID) -> Hypothesis:
        """Get a hypothesis by ID. Raises KeyError if not found."""
        if hypothesis_id not in self._hypotheses:
            raise KeyError(f"Hypothesis {hypothesis_id} not found")
        return self._hypotheses[hypothesis_id]

    def add_evidence(self, hypothesis_id: UUID, evidence: Evidence) -> None:
        """Add evidence to an existing hypothesis. Saves immediately."""
        hypo = self.get(hypothesis_id)
        hypo.evidence.append(evidence)
        self.save()

    def update_status(
        self, hypothesis_id: UUID, status: HypothesisStatus
    ) -> None:
        """Update hypothesis status. Saves immediately."""
        hypo = self.get(hypothesis_id)
        hypo.status = status
        self.save()

    def list_all(self) -> list[Hypothesis]:
        """List all hypotheses."""
        return list(self._hypotheses.values())

    def list_active(self) -> list[Hypothesis]:
        """List hypotheses that are proposed or supported."""
        return [
            h
            for h in self._hypotheses.values()
            if h.status
            in (HypothesisStatus.PROPOSED, HypothesisStatus.SUPPORTED)
        ]

    def summary(self) -> str:
        """Human-readable summary of all hypotheses."""
        if not self._hypotheses:
            return "No hypotheses registered."
        lines = []
        for h in self._hypotheses.values():
            n_for = sum(1 for e in h.evidence if e.supports)
            n_against = sum(1 for e in h.evidence if not e.supports)
            lines.append(
                f"[{h.status.value}] {h.statement} "
                f"(+{n_for}/-{n_against}, conf={h.confidence:.2f})"
            )
        return "\n".join(lines)
