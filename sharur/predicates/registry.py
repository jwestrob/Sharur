"""
Predicate registry for managing predicate definitions.

Predicates are boolean properties computed over proteins. Each predicate
has a definition including an evaluation query (SQL) or Python callable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Callable, Literal, Optional, Union

# Global registry instance
_registry: Optional["PredicateRegistry"] = None


@dataclass
class PredicateDefinition:
    """
    Definition of a predicate.

    Predicates can be evaluated via:
    1. SQL query (eval_query) - returns protein_ids where predicate is true
    2. Python callable (eval_func) - takes protein dict, returns bool
    """

    predicate_id: str
    name: str
    description: str
    category: Literal["size", "annotation", "composition", "context", "custom"]
    entity_type: Literal["protein", "contig", "bin"] = "protein"
    eval_query: Optional[str] = None  # SQL returning protein_ids
    eval_func: Optional[Callable[[dict], bool]] = None  # Python evaluator
    window_bp: Optional[int] = None  # For context predicates
    references: list[str] = field(default_factory=list)
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    def __post_init__(self):
        if self.eval_query is None and self.eval_func is None:
            raise ValueError("Predicate must have either eval_query or eval_func")


class PredicateRegistry:
    """Registry for predicate definitions."""

    def __init__(self):
        self._predicates: dict[str, PredicateDefinition] = {}

    def register(self, predicate: PredicateDefinition) -> None:
        """Register a predicate definition."""
        self._predicates[predicate.predicate_id] = predicate

    def get(self, predicate_id: str) -> Optional[PredicateDefinition]:
        """Get a predicate by ID."""
        return self._predicates.get(predicate_id)

    def list_predicates(
        self, category: Optional[str] = None
    ) -> list[PredicateDefinition]:
        """List all predicates, optionally filtered by category."""
        preds = list(self._predicates.values())
        if category:
            preds = [p for p in preds if p.category == category]
        return sorted(preds, key=lambda p: (p.category, p.predicate_id))

    def list_categories(self) -> list[str]:
        """List unique predicate categories."""
        return sorted(set(p.category for p in self._predicates.values()))

    def exists(self, predicate_id: str) -> bool:
        """Check if predicate exists."""
        return predicate_id in self._predicates

    def __len__(self) -> int:
        return len(self._predicates)

    def __contains__(self, predicate_id: str) -> bool:
        return predicate_id in self._predicates


def get_registry() -> PredicateRegistry:
    """Get the global predicate registry, initializing defaults if needed."""
    global _registry
    if _registry is None:
        _registry = PredicateRegistry()
        # Register defaults
        from sharur.predicates.defaults import register_defaults

        register_defaults(_registry)
    return _registry


def reset_registry() -> None:
    """Reset the global registry (for testing)."""
    global _registry
    _registry = None


__all__ = [
    "PredicateDefinition",
    "PredicateRegistry",
    "get_registry",
    "reset_registry",
]
