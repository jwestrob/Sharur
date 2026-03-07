"""
Predicate V2: Typed Semantic Atom System.

Replaces flat boolean predicates with typed semantic atoms, declarative
composite rules, and an unmapped-accession review queue. Designed to run
in shadow mode alongside V1 for validation before switchover.

Public API:
    - generate_and_persist_v2: Full pipeline: generate + persist to DuckDB
    - shadow_diff: Compare V1 and V2 outputs on the same data
    - build_review_queue: Surface unmapped accessions for curation
    - semantic_state_to_predicates: V2 state -> V1 flat predicate list
"""

from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
    create_v2_tables,
)
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.composites import evaluate_composites, load_composites
from sharur.predicates_v2.compat import semantic_state_to_predicates
from sharur.predicates_v2.review_queue import build_review_queue, format_review_queue_tsv
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.predicates_v2.shadow import shadow_diff, run_shadow_diff_on_store

__all__ = [
    # Data model
    "ClaimRelation",
    "SemanticAtom",
    "SemanticFacet",
    "SemanticState",
    # Pipeline
    "AtomGenerator",
    "aggregate_atoms",
    "evaluate_composites",
    "load_composites",
    "create_v2_tables",
    "generate_and_persist_v2",
    # Compatibility
    "semantic_state_to_predicates",
    # Review
    "build_review_queue",
    "format_review_queue_tsv",
    # Shadow diff
    "shadow_diff",
    "run_shadow_diff_on_store",
]
