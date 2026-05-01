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

from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.compat import semantic_state_to_predicates
from sharur.predicates_v2.composites import (
    evaluate_composites,
    explain_composites,
    load_composites,
)
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
    create_v2_tables,
)
from sharur.predicates_v2.persistence import (
    generate_and_persist_v2,
    materialize_legacy_predicates_from_v2,
    materialize_semantic_terms_from_v2,
    refresh_composite_predicates_from_v2,
)
from sharur.predicates_v2.review_queue import build_review_queue, format_review_queue_tsv
from sharur.predicates_v2.shadow import (
    evaluate_shadow_gate,
    run_shadow_diff_on_store,
    shadow_diff,
)
from sharur.predicates_v2.validated_systems import (
    fetch_validated_system_annotations,
    materialize_system_proteins,
)


__all__ = [
    "AtomGenerator",
    "ClaimRelation",
    "SemanticAtom",
    "SemanticFacet",
    "SemanticState",
    "aggregate_atoms",
    "build_review_queue",
    "create_v2_tables",
    "evaluate_composites",
    "evaluate_shadow_gate",
    "explain_composites",
    "fetch_validated_system_annotations",
    "format_review_queue_tsv",
    "generate_and_persist_v2",
    "load_composites",
    "materialize_legacy_predicates_from_v2",
    "materialize_semantic_terms_from_v2",
    "materialize_system_proteins",
    "refresh_composite_predicates_from_v2",
    "run_shadow_diff_on_store",
    "semantic_state_to_predicates",
    "shadow_diff",
]
