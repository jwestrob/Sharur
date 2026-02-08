"""
Predicate system for Bennu.

Predicates are boolean properties computed over proteins that enable
semantic filtering (e.g., "giant unannotated proteins").

Key components:
- vocabulary: Comprehensive predicate definitions with hierarchy
- mappings: PFAM/KEGG/CAZy annotation to predicate mappings
- generator: Main interface for computing predicates from annotations
- topology: Transmembrane topology prediction (requires pyTMHMM)
- registry: Legacy predicate registration system
- evaluator: Legacy SQL-based predicate evaluation

For new code, prefer using the generator:
    from bennu.predicates.generator import PredicateGenerator
    gen = PredicateGenerator()
    predicates = gen.generate_for_protein(protein, annotations)

For topology predictions (requires bennu[topology]):
    from bennu.predicates.topology import predict_topology, get_topology_predicates
    prediction = predict_topology(sequence)
    if prediction:
        predicates = get_topology_predicates(prediction)
"""

from bennu.predicates.registry import (
    PredicateDefinition,
    PredicateRegistry,
    get_registry,
)
from bennu.predicates.evaluator import (
    evaluate_predicate,
    compute_all_predicates,
    compute_predicates_for_protein,
    persist_predicates,
)
from bennu.predicates.vocabulary import (
    PredicateVocab,
    ALL_PREDICATES,
    PREDICATE_BY_ID,
    get_predicate,
    list_predicates as list_vocab_predicates,
    list_categories,
    get_hierarchy,
)
from bennu.predicates.generator import (
    AnnotationRecord,
    ProteinRecord,
    PredicateGenerator,
    generate_predicates_for_proteins,
    persist_generated_predicates,
)
from bennu.predicates.topology import (
    is_available as topology_is_available,
    TopologyPrediction,
    predict_topology,
    get_topology_predicates,
    predict_topology_batch,
)

__all__ = [
    # Registry (legacy)
    "PredicateDefinition",
    "PredicateRegistry",
    "get_registry",
    # Evaluator (legacy)
    "evaluate_predicate",
    "compute_all_predicates",
    "compute_predicates_for_protein",
    "persist_predicates",
    # Vocabulary
    "PredicateVocab",
    "ALL_PREDICATES",
    "PREDICATE_BY_ID",
    "get_predicate",
    "list_vocab_predicates",
    "list_categories",
    "get_hierarchy",
    # Generator (preferred)
    "AnnotationRecord",
    "ProteinRecord",
    "PredicateGenerator",
    "generate_predicates_for_proteins",
    "persist_generated_predicates",
    # Topology (optional - requires pyTMHMM)
    "topology_is_available",
    "TopologyPrediction",
    "predict_topology",
    "get_topology_predicates",
    "predict_topology_batch",
]
