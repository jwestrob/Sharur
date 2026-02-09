"""
Sharur: Agent-driven metagenomic discovery system.

Sharur provides an LLM-powered interface for exploring metagenomic datasets.
Ask natural language questions and get structured, reproducible answers.
"""

__version__ = "0.1.0"

from sharur.core.models import (
    Annotation,
    AnnotationSource,
    AnomalyHit,
    ComparisonResult,
    GenomicRegion,
    Hypothesis,
    Locus,
    ProteinHit,
    ProteinSet,
)
from sharur.core import (
    ExplorationSession,
    Evidence,
    HypothesisStatus,
    WorkingSet,
    AnomalyReport,
    AnomalySignal,
    AnomalyType,
    Protein,
    Contig,
    Bin,
    FocusEntity,
    ProvenanceEntry,
    Strand,
)
# Backwards compatibility
from sharur.core.session import ExplorationSession as SessionState

__all__ = [
    # Version
    "__version__",
    # Core models
    "Annotation",
    "AnnotationSource",
    "AnomalyHit",
    "ComparisonResult",
    "GenomicRegion",
    "Hypothesis",
    "Locus",
    "ProteinHit",
    "ProteinSet",
    # New core types
    "ExplorationSession",
    "WorkingSet",
    "Evidence",
    "HypothesisStatus",
    "AnomalyReport",
    "AnomalySignal",
    "AnomalyType",
    "Protein",
    "Contig",
    "Bin",
    "FocusEntity",
    "ProvenanceEntry",
    "Strand",
    # Session
    "SessionState",
]
