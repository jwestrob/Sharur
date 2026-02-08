"""
Bennu: Agent-driven metagenomic discovery system.

Bennu provides an LLM-powered interface for exploring metagenomic datasets.
Ask natural language questions and get structured, reproducible answers.
"""

__version__ = "0.1.0"

from bennu.core.models import (
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
from bennu.core import (
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
from bennu.core.session import ExplorationSession as SessionState
from bennu.agent.orchestrator import BennuAgent

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
    # Agent
    "BennuAgent",
]
