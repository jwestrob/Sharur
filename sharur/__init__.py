"""
Sharur: Agent-driven metagenomic discovery system.

Sharur provides an LLM-powered interface for exploring metagenomic datasets.
Ask natural language questions and get structured, reproducible answers.
"""

__version__ = "0.1.0"

from sharur.core import (
    AnomalyReport,
    AnomalySignal,
    AnomalyType,
    AssemblyEvidenceRecord,
    Bin,
    CaseEntityType,
    CaseRecord,
    ClaimValidationError,
    ContextComparison,
    ContextFeature,
    ContextWindow,
    Contig,
    Evidence,
    EvidenceLevel,
    ExplorationSession,
    FindingDraft,
    FocusEntity,
    Hypothesis,
    HypothesisStatus,
    Protein,
    ProvenanceEntry,
    ScientificClaim,
    Strand,
    VerificationRecord,
    WorkingSet,
)
from sharur.core.models import (
    Annotation,
    AnnotationSource,
    AnomalyHit,
    ComparisonResult,
    GenomicRegion,
    LegacyHypothesis,
    Locus,
    ProteinHit,
    ProteinSet,
)

# Backwards compatibility
from sharur.core.session import ExplorationSession as SessionState


__all__ = [
    # Version
    "__version__",
    # Core models
    "Annotation",
    "AnnotationSource",
    "AssemblyEvidenceRecord",
    "AnomalyHit",
    "ComparisonResult",
    "CaseEntityType",
    "CaseRecord",
    "ClaimValidationError",
    "ContextComparison",
    "ContextFeature",
    "ContextWindow",
    "GenomicRegion",
    "Hypothesis",
    "LegacyHypothesis",
    "Locus",
    "ProteinHit",
    "ProteinSet",
    # New core types
    "ExplorationSession",
    "WorkingSet",
    "Evidence",
    "EvidenceLevel",
    "HypothesisStatus",
    "AnomalyReport",
    "AnomalySignal",
    "AnomalyType",
    "Protein",
    "Contig",
    "Bin",
    "FocusEntity",
    "FindingDraft",
    "ProvenanceEntry",
    "Strand",
    "ScientificClaim",
    "VerificationRecord",
    # Session
    "SessionState",
]
