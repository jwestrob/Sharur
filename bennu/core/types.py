"""
Canonical Pydantic models for Bennu (per BENNU_PROJECT_SEED.md).

These types are the single source of truth for domain entities and are
imported by every other module. Keep them minimal and free of runtime
dependencies beyond the standard library and Pydantic.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Literal, Optional
from uuid import UUID, uuid4

from pydantic import BaseModel, Field


# --------------------------------------------------------------------------- #
# Genomic core types
# --------------------------------------------------------------------------- #


class Strand(str, Enum):
    """Strand orientation on a contig."""

    FORWARD = "+"
    REVERSE = "-"
    UNKNOWN = "?"


class Protein(BaseModel):
    """A protein-coding gene from a metagenome."""

    protein_id: str = Field(
        ..., description="Unique identifier, e.g., 'bin_001_contig_042_gene_007'"
    )
    contig_id: str
    bin_id: Optional[str] = None
    start: int = Field(..., ge=0, description="Start coordinate (0-indexed, inclusive)")
    end: int = Field(..., gt=0, description="End coordinate (0-indexed, exclusive)")
    strand: Strand
    sequence: Optional[str] = None  # Amino acid sequence, loaded on demand

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    @property
    def length_aa(self) -> int:
        return self.length_bp // 3


class Contig(BaseModel):
    """A contiguous sequence from assembly."""

    contig_id: str
    bin_id: Optional[str] = None
    length: int
    gc_content: float = Field(..., ge=0, le=1)
    is_circular: bool = False
    taxonomy: Optional[str] = None  # GTDB taxonomy string


class Bin(BaseModel):
    """A metagenome-assembled genome (MAG)."""

    bin_id: str
    completeness: float = Field(..., ge=0, le=100)
    contamination: float = Field(..., ge=0, le=100)
    taxonomy: Optional[str] = None
    n_contigs: int
    total_length: int


class Annotation(BaseModel):
    """A functional annotation on a protein."""

    protein_id: str
    source: Literal["pfam", "kegg", "cazy", "bgc", "cog", "custom"]
    accession: str  # e.g., "PF00142", "K00001"
    name: Optional[str] = None
    description: Optional[str] = None
    evalue: Optional[float] = None
    score: Optional[float] = None
    start_aa: Optional[int] = None  # Domain coordinates within protein
    end_aa: Optional[int] = None


# --------------------------------------------------------------------------- #
# Loci
# --------------------------------------------------------------------------- #


class LocusType(str, Enum):
    BGC = "bgc"  # Biosynthetic gene cluster
    PROPHAGE = "prophage"  # Integrated phage
    CRISPR = "crispr"  # CRISPR-Cas system
    OPERON = "operon"  # Co-transcribed genes
    METABOLIC = "metabolic"  # Metabolic pathway cluster
    CUSTOM = "custom"  # User-defined


class Locus(BaseModel):
    """A genomic region containing multiple genes."""

    locus_id: str
    locus_type: LocusType
    contig_id: str
    start: int
    end: int
    proteins: list[str]  # Protein IDs in order
    confidence: float = Field(..., ge=0, le=1)
    metadata: dict = Field(default_factory=dict)  # Type-specific info

    @property
    def n_genes(self) -> int:
        return len(self.proteins)


# --------------------------------------------------------------------------- #
# Anomaly detection
# --------------------------------------------------------------------------- #


class AnomalyType(str, Enum):
    STATISTICAL = "statistical"  # Embedding outlier, unusual length, etc.
    PATTERN = "pattern"  # Rule violation (incomplete system, etc.)
    NOVEL = "novel"  # Unannotated but interesting


class AnomalySignal(BaseModel):
    """A single anomaly signal for a protein."""

    protein_id: str
    signal_type: AnomalyType
    signal_name: str  # e.g., "isolation_score", "incomplete_nitrogenase"
    score: float  # Higher = more anomalous
    details: Optional[str] = None


class AnomalyReport(BaseModel):
    """Aggregated anomaly information for a protein."""

    protein_id: str
    signals: list[AnomalySignal]
    composite_score: float  # Weighted combination
    rank: int  # Rank within dataset


# --------------------------------------------------------------------------- #
# Session state
# --------------------------------------------------------------------------- #


class WorkingSet(BaseModel):
    """A named collection of proteins for iterative analysis."""

    set_id: UUID = Field(default_factory=uuid4)
    name: str
    description: Optional[str] = None
    protein_ids: list[str]
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    def __len__(self) -> int:
        return len(self.protein_ids)


class HypothesisStatus(str, Enum):
    PROPOSED = "proposed"
    SUPPORTED = "supported"
    REFUTED = "refuted"
    UNCERTAIN = "uncertain"


class Evidence(BaseModel):
    """A piece of evidence for/against a hypothesis."""

    query: str
    result_summary: str
    supports: bool  # True if supports, False if refutes
    confidence: float = Field(..., ge=0, le=1)
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class Hypothesis(BaseModel):
    """A tracked scientific hypothesis."""

    hypothesis_id: UUID = Field(default_factory=uuid4)
    statement: str
    status: HypothesisStatus = HypothesisStatus.PROPOSED
    evidence: list[Evidence] = Field(default_factory=list)
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    @property
    def confidence(self) -> float:
        if not self.evidence:
            return 0.5
        supporting = sum(e.confidence for e in self.evidence if e.supports)
        refuting = sum(e.confidence for e in self.evidence if not e.supports)
        total = supporting + refuting
        return supporting / total if total > 0 else 0.5


class FocusEntity(BaseModel):
    """An entity in the focus stack for pronoun resolution."""

    entity_type: Literal["protein", "locus", "bin", "contig", "set"]
    entity_id: str
    mentioned_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class ProvenanceEntry(BaseModel):
    """An entry in the provenance log."""

    entry_id: UUID = Field(default_factory=uuid4)
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    query: str
    tool_calls: list[dict]  # Tool name + parameters
    results_summary: str
    duration_ms: int
    error: Optional[str] = None


__all__ = [
    # Core genome types
    "Strand",
    "Protein",
    "Contig",
    "Bin",
    "Annotation",
    # Loci
    "LocusType",
    "Locus",
    # Anomalies
    "AnomalyType",
    "AnomalySignal",
    "AnomalyReport",
    # Session
    "WorkingSet",
    "HypothesisStatus",
    "Evidence",
    "Hypothesis",
    "FocusEntity",
    "ProvenanceEntry",
]
