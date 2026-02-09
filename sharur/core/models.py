"""
Core data models for Sharur.

All data flows through Pydantic models. No raw dicts.
Every other module imports from here.
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field


# =============================================================================
# ENUMS
# =============================================================================


class AnnotationSource(str, Enum):
    """Known annotation database sources."""

    PFAM = "pfam"
    KEGG = "kegg"
    COG = "cog"
    CAZY = "cazy"
    TIGR = "tigr"
    EGGNOG = "eggnog"
    INTERPRO = "interpro"
    CUSTOM = "custom"


class LocusType(str, Enum):
    """Types of genomic loci that can be detected."""

    PROPHAGE = "prophage"
    BGC = "bgc"  # Biosynthetic gene cluster
    CRISPR = "crispr"
    OPERON = "operon"
    ISLAND = "island"  # Genomic island
    TRANSPOSON = "transposon"
    CUSTOM = "custom"


class AnomalyMetric(str, Enum):
    """Pre-computed anomaly detection metrics."""

    ISOLATION_SCORE = "isolation_score"
    GC_DEVIATION = "gc_deviation"
    LENGTH_ZSCORE = "length_zscore"
    ARCHITECTURE_RARITY = "architecture_rarity"
    TAX_INCONGRUENCE = "tax_incongruence"


# =============================================================================
# ANNOTATION MODELS
# =============================================================================


class Annotation(BaseModel):
    """A single annotation on a protein."""

    model_config = ConfigDict(frozen=True)

    source: AnnotationSource
    accession: str  # e.g., "PF00142", "K00370"
    description: Optional[str] = None  # e.g., "Nitrogenase iron protein"
    evalue: Optional[float] = None
    score: Optional[float] = None
    domain_start: Optional[int] = None  # Position within protein
    domain_end: Optional[int] = None

    @property
    def confidence(self) -> Literal["high", "medium", "low"]:
        """Categorize annotation confidence based on e-value."""
        if self.evalue is None:
            return "low"
        if self.evalue < 1e-20:
            return "high"
        if self.evalue < 1e-5:
            return "medium"
        return "low"

    def __str__(self) -> str:
        desc = self.description or self.accession
        if self.evalue:
            return f"{desc} (e={self.evalue:.1e})"
        return desc


# =============================================================================
# PROTEIN MODELS
# =============================================================================


class ProteinHit(BaseModel):
    """
    A protein with its metadata and annotations.

    This is the primary entity in Sharur. Proteins are returned from
    searches, stored in working sets, and linked to genomic regions.
    """

    model_config = ConfigDict(frozen=False)  # Mutable for adding computed fields

    # Core identifiers
    id: str
    contig_id: str
    bin_id: Optional[str] = None

    # Genomic coordinates
    start: int
    end: int
    strand: Literal["+", "-"]
    gene_index: int  # Position in contig (0-indexed)

    # Sequence (optional - can be large)
    sequence: Optional[str] = None
    sequence_length: int

    # Annotations
    annotations: List[Annotation] = Field(default_factory=list)
    annotation_text: Optional[str] = None  # Primary annotation string

    # Taxonomy (from parent bin)
    taxonomy: Optional[str] = None  # GTDB format: d__Archaea;p__...

    # Quality indicators (from parent bin)
    bin_completeness: Optional[float] = None
    bin_contamination: Optional[float] = None

    # Computed fields (added during queries)
    similarity_score: Optional[float] = None  # For similarity searches
    anomaly_scores: Dict[str, float] = Field(default_factory=dict)

    @property
    def length_bp(self) -> int:
        """Length in base pairs."""
        return abs(self.end - self.start)

    @property
    def domains(self) -> List[str]:
        """List of PFAM domain accessions."""
        return [
            a.accession for a in self.annotations if a.source == AnnotationSource.PFAM
        ]

    @property
    def all_accessions(self) -> List[str]:
        """All annotation accessions regardless of source."""
        return [a.accession for a in self.annotations]

    @property
    def best_annotation(self) -> Optional[Annotation]:
        """Annotation with lowest e-value."""
        valid = [a for a in self.annotations if a.evalue is not None]
        if not valid:
            return self.annotations[0] if self.annotations else None
        return min(valid, key=lambda a: a.evalue)

    @property
    def is_hypothetical(self) -> bool:
        """Check if protein lacks meaningful annotation."""
        if not self.annotations and not self.annotation_text:
            return True
        if self.annotation_text:
            lower = self.annotation_text.lower()
            return "hypothetical" in lower or "unknown" in lower
        return False

    @property
    def reliability_warnings(self) -> List[str]:
        """List of reliability concerns for this protein."""
        warnings = []
        if self.bin_contamination and self.bin_contamination > 5:
            warnings.append(
                f"High bin contamination ({self.bin_contamination:.1f}%)"
            )
        if self.bin_completeness and self.bin_completeness < 50:
            warnings.append(f"Low bin completeness ({self.bin_completeness:.1f}%)")
        best = self.best_annotation
        if best and best.evalue and best.evalue > 1e-5:
            warnings.append(f"Weak annotation (e={best.evalue:.1e})")
        if not self.annotations:
            warnings.append("No annotations")
        return warnings

    @property
    def taxonomy_phylum(self) -> Optional[str]:
        """Extract phylum from GTDB taxonomy string."""
        if not self.taxonomy:
            return None
        parts = self.taxonomy.split(";")
        for part in parts:
            if part.startswith("p__"):
                return part.replace("p__", "")
        return None

    def __str__(self) -> str:
        annot = self.annotation_text or "hypothetical"
        return f"{self.id}: {annot} ({self.sequence_length} aa)"

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ProteinHit):
            return self.id == other.id
        return False


class ProteinSet(BaseModel):
    """
    Collection of proteins with summary statistics.

    Used for search results, working sets, and any grouped protein collection.
    Automatically generates summary statistics on creation.
    """

    model_config = ConfigDict(frozen=False)

    proteins: List[ProteinHit]

    # Auto-generated metadata
    count: int = 0
    summary: str = ""
    notable: List[str] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)

    # Optional name for working sets
    set_id: Optional[str] = None
    set_name: Optional[str] = None
    created_at: Optional[datetime] = None

    def model_post_init(self, __context: Any) -> None:
        """Compute summary statistics after initialization."""
        self.count = len(self.proteins)
        if not self.summary:
            self.summary = self._generate_summary()
        if not self.notable:
            self.notable = self._find_notable()
        if not self.warnings:
            self.warnings = self._collect_warnings()

    def _generate_summary(self) -> str:
        """Generate natural language summary."""
        if not self.proteins:
            return "Empty protein set"

        # Count by taxonomy
        tax_counts: Dict[str, int] = {}
        for p in self.proteins:
            phylum = p.taxonomy_phylum or "Unknown"
            tax_counts[phylum] = tax_counts.get(phylum, 0) + 1

        # Count by bin
        bin_count = len(set(p.bin_id for p in self.proteins if p.bin_id))

        # Most common annotation
        annot_counts: Dict[str, int] = {}
        for p in self.proteins:
            if p.annotation_text and "hypothetical" not in p.annotation_text.lower():
                annot_counts[p.annotation_text] = (
                    annot_counts.get(p.annotation_text, 0) + 1
                )

        top_annot = (
            max(annot_counts.items(), key=lambda x: x[1])[0]
            if annot_counts
            else "various functions"
        )

        # Build summary
        if len(tax_counts) == 1:
            tax_str = list(tax_counts.keys())[0]
        else:
            top_tax = sorted(tax_counts.items(), key=lambda x: -x[1])[:2]
            tax_str = ", ".join([f"{t[0]} ({t[1]})" for t in top_tax])

        return f"{self.count} proteins across {bin_count} bins; {tax_str}; primarily {top_annot}"

    def _find_notable(self) -> List[str]:
        """Identify notable patterns."""
        notable = []

        if not self.proteins:
            return notable

        # Length outliers
        lengths = [p.sequence_length for p in self.proteins]
        if lengths:
            avg_len = sum(lengths) / len(lengths)
            long_proteins = [p for p in self.proteins if p.sequence_length > avg_len * 2]
            if long_proteins:
                notable.append(
                    f"{len(long_proteins)} unusually long proteins (>2x average)"
                )

        # Hypothetical proportion
        hypo_count = sum(1 for p in self.proteins if p.is_hypothetical)
        if self.count > 0 and hypo_count > self.count * 0.5:
            notable.append(f"{hypo_count}/{self.count} are hypothetical")

        # Bin clustering
        bin_counts: Dict[str, int] = {}
        for p in self.proteins:
            if p.bin_id:
                bin_counts[p.bin_id] = bin_counts.get(p.bin_id, 0) + 1
        if bin_counts:
            top_bin = max(bin_counts.items(), key=lambda x: x[1])
            if top_bin[1] > self.count * 0.5:
                notable.append(
                    f"Concentrated in {top_bin[0]} ({top_bin[1]} proteins)"
                )

        return notable

    def _collect_warnings(self) -> List[str]:
        """Collect reliability warnings across all proteins."""
        warnings = []

        # High contamination bins
        high_contam = sum(
            1
            for p in self.proteins
            if p.bin_contamination and p.bin_contamination > 5
        )
        if high_contam > 0:
            warnings.append(
                f"{high_contam} proteins from high-contamination bins (>5%)"
            )

        # Low completeness bins
        low_complete = sum(
            1
            for p in self.proteins
            if p.bin_completeness and p.bin_completeness < 50
        )
        if low_complete > 0:
            warnings.append(
                f"{low_complete} proteins from low-completeness bins (<50%)"
            )

        return warnings

    def to_id_list(self) -> List[str]:
        """Extract just protein IDs."""
        return [p.id for p in self.proteins]

    def filter(self, predicate: callable) -> "ProteinSet":
        """Return a new ProteinSet with filtered proteins."""
        filtered = [p for p in self.proteins if predicate(p)]
        return ProteinSet(proteins=filtered)

    def __len__(self) -> int:
        return self.count

    def __iter__(self):
        return iter(self.proteins)

    def __getitem__(self, idx: int) -> ProteinHit:
        return self.proteins[idx]


# =============================================================================
# GENOMIC CONTEXT MODELS
# =============================================================================


class GenomicRegion(BaseModel):
    """A contiguous genomic region with its genes."""

    model_config = ConfigDict(frozen=False)

    contig_id: str
    start: int
    end: int
    strand: Optional[Literal["+", "-"]] = None  # None if mixed

    proteins: List[ProteinHit]

    # Context
    anchor_protein_id: Optional[str] = None  # The protein this region is centered on
    window_type: Literal["bp", "genes"] = "genes"
    window_size: int = 10

    # Computed properties
    gc_content: Optional[float] = None
    intergenic_lengths: List[int] = Field(default_factory=list)

    # Visualization
    ascii_diagram: Optional[str] = None

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    @property
    def gene_count(self) -> int:
        return len(self.proteins)

    @property
    def has_strand_break(self) -> bool:
        """Check if region has genes on both strands."""
        strands = set(p.strand for p in self.proteins)
        return len(strands) > 1

    def get_anchor(self) -> Optional[ProteinHit]:
        """Get the anchor protein if set."""
        if not self.anchor_protein_id:
            return None
        for p in self.proteins:
            if p.id == self.anchor_protein_id:
                return p
        return None


class Locus(BaseModel):
    """A detected functional gene cluster."""

    model_config = ConfigDict(frozen=False)

    locus_id: str
    locus_type: LocusType

    contig_id: str
    start: int
    end: int

    proteins: List[ProteinHit]

    # Detection metadata
    confidence: float  # 0.0 to 1.0
    detection_method: str  # "domain_pattern", "hmm", "ml", etc.
    evidence: List[str] = Field(default_factory=list)

    # Optional detailed annotations
    subtype: Optional[str] = None  # e.g., "Mu-like" for prophages
    metadata: Dict[str, Any] = Field(default_factory=dict)

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    @property
    def gene_count(self) -> int:
        return len(self.proteins)


# =============================================================================
# ANOMALY MODELS
# =============================================================================


class AnomalyHit(BaseModel):
    """A detected statistical anomaly."""

    model_config = ConfigDict(frozen=True)

    protein_id: str
    protein: Optional[ProteinHit] = None

    metric: AnomalyMetric
    score: float  # Raw metric value
    zscore: float  # Standard deviations from mean
    context: str  # "global", "within_bin", "within_phylum"

    interpretation: str  # Human-readable explanation

    @property
    def severity(self) -> Literal["extreme", "high", "moderate"]:
        """Categorize anomaly severity."""
        z = abs(self.zscore)
        if z > 4:
            return "extreme"
        if z > 3:
            return "high"
        return "moderate"


# =============================================================================
# COMPARISON MODELS
# =============================================================================


class ComparisonResult(BaseModel):
    """Cross-genome comparison output."""

    model_config = ConfigDict(frozen=False)

    feature: str  # What was compared
    groupby: str  # How it was grouped

    # Flexible data structure
    data: Dict[str, Any]

    # Summary
    summary: str
    visualization_hint: Literal["heatmap", "bar", "scatter", "table"] = "table"

    # Statistical notes
    total_groups: int = 0
    groups_with_feature: int = 0


# =============================================================================
# HYPOTHESIS TRACKING
# =============================================================================


class Hypothesis(BaseModel):
    """A tracked scientific hypothesis."""

    model_config = ConfigDict(frozen=False)

    hypothesis_id: str
    statement: str  # "MAG_023 is a predator"

    evidence_for: List[str] = Field(default_factory=list)  # Step IDs supporting
    evidence_against: List[str] = Field(default_factory=list)  # Step IDs contradicting

    status: Literal["proposed", "supported", "weakened", "refuted"] = "proposed"
    confidence: float = 0.5  # 0.0 to 1.0

    created_at: datetime
    updated_at: datetime

    notes: str = ""

    def add_evidence(self, step_id: str, supports: bool) -> None:
        """Add evidence for or against the hypothesis."""
        if supports:
            self.evidence_for.append(step_id)
        else:
            self.evidence_against.append(step_id)

        # Update confidence
        total = len(self.evidence_for) + len(self.evidence_against)
        if total > 0:
            self.confidence = len(self.evidence_for) / total

        # Update status
        if self.confidence > 0.7:
            self.status = "supported"
        elif self.confidence < 0.3:
            self.status = "weakened"

        self.updated_at = datetime.now()


# =============================================================================
# EXPLORATION HISTORY
# =============================================================================


class ExplorationStep(BaseModel):
    """A single step in the exploration history."""

    model_config = ConfigDict(frozen=False)

    step_id: str
    timestamp: datetime

    # User input
    user_query: str
    resolved_references: Dict[str, str] = Field(
        default_factory=dict
    )  # "those" -> "ProteinSet[P001,P002]"

    # Execution
    tool_calls: List[Dict[str, Any]] = Field(
        default_factory=list
    )  # [{tool: "find_proteins", params: {...}, duration_ms: 123}]

    # Results
    result_summary: str
    entities_found: List[str] = Field(default_factory=list)  # IDs added to focus

    # Agent response
    response_text: str
    suggested_followups: List[str] = Field(default_factory=list)
    warnings_surfaced: List[str] = Field(default_factory=list)


# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    # Enums
    "AnnotationSource",
    "LocusType",
    "AnomalyMetric",
    # Core models
    "Annotation",
    "ProteinHit",
    "ProteinSet",
    "GenomicRegion",
    "Locus",
    "AnomalyHit",
    "ComparisonResult",
    "Hypothesis",
    "ExplorationStep",
]
