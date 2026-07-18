"""Typed records for biological case inspection and context comparison.

The operator layer owns database access.  These models are deliberately
storage-agnostic so cases and comparisons can be serialized into durable
evidence bundles without keeping a DuckDB connection alive.
"""

from __future__ import annotations

import json
import re
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


CASE_SCHEMA_VERSION = "1.0"
DEFAULT_CONTEXT_ORFS = 10


class CaseEntityType(str, Enum):
    """Entity types accepted by :meth:`Sharur.inspect`."""

    protein = "protein"
    system = "system"
    locus = "locus"
    contig = "contig"
    bin = "bin"


class EvidenceLevel(str, Enum):
    """Claim/evidence boundary used throughout case records."""

    observed = "observed"
    caller_named = "caller_named"
    inferred = "inferred"
    unverified = "unverified"
    experimental = "experimental"


class EntityReference(BaseModel):
    """Stable reference to the entity anchoring a case."""

    model_config = ConfigDict(frozen=True)

    entity_id: str
    entity_type: CaseEntityType
    source_table: str
    name: str | None = None
    subtype: str | None = None


class AnnotationEvidence(BaseModel):
    """Raw per-protein annotation evidence.

    These rows are observations.  Purpose-built caller output is represented
    separately by :class:`NamedCallEvidence`.
    """

    model_config = ConfigDict(frozen=True)

    protein_id: str
    source: str
    accession: str
    name: str | None = None
    description: str | None = None
    evalue: float | None = None
    score: float | None = None
    start_aa: int | None = None
    end_aa: int | None = None
    evidence_level: EvidenceLevel = EvidenceLevel.observed


class NamedCallEvidence(BaseModel):
    """One system/locus name emitted by a structured caller resource."""

    model_config = ConfigDict(frozen=True)

    call_id: str
    call_type: str
    call_subtype: str | None = None
    source_table: str
    system_source: str | None = None
    protein_id: str | None = None
    profile_name: str | None = None
    position: int | None = None
    score: float | None = None
    evidence_level: EvidenceLevel = EvidenceLevel.caller_named


class ProteinContextRecord(BaseModel):
    """One protein in the requested case neighborhood."""

    protein_id: str
    bin_id: str | None = None
    contig_id: str
    start: int
    end: int
    strand: str
    gene_index: int | None = None
    sequence_length: int | None = None
    sequence: str | None = None
    relative_orf: int | None = None
    region_role: Literal["upstream", "anchor", "downstream", "context"] = "context"
    is_component: bool = False
    predicates: list[str] = Field(default_factory=list)
    annotations: list[AnnotationEvidence] = Field(default_factory=list)
    named_calls: list[NamedCallEvidence] = Field(default_factory=list)


class ContextWindow(BaseModel):
    """Requested and realized ORF context around an entity."""

    default_orfs: int = DEFAULT_CONTEXT_ORFS
    upstream_orfs: int = DEFAULT_CONTEXT_ORFS
    downstream_orfs: int = DEFAULT_CONTEXT_ORFS
    orientation: Literal["+", "-", "mixed", "unknown"] = "unknown"
    orientation_basis: str = "anchor component strands"
    requested_min_gene_index: int | None = None
    requested_max_gene_index: int | None = None
    realized_upstream_orfs: int = 0
    realized_downstream_orfs: int = 0
    complete_upstream: bool = False
    complete_downstream: bool = False

    @model_validator(mode="after")
    def _non_negative(self) -> ContextWindow:
        for field_name in ("default_orfs", "upstream_orfs", "downstream_orfs"):
            if getattr(self, field_name) < 0:
                raise ValueError(f"{field_name} must be non-negative")
        return self


class AssemblyEvidenceRecord(BaseModel):
    """Optional contig-level assembly/host-assignment evidence."""

    model_config = ConfigDict(allow_inf_nan=False)

    bin_id: str
    contig_id: str
    coverage_mean: float | None = Field(default=None, ge=0)
    coverage_sd: float | None = Field(default=None, ge=0)
    coverage_cv: float | None = Field(default=None, ge=0)
    coverage_ratio_to_bin_median: float | None = Field(default=None, ge=0)
    mapped_reads: int | None = Field(default=None, ge=0)
    proper_pair_fraction: float | None = Field(default=None, ge=0, le=1)
    insert_size_median: float | None = Field(default=None, ge=0)
    insert_size_mad: float | None = Field(default=None, ge=0)
    snv_count: int | None = Field(default=None, ge=0)
    snv_density_per_kb: float | None = Field(default=None, ge=0)
    assembly_graph_degree: int | None = Field(default=None, ge=0)
    assembly_graph_component: str | None = None
    taxonomy: str | None = None
    taxonomy_method: str | None = None
    taxonomy_confidence: float | None = None
    taxonomy_congruent: bool | None = None
    gc_zscore: float | None = None
    tetranucleotide_distance: float | None = Field(default=None, ge=0, le=2)
    tetranucleotide_percentile: float | None = Field(
        default=None,
        ge=0,
        le=100,
    )
    source: str | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)

    @property
    def available_fields(self) -> list[str]:
        """Return populated analytical fields, excluding identifiers/metadata."""
        excluded = {"bin_id", "contig_id", "source", "metadata"}
        return [
            name
            for name, value in self.model_dump().items()
            if name not in excluded and value is not None
        ]


class CaseRecord(BaseModel):
    """Serializable, self-contained biological evidence case."""

    schema_version: str = CASE_SCHEMA_VERSION
    entity: EntityReference
    bin: dict[str, Any] | None = None
    contig: dict[str, Any] | None = None
    components: list[str] = Field(default_factory=list)
    context_window: ContextWindow | None = None
    proteins: list[ProteinContextRecord] = Field(default_factory=list)
    named_calls: list[NamedCallEvidence] = Field(default_factory=list)
    assembly_evidence: AssemblyEvidenceRecord | None = None
    assembly_evidence_state: Literal["available", "unavailable", "failed"] = "unavailable"
    limitations: list[str] = Field(default_factory=list)
    trace: dict[str, Any] = Field(default_factory=dict)

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.model_dump(mode="json"), indent=indent, default=str)


class ContextFeature(BaseModel):
    """One testable feature in a genomic-context comparison."""

    feature_id: str | None = None
    kind: Literal[
        "annotation",
        "annotation_name",
        "predicate",
        "called_system",
        "other_called_system",
        "locus_type",
    ]
    value: str | None = None
    source: str | None = None
    side: Literal["either", "upstream", "downstream", "anchor"] = "either"
    max_orfs: int | None = None
    include_anchor_call: bool = False

    @model_validator(mode="after")
    def _validate_feature(self) -> ContextFeature:
        if self.kind not in {"other_called_system"} and not self.value:
            raise ValueError(f"{self.kind} requires a value")
        if self.max_orfs is not None and self.max_orfs < 0:
            raise ValueError("max_orfs must be non-negative")
        if self.feature_id is None:
            parts = [value for value in (self.kind, self.source, self.value, self.side) if value]
            if self.max_orfs is not None:
                parts.append(f"max_orfs_{self.max_orfs}")
            if self.include_anchor_call:
                parts.append("include_anchor_call")
            raw = ":".join(parts)
            self.feature_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("_")
        return self

    @classmethod
    def parse(cls, value: ContextFeature | str | dict[str, Any]) -> ContextFeature:
        """Parse a compact feature string or validate an explicit record.

        Compact forms:

        - ``pfam:PF00589`` or ``annotation:PF00589``
        - ``name:integrase``
        - ``predicate:mobile_element``
        - ``system:RM_Type_III``
        - ``other_called_system``
        - ``locus:crispr``
        """
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            return cls(**value)
        if not isinstance(value, str) or not value.strip():
            raise ValueError("Context features must be non-empty strings or mappings")

        text = value.strip()
        if text == "other_called_system":
            return cls(kind="other_called_system")
        prefix, separator, payload = text.partition(":")
        if not separator or not payload:
            raise ValueError(
                f"Invalid context feature {text!r}; use source:accession, "
                "name:TEXT, predicate:ID, system:TYPE, locus:TYPE, or "
                "other_called_system"
            )
        normalized = prefix.lower()
        if normalized in {"annotation", "accession"}:
            return cls(kind="annotation", value=payload)
        if normalized == "name":
            return cls(kind="annotation_name", value=payload)
        if normalized == "predicate":
            return cls(kind="predicate", value=payload)
        if normalized in {"system", "called_system"}:
            return cls(kind="called_system", value=payload)
        if normalized in {"locus", "locus_type"}:
            return cls(kind="locus_type", value=payload)
        return cls(kind="annotation", source=prefix, value=payload)


class ContextGroupSummary(BaseModel):
    """Presence/absence counts for one comparison group."""

    group: Literal["foreground", "background"]
    positive: int
    total: int
    entity_ids: list[str] = Field(default_factory=list)
    unit_keys: list[str] = Field(default_factory=list)


class ContextStatistic(BaseModel):
    """Exact 2x2 comparison result for one feature or composite."""

    feature_id: str
    foreground_positive: int
    foreground_total: int
    background_positive: int
    background_total: int
    odds_ratio: float | None = None
    odds_ratio_ci95: tuple[float | None, float | None] | None = None
    fisher_p: float
    fisher_q: float | None = None
    alternative: Literal["greater", "less", "two-sided"] = "greater"


class ContextComparison(BaseModel):
    """Serializable output from :meth:`BiologicalCase.compare_context`."""

    schema_version: str = CASE_SCHEMA_VERSION
    case_entity_id: str
    case_entity_type: CaseEntityType
    foreground_definition: dict[str, Any]
    background_definition: dict[str, Any]
    features: list[ContextFeature]
    combine: Literal["all", "any"] = "all"
    context_window: ContextWindow
    require_full_context: bool = True
    deduplicate_by: Literal["entity", "replicon", "bin"] = "replicon"
    foreground: ContextGroupSummary
    background: ContextGroupSummary
    composite: ContextStatistic
    per_feature: list[ContextStatistic] = Field(default_factory=list)
    entity_matrix: list[dict[str, Any]] = Field(default_factory=list)
    recipe: dict[str, Any] = Field(default_factory=dict)
    limitations: list[str] = Field(default_factory=list)

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.model_dump(mode="json"), indent=indent, default=str)

    def to_markdown(self) -> str:
        """Render the comparison without hiding denominators or test choices."""
        lines = [
            f"# Context comparison: {self.case_entity_id}",
            "",
            f"- Context: {self.context_window.upstream_orfs} upstream / "
            f"{self.context_window.downstream_orfs} downstream ORFs",
            f"- Composite: `{self.combine}` across {len(self.features)} feature(s)",
            f"- Unit: `{self.deduplicate_by}`",
            f"- Full context required: `{self.require_full_context}`",
            "",
            "| Result | Foreground | Background | Odds ratio (95% CI) | Fisher p | BH q |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        statistics = [self.composite, *self.per_feature]
        for statistic in statistics:
            interval = statistic.odds_ratio_ci95 or (None, None)
            odds_ratio = "NA" if statistic.odds_ratio is None else f"{statistic.odds_ratio:.4g}"
            interval_text = (
                "NA"
                if interval[0] is None or interval[1] is None
                else f"{interval[0]:.4g} to {interval[1]:.4g}"
            )
            q_value = "NA" if statistic.fisher_q is None else f"{statistic.fisher_q:.4g}"
            lines.append(
                f"| `{statistic.feature_id}` | "
                f"{statistic.foreground_positive}/{statistic.foreground_total} | "
                f"{statistic.background_positive}/{statistic.background_total} | "
                f"{odds_ratio} ({interval_text}) | "
                f"{statistic.fisher_p:.4g} | {q_value} |"
            )
        if self.limitations:
            lines.extend(["", "## Limitations", ""])
            lines.extend(f"- {limitation}" for limitation in self.limitations)
        return "\n".join(lines)


__all__ = [
    "CASE_SCHEMA_VERSION",
    "DEFAULT_CONTEXT_ORFS",
    "AnnotationEvidence",
    "AssemblyEvidenceRecord",
    "CaseEntityType",
    "CaseRecord",
    "ContextComparison",
    "ContextFeature",
    "ContextGroupSummary",
    "ContextStatistic",
    "ContextWindow",
    "EntityReference",
    "EvidenceLevel",
    "NamedCallEvidence",
    "ProteinContextRecord",
]
