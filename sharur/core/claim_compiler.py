"""Executable claim-boundary validation for Sharur findings.

The compiler does not infer biology.  It checks that the author's declared
claim level is supported by the attached case, that caller-emitted names are
not silently promoted from raw observations, and that quantitative statements
have replayable verification records.
"""

from __future__ import annotations

import json
import math
import re
from typing import TYPE_CHECKING, Any, Literal

from pydantic import BaseModel, Field, PrivateAttr

from sharur.core.case_models import ContextComparison, EvidenceLevel


if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path

    from sharur.operators.cases import BiologicalCase


_NUMBER_RE = re.compile(
    r"(?<![A-Za-z0-9_.])"
    r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?%?"
    r"(?![A-Za-z0-9_.])"
)
_PRIORITY_RE = re.compile(
    r"\b(first(?:-ever)?|unprecedented|to our knowledge|novel occurrence)\b",
    re.IGNORECASE,
)
_DEFINITIVE_RE = re.compile(
    r"\b(confirms?|proves?|demonstrates?)\b",
    re.IGNORECASE,
)


class ClaimValidationError(ValueError):
    """Raised when a finding draft crosses an unsupported claim boundary."""


class VerificationRecord(BaseModel):
    """One executable or concretely replayable verification."""

    claim: str
    query: str
    expected: Any


class ScientificClaim(BaseModel):
    """One explicitly typed scientific statement."""

    text: str = Field(min_length=1)
    level: EvidenceLevel
    evidence_refs: list[str] = Field(default_factory=list)
    verification: list[VerificationRecord] = Field(default_factory=list)
    named_call_ids: list[str] = Field(default_factory=list)
    notes: str | None = None


class ClaimIssue(BaseModel):
    """One claim-compiler issue."""

    severity: Literal["error", "warning"]
    code: str
    message: str
    claim_index: int | None = None


class ClaimValidationReport(BaseModel):
    """Complete validation result for a draft."""

    valid: bool
    issues: list[ClaimIssue] = Field(default_factory=list)

    @property
    def errors(self) -> list[ClaimIssue]:
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ClaimIssue]:
        return [issue for issue in self.issues if issue.severity == "warning"]


def _sql_literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _sql_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _case_evidence_refs(case: BiologicalCase) -> set[str]:
    refs: set[str] = {f"entity:{case.record.entity.entity_id}"}
    if case.record.bin and case.record.bin.get("bin_id"):
        refs.add(f"bin:{case.record.bin['bin_id']}")
    if case.record.contig and case.record.contig.get("contig_id"):
        bin_id = (
            str(case.record.bin.get("bin_id"))
            if case.record.bin and case.record.bin.get("bin_id")
            else ""
        )
        refs.add(f"contig:{bin_id}:{case.record.contig['contig_id']}")
    for protein in case.record.proteins:
        refs.add(f"protein:{protein.protein_id}")
        for annotation in protein.annotations:
            refs.add(f"annotation:{protein.protein_id}:{annotation.source}:{annotation.accession}")
        for call in protein.named_calls:
            refs.add(f"call:{call.call_id}")
    for call in case.record.named_calls:
        refs.add(f"call:{call.call_id}")
    return refs


def _auto_verification(
    claim: ScientificClaim,
    case: BiologicalCase,
    comparison: ContextComparison | None,
) -> list[VerificationRecord]:
    """Generate deterministic verification from exact structured references."""
    generated: list[VerificationRecord] = []
    call_by_id = {
        call.call_id: call
        for call in [
            *case.record.named_calls,
            *(
                nested_call
                for protein in case.record.proteins
                for nested_call in protein.named_calls
            ),
        ]
    }
    annotation_by_ref = {
        (f"annotation:{protein.protein_id}:{annotation.source}:{annotation.accession}"): annotation
        for protein in case.record.proteins
        for annotation in protein.annotations
    }
    for ref in claim.evidence_refs:
        if ref.startswith("protein:"):
            protein_id = ref.split(":", 1)[1]
            generated.append(
                VerificationRecord(
                    claim=f"Protein {protein_id} exists in the dataset",
                    query=(
                        "SELECT COUNT(*) FROM proteins WHERE protein_id = "
                        f"{_sql_literal(protein_id)}"
                    ),
                    expected=1,
                )
            )
        elif ref in annotation_by_ref:
            _, protein_id, source, accession = ref.split(":", 3)
            generated.append(
                VerificationRecord(
                    claim=(f"{protein_id} carries observed {source}:{accession}"),
                    query=(
                        "SELECT COUNT(DISTINCT protein_id) FROM annotations "
                        f"WHERE protein_id = {_sql_literal(protein_id)} "
                        f"AND source = {_sql_literal(source)} "
                        f"AND accession = {_sql_literal(accession)}"
                    ),
                    expected=1,
                )
            )
        elif ref.startswith("call:"):
            call_id = ref.split(":", 1)[1]
            call = call_by_id.get(call_id)
            if call is not None:
                columns = {
                    str(row[0])
                    for row in case.store.execute(
                        """
                        SELECT column_name
                        FROM information_schema.columns
                        WHERE table_schema = 'main' AND table_name = ?
                        """,
                        [call.source_table],
                    )
                }
                if {"system_id", "system_type"} <= columns:
                    id_column = "system_id"
                    type_column = "system_type"
                elif {"locus_id", "locus_type"} <= columns:
                    id_column = "locus_id"
                    type_column = "locus_type"
                else:
                    continue
                generated.append(
                    VerificationRecord(
                        claim=(f"Structured caller emitted {call.call_type} as {call_id}"),
                        query=(
                            f"SELECT {_sql_identifier(type_column)} FROM "
                            f"{_sql_identifier(call.source_table)} WHERE "
                            f"{_sql_identifier(id_column)} = {_sql_literal(call_id)}"
                        ),
                        expected=call.call_type,
                    )
                )
        elif ref == "comparison:composite" and comparison is not None:
            generated.append(
                VerificationRecord(
                    claim="Context comparison reproduces the composite table",
                    query=("Re-run Sharur compare_context with evidence.comparison.recipe"),
                    expected={
                        "foreground_positive": comparison.composite.foreground_positive,
                        "foreground_total": comparison.composite.foreground_total,
                        "background_positive": comparison.composite.background_positive,
                        "background_total": comparison.composite.background_total,
                        "fisher_p": comparison.composite.fisher_p,
                    },
                )
            )
    deduplicated: list[VerificationRecord] = []
    seen = set()
    for verification in generated:
        key = (verification.query, json.dumps(verification.expected, sort_keys=True))
        if key not in seen:
            deduplicated.append(verification)
            seen.add(key)
    return deduplicated


def _expected_numeric_values(value: Any) -> list[float]:
    """Flatten finite numeric expected outputs, excluding booleans."""
    if isinstance(value, bool) or value is None:
        return []
    if isinstance(value, (int, float)):
        numeric = float(value)
        return [numeric] if math.isfinite(numeric) else []
    if isinstance(value, dict):
        return [numeric for item in value.values() for numeric in _expected_numeric_values(item)]
    if isinstance(value, (list, tuple)):
        return [numeric for item in value for numeric in _expected_numeric_values(item)]
    return []


def _numeric_token_candidates(token: str) -> list[tuple[float, float]]:
    """Return raw/fractional interpretations plus rounding tolerance."""
    percentage = token.endswith("%")
    text = token[:-1] if percentage else token
    numeric = float(text)
    mantissa, separator, exponent_text = text.lower().partition("e")
    exponent = int(exponent_text) if separator else 0
    decimals = len(mantissa.partition(".")[2]) if "." in mantissa else 0 if separator else None
    tolerance = (
        0.5 * (10.0 ** (exponent - decimals))
        if decimals is not None
        else 1e-12 * max(1.0, abs(numeric))
    )
    candidates = [(numeric, tolerance)]
    if percentage:
        candidates.append((numeric / 100.0, tolerance / 100.0))
    return candidates


def _token_matches_value(token: str, value: float) -> bool:
    return any(
        math.isclose(
            value,
            candidate,
            rel_tol=1e-12,
            abs_tol=max(tolerance, 1e-15),
        )
        for candidate, tolerance in _numeric_token_candidates(token)
    )


def _match_numeric_tokens(
    text: str,
    verification: Sequence[VerificationRecord],
) -> tuple[list[str], list[float]]:
    """Match every textual number to a distinct expected query output."""
    remaining = [
        numeric for item in verification for numeric in _expected_numeric_values(item.expected)
    ]
    unmatched: list[str] = []
    matched: list[float] = []
    for token in _NUMBER_RE.findall(text):
        match_index = next(
            (index for index, value in enumerate(remaining) if _token_matches_value(token, value)),
            None,
        )
        if match_index is None:
            unmatched.append(token)
        else:
            matched.append(remaining.pop(match_index))
    return unmatched, matched


class FindingDraft(BaseModel):
    """Claim-checked finding draft attached to a live biological case."""

    title: str = Field(min_length=1)
    category: str = Field(min_length=1)
    description: str = Field(min_length=1)
    claims: list[ScientificClaim] = Field(min_length=1)
    novelty: int = Field(default=0, ge=0, le=3)
    falsification: str | None = None
    literature_status: Literal["unchecked", "searched", "priority_supported"] = "unchecked"
    literature_references: list[str] = Field(default_factory=list)
    comparison: ContextComparison | None = None
    _case: BiologicalCase = PrivateAttr()

    @classmethod
    def from_case(
        cls,
        case: BiologicalCase,
        *,
        title: str,
        category: str,
        description: str,
        claims: Sequence[Mapping[str, Any]],
        novelty: int = 0,
        falsification: str | None = None,
        literature_status: str = "unchecked",
        literature_references: Sequence[str] | None = None,
        comparison: ContextComparison | None = None,
    ) -> FindingDraft:
        draft = cls(
            title=title,
            category=category,
            description=description,
            claims=[ScientificClaim(**dict(claim)) for claim in claims],
            novelty=novelty,
            falsification=falsification,
            literature_status=literature_status,
            literature_references=list(literature_references or []),
            comparison=comparison,
        )
        draft._case = case
        return draft

    @property
    def case(self) -> BiologicalCase:
        return self._case

    def validate_claims(self) -> ClaimValidationReport:
        """Validate provenance, naming, quantitation, and priority boundaries."""
        issues: list[ClaimIssue] = []
        has_replayable_verification = False
        supported_summary_numbers: list[float] = []
        allowed_refs = _case_evidence_refs(self.case)
        if self.comparison is not None:
            allowed_refs.add("comparison:composite")
            allowed_refs.update(
                f"comparison:{statistic.feature_id}" for statistic in self.comparison.per_feature
            )
        all_named_calls = [
            *self.case.record.named_calls,
            *(call for protein in self.case.record.proteins for call in protein.named_calls),
        ]
        named_calls = {call.call_id: call.call_type for call in all_named_calls}
        named_types = sorted(set(named_calls.values()), key=len, reverse=True)

        if self.novelty >= 2 and not (self.falsification or "").strip():
            issues.append(
                ClaimIssue(
                    severity="error",
                    code="missing_falsification",
                    message="Findings with novelty >= 2 require a testable falsification.",
                )
            )

        full_text = " ".join([self.title, self.description, *(claim.text for claim in self.claims)])
        if _PRIORITY_RE.search(full_text):
            if self.literature_status == "unchecked":
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="literature_unchecked",
                        message=(
                            "Priority language is blocked until a literature audit "
                            "is explicitly attached."
                        ),
                    )
                )
            elif not self.literature_references:
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="literature_references_missing",
                        message=("Literature status is checked but no references were attached."),
                    )
                )

        for index, claim in enumerate(self.claims):
            if not claim.text.strip():
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="empty_claim",
                        message="Claims must contain text.",
                        claim_index=index,
                    )
                )
            unknown_refs = sorted(
                ref
                for ref in set(claim.evidence_refs) - allowed_refs
                if not (
                    claim.level == EvidenceLevel.experimental and ref.startswith("experimental:")
                )
            )
            if unknown_refs:
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="unknown_evidence_reference",
                        message=f"Unknown evidence references: {', '.join(unknown_refs)}",
                        claim_index=index,
                    )
                )
            if claim.level != EvidenceLevel.unverified and not claim.evidence_refs:
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="missing_evidence_reference",
                        message=(
                            f"{claim.level.value} claims require explicit evidence references."
                        ),
                        claim_index=index,
                    )
                )
            if claim.level == EvidenceLevel.caller_named:
                declared_ids = set(claim.named_call_ids)
                referenced_ids = {
                    ref.split(":", 1)[1] for ref in claim.evidence_refs if ref.startswith("call:")
                }
                call_ids = declared_ids | referenced_ids
                if not call_ids:
                    issues.append(
                        ClaimIssue(
                            severity="error",
                            code="missing_caller_call",
                            message=(
                                "caller_named claims require at least one structured call ID."
                            ),
                            claim_index=index,
                        )
                    )
                unsupported = sorted(call_ids - named_calls.keys())
                if unsupported:
                    issues.append(
                        ClaimIssue(
                            severity="error",
                            code="unsupported_caller_name",
                            message=(
                                "Claim cites call IDs absent from the case's structured "
                                f"caller evidence: {', '.join(unsupported)}"
                            ),
                            claim_index=index,
                        )
                    )
            if claim.level == EvidenceLevel.observed:
                promoted = [
                    named_type
                    for named_type in named_types
                    if named_type.lower() in claim.text.lower()
                ]
                if promoted:
                    issues.append(
                        ClaimIssue(
                            severity="error",
                            code="named_claim_marked_observed",
                            message=(
                                "A caller-emitted name appears in a claim declared as "
                                "observed; declare it caller_named or weaken the wording."
                            ),
                            claim_index=index,
                        )
                    )
            if claim.level == EvidenceLevel.experimental and not any(
                ref.startswith("experimental:") for ref in claim.evidence_refs
            ):
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="experimental_evidence_missing",
                        message=(
                            "Experimental claims require an experimental:* evidence reference."
                        ),
                        claim_index=index,
                    )
                )
            if claim.level != EvidenceLevel.experimental and _DEFINITIVE_RE.search(claim.text):
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="overstated_language",
                        message=(
                            "Definitive language is reserved for direct experimental evidence."
                        ),
                        claim_index=index,
                    )
                )

            verification = list(claim.verification)
            if not verification:
                verification = _auto_verification(
                    claim,
                    self.case,
                    self.comparison,
                )
            has_replayable_verification = has_replayable_verification or bool(verification)
            if claim.level != EvidenceLevel.unverified and not verification:
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="claim_not_replayable",
                        message=(
                            f"{claim.level.value} claim has no executable or "
                            "concretely replayable verification."
                        ),
                        claim_index=index,
                    )
                )
            numeric_tokens = _NUMBER_RE.findall(claim.text)
            if numeric_tokens and not verification:
                issues.append(
                    ClaimIssue(
                        severity="error",
                        code="quantitative_claim_unverified",
                        message=("Quantitative claim contains numbers but no verification."),
                        claim_index=index,
                    )
                )
            elif numeric_tokens:
                unmatched_numbers, matched_numbers = _match_numeric_tokens(
                    claim.text,
                    verification,
                )
                supported_summary_numbers.extend(matched_numbers)
                if unmatched_numbers:
                    issues.append(
                        ClaimIssue(
                            severity="error",
                            code="numeric_decomposition_unverified",
                            message=(
                                "Claim numbers lack matching expected query outputs: "
                                + ", ".join(unmatched_numbers)
                            ),
                            claim_index=index,
                        )
                    )
            if claim.level == EvidenceLevel.inferred and not self.falsification:
                issues.append(
                    ClaimIssue(
                        severity="warning",
                        code="inference_without_falsification",
                        message=("Inferred claim has no finding-level falsification condition."),
                        claim_index=index,
                    )
                )

        summary_numbers = _NUMBER_RE.findall(f"{self.title} {self.description}")
        unsupported_summary_numbers = [
            token
            for token in summary_numbers
            if not any(_token_matches_value(token, value) for value in supported_summary_numbers)
        ]
        if unsupported_summary_numbers:
            issues.append(
                ClaimIssue(
                    severity="error",
                    code="summary_number_not_claimed",
                    message=(
                        "Numbers in the title/description must also occur in a "
                        "verified claim: " + ", ".join(unsupported_summary_numbers)
                    ),
                )
            )

        if not has_replayable_verification:
            issues.append(
                ClaimIssue(
                    severity="error",
                    code="finding_without_verification",
                    message=(
                        "A canonical finding requires at least one executable "
                        "or concretely replayable verification."
                    ),
                )
            )

        return ClaimValidationReport(
            valid=not any(issue.severity == "error" for issue in issues),
            issues=issues,
        )

    def compile(self, *, strict: bool = True) -> dict[str, Any]:
        """Compile a canonical finding-shaped record without writing it."""
        report = self.validate_claims()
        if strict and not report.valid:
            raise ClaimValidationError(
                "Finding draft failed claim validation: "
                + "; ".join(issue.message for issue in report.errors)
            )

        compiled_claims = []
        verification: list[dict[str, Any]] = []
        for claim in self.claims:
            claim_verification = list(claim.verification) or _auto_verification(
                claim,
                self.case,
                self.comparison,
            )
            compiled_claims.append(
                {
                    **claim.model_dump(mode="json", exclude={"verification"}),
                    "verification": [item.model_dump(mode="json") for item in claim_verification],
                }
            )
            verification.extend(item.model_dump(mode="json") for item in claim_verification)

        evidence: dict[str, Any] = {
            "case_schema_version": self.case.record.schema_version,
            "entity": self.case.record.entity.model_dump(mode="json"),
            "claims": compiled_claims,
            "assembly_evidence_state": self.case.record.assembly_evidence_state,
            "claim_validation": report.model_dump(mode="json"),
        }
        if self.comparison is not None:
            evidence["comparison"] = self.comparison.model_dump(mode="json")

        finding: dict[str, Any] = {
            "title": self.title,
            "category": self.category,
            "description": self.description,
            "evidence": evidence,
            "verification": verification,
            "protein_ids": list(self.case.record.components),
            "contigs": (
                [str(self.case.record.contig["contig_id"])]
                if self.case.record.contig and self.case.record.contig.get("contig_id")
                else []
            ),
            "novelty": self.novelty,
            "provenance": {
                "query_used": "Sharur inspect + claim compiler",
                "dataset": self.case.record.trace.get("database"),
                "case_entity_id": self.case.record.entity.entity_id,
            },
        }
        if self.falsification:
            finding["falsification"] = self.falsification
        if self.literature_status != "unchecked" or self.literature_references:
            finding["provenance"]["literature_status"] = self.literature_status
            finding["provenance"]["literature_references"] = list(self.literature_references)
        return finding

    def to_json(self, *, indent: int | None = 2, strict: bool = True) -> str:
        return json.dumps(self.compile(strict=strict), indent=indent, default=str)

    def publish_bundle(
        self,
        output_dir: str | Path,
        *,
        include_sequences: bool = True,
        include_plot: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """Validate and export this finding with its complete evidence case."""
        _ = self.compile(strict=True)
        return self.case.export(
            output_dir,
            comparison=self.comparison,
            finding=self,
            include_sequences=include_sequences,
            include_plot=include_plot,
            overwrite=overwrite,
        )


__all__ = [
    "ClaimIssue",
    "ClaimValidationError",
    "ClaimValidationReport",
    "FindingDraft",
    "ScientificClaim",
    "VerificationRecord",
]
