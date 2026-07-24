"""Versioned configuration models for the scientific review DAG."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator


REVIEW_POLICY_SCHEMA = "sharur-review-policy/1.0"
DEFAULT_POLICY_PATH = Path(__file__).with_name("default_policy.yaml")


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ExecutionProfile(StrictModel):
    """Semantic queue capability resolved to an exact execution identity."""

    provider: str = Field(min_length=1)
    model: str = Field(min_length=1)
    reasoning_effort: str = Field(min_length=1)
    capability: str = Field(min_length=1)
    max_pending_tasks: int = Field(default=1_000, ge=1)
    resource_request: dict[str, int] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)


class ReviewTier(StrictModel):
    """Quorum and routing rules for one review tier."""

    profiles: list[str] = Field(min_length=1)
    minimum_reviews: int = Field(default=1, ge=1)
    minimum_promote_reviews: int = Field(default=1, ge=1)
    minimum_confidence: float = Field(default=0.0, ge=0.0, le=1.0)
    require_distinct_providers: bool = False
    require_distinct_reviewers: bool = False
    blind_to_prior_scores: bool = True
    blind_to_other_reviews: bool = True
    next_tier: str | None = None

    @model_validator(mode="after")
    def validate_quorum(self) -> ReviewTier:
        if self.minimum_promote_reviews > self.minimum_reviews:
            raise ValueError(
                "minimum_promote_reviews cannot exceed minimum_reviews"
            )
        return self


class AuditPolicy(StrictModel):
    """Deterministic negative-control and rejected-case audit rates."""

    sampling_seed: str = Field(min_length=1)
    clear_unit_rate: float = Field(default=0.01, ge=0.0, le=1.0)
    rejected_cluster_rate: float = Field(default=0.05, ge=0.0, le=1.0)
    clear_profile: str
    rejected_profile: str
    stratification_fields: list[str] = Field(default_factory=list)


class RoutingPolicy(StrictModel):
    initial_tier: str
    independent_tier: str
    adjudication_tier: str
    canonical_tier: str
    materialization_profile: str


class ReviewLimits(StrictModel):
    max_evidence_tasks_per_review: int = Field(default=8, ge=0, le=100)
    event_batch_size: int = Field(default=500, ge=1, le=1_000)


class ReviewPolicy(StrictModel):
    """Complete, content-hashed routing contract."""

    schema_version: str = REVIEW_POLICY_SCHEMA
    name: str = Field(min_length=1)
    version: str = Field(min_length=1)
    controller_id: str = Field(min_length=1)
    rubric_version: str = Field(min_length=1)
    profiles: dict[str, ExecutionProfile]
    tiers: dict[str, ReviewTier]
    routing: RoutingPolicy
    audit: AuditPolicy
    limits: ReviewLimits = Field(default_factory=ReviewLimits)

    @model_validator(mode="after")
    def validate_references(self) -> ReviewPolicy:
        if self.schema_version != REVIEW_POLICY_SCHEMA:
            raise ValueError(
                f"Unsupported review policy schema {self.schema_version!r}"
            )
        profile_names = set(self.profiles)
        tier_names = set(self.tiers)
        referenced_profiles = {
            profile
            for tier in self.tiers.values()
            for profile in tier.profiles
        }
        referenced_profiles.update(
            {
                self.audit.clear_profile,
                self.audit.rejected_profile,
                self.routing.materialization_profile,
            }
        )
        missing_profiles = sorted(referenced_profiles - profile_names)
        if missing_profiles:
            raise ValueError(
                "Review policy references missing profiles: "
                + ", ".join(missing_profiles)
            )
        referenced_tiers = {
            self.routing.initial_tier,
            self.routing.independent_tier,
            self.routing.adjudication_tier,
            self.routing.canonical_tier,
        }
        referenced_tiers.update(
            tier.next_tier for tier in self.tiers.values() if tier.next_tier
        )
        missing_tiers = sorted(referenced_tiers - tier_names)
        if missing_tiers:
            raise ValueError(
                "Review policy references missing tiers: "
                + ", ".join(missing_tiers)
            )
        for name, tier in self.tiers.items():
            providers = {
                self.profiles[profile].provider for profile in tier.profiles
            }
            if tier.require_distinct_providers and len(providers) < tier.minimum_reviews:
                raise ValueError(
                    f"Tier {name!r} requires {tier.minimum_reviews} distinct "
                    "providers but its profiles cannot supply them"
                )
        return self

    @property
    def policy_hash(self) -> str:
        payload = self.model_dump(mode="json", exclude_none=False)
        return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()

    def profile(self, name: str) -> ExecutionProfile:
        try:
            return self.profiles[name]
        except KeyError as exc:
            raise KeyError(f"Unknown execution profile {name!r}") from exc

    def tier(self, name: str) -> ReviewTier:
        try:
            return self.tiers[name]
        except KeyError as exc:
            raise KeyError(f"Unknown review tier {name!r}") from exc


def load_review_policy(path: str | Path | None = None) -> ReviewPolicy:
    """Load and validate a review policy, defaulting to the packaged policy."""

    resolved = (
        Path(path).expanduser().resolve()
        if path is not None
        else DEFAULT_POLICY_PATH
    )
    try:
        raw = yaml.safe_load(resolved.read_text())
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"Could not load review policy {resolved}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"Review policy must be a YAML object: {resolved}")
    return ReviewPolicy.model_validate(raw)


__all__ = [
    "DEFAULT_POLICY_PATH",
    "REVIEW_POLICY_SCHEMA",
    "AuditPolicy",
    "ExecutionProfile",
    "ReviewPolicy",
    "ReviewTier",
    "load_review_policy",
]
