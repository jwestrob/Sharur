"""Unit tests for the worker executors.

These cover the pure logic that sits between the Ops/Query contracts and the
model CLIs — the parts most likely to silently corrupt a campaign: output
parsing, the sequence-free guard, reduction-key stability, and rate-limit
classification. Live-service paths are exercised separately by
`sharur-worker atlas-scan --dry-run`.
"""

from __future__ import annotations

import json

import pytest

from sharur.workers import model_cli
from sharur.workers.atlas_scan import (
    SCAN_OUTPUT_SCHEMA,
    SCAN_SYSTEM_PROMPT,
    _assert_sequence_free,
    _canonical,
    _stable_key,
)


class TestExtractJson:
    def test_plain_object(self):
        assert model_cli._extract_json('{"candidates": []}') == {"candidates": []}

    def test_object_wrapped_in_prose(self):
        text = 'Here is my analysis.\n{"candidates": [{"a": 1}]}\nHope that helps.'
        assert model_cli._extract_json(text) == {"candidates": [{"a": 1}]}

    def test_takes_last_balanced_object(self):
        text = '{"draft": true}\nActually, revised:\n{"candidates": []}'
        assert model_cli._extract_json(text) == {"candidates": []}

    def test_nested_braces(self):
        payload = {"candidates": [{"signature": {"acc": {"pfam": "PF00485"}}}]}
        assert model_cli._extract_json(json.dumps(payload)) == payload

    def test_empty_stdout_raises(self):
        with pytest.raises(model_cli.ModelError, match="empty stdout"):
            model_cli._extract_json("   ")

    def test_no_json_raises(self):
        with pytest.raises(model_cli.ModelError, match="no JSON object"):
            model_cli._extract_json("I could not complete this task.")

    def test_non_object_json_raises(self):
        with pytest.raises(model_cli.ModelError):
            model_cli._extract_json("[1, 2, 3]")


class TestRateLimitClassification:
    @pytest.mark.parametrize(
        "stderr",
        [
            "Error: rate limit exceeded",
            "HTTP 429 Too Many Requests",
            "usage limit reached for this window",
            "RESOURCE_EXHAUSTED",
        ],
    )
    def test_detects_rate_limits(self, stderr):
        assert model_cli._looks_rate_limited(stderr, "")

    def test_ignores_unrelated_failures(self):
        assert not model_cli._looks_rate_limited("connection refused", "")

    def test_checks_stdout_too(self):
        assert model_cli._looks_rate_limited("", "quota exhausted")


class TestSequenceGuard:
    """The hard rule: raw biological sequence must never reach a model."""

    def test_clean_payload_passes(self):
        payload = {
            "contigs": [
                {"contig_id": "c1", "proteins": [{"protein_id": "p1", "length_aa": 300}]}
            ]
        }
        _assert_sequence_free(payload, "bin1")

    def test_sequence_key_raises(self):
        payload = {
            "contigs": [
                {"contig_id": "c1", "proteins": [{"protein_id": "p1", "sequence": "MKV"}]}
            ]
        }
        with pytest.raises(RuntimeError, match="carried sequence data"):
            _assert_sequence_free(payload, "bin1")

    def test_sequence_length_is_not_a_sequence(self):
        payload = {
            "contigs": [
                {"contig_id": "c1", "proteins": [{"protein_id": "p1", "sequence_length": 300}]}
            ]
        }
        _assert_sequence_free(payload, "bin1")

    def test_handles_empty_and_missing(self):
        _assert_sequence_free({}, "bin1")
        _assert_sequence_free({"contigs": []}, "bin1")
        _assert_sequence_free({"contigs": [{"contig_id": "c1"}]}, "bin1")


class TestReductionKeyStability:
    """Signatures are the exact reduction key; identical findings must collide."""

    def test_canonical_is_key_order_independent(self):
        assert _canonical({"b": 1, "a": 2}) == _canonical({"a": 2, "b": 1})

    def test_stable_key_is_deterministic(self):
        a = _stable_key("task-1", "0", _canonical({"x": 1}))
        b = _stable_key("task-1", "0", _canonical({"x": 1}))
        assert a == b and len(a) == 32

    def test_stable_key_separates_distinct_inputs(self):
        assert _stable_key("task-1", "0", "{}") != _stable_key("task-1", "1", "{}")

    def test_key_ignores_signature_field_order(self):
        a = _stable_key("t", "0", _canonical({"acc": "PF00485", "n": 3}))
        b = _stable_key("t", "0", _canonical({"n": 3, "acc": "PF00485"}))
        assert a == b


class TestScanContract:
    def test_schema_forbids_extra_fields(self):
        assert SCAN_OUTPUT_SCHEMA["additionalProperties"] is False
        item = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]
        assert item["additionalProperties"] is False

    def test_schema_requires_reduction_fields(self):
        required = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]["required"]
        for field in ("candidate_type", "signature", "evidence", "subject_refs"):
            assert field in required

    def test_schema_is_json_serialisable(self):
        json.dumps(SCAN_OUTPUT_SCHEMA)

    def test_prompt_carries_the_known_confounds(self):
        """Scanners meet these at corpus scale with no human in the loop."""
        for marker in (
            "Complex I",
            "NiFeSe_Hases",
            "hyddb_needs_curation",
            "PRK",
            "Ald_Xan_dh_C",
        ):
            assert marker in SCAN_SYSTEM_PROMPT, f"prompt lost the {marker} guard"

    def test_prompt_states_the_sequence_rule(self):
        assert "no biological sequences" in SCAN_SYSTEM_PROMPT

    def test_prompt_separates_observed_from_named(self):
        assert "OBSERVED" in SCAN_SYSTEM_PROMPT and "NAMED" in SCAN_SYSTEM_PROMPT


class TestModelRunProvenance:
    def test_provenance_carries_audit_fields(self):
        run = model_cli.ModelRun(
            provider="openai",
            model="gpt-5.6-terra",
            reasoning_effort="medium",
            record={"candidates": []},
            prompt_sha256="a" * 64,
            result_sha256="b" * 64,
            exit_status=0,
            stderr="",
            usage={"input_tokens": 10},
        )
        prov = run.provenance()
        for field in (
            "provider",
            "model",
            "reasoning_effort",
            "prompt_sha256",
            "result_sha256",
            "exit_status",
            "usage",
        ):
            assert field in prov

    def test_stderr_is_truncated(self):
        run = model_cli.ModelRun(
            provider="openai",
            model="m",
            reasoning_effort="medium",
            record={},
            prompt_sha256="",
            result_sha256="",
            exit_status=1,
            stderr="x" * 10_000,
        )
        assert len(run.provenance()["stderr_tail"]) == 2000


class TestProviderRegistry:
    def test_both_policy_providers_have_drivers(self):
        from sharur.review.models import load_review_policy

        policy = load_review_policy()
        for profile in policy.profiles.values():
            assert profile.provider in model_cli.PROVIDERS, (
                f"policy profile uses provider {profile.provider!r} with no driver"
            )

    def test_unknown_provider_raises(self):
        with pytest.raises(model_cli.ModelError, match="no driver"):
            model_cli.run_profile(
                provider="nope",
                model="m",
                reasoning_effort="medium",
                system_prompt="",
                payload_text="{}",
                output_schema={},
            )
