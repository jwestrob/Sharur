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
    _parse_json_field,
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


class TestTransientClassification:
    """A DNS blip must not consume a task attempt.

    Biotite resolves against external nameservers with `options timeout:2
    attempts:3`, so lookups fail intermittently under load. Observed in a live
    scan: `failed to lookup address information: Try again,
    url: wss://chatgpt.com/backend-api/codex/responses`.
    """

    @pytest.mark.parametrize(
        "stderr",
        [
            "failed to lookup address information: Try again",
            "failed to connect to websocket: IO error",
            "Temporary failure in name resolution",
            "connection reset by peer",
            "network is unreachable",
        ],
    )
    def test_detects_transport_failures(self, stderr):
        assert model_cli._looks_transient(stderr, "")

    def test_ignores_genuine_model_errors(self):
        assert not model_cli._looks_transient("invalid model name", "")

    def test_transient_raises_its_own_class(self):
        with pytest.raises(model_cli.ModelTransient):
            model_cli._classify_failure(1, "failed to lookup address information", "")

    def test_rate_limit_wins_over_transient(self):
        """A 429 can also mention a timeout; the two are handled differently."""
        with pytest.raises(model_cli.ModelRateLimited):
            model_cli._classify_failure(1, "429 rate limit; connection timed out", "")

    def test_unclassified_failure_is_a_model_error(self):
        with pytest.raises(model_cli.ModelError, match="exited 1"):
            model_cli._classify_failure(1, "some novel breakage", "")

    def test_transient_is_not_a_model_error_subclass(self):
        """Otherwise the generic handler would swallow it and burn an attempt."""
        assert not issubclass(model_cli.ModelTransient, model_cli.ModelError)
        assert not issubclass(model_cli.ModelRateLimited, model_cli.ModelError)


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
        for field in (
            "candidate_type",
            "signature",
            "evidence_json",
            "subject_refs_json",
        ):
            assert field in required

    def test_schema_is_json_serialisable(self):
        json.dumps(SCAN_OUTPUT_SCHEMA)

    def test_schema_is_openai_strict_compliant(self):
        """OpenAI strict mode rejects open objects and optional properties.

        A live scan failed with `invalid_json_schema` because signature/
        evidence/subject_refs were declared as free-form objects. Walk the
        whole schema so that cannot regress silently.
        """

        def walk(node, path="root"):
            if not isinstance(node, dict):
                return
            if node.get("type") == "object":
                assert node.get("additionalProperties") is False, f"{path} is not closed"
                assert set(node.get("required", [])) == set(node.get("properties", {})), (
                    f"{path}: every property must be required under strict mode"
                )
                for key, child in node.get("properties", {}).items():
                    walk(child, f"{path}.{key}")
            if node.get("type") == "array":
                walk(node.get("items", {}), f"{path}[]")

        walk(SCAN_OUTPUT_SCHEMA)

    def test_free_form_fields_travel_as_strings(self):
        """Evidence and refs stay free-form; only the reduction key is pinned."""
        props = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]["properties"]
        for field in ("evidence_json", "subject_refs_json"):
            assert props[field]["type"] == "string"
        assert props["signature"]["type"] == "object"

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


class TestJsonCarrierFields:
    """signature/evidence/subject_refs arrive as JSON strings under strict mode."""

    def test_parses_a_json_string(self):
        assert _parse_json_field('{"a": 1}', "signature") == {"a": 1}

    def test_passes_a_dict_through(self):
        """The Anthropic driver is not schema-constrained and may send a dict."""
        assert _parse_json_field({"a": 1}, "signature") == {"a": 1}

    def test_malformed_json_returns_none(self):
        assert _parse_json_field("{not json", "signature") is None

    def test_non_object_json_returns_none(self):
        assert _parse_json_field("[1,2,3]", "signature") is None

    def test_empty_returns_none(self):
        assert _parse_json_field("", "signature") is None
        assert _parse_json_field(None, "signature") is None


class TestSignatureNormalisation:
    """The signature is the exact reduction key and must be genome-invariant.

    The first pilot run produced a 1.00x collapse ratio — 398 of 399 clusters
    were singletons — because 91% of signatures carried `contig_id`. Identical
    biology in two genomes therefore hashed differently and cross-genome
    reduction was structurally impossible. Prompting alone cannot guarantee
    this, so it is enforced here.
    """

    @pytest.mark.parametrize(
        "key",
        [
            "contig_id",
            "contig_ids",
            "protein_ids",
            "coordinates",
            "gene_indices",
            "gene_indexes",
            "gene_index_start",
            "locus_end",
            "start",
            "end",
            "position",
            "offset",
        ],
    )
    def test_locator_keys_are_stripped(self, key):
        from sharur.workers.atlas_scan import _split_signature

        inv, loc = _split_signature({"accessions": ["A"], key: "whatever"}, "g1")
        assert key not in inv
        assert key in loc

    def test_invariant_fields_survive(self):
        from sharur.workers.atlas_scan import _split_signature

        sig = {"accessions": ["CofC", "CofD"], "n_genes": 2, "call_type": "defense"}
        inv, loc = _split_signature(sig, "g1")
        assert inv == sig
        assert loc == {}

    def test_values_embedding_the_genome_id_are_stripped(self):
        """Catches locator-ish values under a key the patterns do not match."""
        from sharur.workers.atlas_scan import _split_signature

        inv, loc = _split_signature(
            {"accessions": ["A"], "note": "found in GCA_963850385.1"}, "GCA_963850385.1"
        )
        assert "note" not in inv
        assert "note" in loc

    def test_same_finding_in_two_genomes_collapses(self):
        """The whole point: identical biology must produce one reduction key."""
        from sharur.workers.atlas_scan import _canonical, _split_signature

        a = {
            "accessions": ["CofC", "CofD", "F420_ligase"],
            "n_genes": 3,
            "contig_id": "genome_a__c1",
            "coordinates": [100, 900],
        }
        b = {
            "accessions": ["CofC", "CofD", "F420_ligase"],
            "n_genes": 3,
            "contig_id": "genome_b__c7",
            "coordinates": [50000, 50800],
        }
        ka, _ = _split_signature(a, "genome_a")
        kb, _ = _split_signature(b, "genome_b")
        assert _canonical(ka) == _canonical(kb)

    def test_distinct_findings_still_separate(self):
        from sharur.workers.atlas_scan import _canonical, _split_signature

        a, _ = _split_signature({"accessions": ["CofC"], "n_genes": 1}, "g")
        b, _ = _split_signature({"accessions": ["PRK"], "n_genes": 1}, "g")
        assert _canonical(a) != _canonical(b)


class TestCanonicalSignature:
    """Schema drift is a second, independent reduction killer.

    Stripping locators is not sufficient. Left free-form, the model invents a
    different schema per finding — six competing shapes appeared across 37
    signatures in the first pilot run (`accessions` vs `caller_accessions` vs
    `defense_calls` vs `components` vs `profiles`) — so the same system in two
    genomes still hashed differently. The reduction key is now a fixed
    structure, canonicalised here.
    """

    def test_accession_order_does_not_matter(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a = _canonical_signature({"accessions": ["B", "A"], "n_genes": 2, "system": ""})
        b = _canonical_signature({"accessions": ["A", "B"], "n_genes": 2, "system": ""})
        assert _canonical(a) == _canonical(b)

    def test_duplicate_accessions_collapse(self):
        from sharur.workers.atlas_scan import _canonical_signature

        sig = _canonical_signature(
            {"accessions": ["MzaB", "MzaC", "MzaB"], "n_genes": 3, "system": "Gao_Mza"}
        )
        assert sig["accessions"] == ["MzaB", "MzaC"]

    def test_shape_is_fixed_regardless_of_input(self):
        from sharur.workers.atlas_scan import _canonical_signature

        assert set(_canonical_signature({}).keys()) == {"accessions", "n_genes", "system"}

    def test_missing_fields_get_stable_defaults(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a = _canonical_signature({"accessions": ["A"]})
        b = _canonical_signature({"accessions": ["A"], "n_genes": 0, "system": ""})
        assert _canonical(a) == _canonical(b)

    def test_scalar_accession_is_tolerated(self):
        from sharur.workers.atlas_scan import _canonical_signature

        assert _canonical_signature({"accessions": "A"})["accessions"] == ["A"]

    def test_signature_is_a_structured_schema_field_not_a_json_string(self):
        """Only the reduction key needs pinning; evidence/refs stay free-form."""
        props = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]["properties"]
        assert props["signature"]["type"] == "object"
        assert props["signature"]["additionalProperties"] is False
        assert set(props["signature"]["required"]) == {"accessions", "n_genes", "system"}
        assert props["evidence_json"]["type"] == "string"
        assert props["subject_refs_json"]["type"] == "string"
