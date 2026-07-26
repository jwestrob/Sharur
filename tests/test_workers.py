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
from pathlib import Path
from sharur.workers.atlas_scan import (
    SCAN_OUTPUT_SCHEMA,
    SCAN_SYSTEM_PROMPT,
    _assert_sequence_free,
    _canonical,
    _parse_json_field,
    _conflict_detail,
    _is_lease_failure,
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
    """The reduction key must be the coarse identity of a finding.

    Two rounds of measurement on live output drove this shape:

    * Free-form signatures produced six competing schemas and a 1.00x collapse.
    * Pinning the schema fixed the form but not the granularity — keying on the
      full accession set gave 1.12x with 94% singletons, because two genuine
      Gabija systems differ by a domain or two and never collide. Sharing fell
      from 16% of clusters at five accessions to 2% at ten.
    * Keying on the normalised system name gave 2.06x, with real clusters
      (surface-polysaccharide across 17 genomes, mokosh-type-i across 17).

    The accession set is preserved in reduction_features so the variation is
    still visible to reviewers without fragmenting the key.
    """

    def test_same_system_different_gene_complement_collapses(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature(
            {"accessions": ["GajA", "GajB"], "n_genes": 2, "system": "Gabija"},
            "atlas-scan:defense-locus",
        )
        b, _ = _canonical_signature(
            {"accessions": ["GajA", "GajB", "UvrD"], "n_genes": 3, "system": "gabija-system"},
            "atlas-scan:defense-locus",
        )
        assert _canonical(a) == _canonical(b)

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("surface-polysaccharide-biosynthesis", "surface-polysaccharide"),
            ("surface-polysaccharide-locus", "surface-polysaccharide"),
            ("Mokosh type I", "mokosh-type-i"),
            ("Gabija", "gabija"),
            ("CBASS_operon", "cbass"),
        ],
    )
    def test_system_label_drift_is_normalised(self, raw, expected):
        from sharur.workers.atlas_scan import _normalise_system

        assert _normalise_system(raw) == expected

    def test_different_systems_stay_separate(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature({"system": "gabija"}, "atlas-scan:defense-locus")
        b, _ = _canonical_signature({"system": "hachiman"}, "atlas-scan:defense-locus")
        assert _canonical(a) != _canonical(b)

    def test_same_system_different_class_stays_separate(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature({"system": "x"}, "atlas-scan:defense-locus")
        b, _ = _canonical_signature({"system": "x"}, "atlas-scan:mobile-element-cargo")
        assert _canonical(a) != _canonical(b)

    def test_unnamed_findings_fall_back_to_accessions(self):
        """A coarse "unnamed" key would wrongly merge unrelated findings."""
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature(
            {"accessions": ["A", "B"], "system": ""}, "atlas-scan:novel-gene-cluster"
        )
        b, _ = _canonical_signature(
            {"accessions": ["C", "D"], "system": ""}, "atlas-scan:novel-gene-cluster"
        )
        assert _canonical(a) != _canonical(b)
        assert "accessions" in a

    def test_accession_order_and_duplicates_do_not_matter(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, fa = _canonical_signature(
            {"accessions": ["B", "A", "A"], "system": ""}, "atlas-scan:x"
        )
        b, fb = _canonical_signature({"accessions": ["A", "B"], "system": ""}, "atlas-scan:x")
        assert _canonical(a) == _canonical(b)
        assert fa["accessions"] == fb["accessions"] == ["A", "B"]

    def test_features_preserve_what_the_key_drops(self):
        from sharur.workers.atlas_scan import _canonical_signature

        key, feat = _canonical_signature(
            {"accessions": ["GajA", "GajB"], "n_genes": 2, "system": "Gabija"},
            "atlas-scan:defense-locus",
        )
        assert "accessions" not in key
        assert feat["accessions"] == ["GajA", "GajB"]
        assert feat["n_genes"] == 2
        assert "subtype" in feat

    def test_scalar_accession_is_tolerated(self):
        from sharur.workers.atlas_scan import _canonical_signature

        _, feat = _canonical_signature({"accessions": "A", "system": ""}, "atlas-scan:x")
        assert feat["accessions"] == ["A"]

    def test_signature_is_a_structured_schema_field_not_a_json_string(self):
        """Only the reduction key needs pinning; evidence/refs stay free-form."""
        props = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]["properties"]
        assert props["signature"]["type"] == "object"
        assert props["signature"]["additionalProperties"] is False
        assert props["evidence_json"]["type"] == "string"
        assert props["subject_refs_json"]["type"] == "string"



class TestBaseTypeAndSubtype:
    """`system` groups; `subtype` refines. Only `system` enters the key.

    Free-form labels fragmented the groups they should form: five defence
    findings on live output arrived as five singleton labels
    (defense-island-cbass-pmt-cargo, -mokosh-signaling-array, ...) rather than
    one group of eight, costing ~25% of achievable collapse. Detail past the
    base type is demoted to `subtype`, which rides in reduction_features.
    """

    @pytest.mark.parametrize(
        "emitted,system,subtype",
        [
            ("defense-island-cbass-pmt-cargo", "defense-island", "cbass-pmt-cargo"),
            ("secretion-system-t2ss", "secretion-system", "t2ss"),
            ("glycan-locus", "glycan-locus", ""),
            ("defense-island", "defense-island", ""),
        ],
    )
    def test_compound_labels_are_demoted(self, emitted, system, subtype):
        from sharur.workers.atlas_scan import _canonical_signature

        key, feat = _canonical_signature({"system": emitted}, "atlas-scan:x")
        assert key["system"] == system
        assert feat["subtype"] == subtype

    def test_variants_of_one_module_share_a_key(self):
        """The whole point: eight defence islands must form one group."""
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature(
            {"system": "defense-island-cbass-pmt-cargo"}, "atlas-scan:defense-locus"
        )
        b, _ = _canonical_signature(
            {"system": "defense-island-mokosh-signaling-array"}, "atlas-scan:defense-locus"
        )
        assert _canonical(a) == _canonical(b)

    def test_subtype_never_enters_the_key(self):
        from sharur.workers.atlas_scan import _canonical, _canonical_signature

        a, _ = _canonical_signature(
            {"system": "glycan-locus", "subtype": "pul"}, "atlas-scan:x"
        )
        b, _ = _canonical_signature(
            {"system": "glycan-locus", "subtype": "capsule-eps"}, "atlas-scan:x"
        )
        assert _canonical(a) == _canonical(b)

    def test_subtype_is_preserved_as_a_feature(self):
        from sharur.workers.atlas_scan import _canonical_signature

        _, feat = _canonical_signature(
            {"system": "glycan-locus", "subtype": "PUL"}, "atlas-scan:x"
        )
        assert feat["subtype"] == "pul"

    def test_base_types_survive_drift_stripping(self):
        """`-locus`/`-system`/`-island` are drift tokens AND real base types."""
        from sharur.workers.atlas_scan import _BASE_MODULE_TYPES, _normalise_system

        for base in _BASE_MODULE_TYPES:
            assert _normalise_system(base) == base, base

    def test_novel_labels_are_left_intact(self):
        from sharur.workers.atlas_scan import _canonical_signature

        key, _ = _canonical_signature({"system": "novel-weird-thing"}, "atlas-scan:x")
        assert key["system"] == "novel-weird-thing"

    def test_schema_requires_both_fields(self):
        props = SCAN_OUTPUT_SCHEMA["properties"]["candidates"]["items"]["properties"]
        assert set(props["signature"]["required"]) == {
            "accessions",
            "n_genes",
            "system",
            "subtype",
        }


class TestUnboundedFramesWithStallDetection:
    """A working call runs as long as it wants; a silent one is killed.

    A hard wall-clock cap discards a hard genome's work for being slow, which
    is the wrong failure to optimise for. `codex exec --json` streams events
    throughout, so silence — not duration — is what indicates a dead call.
    """

    def test_an_active_process_outlives_the_stall_window(self):
        import time

        from sharur.workers.model_cli import _run

        started = time.monotonic()
        code, out, _ = _run(
            ["bash", "-c", "for i in 1 2 3 4; do echo tick; sleep 1; done"], "", 2
        )
        assert code == 0
        assert time.monotonic() - started > 2, "should not have been killed at the stall window"
        assert out.count("tick") == 4

    def test_a_silent_process_is_killed(self):
        from sharur.workers.model_cli import _run

        code, _, err = _run(["bash", "-c", "sleep 30"], "", 2)
        assert code != 0
        assert "stalled connection" in err

    def test_a_stall_is_classified_transient_not_a_defect(self):
        """So it retries in place rather than consuming a task attempt."""
        from sharur.workers.model_cli import _looks_transient

        assert _looks_transient(
            "[sharur] terminated after 900s with no output; "
            "treating as a stalled connection",
            "",
        )

    def test_stdout_and_stderr_are_both_captured(self):
        from sharur.workers.model_cli import _run

        code, out, err = _run(["bash", "-c", "echo O; echo E >&2"], "", 10)
        assert code == 0 and "O" in out and "E" in err

    def test_stdin_reaches_the_process(self):
        from sharur.workers.model_cli import _run

        code, out, _ = _run(["cat"], "payload-here", 10)
        assert code == 0 and "payload-here" in out

    def test_keepalive_exists_and_is_a_context_manager(self):
        from sharur.workers.atlas_scan import AtlasScanWorker

        assert hasattr(AtlasScanWorker, "_lease_keepalive")
        assert hasattr(AtlasScanWorker._lease_keepalive, "__wrapped__")


class TestCandidateIdempotencyKey:
    """The key must separate distinct records that share a collapsed signature.

    Regression for the livelock that stalled the Dormibacteria pilot at 99/944.
    Reduction deliberately strips locators from the signature so the same system
    collapses across genomes; that also makes two unrelated findings in one
    frame share a signature. Keying only on the signature made them collide,
    the Ops service rejected the second with HTTP 409, and the worker read every
    409 as a lost lease and abandoned the genome — forever, on every retry.
    """

    @staticmethod
    def _key(cand, *, task="t1", frame="0"):
        return _stable_key(
            task,
            frame,
            _canonical(
                {
                    "candidate_type": cand["candidate_type"],
                    "signature": cand["signature"],
                    "subject_refs": cand["subject_refs"],
                }
            ),
        )

    @staticmethod
    def _cand(**over):
        base = {
            "candidate_type": "atlas-scan:co-located-pathway",
            "signature": {"system": "glycan-locus", "class": "co-located-pathway"},
            "evidence": {"note": "six glycosyltransferases and a Wzy polymerase"},
            "subject_refs": {"protein_ids": ["p1", "p2"]},
            "reason_codes": ["colocated"],
            "uncertainty": {},
            "reduction_features": {},
        }
        base.update(over)
        return base

    def test_same_signature_different_locus_gets_distinct_keys(self):
        a = self._cand()
        b = self._cand(
            evidence={"note": "a separate locus on another contig"},
            subject_refs={"protein_ids": ["p9", "p10"]},
        )
        assert a["signature"] == b["signature"], "precondition: signatures collapse"
        assert self._key(a) != self._key(b)

    def test_reworded_evidence_for_the_same_locus_keeps_one_key(self):
        """Replay after a resume rewrites prose; identity must not move.

        Keying on evidence made a replayed frame insert a duplicate row, which
        then inflated the persisted count that finalize must match -- so the
        genome could never complete.
        """
        a = self._cand(evidence={"note": "six glycosyltransferases and a Wzy polymerase"})
        b = self._cand(evidence={"note": "a Wzy-dependent locus with six GTs"})
        assert self._key(a) == self._key(b)

    def test_reason_code_drift_does_not_move_identity(self):
        a = self._cand(reason_codes=["colocated"])
        b = self._cand(reason_codes=["COLOCALIZED_CARBOHYDRATE_ENZYMES", "UNVERIFIED"])
        assert self._key(a) == self._key(b)

    def test_identical_replay_collapses_onto_one_key(self):
        # A resume replays the frame; byte-identical output must be idempotent.
        assert self._key(self._cand()) == self._key(self._cand())

    def test_distinct_frames_stay_distinct(self):
        c = self._cand()
        assert self._key(c, frame="0") != self._key(c, frame="1")

    def test_distinct_tasks_stay_distinct(self):
        c = self._cand()
        assert self._key(c, task="t1") != self._key(c, task="t2")


class TestConflictDisambiguation:
    """A 409 means either a dead lease or a rejected record. Only the first is fatal."""

    @staticmethod
    def _exc(detail):
        import requests

        resp = requests.Response()
        resp.status_code = 409
        resp._content = json.dumps({"detail": detail}).encode()
        resp.headers["Content-Type"] = "application/json"
        return requests.HTTPError(response=resp)

    def test_lease_message_is_fatal(self):
        exc = self._exc("Candidate task has no active lease owned by the submitting agent")
        assert _is_lease_failure(exc) is True

    def test_fence_message_is_fatal(self):
        exc = self._exc("Task abc has no live attempt 3 lease owned by terra-scan-02")
        assert _is_lease_failure(exc) is True

    def test_record_rejection_is_not_fatal(self):
        exc = self._exc(
            "Candidate-occurrence idempotency key conflicts with an existing record"
        )
        assert _is_lease_failure(exc) is False

    def test_detail_is_extracted_for_logging(self):
        exc = self._exc("some rejection reason")
        assert "some rejection reason" in _conflict_detail(exc)


class TestCheckpointResumeGuard:
    """A checkpoint whose frame list disagrees with its cursor must be discarded.

    Regression for genomes that died at finalize with "Coverage repeats
    frame_id". The frame list is restored and appended to on every attempt but
    only reaches the store on a successful checkpoint write, so an attempt that
    walked frames and then died left the persisted cursor behind the persisted
    list. The next resume re-walked frames the list already held, appended them
    a second time, and the manifest builder rejected the duplicate -- after the
    genome had been paid for in full. Observed live: "resuming from frame 9"
    followed by frames 3, 4, 5.
    """

    @staticmethod
    def _frames(n, *, start=0):
        return [{"frame_id": f"genome-frame-{i:04d}", "frame_index": i} for i in range(start, start + n)]

    def _resume_with(self, monkeypatch, *, frames, cursor):
        from sharur.workers.atlas_scan import AtlasScanWorker

        worker = AtlasScanWorker.__new__(AtlasScanWorker)

        class _Ops:
            def get_task_checkpoint(self, task_id, key):
                return {"cursor": cursor, "payload": {"frames": frames}}

        worker.ops = _Ops()
        return worker._resume("t1", "atlas_progress")

    def test_healthy_checkpoint_resumes(self, monkeypatch):
        cur, frames, cands, state = self._resume_with(
            monkeypatch, frames=self._frames(3), cursor="abc"
        )
        assert cur == "abc" and len(frames) == 3

    def test_duplicate_frame_id_discards(self, monkeypatch):
        bad = self._frames(3) + [{"frame_id": "genome-frame-0000", "frame_index": 3}]
        cur, frames, _, _ = self._resume_with(monkeypatch, frames=bad, cursor="abc")
        assert cur is None and frames == []

    def test_missing_cursor_with_frames_discards(self, monkeypatch):
        cur, frames, _, _ = self._resume_with(
            monkeypatch, frames=self._frames(3), cursor=None
        )
        assert cur is None and frames == []

    def test_non_contiguous_indices_discard(self, monkeypatch):
        # The live failure: list claims 9 frames, cursor sits at frame 3.
        cur, frames, _, _ = self._resume_with(
            monkeypatch, frames=self._frames(3, start=6), cursor="abc"
        )
        assert cur is None and frames == []

    def test_empty_checkpoint_is_not_a_discard(self, monkeypatch):
        cur, frames, _, _ = self._resume_with(monkeypatch, frames=[], cursor=None)
        assert cur is None and frames == []


class TestCursorAlignmentGuard:
    """Frame-list contiguity does not prove the cursor resumes where the list ends."""

    @staticmethod
    def _cursor(frame_index):
        import base64

        raw = json.dumps({"frame_index": frame_index, "after_contig_id": "c"}).encode()
        return base64.urlsafe_b64encode(raw).decode().rstrip("=")

    def _resume_with(self, frames, cursor):
        from sharur.workers.atlas_scan import AtlasScanWorker

        worker = AtlasScanWorker.__new__(AtlasScanWorker)

        class _Ops:
            def get_task_checkpoint(self, task_id, key):
                return {"cursor": cursor, "payload": {"frames": frames}}

        worker.ops = _Ops()
        return worker._resume("t1", "atlas_progress")

    @staticmethod
    def _frames(n):
        return [
            {"frame_id": f"genome-frame-{i:04d}", "frame_index": i} for i in range(n)
        ]

    def test_contiguous_list_with_lagging_cursor_is_discarded(self):
        """Receipts 0..8 are contiguous while the cursor still points at 3."""
        cur, frames, _, _ = self._resume_with(self._frames(9), self._cursor(3))
        assert cur is None and frames == []

    def test_aligned_cursor_resumes(self):
        cur, frames, _, _ = self._resume_with(self._frames(9), self._cursor(9))
        assert cur is not None and len(frames) == 9

    def test_undecodable_cursor_is_not_treated_as_mismatch(self):
        cur, frames, _, _ = self._resume_with(self._frames(3), "not-base64-at-all")
        assert cur == "not-base64-at-all" and len(frames) == 3


class TestAnthropicStreamParsing:
    """The driver must read the terminal `result` event, not guess at stdout."""

    def test_extracts_result_and_usage(self):
        stream = "\n".join([
            json.dumps({"type": "system", "subtype": "init"}),
            json.dumps({"type": "stream_event", "event": {"type": "content_block_delta"}}),
            json.dumps({
                "type": "result",
                "subtype": "success",
                "is_error": False,
                "result": '{"candidates": [], "notes": "ok"}',
                "usage": {"input_tokens": 10, "output_tokens": 3},
            }),
        ])
        usage, body = model_cli._parse_anthropic_stream(stream)
        assert usage["input_tokens"] == 10
        assert json.loads(body)["notes"] == "ok"

    def test_error_result_raises_rather_than_returning_empty(self):
        stream = json.dumps({
            "type": "result",
            "subtype": "error_during_execution",
            "is_error": True,
            "result": "tool use blocked",
        })
        with pytest.raises(model_cli.ModelError, match="turn failed"):
            model_cli._parse_anthropic_stream(stream)

    def test_falls_back_to_raw_stdout_when_no_result_event(self):
        usage, body = model_cli._parse_anthropic_stream('{"candidates": []}')
        assert usage == {}
        assert json.loads(body) == {"candidates": []}

    def test_ignores_non_json_noise(self):
        stream = "warning: something\n" + json.dumps({
            "type": "result", "subtype": "success", "is_error": False,
            "result": '{"candidates": []}', "usage": {},
        })
        _, body = model_cli._parse_anthropic_stream(stream)
        assert json.loads(body) == {"candidates": []}


class TestCheckpointSpill:
    """Bulk resume state must not ride inline in the checkpoint payload.

    The store caps inline JSON at 256 KB. Frames and contig_state both grow
    with contig count, so a fragmented assembly (215 contigs per frame) crossed
    the cap at frame 4: the checkpoint was refused with a 409 that the worker
    read as lease loss, and the genome failed identically on every retry.
    """

    def test_spill_roundtrip(self, tmp_path):
        from sharur.workers.atlas_scan import _read_spill, _write_spill

        frames = [{"frame_id": f"f{i}", "frame_index": i} for i in range(5)]
        state = {"c1": {"contig_id": "c1", "protein_count": 7}}
        p = tmp_path / "spill" / "t1.json"
        _write_spill(p, frames, state)
        got_frames, got_state = _read_spill(p)
        assert got_frames == frames and got_state == state

    def test_missing_spill_reads_as_empty(self, tmp_path):
        from sharur.workers.atlas_scan import _read_spill

        assert _read_spill(tmp_path / "absent.json") == ([], {})

    def test_corrupt_spill_reads_as_empty(self, tmp_path):
        from sharur.workers.atlas_scan import _read_spill

        p = tmp_path / "bad.json"
        p.write_text("{truncated")
        assert _read_spill(p) == ([], {})

    def test_write_is_atomic_leaving_no_temp(self, tmp_path):
        from sharur.workers.atlas_scan import _write_spill

        p = tmp_path / "t.json"
        _write_spill(p, [], {})
        assert p.exists()
        assert not list(tmp_path.glob("*.tmp"))

    def test_path_is_none_without_manifest(self):
        from sharur.workers.atlas_scan import _spill_path

        assert _spill_path({}, "t1") is None

    def test_path_sits_beside_the_manifest_tree(self):
        from sharur.workers.atlas_scan import _spill_path

        p = _spill_path({"coverage_manifest": "/data/atlas/coverage/x.json"}, "t9")
        assert p == Path("/data/atlas/spill/t9.json")


class TestRateLimitPrecision:
    """`429` as a bare substring turned DNS blips into rate-limit backoffs.

    Rate-limit classification runs before transient, and the resulting backoff
    doubles and persists across tasks (60s -> 1800s). On a host where DNS
    failures are endemic, one incidental "429" in a verbose stderr -- a request
    id, a token count, a byte offset -- could idle a worker for half an hour
    over a hiccup it should have retried in seconds.
    """

    @pytest.mark.parametrize(
        "blob",
        [
            "req_id=a429f2 failed to lookup address information: Try again",
            "tokens_used=4290\nconnection reset by peer",
            "wrote 1429 bytes then the stream closed",
        ],
    )
    def test_incidental_429_is_not_a_rate_limit(self, blob):
        assert model_cli._looks_rate_limited(blob, "") is False
        assert model_cli._looks_transient(blob, "") is True

    @pytest.mark.parametrize(
        "blob",
        [
            "HTTP status 429 returned by endpoint",
            "429 too many requests",
            "error code: 429",
            "Error: rate limit exceeded",
            "you have hit your usage limit",
            "RESOURCE_EXHAUSTED",
        ],
    )
    def test_genuine_rate_limits_still_classify(self, blob):
        assert model_cli._looks_rate_limited(blob, "") is True

    def test_dns_failure_raises_transient_not_rate_limited(self):
        blob = "req_id=a429f2 failed to connect to websocket: failed to lookup address information"
        with pytest.raises(model_cli.ModelTransient):
            model_cli._classify_failure(1, blob, "")
