"""Claim-boundary enforcement and compact evidence bundles."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys

from sharur.core.claim_compiler import ClaimValidationError
from sharur.operators import Sharur


def _case(case_database):
    return Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
        window=2,
    )


def test_observed_and_caller_named_claims_require_their_actual_provenance(
    case_database,
):
    case = _case(case_database)
    observed = case.propose_finding(
        title="Observed domain",
        category="test",
        description="A raw domain observation.",
        claims=[
            {
                "text": "The flanking protein carries PFTEST.",
                "level": "observed",
                "evidence_refs": [
                    "annotation:fg_plus_p1:pfam:PFTEST",
                ],
            }
        ],
    )
    assert observed.validate_claims().valid
    verification = observed.compile()["verification"]
    assert verification[0]["expected"] == 1
    assert "annotations" in verification[0]["query"]

    wrongly_typed = case.propose_finding(
        title="Wrong boundary",
        category="test",
        description="A caller name mislabeled as a raw observation.",
        claims=[
            {
                "text": "Target is present.",
                "level": "observed",
                "evidence_refs": ["protein:fg_plus_p3"],
            }
        ],
    )
    assert not wrongly_typed.validate_claims().valid
    assert {issue.code for issue in wrongly_typed.validate_claims().errors} == {
        "named_claim_marked_observed"
    }

    caller_named = case.propose_finding(
        title="Caller output",
        category="test",
        description="A name emitted by the structured caller.",
        claims=[
            {
                "text": "The structured caller emitted Target.",
                "level": "caller_named",
                "evidence_refs": ["call:target_plus"],
            }
        ],
    )
    assert caller_named.validate_claims().valid
    assert '"defense_systems"' in caller_named.compile()["verification"][0]["query"]


def test_numeric_decomposition_and_priority_language_fail_closed(case_database):
    case = _case(case_database)
    numeric = case.propose_finding(
        title="Unsupported decomposition",
        category="test",
        description="The counts require distinct verification outputs.",
        claims=[
            {
                "text": "There are 2 positive cases among 3 controls.",
                "level": "observed",
                "evidence_refs": ["protein:fg_plus_p1"],
            }
        ],
    )
    report = numeric.validate_claims()
    assert not report.valid
    assert "numeric_decomposition_unverified" in {issue.code for issue in report.errors}

    priority = case.propose_finding(
        title="First occurrence in the clade",
        category="test",
        description="Priority language requires a literature audit.",
        claims=[
            {
                "text": "The structured caller emitted Target.",
                "level": "caller_named",
                "evidence_refs": ["call:target_plus"],
            }
        ],
    )
    assert "literature_unchecked" in {issue.code for issue in priority.validate_claims().errors}

    try:
        priority.compile()
    except ClaimValidationError:
        pass
    else:  # pragma: no cover - explicit assertion message is clearer
        raise AssertionError("strict compilation accepted unchecked priority language")

    hypothesis_only = case.propose_finding(
        title="Unverified idea",
        category="test",
        description="This belongs in the hypothesis registry until tested.",
        claims=[
            {
                "text": "A mechanism might link these genes.",
                "level": "unverified",
            }
        ],
    )
    assert "finding_without_verification" in {
        issue.code for issue in hypothesis_only.validate_claims().errors
    }


def test_numeric_values_and_summary_text_must_match_expected_outputs(
    case_database,
):
    case = _case(case_database)
    wrong_value = case.propose_finding(
        title="Unsupported count",
        category="test",
        description="A numeric value cannot borrow an unrelated existence check.",
        claims=[
            {
                "text": "There are 9 matching proteins.",
                "level": "observed",
                "evidence_refs": ["protein:fg_plus_p1"],
            }
        ],
    )
    assert "numeric_decomposition_unverified" in {
        issue.code for issue in wrong_value.validate_claims().errors
    }

    unclaimed_summary = case.propose_finding(
        title="4 supported contexts",
        category="test",
        description="The summary count is not present in a verified claim.",
        claims=[
            {
                "text": "The structured caller emitted Target.",
                "level": "caller_named",
                "evidence_refs": ["call:target_plus"],
            }
        ],
    )
    assert "summary_number_not_claimed" in {
        issue.code for issue in unclaimed_summary.validate_claims().errors
    }


def test_comparison_expected_object_can_verify_multiple_numeric_values(
    case_database,
):
    case = _case(case_database)
    comparison = case.compare_context(
        features=["pfam:PFTEST"],
        window=2,
        min_components=2,
        deduplicate_by="entity",
    )
    draft = case.propose_finding(
        title="Context association",
        category="test",
        description="An explicitly qualified cohort comparison.",
        novelty=2,
        falsification="Wrong if the caller cohort or feature mapping is invalid.",
        comparison=comparison,
        claims=[
            {
                "text": (
                    "The feature occurs in 2 of 2 foreground entities and "
                    "2 of 3 background entities."
                ),
                "level": "inferred",
                "evidence_refs": ["comparison:composite"],
            }
        ],
    )
    assert draft.validate_claims().valid


def test_bundle_is_compact_hashed_and_replays_the_comparison(
    case_database,
    tmp_path,
):
    case = _case(case_database)
    comparison = case.compare_context(
        features=["pfam:PFTEST"],
        window=2,
        min_components=2,
        deduplicate_by="replicon",
    )
    draft = case.propose_finding(
        title="Caller-linked context",
        category="test",
        description="A compact evidence-bundle fixture.",
        comparison=comparison,
        claims=[
            {
                "text": "The structured caller emitted Target.",
                "level": "caller_named",
                "evidence_refs": ["call:target_plus"],
            }
        ],
    )
    output = draft.publish_bundle(tmp_path / "bundle")

    expected_files = {
        "README.md",
        "case.json",
        "case.md",
        "comparison.json",
        "comparison_recipe.json",
        "verify_comparison.py",
        "finding.json",
        "verification.sql",
        "components.faa",
        "manifest.json",
    }
    assert {path.name for path in output.iterdir()} == expected_files
    assert not any(path.suffix == ".duckdb" for path in output.iterdir())

    manifest = json.loads((output / "manifest.json").read_text())
    for entry in manifest["files"]:
        digest = hashlib.sha256((output / entry["path"]).read_bytes()).hexdigest()
        assert digest == entry["sha256"]

    replay = subprocess.run(
        [
            sys.executable,
            str(output / "verify_comparison.py"),
            "--db",
            str(case_database),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    replayed = json.loads(replay.stdout)
    assert replayed["foreground_total"] == 2
    assert replayed["background_total"] == 2
