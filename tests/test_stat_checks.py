"""Statistical checks and pre-registration records (sharur/stat_checks.py, sharur/preregistration.py).

The worked example is a synthetic study, "study-x" (tests/data/prereg_worked_example.json, made by
tests/data/make_prereg_worked_example.py; every label and number is invented). It pre-registers ten
predictions (H1-H10), a Benjamini-Hochberg family over them and the verdict rule "q < 0.05, CI
excludes 0 in the predicted direction, placebo p <= 0.1", with a calibration rule set on H1-H3.
Its recorded q values come from statsmodels and its verdicts from plain decision logic, so the tests
check Sharur's rule engine and BH/Holm against independent implementations; every verdict rule fires
at least once.
"""

from __future__ import annotations

import copy
import hashlib
import json
import stat
import sys
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from sharur.preregistration import (
    PreregError,
    Preregistration,
    load_preregistration,
    sha256_file,
    verify_preregistration,
    write_preregistration,
)
from sharur.stat_checks import (
    RuleError,
    StatisticalClaim,
    apply_rules,
    bh,
    check_claim,
    compare_with_tolerance,
    evaluate_rule,
    holm,
    statistical_check,
)


FIXTURE = json.loads((Path(__file__).parent / "data" / "prereg_worked_example.json").read_text())
CASES = {c["case"]: c for c in FIXTURE["cases"]}
COMPUTED_AT = "2026-01-15T14:00:00Z"  # after the pre-registration's written_at
PLACEBO_FAILS = ("NOT SUPPORTED (Delta passes; placebo control fails: {placebos_at_or_beyond} of 10 "
                 "placebo means at or beyond the target)")
RAW_ONLY = "NOT SUPPORTED (raw change at the target sites, matched controls change as much)"


def verdict_rules(direction: str) -> list[dict]:
    """The pre-registered verdict rule, written as ordered expressions."""
    if direction == ">":
        return [{"label": "SUPPORTED", "when": "q_BH < 0.05 and lo > 0 and p_placebo <= 0.1"},
                {"label": PLACEBO_FAILS, "when": "q_BH < 0.05 and lo > 0"},
                {"label": RAW_ONLY, "when": "p_raw < 0.05"},
                {"label": "NOT SUPPORTED", "when": "true"}]
    return [{"label": "SUPPORTED", "when": "q_BH < 0.05 and hi < 0 and p_placebo <= 0.1"},
            {"label": PLACEBO_FAILS, "when": "q_BH < 0.05 and hi < 0"},
            {"label": "NOT SUPPORTED", "when": "true"}]


CALIBRATION = [{"label": "RECOVERED", "when": "q_BH < 0.05 and lo > 0 and p_placebo <= 0.1"},
               {"label": "DIRECTION ONLY", "when": "mean > 0 and n_pred >= 7"},
               {"label": "NOT RECOVERED", "when": "true"}]


def study_prereg() -> Preregistration:
    preds = []
    for cid, c in CASES.items():
        preds.append({
            "id": cid, "text": c["label"],
            "test": ("one-sided Wilcoxon signed-rank on per-contrast Delta (target minus matched control "
                     "sites, loss minus sister); mean Delta with a 10,000-resample bootstrap 95% CI"),
            "statistic": "mean_delta", "direction": "greater" if c["direction"] == ">" else "less",
            "verdict_rules": verdict_rules(c["direction"]),
            "extra_rules": {"calibration": CALIBRATION} if c["calibration"] else {}})
    return Preregistration.model_validate({
        "id": "study-x", "title": "study-x: contact-site departures after factor loss (synthetic)",
        "dataset": FIXTURE["dataset"],
        "document": {"path": "analyses/study-x/" + FIXTURE["document"]["name"],
                     "sha256": FIXTURE["document"]["sha256"], "bytes": FIXTURE["document"]["bytes"]},
        "design_inputs": ["presence calls", "species tree", "loss events and sister sets", "placebo factors"],
        "predictions": preds,
        "correction": {"method": "BH", "family": list(CASES), "alpha": 0.05}})


def claim_for(cid: str, prereg: Preregistration, **overrides) -> StatisticalClaim:
    c = CASES[cid]
    label, _, calibration = c["verdict"].partition("; calibration ")
    pred = prereg.prediction(cid)
    values = {k: c[k] for k in ("mean", "lo", "hi", "q_BH", "p_raw", "p_placebo", "placebos_at_or_beyond",
                                "n_pred", "p_random")}
    data = {
        "prediction": cid,
        "preregistration": {"id": prereg.id, "content_sha256": prereg.content_sha256},
        "dataset_id": prereg.dataset.dataset_id,
        "statistic": {"name": "W+", "test": "Wilcoxon signed-rank", "p": c["p"], "alternative": pred.direction},
        "effect": {"name": "mean_delta", "value": c["mean"], "interval": [c["lo"], c["hi"]],
                   "method": "bootstrap over contrasts, 10,000 resamples"},
        "n_events": c["events"], "event_unit": "internal loss edge", "n_used": c["contrasts"],
        "multiple_testing": {"method": "BH", "family": "H1-H10", "family_size": 10, "q": c["q_BH"]},
        "controls": [
            {"name": "placebos", "kind": "placebo", "n": 10, "observed": c["placebo_max"],
             "at_or_beyond": c["placebos_at_or_beyond"], "p_empirical": c["p_placebo"],
             "rule": "p_placebo <= 0.1", "passed": c["p_placebo"] <= 0.1},
            {"name": "random edges", "kind": "random_edge", "n": 200, "p_empirical": c["p_random"]},
        ],
        "values": values,
        "verdict": {"label": label, "rules": [r.model_dump() for r in pred.verdict_rules],
                    "extra": {"calibration": calibration} if calibration else {},
                    "extra_rules": {k: [r.model_dump() for r in v] for k, v in pred.extra_rules.items()}},
        "computed_at": COMPUTED_AT,
    }
    data.update(overrides)
    return StatisticalClaim.model_validate(data)


@pytest.fixture
def prereg(tmp_path: Path) -> tuple[Path, Preregistration]:
    return write_preregistration(study_prereg(), tmp_path / "prereg", written_at=FIXTURE["written_at"])


# --------------------------------------------------------------------------- the worked example
def test_worked_example_verdicts_follow_from_the_preregistered_rules(prereg):
    _, record = prereg
    supported, texts = set(), set()
    for cid, case in CASES.items():
        claim = claim_for(cid, record)
        assert claim.evaluate().text() == case["verdict"], cid
        assert check_claim(claim, prereg=record) == [], cid
        texts.add(case["verdict"])
        if claim.verdict.label == "SUPPORTED":
            supported.add(cid)
    assert supported == {"H2", "H4", "H8"}
    # every rule fires: both placebo-failure labels (with their counts), the raw-change label, the
    # plain negative in both directions, and all three calibration outcomes
    assert claim_for("H5", record).evaluate().label == PLACEBO_FAILS.format(placebos_at_or_beyond=2)
    assert claim_for("H9", record).evaluate().label == PLACEBO_FAILS.format(placebos_at_or_beyond=1)
    assert {RAW_ONLY, "SUPPORTED; calibration RECOVERED", "NOT SUPPORTED; calibration DIRECTION ONLY",
            "NOT SUPPORTED; calibration NOT RECOVERED"} <= texts
    assert [c for c, case in CASES.items() if case["verdict"] == "NOT SUPPORTED"] == ["H7", "H10"]


def test_recorded_q_values_recompute_within_tolerance():
    order = list(CASES)
    q = bh([CASES[c]["p"] for c in order])
    h = holm([CASES[c]["p"] for c in order])
    for cid, qq, hh in zip(order, q, h, strict=True):
        assert compare_with_tolerance(qq, CASES[cid]["q_BH"], {"*": {"rel": 0.001}}) == [], cid
        assert compare_with_tolerance(hh, CASES[cid]["p_Holm"], {"*": {"rel": 0.001}}) == [], cid


def test_preregistration_is_frozen_hashed_and_tamper_evident(prereg, tmp_path):
    path, record = prereg
    assert path.name == "study-x_20260115T100000Z.prereg.json"
    assert not path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    loaded = load_preregistration(path)
    assert loaded.content_sha256 == record.content_sha256 == loaded.compute_hash()
    assert verify_preregistration(path, check_document=False).valid
    missing_doc = verify_preregistration(path)
    assert not missing_doc.valid and "document not found" in missing_doc.issues[0]
    with pytest.raises(PreregError, match="never replaced"):
        write_preregistration(study_prereg(), path.parent, written_at=FIXTURE["written_at"])
    with pytest.raises(PreregError, match="already carries"):
        write_preregistration(loaded, tmp_path / "again")
    edited = json.loads(path.read_text())
    edited["predictions"][3]["verdict_rules"][0]["when"] = "q_BH < 0.1 and lo > 0"
    other = tmp_path / "edited.prereg.json"
    other.write_text(json.dumps(edited))
    result = verify_preregistration(other, check_document=False)
    assert not result.valid and "edited after it was written" in result.issues[0]


def test_document_hash_binds_the_full_text(tmp_path):
    doc = tmp_path / FIXTURE["document"]["name"]
    doc.write_text(FIXTURE["document"]["text"])
    assert sha256_file(doc) == FIXTURE["document"]["sha256"] == hashlib.sha256(doc.read_bytes()).hexdigest()
    data = study_prereg().payload()
    data["document"] = {"path": str(doc), "sha256": sha256_file(doc)}
    path, _ = write_preregistration(Preregistration.model_validate(data), tmp_path / "p")
    assert verify_preregistration(path).valid
    doc.write_text(FIXTURE["document"]["text"] + "Amendment after the analysis.\n")
    assert "changed since it was registered" in verify_preregistration(path).issues[0]


def test_check_claim_names_each_inconsistency(prereg, tmp_path):
    _, record = prereg
    base = claim_for("H4", record)
    bad_ci = base.model_copy(update={"effect": base.effect.model_copy(update={"interval": (0.2, 0.3)})})
    assert any("outside its interval" in i for i in check_claim(bad_ci, prereg=record))
    wrong_label = base.model_copy(update={"verdict": base.verdict.model_copy(update={"label": "NOT SUPPORTED"})})
    assert any("the rules give 'SUPPORTED'" in i for i in check_claim(wrong_label, prereg=record))
    early = base.model_copy(update={"computed_at": "2026-01-01T00:00:00Z"})
    assert any("not after the pre-registration" in i for i in check_claim(early, prereg=record))
    other_rules = base.model_copy(update={"verdict": base.verdict.model_copy(
        update={"rules": base.verdict.rules[:1]})})
    assert any("differ from those pre-registered" in i for i in check_claim(other_rules, prereg=record))
    other_ds = base.model_copy(update={"dataset_id": "0" * 64})
    assert any("differs from the pre-registered" in i for i in check_claim(other_ds, prereg=record))
    assert any("dataset seal now names" in i for i in check_claim(base, prereg=record, seal_dataset_id="1" * 64))
    lied = copy.deepcopy(base.model_dump())
    lied["controls"][0]["passed"] = False
    assert any("recorded passed=False" in i for i in check_claim(lied, prereg=record))
    unknown = base.model_copy(update={"prediction": "H99"})
    assert any("H99 is not in" in i for i in check_claim(unknown, prereg=record))
    inp = tmp_path / "contrasts.tsv"
    inp.write_text("a\t1\n")
    with_input = StatisticalClaim.model_validate(
        {**base.model_dump(), "reproduce": {"inputs": [{"path": str(inp), "sha256": sha256_file(inp)}]}})
    assert check_claim(with_input, prereg=record) == []
    inp.write_text("a\t2\n")
    assert any("input changed" in i for i in check_claim(with_input, prereg=record))


def test_rule_language_is_small_and_safe():
    assert evaluate_rule("q < 0.05 and (lo > 0 or hi < 0)", {"q": 0.01, "lo": -1, "hi": -0.2})
    assert evaluate_rule("abs(mean - 0.1) <= 0.03 and not flag", {"mean": 0.12, "flag": False})
    assert not evaluate_rule("p_raw < 0.05", {"p_raw": None})
    assert apply_rules([{"label": "{k} of 10", "when": "k > 0"}, {"label": "none", "when": "true"}],
                       {"k": 3.0}) == "3 of 10"
    for bad in ("__import__('os')", "x.real > 0", "[1][0] == 1", "q <"):
        with pytest.raises(RuleError):
            evaluate_rule(bad, {"x": 1, "q": 1})
    with pytest.raises(RuleError, match="unknown value"):
        evaluate_rule("missing > 0", {})


def test_tolerance_comparison():
    expected = {"mean": 0.118, "q_BH": 0.00525, "n_used": 26}
    tol = {"mean": {"abs": 0.005}, "*": {"rel": 0.05}}
    assert compare_with_tolerance({"mean": 0.1221, "q_BH": 0.0054, "n_used": 26, "extra": 1}, expected, tol) == []
    misses = compare_with_tolerance({"mean": 0.131, "q_BH": 0.0060, "n_used": 25}, expected, tol)
    assert [m.split()[0] for m in misses] == ["mean", "q_BH", "n_used"]
    assert compare_with_tolerance([[0.5, 3]], [0.51, 3], {"*": {"abs": 0.02}}) == []
    assert "q_BH missing" in compare_with_tolerance({"mean": 0.118}, expected, tol)


# --------------------------------------------------------------------------- re-verification
def _repo(tmp_path: Path, dataset_id: str) -> Path:
    root = tmp_path / "repo"
    ds = root / "data" / "ds"
    (ds / "exploration").mkdir(parents=True)
    con = duckdb.connect(str(ds / "sharur.duckdb"))
    con.execute("CREATE TABLE bins (bin_id VARCHAR)")
    con.close()
    (ds / "dataset.seal.json").write_text(json.dumps({"dataset_id": dataset_id}))
    return root


def _replay_values(cid: str, **changes) -> dict:
    c = CASES[cid]
    values = {k: c[k] for k in ("mean", "lo", "hi", "q_BH", "p_raw", "p_placebo", "placebos_at_or_beyond")}
    values.update(changes)
    return values


def test_reverify_recomputes_statistical_checks(tmp_path, prereg):
    from sharur.reverify import reverify

    _, record = prereg
    root = _repo(tmp_path, record.dataset.dataset_id)
    ds = root / "data" / "ds"
    prereg_copy = root / "prereg.json"
    prereg_copy.write_text(Path(prereg[0]).read_text())
    claim = claim_for("H4", record, preregistration={"id": record.id, "content_sha256": record.content_sha256,
                                                      "path": "prereg.json"})
    tol = {"mean": {"abs": 0.005}, "*": {"rel": 0.05}}
    expected = _replay_values("H4")
    # the reproducing script reads an input table; its hash is part of the claim
    table = ds / "exploration" / "h4.json"
    table.write_text(json.dumps(_replay_values("H4", mean=0.1191)))
    script = ds / "exploration" / "h4.py"
    script.write_text("import json, sys\nprint(json.dumps(json.load(open(sys.argv[1]))))\n")
    scripted = StatisticalClaim.model_validate({**claim.model_dump(), "reproduce": {
        "script": "data/ds/exploration/h4.py", "script_sha256": sha256_file(script),
        "inputs": [{"path": "data/ds/exploration/h4.json", "sha256": sha256_file(table)}]}})
    items = {
        "within": statistical_check("H4 within tolerance", f"python: {_replay_values('H4', mean=0.1212)!r}",
                                    claim, expected, tol),
        "drift": statistical_check("H4 drifts", f"python: {_replay_values('H4', mean=0.131)!r}",
                                   claim, expected, tol),
        "flip": statistical_check("H4 verdict flips", f"python: {_replay_values('H4', q_BH=0.08)!r}",
                                  claim, expected, {"*": {"rel": 100.0}}),
        "script": statistical_check("H4 from its script",
                                    f"shell: {sys.executable} data/ds/exploration/h4.py data/ds/exploration/h4.json",
                                    scripted, expected, tol),
        "count": {"claim": "count-style triple unchanged", "query": "SELECT COUNT(*) FROM bins", "expected": 0},
    }
    findings = [{"id": name, "title": name, "verification": [item]} for name, item in items.items()]
    (ds / "exploration" / "findings.jsonl").write_text("".join(json.dumps(f) + "\n" for f in findings))

    report = reverify([ds], workers=2)
    by = {r.triple.finding_id: r for r in report.results}
    assert by["within"].replay.status == "OK" and by["within"].flags == []
    assert by["drift"].replay.status == "FAIL"
    assert any(f.startswith("TOLERANCE:mean") for f in by["drift"].flags)
    assert by["flip"].replay.status == "FAIL"
    assert any(f.startswith("VERDICT_CHANGED:SUPPORTED -> NOT SUPPORTED (raw change") for f in by["flip"].flags)
    assert by["script"].replay.status == "OK"
    assert by["count"].replay.status == "OK" and by["count"].flags == []
    assert "statistical check: TOLERANCE:mean" in report.to_markdown()

    # an edited input, a new dataset seal and an edited pre-registration are all flagged
    table.write_text(json.dumps(_replay_values("H4", mean=0.2)))
    (ds / "dataset.seal.json").write_text(json.dumps({"dataset_id": "f" * 64}))
    edited = json.loads(prereg_copy.read_text())
    edited["title"] = "edited later"
    prereg_copy.chmod(0o644)
    prereg_copy.write_text(json.dumps(edited))
    again = {r.triple.finding_id: r for r in reverify([ds], ids=["script"], workers=1).results}["script"]
    assert again.replay.status == "FAIL"
    assert {"INPUT_CHANGED:data/ds/exploration/h4.json", "SEAL_CHANGED", "PREREG_MISMATCH"} <= set(again.flags)
    # a pre-registration that is gone
    prereg_copy.unlink()
    gone = {r.triple.finding_id: r for r in reverify([ds], ids=["within"], workers=1).results}["within"]
    assert "PREREG_MISSING" in gone.flags


def test_strict_writer_keeps_statistical_checks(tmp_path, prereg):
    from sharur.core.analysis_record_io import append_finding_record

    _, record = prereg
    item = statistical_check("H4 contact sites", "python: {'mean': 0.118}", claim_for("H4", record),
                             {"mean": 0.118}, {"mean": {"abs": 0.005}})
    path = tmp_path / "findings.jsonl"
    append_finding_record(path, {"title": "Loss of factor B changes contact sites on partner B (synthetic)",
                                 "category": "general",
                                 "description": "Pre-registered H4 passes its verdict rule.",
                                 "evidence": {"case": "H4"}, "verification": [item]}, phase="exploration")
    stored = json.loads(path.read_text())
    kept = stored["verification"][0]
    assert kept["check"] == "statistical" and kept["tolerance"] == {"mean": {"abs": 0.005}}
    assert StatisticalClaim.model_validate(kept["statistical_claim"]).verdict.label == "SUPPORTED"


def test_stats_cli(tmp_path):
    from sharur.cli import app

    runner = CliRunner()
    doc = tmp_path / "study_plan.txt"
    doc.write_text(FIXTURE["document"]["text"])
    spec = study_prereg().payload()
    spec.pop("document")
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    out = runner.invoke(app, ["stats", "preregister", str(tmp_path / "spec.json"), "--out-dir",
                              str(tmp_path / "prereg"), "--document", str(doc)])
    assert out.exit_code == 0, out.output
    made = json.loads(out.stdout)
    assert len(made["content_sha256"]) == 64 and made["predictions"][:2] == ["H1", "H2"]
    ok = runner.invoke(app, ["stats", "verify-prereg", made["path"]])
    assert ok.exit_code == 0 and "verified" in ok.stdout
    record = load_preregistration(made["path"])
    good = statistical_check("H8", "python: {}", claim_for("H8", record, computed_at="2099-01-01T00:00:00Z"))
    bad = copy.deepcopy(good)
    bad["statistical_claim"]["verdict"]["label"] = "NOT SUPPORTED"
    findings = tmp_path / "findings.jsonl"
    findings.write_text(json.dumps({"id": "F1", "verification": [good, bad]}) + "\n")
    res = runner.invoke(app, ["stats", "check", str(findings), "--root", str(tmp_path), "--strict"])
    assert res.exit_code == 1
    assert "2 found, 1 consistent, 1 with issues" in res.stdout and "the rules give 'SUPPORTED'" in res.stdout
