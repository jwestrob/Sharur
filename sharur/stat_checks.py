"""Statistical checks: verification records that hold a statistical claim.

A count-style verification triple is ``{claim, query, expected}`` and replays to an exact
number. A statistical check is the same triple with three more keys, so every existing
reader (normalizer, strict writer, ``sharur reverify``) accepts it unchanged::

    {
      "claim": "H4: contact sites on partner B change more where factor B was lost",
      "query": "shell: python3 analyses/study-x/outcomes.py --case H4 --json",
      "expected": {"mean": 0.118, "lo": 0.061, "hi": 0.175, "q_BH": 0.00525, "p_placebo": 0.09091},
      "check": "statistical",
      "tolerance": {"mean": {"abs": 0.005}, "*": {"rel": 0.05}},
      "statistical_claim": { ... StatisticalClaim ... }
    }

``statistical_claim`` records the statistic, the effect size with its interval, the number
of independent events, the multiple-testing correction, the control results the claim must
beat (placebos, label shuffles, random edges, matched sets: each with the rule it has to
pass), the verdict rules and recorded verdict, the pre-registration it implements
(``{id, content_sha256}``), the dataset seal, and the script, command and input hashes that
reproduce it.

Re-verification (``sharur reverify``) replays ``query`` like any triple, compares the observed
values with ``expected`` under ``tolerance`` (``{"name": {"abs": a, "rel": r}}``; ``"*"`` is
the default for float values; integers and values without a tolerance keep the precision rule
of count checks), then
re-applies the verdict rules to the observed values. The triple fails when a value leaves its
tolerance or the verdict changes, and it is flagged when an input file's hash, the dataset
seal or the pre-registration no longer match the record. ``sharur stats check`` runs the
static part (rules, intervals, controls, pre-registration linkage) without replaying.

Verdict rules are expressions over named values: numbers, comparisons, ``and``/``or``/``not``,
``+ - * /``, ``abs``/``min``/``max``, and ``true``/``false``. Comparisons with a missing (null)
value are false.
"""

from __future__ import annotations

import ast
import json
import math
import operator
from collections.abc import Iterable, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from sharur.preregistration import (
    PreregError,
    Preregistration,
    VerdictRule,
    load_preregistration,
    parse_time,
    sha256_file,
)


CHECK_KIND = "statistical"


class RuleError(ValueError):
    """A verdict rule is not a valid expression or names an unknown value."""


# --------------------------------------------------------------------------- #
# Verdict-rule expressions
# --------------------------------------------------------------------------- #
_CMP = {ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt, ast.GtE: operator.ge,
        ast.Eq: operator.eq, ast.NotEq: operator.ne}
_BIN = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
_FUNCS = {"abs": abs, "min": min, "max": max}
_CONST = {"true": True, "false": False, "True": True, "False": False, "null": None, "None": None}


def _eval(node: ast.AST, values: dict[str, Any]) -> Any:
    if isinstance(node, ast.Expression):
        return _eval(node.body, values)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, str, bool, type(None))):
        return node.value
    if isinstance(node, ast.Name):
        if node.id in _CONST:
            return _CONST[node.id]
        if node.id not in values:
            raise RuleError(f"unknown value {node.id!r}")
        return values[node.id]
    if isinstance(node, ast.BoolOp):
        results = (bool(_eval(v, values)) for v in node.values)
        return all(results) if isinstance(node.op, ast.And) else any(results)
    if isinstance(node, ast.UnaryOp):
        operand = _eval(node.operand, values)
        if isinstance(node.op, ast.Not):
            return not operand
        if operand is None:
            return None
        return -operand if isinstance(node.op, ast.USub) else +operand
    if isinstance(node, ast.BinOp) and type(node.op) in _BIN:
        left, right = _eval(node.left, values), _eval(node.right, values)
        if left is None or right is None:
            return None
        return _BIN[type(node.op)](left, right)
    if isinstance(node, ast.Compare):
        left = _eval(node.left, values)
        for op, comparator in zip(node.ops, node.comparators):
            if type(op) not in _CMP:
                raise RuleError(f"operator {type(op).__name__} is not allowed")
            right = _eval(comparator, values)
            if left is None or right is None:
                return False
            if (isinstance(left, float) and math.isnan(left)) or (isinstance(right, float) and math.isnan(right)):
                return False
            if not _CMP[type(op)](left, right):
                return False
            left = right
        return True
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in _FUNCS
            and not node.keywords):
        args = [_eval(a, values) for a in node.args]
        return None if any(a is None for a in args) else _FUNCS[node.func.id](*args)
    raise RuleError(f"{type(node).__name__} is not allowed in a verdict rule")


def evaluate_rule(expression: str, values: dict[str, Any]) -> bool:
    """Truth of one rule expression over named values."""
    try:
        tree = ast.parse(expression.strip(), mode="eval")
    except SyntaxError as exc:
        raise RuleError(f"not an expression: {expression!r}") from exc
    return bool(_eval(tree, values))


def rule_names(expression: str) -> set[str]:
    tree = ast.parse(expression.strip(), mode="eval")
    return {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and n.id not in _CONST
            and n.id not in _FUNCS}


class _FormatValues(dict):
    def __missing__(self, key: str) -> str:
        return "{" + key + "}"


def apply_rules(rules: Sequence[VerdictRule | dict[str, str]], values: dict[str, Any]) -> str | None:
    """Label of the first rule whose condition holds (labels may name values in braces)."""
    for rule in rules:
        rule = rule if isinstance(rule, VerdictRule) else VerdictRule.model_validate(rule)
        if evaluate_rule(rule.when, values):
            shown = {k: _plain_number(v) for k, v in values.items()}
            return rule.label.format_map(_FormatValues(shown))
    return None


def _plain_number(value: Any) -> Any:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


# --------------------------------------------------------------------------- #
# Multiple testing (for recomputing recorded q values)
# --------------------------------------------------------------------------- #
def bh(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p values (step-up, monotone, capped at 1)."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    out = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * n / rank)
        out[i] = min(1.0, running)
    return out


def holm(pvalues: Sequence[float]) -> list[float]:
    """Holm step-down adjusted p values."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    out = [1.0] * n
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, pvalues[i] * (n - rank)))
        out[i] = running
    return out


# --------------------------------------------------------------------------- #
# The claim
# --------------------------------------------------------------------------- #
Number = Optional[float]


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Statistic(_Model):
    name: str
    test: str
    value: Number = None
    p: Number = Field(default=None, ge=0, le=1)
    alternative: Literal["greater", "less", "two-sided"] | None = None


class Effect(_Model):
    name: str
    value: float
    interval: tuple[float, float] | None = None
    level: float = Field(default=0.95, gt=0, lt=1)
    method: str | None = None


class MultipleTesting(_Model):
    method: Literal["BH", "Holm", "Bonferroni", "none"]
    family: str | None = None
    family_size: int | None = Field(default=None, ge=1)
    q: Number = Field(default=None, ge=0, le=1)


class ControlResult(_Model):
    """One control the claim must beat: what was run, what it gave, and the rule it must pass."""

    name: str
    kind: Literal["placebo", "shuffle", "random_edge", "matched", "planted", "calibration", "other"]
    n: int | None = Field(default=None, ge=0, description="Placebos, permutations or draws run.")
    observed: Number = Field(default=None, description="The control's own value (e.g. the largest placebo mean).")
    at_or_beyond: int | None = Field(default=None, ge=0, description="Controls at or beyond the target.")
    p_empirical: Number = Field(default=None, ge=0, le=1)
    rule: str | None = Field(default=None, description="Expression the claim must satisfy; None = reported only.")
    passed: bool | None = None
    note: str | None = None


class InputFile(_Model):
    path: str
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class Reproduce(_Model):
    command: str | None = Field(default=None, description="Replay form; the triple's `query` when omitted.")
    script: str | None = None
    script_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    inputs: list[InputFile] = Field(default_factory=list)
    outputs: list[str] = Field(default_factory=list)


class PreregRef(_Model):
    id: str
    content_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    path: str | None = None


class Verdict(_Model):
    label: str
    rules: list[VerdictRule] = Field(default_factory=list)
    extra: dict[str, str] = Field(default_factory=dict, description="Labels of further named rule sets.")
    extra_rules: dict[str, list[VerdictRule]] = Field(default_factory=dict)

    def text(self) -> str:
        return "; ".join([self.label, *[f"{k} {v}" for k, v in self.extra.items()]])


class StatisticalClaim(_Model):
    prediction: str | None = Field(default=None, description="Prediction ID in the pre-registration.")
    preregistration: PreregRef | None = None
    dataset_id: str | None = None
    statistic: Statistic
    effect: Effect | None = None
    n_events: int = Field(ge=0, description="Independent events (losses, contrasts, clades) behind the test.")
    event_unit: str = "event"
    n_used: int | None = Field(default=None, ge=0, description="Events with data on both sides.")
    multiple_testing: MultipleTesting | None = None
    controls: list[ControlResult] = Field(default_factory=list)
    values: dict[str, Any] = Field(default_factory=dict, description="Named values the rules read.")
    verdict: Verdict | None = None
    reproduce: Reproduce = Field(default_factory=Reproduce)
    tolerance: dict[str, dict[str, float]] = Field(default_factory=dict)
    computed_at: str | None = None
    notes: str | None = None

    def namespace(self, observed: dict[str, Any] | None = None) -> dict[str, Any]:
        """Values the rules read: ``values``, then fields of the claim, then ``observed``."""
        ns: dict[str, Any] = {"n_events": self.n_events, "n_used": self.n_used,
                              "statistic": self.statistic.value, "p": self.statistic.p}
        if self.effect is not None:
            ns["effect"] = self.effect.value
            if self.effect.interval is not None:
                ns["effect_lo"], ns["effect_hi"] = self.effect.interval
        if self.multiple_testing is not None:
            ns["q"] = self.multiple_testing.q
        for control in self.controls:
            key = "".join(ch if ch.isalnum() else "_" for ch in control.name.lower())
            ns[f"{key}_p"] = control.p_empirical
            ns[f"{key}_at_or_beyond"] = control.at_or_beyond
        ns.update(self.values)
        if observed:
            ns.update(observed)
        return ns

    def evaluate(self, observed: dict[str, Any] | None = None) -> Verdict | None:
        """The verdict the recorded rules give for these values (observed values override)."""
        if self.verdict is None or not self.verdict.rules:
            return None
        ns = self.namespace(observed)
        label = apply_rules(self.verdict.rules, ns) or "NO RULE MATCHED"
        extra = {name: apply_rules(rules, ns) or "NO RULE MATCHED"
                 for name, rules in self.verdict.extra_rules.items()}
        return self.verdict.model_copy(update={"label": label, "extra": extra})


def statistical_check(claim_text: str, query: str, claim: StatisticalClaim | dict[str, Any],
                      expected: dict[str, Any] | None = None,
                      tolerance: dict[str, dict[str, float]] | None = None) -> dict[str, Any]:
    """A verification record for ``findings.jsonl`` holding a statistical claim.

    ``expected`` defaults to the claim's ``values``; ``tolerance`` to the claim's.
    """
    model = claim if isinstance(claim, StatisticalClaim) else StatisticalClaim.model_validate(claim)
    return {"claim": claim_text, "query": query,
            "expected": expected if expected is not None else dict(model.values),
            "check": CHECK_KIND, "tolerance": tolerance if tolerance is not None else dict(model.tolerance),
            "statistical_claim": model.model_dump(mode="json", exclude_none=True)}


# --------------------------------------------------------------------------- #
# Static consistency
# --------------------------------------------------------------------------- #
def _resolve(path: str, root: Path | None) -> Path:
    p = Path(path).expanduser()
    return p if p.is_absolute() or root is None else root / p


def check_claim(claim: StatisticalClaim | dict[str, Any], *, root: Path | None = None,
                prereg: Preregistration | None = None, seal_dataset_id: str | None = None) -> list[str]:
    """Problems with a claim as recorded (an empty list means it is consistent).

    Checks the interval, the controls against their rules, the recorded verdict against its
    rules, the pre-registration (hash, timing, prediction, rules, dataset), the dataset seal,
    and the hashes of the script and input files.
    """
    c = claim if isinstance(claim, StatisticalClaim) else StatisticalClaim.model_validate(claim)
    issues: list[str] = []
    ns = c.namespace()
    if c.effect is not None and c.effect.interval is not None:
        lo, hi = c.effect.interval
        if not lo <= c.effect.value <= hi:
            issues.append(f"effect {c.effect.value} lies outside its interval [{lo}, {hi}]")
    if c.n_used is not None and c.n_used > c.n_events:
        issues.append(f"n_used {c.n_used} exceeds n_events {c.n_events}")
    for control in c.controls:
        if control.rule is None:
            continue
        try:
            ok = evaluate_rule(control.rule, ns)
        except RuleError as exc:
            issues.append(f"control {control.name}: {exc}")
            continue
        if control.passed is not None and control.passed != ok:
            issues.append(f"control {control.name}: recorded passed={control.passed}, rule {control.rule!r} gives {ok}")
    if c.verdict is not None and c.verdict.rules:
        try:
            got = c.evaluate()
        except RuleError as exc:
            issues.append(f"verdict rules: {exc}")
            got = None
        if got is not None and got.text() != c.verdict.text():
            issues.append(f"recorded verdict {c.verdict.text()!r}; the rules give {got.text()!r}")
    if c.preregistration is not None:
        if prereg is None and c.preregistration.path:
            try:
                prereg = load_preregistration(_resolve(c.preregistration.path, root))
            except PreregError as exc:
                issues.append(f"pre-registration: {exc}")
        if prereg is not None:
            issues += check_against_preregistration(c, prereg)
    if seal_dataset_id is not None and c.dataset_id is not None and seal_dataset_id != c.dataset_id:
        issues.append(f"dataset seal now names {seal_dataset_id}; the claim was computed on {c.dataset_id}")
    issues += check_inputs(c, root)
    return issues


def check_inputs(claim: StatisticalClaim, root: Path | None = None) -> list[str]:
    issues = []
    files = list(claim.reproduce.inputs)
    if claim.reproduce.script and claim.reproduce.script_sha256:
        files.append(InputFile(path=claim.reproduce.script, sha256=claim.reproduce.script_sha256))
    for item in files:
        p = _resolve(item.path, root)
        if not p.is_file():
            issues.append(f"input missing: {item.path}")
        elif sha256_file(p) != item.sha256:
            issues.append(f"input changed: {item.path}")
    return issues


def check_against_preregistration(claim: StatisticalClaim, prereg: Preregistration) -> list[str]:
    """The claim implements the pre-registration it cites, unchanged and written before it."""
    issues = []
    ref = claim.preregistration
    recomputed = prereg.compute_hash()
    if prereg.content_sha256 != recomputed:
        issues.append("pre-registration was edited after it was written (content_sha256 mismatch)")
    if ref is not None and ref.content_sha256 != recomputed:
        issues.append("claim cites a different pre-registration hash than the record holds")
    if ref is not None and ref.id != prereg.id:
        issues.append(f"claim cites pre-registration {ref.id}, record is {prereg.id}")
    if prereg.written_at and claim.computed_at:
        if parse_time(claim.computed_at) <= parse_time(prereg.written_at):
            issues.append(f"claim computed at {claim.computed_at}, not after the pre-registration "
                          f"({prereg.written_at})")
    elif not claim.computed_at:
        issues.append("claim has no computed_at, so its order relative to the pre-registration is unknown")
    if claim.dataset_id and claim.dataset_id != prereg.dataset.dataset_id:
        issues.append(f"claim dataset {claim.dataset_id} differs from the pre-registered {prereg.dataset.dataset_id}")
    if claim.prediction is None:
        issues.append("claim cites a pre-registration but names no prediction")
        return issues
    try:
        pred = prereg.prediction(claim.prediction)
    except KeyError:
        issues.append(f"prediction {claim.prediction} is not in pre-registration {prereg.id}")
        return issues
    if claim.verdict is not None:
        if [r.model_dump() for r in claim.verdict.rules] != [r.model_dump() for r in pred.verdict_rules]:
            issues.append(f"verdict rules differ from those pre-registered for {pred.id}")
        mine = {k: [r.model_dump() for r in v] for k, v in claim.verdict.extra_rules.items()}
        theirs = {k: [r.model_dump() for r in v] for k, v in pred.extra_rules.items()}
        if mine != theirs:
            issues.append(f"extra rule sets differ from those pre-registered for {pred.id}")
    if pred.direction and claim.statistic.alternative and pred.direction != claim.statistic.alternative:
        issues.append(f"test alternative {claim.statistic.alternative} differs from the pre-registered "
                      f"direction {pred.direction}")
    if (prereg.correction.method != "none" and pred.inferential and claim.multiple_testing is not None
            and claim.multiple_testing.method != prereg.correction.method):
        issues.append(f"correction {claim.multiple_testing.method} differs from the pre-registered "
                      f"{prereg.correction.method}")
    return issues


# --------------------------------------------------------------------------- #
# Tolerance comparison and reverify hook
# --------------------------------------------------------------------------- #
def _within(observed: Any, expected: Any, spec: dict[str, float] | None) -> bool:
    from sharur.reverify import match  # noqa: PLC0415

    if (spec and isinstance(expected, (int, float)) and not isinstance(expected, bool)
            and isinstance(observed, (int, float)) and not isinstance(observed, bool)):
        if any(isinstance(v, float) and math.isnan(v) for v in (observed, expected)):
            return False
        allowed = max(float(spec.get("abs", 0.0)), float(spec.get("rel", 0.0)) * abs(float(expected)))
        return abs(float(observed) - float(expected)) <= allowed
    return match(observed, expected) == "OK"


def compare_with_tolerance(observed: Any, expected: Any,
                           tolerance: dict[str, dict[str, float]] | None) -> list[str]:
    """Values outside tolerance, as ``name observed vs expected`` strings (empty = all within).

    Dicts compare key by key (observed may carry extra keys); lists position by position;
    scalars as one value. A value's own tolerance applies first; the ``"*"`` default applies to
    floats only, so recorded counts (integers) stay exact unless given their own tolerance.
    """
    tol = tolerance or {}
    star = tol.get("*")

    def default(exp: Any) -> dict[str, float] | None:
        return star if isinstance(exp, float) else None

    if isinstance(expected, dict):
        if isinstance(observed, list) and len(observed) == 1 and isinstance(observed[0], dict):
            observed = observed[0]
        if not isinstance(observed, dict):
            return [f"observed {type(observed).__name__}, expected an object"]
        out = []
        for key, exp in expected.items():
            if key not in observed:
                out.append(f"{key} missing")
            elif not _within(observed[key], exp, tol.get(key, default(exp))):
                out.append(f"{key} {observed[key]!r} vs {exp!r}")
        return out
    if isinstance(expected, list):
        if isinstance(observed, list) and len(observed) == 1 and isinstance(observed[0], list):
            observed = observed[0]
        if not isinstance(observed, list) or len(observed) != len(expected):
            return [f"observed {observed!r} vs {expected!r}"]
        return [f"[{i}] {o!r} vs {e!r}" for i, (o, e) in enumerate(zip(observed, expected))
                if not _within(o, e, default(e))]
    return [] if _within(observed, expected, default(expected)) else [f"{observed!r} vs {expected!r}"]


def _claim_of(triple: Any) -> StatisticalClaim | None:
    raw = getattr(triple, "stat", None)
    if not raw:
        return None
    return StatisticalClaim.model_validate(raw)


def static_flags(triple: Any, root: Path | None) -> list[str]:
    """Flags for a statistical triple before any replay: inputs, seal, pre-registration."""
    try:
        claim = _claim_of(triple)
    except ValueError as exc:
        return [f"STAT_CLAIM_INVALID:{str(exc).splitlines()[0][:120]}"]
    if claim is None:
        return []
    flags = []
    for issue in check_inputs(claim, root):
        flags.append("INPUT_CHANGED:" + issue.split(": ", 1)[-1])
    default_db = getattr(triple, "default_db", None)
    if claim.dataset_id and default_db:
        seal = Path(default_db).parent / "dataset.seal.json"
        if seal.is_file():
            try:
                current = json.loads(seal.read_text()).get("dataset_id")
            except (OSError, json.JSONDecodeError):
                current = None
            if current and current != claim.dataset_id:
                flags.append("SEAL_CHANGED")
    if claim.preregistration is not None and claim.preregistration.path:
        try:
            prereg = load_preregistration(_resolve(claim.preregistration.path, root))
        except PreregError:
            flags.append("PREREG_MISSING")
        else:
            if check_against_preregistration(claim, prereg):
                flags.append("PREREG_MISMATCH")
    return flags


def rejudge(triple: Any, outcome: Any, root: Path | None = None) -> tuple[Any, list[str]]:
    """Re-judge a replayed statistical triple: tolerance, then the verdict rules.

    ``outcome`` is a :class:`sharur.reverify.Outcome`; returns a new one and the flags to add.
    """
    tolerance = getattr(triple, "tolerance", None)
    try:
        claim = _claim_of(triple)
    except ValueError:
        claim = None
    if outcome.status == "ERROR" or (tolerance is None and claim is None):
        return outcome, []
    if tolerance is None and claim is not None:
        tolerance = claim.tolerance or None
    flags: list[str] = []
    status = outcome.status
    if tolerance is not None:
        misses = compare_with_tolerance(outcome.observed, triple.expected, tolerance)
        status = "FAIL" if misses else "OK"
        flags += [f"TOLERANCE:{m}"[:160] for m in misses[:6]]
    if claim is not None and claim.verdict is not None and claim.verdict.rules:
        observed = outcome.observed
        if isinstance(observed, list) and len(observed) == 1 and isinstance(observed[0], dict):
            observed = observed[0]
        try:
            replayed = claim.evaluate(observed if isinstance(observed, dict) else None)
        except RuleError as exc:
            flags.append(f"VERDICT_RULE_ERROR:{exc}"[:160])
            status = "FAIL"
        else:
            if replayed is not None and replayed.text() != claim.verdict.text():
                flags.append(f"VERDICT_CHANGED:{claim.verdict.text()} -> {replayed.text()}"[:200])
                status = "FAIL"
    return replace(outcome, status=status), flags


# --------------------------------------------------------------------------- #
# Findings files
# --------------------------------------------------------------------------- #
def statistical_items(findings: Iterable[Path]) -> list[dict[str, Any]]:
    """Every statistical verification item of the findings files, with its location."""
    out = []
    for path in findings:
        with Path(path).open() as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                record = json.loads(line)
                for index, item in enumerate(record.get("verification") or [], 1):
                    if isinstance(item, dict) and (item.get("check") == CHECK_KIND or "statistical_claim" in item):
                        out.append({"source": str(path), "line": line_no, "index": index,
                                    "finding_id": str(record.get("id") or ""), "item": item})
    return out


def check_findings(findings: Sequence[Path], *, root: Path | None = None) -> list[dict[str, Any]]:
    """Static checks of every statistical item: one row per item with its issues."""
    rows = []
    for entry in statistical_items(findings):
        item = entry["item"]
        issues: list[str] = []
        raw = item.get("statistical_claim")
        if not raw:
            issues.append("check 'statistical' without a statistical_claim")
        else:
            try:
                claim = StatisticalClaim.model_validate(raw)
            except ValueError as exc:
                issues.append(f"invalid statistical_claim: {str(exc).splitlines()[0]}")
            else:
                seal_id = None
                ds = _dataset_seal_for(Path(entry["source"]))
                if ds is not None:
                    seal_id = ds
                issues += check_claim(claim, root=root, seal_dataset_id=seal_id)
                expected = item.get("expected")
                if isinstance(expected, dict):
                    ns = claim.namespace()
                    for key, value in expected.items():
                        if key in ns and ns[key] is not None and not _within(ns[key], value, None):
                            issues.append(f"expected {key}={value!r} differs from the claim's value {ns[key]!r}")
        rows.append({"finding_id": entry["finding_id"], "source": entry["source"], "line": entry["line"],
                     "index": entry["index"], "claim": item.get("claim", ""), "issues": issues})
    return rows


def _dataset_seal_for(findings_path: Path) -> str | None:
    """dataset_id of ``data/<dataset>/dataset.seal.json`` for a findings file under it."""
    for parent in findings_path.resolve().parents:
        seal = parent / "dataset.seal.json"
        if seal.is_file() and (parent / "sharur.duckdb").exists():
            try:
                return json.loads(seal.read_text()).get("dataset_id")
            except (OSError, json.JSONDecodeError):
                return None
        if parent.name == "data":
            break
    return None


__all__ = [
    "CHECK_KIND",
    "ControlResult",
    "Effect",
    "InputFile",
    "MultipleTesting",
    "PreregRef",
    "Reproduce",
    "RuleError",
    "Statistic",
    "StatisticalClaim",
    "Verdict",
    "apply_rules",
    "bh",
    "check_against_preregistration",
    "check_claim",
    "check_findings",
    "check_inputs",
    "compare_with_tolerance",
    "evaluate_rule",
    "holm",
    "rejudge",
    "rule_names",
    "static_flags",
    "statistical_check",
    "statistical_items",
]
