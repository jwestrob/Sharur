#!/usr/bin/env python3
"""Verify the numeric claims in a findings.jsonl by re-executing verification queries.

Per docs/findings_spec.md, every specific number in a finding must carry a
verification query. This script parses those `verification` arrays, executes
the queries against the dataset's sharur.duckdb, and reports pass/fail.

Loophole closed: queries that are just a string reference to a methods section
(e.g. "See exploration/foo.md methods") are NOT executable. They are reported
as STRING_REF and treated as FAIL in --strict mode (the default).

Usage:
    python scripts/verify_claims.py --dataset data/srvp_bacteria_pb
    python scripts/verify_claims.py --dataset data/coronamine_boiler_100nm --no-strict

Exit codes:
    0 — all executable verifications passed (strict mode counts STRING_REF as FAIL)
    1 — at least one verification failed or could not be executed in strict mode
    2 — internal error (missing dataset, malformed findings, etc.)
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import duckdb

# Verification verdicts
PASS = "PASS"
FAIL = "FAIL"
STRING_REF = "STRING_REF"          # query is "See foo.md methods" — un-executable
MANUAL = "MANUAL"                  # query needs human eyeballs (Python one-liner we won't exec)
SHELL_SKIPPED = "SHELL_SKIPPED"    # shell command but --allow-shell not set
ERROR = "ERROR"                    # query ran but raised
MISSING = "MISSING"                # no verification array at all


# Strings that, when leading the query, indicate it is a non-executable
# reference to a methods section or external doc rather than a query.
STRING_REF_PREFIXES = (
    "see ", "refer to ", "per ", "as described in ", "methods:",
    "manual:", "n/a", "none", "tbd",
)


@dataclass
class VerificationResult:
    finding_id: str
    claim: str
    query: str
    expected: Any
    actual: Any = None
    verdict: str = ""
    note: str = ""


@dataclass
class FindingSummary:
    finding_id: str
    title: str
    n_checks: int
    pass_count: int
    fail_count: int
    string_ref_count: int
    manual_count: int
    error_count: int
    results: list[VerificationResult] = field(default_factory=list)


def classify_query(query: str) -> str:
    """Decide what kind of query this is.

    Returns one of: 'sql', 'shell', 'python', 'string_ref'.
    """
    q = query.strip()
    if not q:
        return "string_ref"
    lower = q.lower()
    for prefix in STRING_REF_PREFIXES:
        if lower.startswith(prefix):
            return "string_ref"
    # SQL: starts with SELECT/WITH/SHOW/PRAGMA (read-only allowed verbs)
    if re.match(r"^\s*(SELECT|WITH|SHOW|PRAGMA|DESCRIBE)\b", q, re.IGNORECASE):
        return "sql"
    # Shell: leading bang, or common shell commands
    if q.startswith("!") or re.match(r"^\s*(awk|grep|wc|cut|sort|uniq|head|tail|cat|find|jq)\b", q):
        return "shell"
    # Python one-liner: contains "python" or looks like a python expr
    if re.match(r"^\s*python[0-9]?\s+-c", q) or ("(" in q and "len(" in q):
        return "python"
    # Otherwise we don't know — treat as string_ref
    return "string_ref"


def normalise_value(v: Any) -> Any:
    """Best-effort coerce strings like '42' to int/float for comparison."""
    if isinstance(v, str):
        s = v.strip()
        if re.match(r"^-?\d+$", s):
            return int(s)
        if re.match(r"^-?\d+\.\d+$", s):
            return float(s)
    return v


def _flatten_sql_result(v: Any) -> Any:
    """Normalize SQL fetchall shapes so list/tuple/row containers compare alike.

    DuckDB returns rows as tuples. When a verification expects [251, 208] but
    the actual is [(251, 208)] (one row, two cols), or (251, 208) (the run_sql
    scalar-row path), they should compare equal — the data is the same.
    """
    if isinstance(v, list):
        # Single-row result that we unpack
        if len(v) == 1 and isinstance(v[0], (list, tuple)):
            return _flatten_sql_result(list(v[0]))
        return [_flatten_sql_result(x) for x in v]
    if isinstance(v, tuple):
        return [_flatten_sql_result(x) for x in v]
    return v


def compare_values(actual: Any, expected: Any) -> bool:
    """Pass if actual matches expected, with light type coercion."""
    a = normalise_value(actual)
    e = normalise_value(expected)
    if a == e:
        return True
    # List/tuple/row alignment
    a_norm = _flatten_sql_result(a)
    e_norm = _flatten_sql_result(e)
    if a_norm == e_norm:
        return True
    # Two-list element-by-element with numeric tolerance
    if (isinstance(a_norm, list) and isinstance(e_norm, list)
            and len(a_norm) == len(e_norm)):
        if all(compare_values(ai, ei) for ai, ei in zip(a_norm, e_norm)):
            return True
    # Numeric tolerance for small floats
    if isinstance(a, (int, float)) and isinstance(e, (int, float)):
        return abs(a - e) < 1e-6 or (e != 0 and abs(a - e) / abs(e) < 1e-2)  # 1% tolerance
    # String case-insensitive
    if isinstance(a, str) and isinstance(e, str):
        return a.strip().lower() == e.strip().lower()
    return False


def run_sql(conn: duckdb.DuckDBPyConnection, query: str) -> Any:
    """Run a SELECT/WITH/etc. query and return a scalar or a list of rows."""
    rows = conn.execute(query).fetchall()
    if not rows:
        return None
    # Single scalar
    if len(rows) == 1 and len(rows[0]) == 1:
        return rows[0][0]
    # Single column list
    if len(rows[0]) == 1:
        return [r[0] for r in rows]
    return rows


def run_shell(query: str, cwd: Path) -> Any:
    """Run a shell command and return stdout stripped.

    Strips a leading '!' if present.
    """
    cmd = query.strip()
    if cmd.startswith("!"):
        cmd = cmd[1:].strip()
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, cwd=cwd, timeout=120,
    )
    if result.returncode != 0:
        raise RuntimeError(f"shell command failed (rc={result.returncode}): {result.stderr.strip()}")
    out = result.stdout.strip()
    # Try numeric coercion
    return normalise_value(out)


def verify_one(
    finding_id: str,
    check: dict,
    conn: duckdb.DuckDBPyConnection,
    cwd: Path,
    allow_shell: bool,
) -> VerificationResult:
    claim = check.get("claim", "")
    query = check.get("query", "")
    expected = check.get("expected")
    res = VerificationResult(finding_id=finding_id, claim=claim, query=query, expected=expected)

    kind = classify_query(query)
    if kind == "string_ref":
        res.verdict = STRING_REF
        res.note = "Query is a string reference, not an executable query."
        return res
    if kind == "python":
        res.verdict = MANUAL
        res.note = "Python query — needs manual execution."
        return res
    if kind == "shell" and not allow_shell:
        res.verdict = SHELL_SKIPPED
        res.note = "Shell query — re-run with --allow-shell to execute."
        return res

    try:
        if kind == "sql":
            actual = run_sql(conn, query)
        elif kind == "shell":
            actual = run_shell(query, cwd)
        else:
            res.verdict = MANUAL
            res.note = f"Unknown query kind: {kind}"
            return res
    except Exception as exc:
        res.verdict = ERROR
        res.note = f"{type(exc).__name__}: {exc}"
        return res

    res.actual = actual
    res.verdict = PASS if compare_values(actual, expected) else FAIL
    if res.verdict == FAIL:
        res.note = f"Expected {expected!r}, got {actual!r}"
    return res


def verify_finding(
    finding: dict,
    conn: duckdb.DuckDBPyConnection,
    cwd: Path,
    allow_shell: bool,
) -> FindingSummary:
    fid = finding.get("id") or finding.get("finding_id") or "?"
    title = finding.get("title", "")
    checks = finding.get("verification") or []

    if not isinstance(checks, list) or not checks:
        return FindingSummary(
            finding_id=fid, title=title,
            n_checks=0, pass_count=0, fail_count=0,
            string_ref_count=0, manual_count=0, error_count=0,
            results=[VerificationResult(
                finding_id=fid, claim="(no verification array)",
                query="", expected=None, verdict=MISSING,
                note="Finding has no verification entries.",
            )],
        )

    results = [verify_one(fid, c, conn, cwd, allow_shell) for c in checks if isinstance(c, dict)]
    return FindingSummary(
        finding_id=fid, title=title,
        n_checks=len(results),
        pass_count=sum(1 for r in results if r.verdict == PASS),
        fail_count=sum(1 for r in results if r.verdict == FAIL),
        string_ref_count=sum(1 for r in results if r.verdict == STRING_REF),
        manual_count=sum(1 for r in results if r.verdict in (MANUAL, SHELL_SKIPPED)),
        error_count=sum(1 for r in results if r.verdict == ERROR),
        results=results,
    )


def iter_findings_files(dataset_dir: Path) -> list[Path]:
    """Find every findings.jsonl under the dataset's standard phase dirs."""
    files = []
    for phase in ("survey", "exploration", "characterization", "defense", "metabolism"):
        fp = dataset_dir / phase / "findings.jsonl"
        if fp.exists():
            files.append(fp)
    return files


def load_findings(path: Path) -> list[dict]:
    findings = []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            findings.append(json.loads(line))
        except json.JSONDecodeError as exc:
            print(f"  WARN: {path}:{i} malformed JSON ({exc})", file=sys.stderr)
    return findings


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dataset", required=True, help="Dataset directory (contains sharur.duckdb).")
    p.add_argument("--duckdb", help="Override path to sharur.duckdb.")
    p.add_argument(
        "--findings", help="Verify a single findings.jsonl instead of all phase files."
    )
    p.add_argument(
        "--strict", action="store_true", default=True,
        help="Treat STRING_REF as FAIL. Default: on.",
    )
    p.add_argument("--no-strict", dest="strict", action="store_false")
    p.add_argument(
        "--allow-shell", action="store_true",
        help="Execute shell-style verification queries via subprocess.",
    )
    p.add_argument("--out", help="Write JSONL verification report here.")
    args = p.parse_args()

    dataset = Path(args.dataset).resolve()
    if not dataset.exists():
        print(f"ERROR: dataset not found: {dataset}", file=sys.stderr)
        return 2

    duckdb_path = Path(args.duckdb) if args.duckdb else dataset / "sharur.duckdb"
    if not duckdb_path.exists():
        print(f"ERROR: sharur.duckdb not found at {duckdb_path}", file=sys.stderr)
        return 2

    conn = duckdb.connect(str(duckdb_path), read_only=True)
    cwd = dataset

    if args.findings:
        files = [Path(args.findings)]
    else:
        files = iter_findings_files(dataset)
        if not files:
            print(f"No findings.jsonl files under {dataset}", file=sys.stderr)
            return 2

    summaries: list[FindingSummary] = []
    for fp in files:
        findings = load_findings(fp)
        try:
            label = str(fp.relative_to(dataset.parent))
        except ValueError:
            # dataset is symlinked; fall back to filename
            label = fp.name
        print(f"[{label}] {len(findings)} findings")
        for f in findings:
            summaries.append(verify_finding(f, conn, cwd, args.allow_shell))

    # Aggregate
    tot_checks = sum(s.n_checks for s in summaries)
    tot_pass = sum(s.pass_count for s in summaries)
    tot_fail = sum(s.fail_count for s in summaries)
    tot_strref = sum(s.string_ref_count for s in summaries)
    tot_manual = sum(s.manual_count for s in summaries)
    tot_err = sum(s.error_count for s in summaries)
    n_missing = sum(1 for s in summaries if s.n_checks == 0)

    print()
    print("=" * 70)
    print(f"  Findings audited     : {len(summaries)}")
    print(f"  Findings w/o verify  : {n_missing}")
    print(f"  Total checks         : {tot_checks}")
    print(f"    PASS               : {tot_pass}")
    print(f"    FAIL               : {tot_fail}")
    print(f"    STRING_REF         : {tot_strref}  {'(counted as FAIL in strict mode)' if args.strict else ''}")
    print(f"    MANUAL / SHELL_SKIP: {tot_manual}")
    print(f"    ERROR              : {tot_err}")
    print("=" * 70)

    if args.out:
        out = Path(args.out)
        with out.open("w") as fh:
            for s in summaries:
                payload = asdict(s)
                # Convert dataclasses in results
                payload["results"] = [asdict(r) for r in s.results]
                fh.write(json.dumps(payload, default=str) + "\n")
        print(f"Wrote report to {out}")

    # Print the worst offenders
    worst = sorted(
        [s for s in summaries if s.fail_count or s.error_count or (args.strict and s.string_ref_count)],
        key=lambda s: -(s.fail_count + s.error_count + (s.string_ref_count if args.strict else 0)),
    )[:10]
    if worst:
        print("\nWorst offenders:")
        for s in worst:
            print(f"  {s.finding_id}: {s.title[:70]}")
            for r in s.results:
                if r.verdict in (FAIL, ERROR) or (args.strict and r.verdict == STRING_REF):
                    print(f"    [{r.verdict}] {r.claim[:60]}  — {r.note}")

    fail_in_strict = tot_fail + tot_err + (tot_strref if args.strict else 0)
    return 0 if fail_in_strict == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
