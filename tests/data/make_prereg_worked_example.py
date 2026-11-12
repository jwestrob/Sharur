#!/usr/bin/env python3
"""Generate tests/data/prereg_worked_example.json: a synthetic pre-registered study ("study-x").

Every label and number is invented. The study tests ten predictions (H1-H10) that a factor's
loss changes the contact sites on its partner; H1-H3 also carry a calibration verdict. The
values are chosen so that every verdict rule fires at least once, in both directions.

The recorded q values come from statsmodels and the recorded verdicts from the plain decision
logic below, so tests/test_stat_checks.py checks Sharur's rule engine and its BH/Holm against
independent implementations.

    python tests/data/make_prereg_worked_example.py
"""

import hashlib
import json
from pathlib import Path

from statsmodels.stats.multitest import multipletests

OUT = Path(__file__).with_name("prereg_worked_example.json")
N_PLACEBOS = 10
DOCUMENT = ("study-x pre-registration (synthetic). For each factor, compare contact-site departures "
            "between loss lineages and their sisters against matched control sites. Verdict: SUPPORTED "
            "when BH q < 0.05, the 95% CI of the mean difference excludes 0 in the predicted direction, "
            "and the target exceeds every placebo (placebo p <= 0.1).\n")

# case, label, direction, events, contrasts, n_pred, mean, lo, hi, p, raw, p_raw,
# placebos at or beyond the target, largest placebo value, random-edge p, calibration
CASES = [
    ("H1", "loss of factor A -> contact sites on partner A", ">", 9, 9, 7,
     0.031, -0.012, 0.074, 0.21, 0.035, 0.12, 0, 0.022, 0.06, True),
    ("H2", "loss of factor B' -> contact sites on partner B'", ">", 24, 22, 18,
     0.094, 0.041, 0.147, 0.0004, 0.101, 0.0007, 0, 0.038, 0.004975, True),
    ("H3", "loss of factor C -> contact sites on partner C", ">", 12, 12, 5,
     -0.008, -0.051, 0.036, 0.47, 0.004, 0.42, 3, 0.044, 0.41, True),
    ("H4", "loss of factor B -> contact sites on partner B", ">", 29, 26, 21,
     0.118, 0.061, 0.175, 0.0011, 0.122, 0.0004, 0, 0.051, 0.004975, False),
    ("H5", "loss of factor D -> contact sites on partner D", ">", 31, 31, 22,
     0.072, 0.033, 0.111, 0.0021, 0.069, 0.01, 2, 0.081, 0.02, False),
    ("H6", "loss of factor E -> contact sites on partner E", ">", 27, 27, 13,
     0.012, -0.031, 0.054, 0.31, 0.058, 0.02, 1, 0.047, 0.11, False),
    ("H7", "loss of factor F -> contact sites on partner F", ">", 19, 18, 8,
     -0.015, -0.049, 0.019, 0.66, 0.003, 0.3, 6, 0.039, 0.62, False),
    ("H8", "loss of factor G -> preferred-codon share", "<", 28, 28, 21,
     -0.061, -0.103, -0.019, 0.0016, -0.058, 0.004, 0, -0.022, 0.004975, False),
    ("H9", "loss of factor H -> preferred-codon share", "<", 22, 21, 16,
     -0.044, -0.079, -0.009, 0.0032, -0.04, 0.03, 1, -0.051, 0.03, False),
    ("H10", "loss of factor J -> preferred-codon share", "<", 17, 17, 9,
     -0.012, -0.054, 0.031, 0.12, -0.01, 0.4, 4, -0.033, 0.27, False),
]

PLACEBO_FAILS = ("NOT SUPPORTED (Delta passes; placebo control fails: {k} of 10 placebo means at or beyond "
                 "the target)")
RAW_ONLY = "NOT SUPPORTED (raw change at the target sites, matched controls change as much)"


def sig4(x: float) -> float:
    return float(f"{x:.4g}")


def verdict(direction, mean, lo, hi, q, p_raw, p_placebo, k, n_pred, calibration):
    ci_ok = lo > 0 if direction == ">" else hi < 0
    delta_ok = q < 0.05 and ci_ok
    if delta_ok and p_placebo <= 0.1:
        text = "SUPPORTED"
    elif delta_ok:
        text = PLACEBO_FAILS.format(k=k)
    elif direction == ">" and p_raw < 0.05:
        text = RAW_ONLY
    else:
        text = "NOT SUPPORTED"
    if calibration:
        if text == "SUPPORTED":
            cal = "RECOVERED"
        elif mean > 0 and n_pred >= 7:
            cal = "DIRECTION ONLY"
        else:
            cal = "NOT RECOVERED"
        text += f"; calibration {cal}"
    return text


pvals = [c[9] for c in CASES]
q_bh = multipletests(pvals, method="fdr_bh")[1]
p_holm = multipletests(pvals, method="holm")[1]
cases = []
for c, q, h in zip(CASES, q_bh, p_holm):
    (case, label, direction, events, contrasts, n_pred, mean, lo, hi, p, raw, p_raw, k, pmax, p_random,
     calibration) = c
    p_placebo = sig4((1 + k) / (N_PLACEBOS + 1))
    cases.append({
        "case": case, "label": label, "direction": direction, "events": events, "contrasts": contrasts,
        "n_pred": n_pred, "mean": mean, "lo": lo, "hi": hi, "p": p, "q_BH": sig4(q), "p_Holm": sig4(h),
        "raw": raw, "p_raw": p_raw, "p_placebo": p_placebo, "placebos_at_or_beyond": k,
        "placebo_max": pmax, "p_random": p_random, "calibration": calibration,
        "verdict": verdict(direction, mean, lo, hi, q, p_raw, p_placebo, k, n_pred, calibration)})

fixture = {
    "source": "synthetic: tests/data/make_prereg_worked_example.py (all labels and numbers invented)",
    "document": {"name": "study_plan.txt", "text": DOCUMENT,
                 "sha256": hashlib.sha256(DOCUMENT.encode()).hexdigest(), "bytes": len(DOCUMENT.encode())},
    "written_at": "2026-01-15T10:00:00Z",
    "dataset": {"dataset_id": hashlib.sha256(b"synthetic dataset for study-x").hexdigest(),
                "seal_strength": "structural", "database": "data/study_x/sharur.duckdb"},
    "cases": cases,
}
OUT.write_text(json.dumps(fixture, indent=1) + "\n")
print(f"wrote {OUT}: " + ", ".join(f"{c['case']} {c['verdict']}" for c in cases))
