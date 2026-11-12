# Statistical checks and pre-registration

A finding's verification triples (`{claim, query, expected}`, see
[Structured findings](../findings_spec.md)) replay counts exactly. A statistical result needs
more: the statistic and its p value, the effect size with its interval, how many independent
events stand behind it, which correction family it belongs to, which controls it had to beat,
the verdict rule that turns those numbers into a call, and the pre-registration that fixed the
rule before the numbers existed. Sharur records all of that in two places:

- a **pre-registration record** (`sharur/preregistration.py`), written before the analysis
  and frozen by its content hash;
- a **statistical check**: a verification triple with three extra keys, which cites the
  record and which `sharur reverify` recomputes within a stated tolerance.

Both work alongside count-style triples: the normalizer, the strict findings writer and
`sharur reverify` accept a statistical check as an ordinary triple.

The examples below come from **study-x**, the synthetic worked example of
`tests/test_stat_checks.py`; every label and number in it is invented.

## 1. Pre-register

Write a spec with the predictions, each with its test, direction and ordered verdict rules:

```json
{
  "id": "study-x",
  "title": "study-x: contact-site departures after factor loss (synthetic)",
  "design_inputs": ["presence calls", "species tree", "loss events and sister sets", "placebo factors"],
  "predictions": [
    {"id": "H4", "text": "loss of factor B -> contact sites on partner B",
     "test": "one-sided Wilcoxon signed-rank on per-contrast Delta (target minus matched control sites, loss minus sister); mean Delta with a 10,000-resample bootstrap 95% CI",
     "statistic": "mean_delta", "direction": "greater",
     "verdict_rules": [
       {"label": "SUPPORTED", "when": "q_BH < 0.05 and lo > 0 and p_placebo <= 0.1"},
       {"label": "NOT SUPPORTED (Delta passes; placebo control fails: {placebos_at_or_beyond} of 10 placebo means at or beyond the target)",
        "when": "q_BH < 0.05 and lo > 0"},
       {"label": "NOT SUPPORTED (raw change at the target sites, matched controls change as much)", "when": "p_raw < 0.05"},
       {"label": "NOT SUPPORTED", "when": "true"}]}
  ],
  "correction": {"method": "BH", "family": ["H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10"], "alpha": 0.05}
}
```

```bash
sharur stats preregister spec.json --db data/DATASET/sharur.duckdb \
    --document analyses/study-x/study_plan.txt --out-dir analyses/study-x/prereg
```

The command fills `dataset` from the seal beside the database (`dataset_id`, strength, seal
time; read only), binds the full pre-registration text by SHA-256 when `--document` is given,
stamps `written_at` in UTC, computes `content_sha256` over the whole record, and writes
`<id>_<UTC time>.prereg.json` read-only. It never replaces an existing record. It prints the
path and the hash.

`sharur stats verify-prereg RECORD` recomputes the hash and the document hash; `--seal` also
checks that the dataset still matches the seal the record names.

Verdict rules are expressions over named values: numbers, comparisons, `and`, `or`, `not`,
`+ - * /`, `abs`, `min`, `max`, `true`, `false`. A comparison with a missing value is false.
Rules apply in order and the first true one gives the label; a label can quote values in
braces. `extra_rules` holds further named rule sets, such as a calibration verdict.

## 2. Record the result as a statistical check

```json
{
  "claim": "H4: contact sites on partner B change more where factor B was lost (mean Delta 0.118 [0.061, 0.175], 26 of 29 loss events with data, BH q 0.00525, 0 of 10 placebos at or beyond the target)",
  "query": "shell: python3 analyses/study-x/outcomes.py --case H4 --json",
  "expected": {"mean": 0.118, "lo": 0.061, "hi": 0.175, "q_BH": 0.00525, "p_raw": 0.0004,
               "p_placebo": 0.09091, "placebos_at_or_beyond": 0},
  "check": "statistical",
  "tolerance": {"mean": {"abs": 0.005}, "*": {"rel": 0.05}},
  "statistical_claim": {
    "prediction": "H4",
    "preregistration": {"id": "study-x",
                        "content_sha256": "80a9210ca401bb223da6cb28ba3898be298e3cf506a336b7f0322064aa26dbef",
                        "path": "analyses/study-x/prereg/study-x_20260115T100000Z.prereg.json"},
    "dataset_id": "0f575530a055df85c285821bed7c12a6a9280cf62da6600b063a3dafa2f1d00f",
    "statistic": {"name": "W+", "test": "Wilcoxon signed-rank", "p": 0.0011, "alternative": "greater"},
    "effect": {"name": "mean_delta", "value": 0.118, "interval": [0.061, 0.175], "level": 0.95,
               "method": "bootstrap over contrasts, 10,000 resamples"},
    "n_events": 29, "event_unit": "internal loss edge", "n_used": 26,
    "multiple_testing": {"method": "BH", "family": "H1-H10", "family_size": 10, "q": 0.00525},
    "controls": [
      {"name": "placebos", "kind": "placebo", "n": 10, "observed": 0.051, "at_or_beyond": 0,
       "p_empirical": 0.09091, "rule": "p_placebo <= 0.1", "passed": true},
      {"name": "random edges", "kind": "random_edge", "n": 200, "p_empirical": 0.004975}
    ],
    "values": {"mean": 0.118, "lo": 0.061, "hi": 0.175, "q_BH": 0.00525, "p_raw": 0.0004,
               "p_placebo": 0.09091, "placebos_at_or_beyond": 0, "n_pred": 21, "p_random": 0.004975},
    "verdict": {"label": "SUPPORTED", "rules": ["... the four pre-registered H4 rules, verbatim ..."]},
    "reproduce": {
      "script": "analyses/study-x/outcomes.py",
      "script_sha256": "6ff235bdc460fce3261d3fa4cb38490e84ae275a471d6eed4ed664f05fa3f19e",
      "inputs": [{"path": "analyses/study-x/tables/contrasts.tsv",
                  "sha256": "022bc3cb96c58346119d44e4a92aca581ba3f6b0395eca257b2595a7aac83690"}],
      "outputs": ["analyses/study-x/tables/outcomes.tsv"]
    },
    "computed_at": "2026-01-15T14:00:00Z"
  }
}
```

`query` is any replayable form (`shell:`, `python:`, SQL); a script that prints the case's
values as JSON fits naturally. Build the record in Python with
`sharur.stat_checks.statistical_check(claim_text, query, claim, expected, tolerance)` and write
it with `append_finding_record` as usual.

What each part records:

| Key | Content |
|---|---|
| `statistic` | test, statistic value, p value, alternative |
| `effect` | effect size, interval, level, interval method |
| `n_events`, `n_used`, `event_unit` | independent events behind the test and how many had data |
| `multiple_testing` | method, family, family size, adjusted q |
| `controls` | each control run (placebo, shuffle, random edge, matched set, planted positive, calibration): its size, its own value, how many controls reached the target, its empirical p, and the rule the claim must pass (`passed` records the outcome) |
| `values` | the named values the verdict rules read |
| `verdict` | the label and the rules that produced it (identical to the pre-registered rules) |
| `reproduce` | replay command, script and input files with SHA-256 |
| `preregistration`, `dataset_id`, `computed_at` | what the claim implements, on which sealed dataset, and when |

## 3. Check and recompute

`sharur stats check FINDINGS` reads every statistical check as recorded and reports, per
check: the effect outside its interval; a control whose recorded `passed` disagrees with its
rule; a recorded verdict the rules do not give; a pre-registration whose hash, prediction,
rules, direction, correction method or dataset differ from the claim, or that was written after
the claim was computed; a dataset seal that now names a different dataset; and any script or
input whose hash changed. `--strict` exits 1 on any issue.

`sharur reverify` replays the `query` of each statistical check like any triple, then:

1. compares each observed value with `expected` under `tolerance` (`{"abs": a, "rel": r}`
   per value; `"*"` is the default for float values; integers and values without a tolerance
   keep the precision rule of count checks);
2. applies the verdict rules to the observed values and compares the label with the
   recorded verdict.

The check fails when a value leaves its tolerance (`TOLERANCE:name observed vs expected`) or
the verdict changes (`VERDICT_CHANGED:old -> new`). Before replaying, it flags
`INPUT_CHANGED:path`, `SEAL_CHANGED` (the seal beside the dataset names another dataset),
`PREREG_MISMATCH` and `PREREG_MISSING`. The flags appear in the markdown, JSON and TSV
reports. Count-style triples in the same file replay exactly as before.

## Worked example: study-x

study-x pre-registers ten predictions (H1–H10) that a factor's loss changes the contact sites
on its partner (H8–H10 predict a decrease), one BH family over them, the verdict rule above,
and a calibration rule set on H1–H3. `tests/data/make_prereg_worked_example.py` generates its
outcome table: q values from statsmodels and verdicts from plain decision logic, so the tests
check Sharur's rule engine and its BH and Holm against independent implementations.

| Case | Direction | BH q | Recorded verdict |
|---|---|---|---|
| H1 | > | 0.3 | NOT SUPPORTED; calibration DIRECTION ONLY |
| H2 | > | 0.004 | SUPPORTED; calibration RECOVERED |
| H3 | > | 0.5222 | NOT SUPPORTED; calibration NOT RECOVERED |
| H4 | > | 0.00525 | SUPPORTED |
| H5 | > | 0.00525 | NOT SUPPORTED (Delta passes; placebo control fails: 2 of 10 placebo means at or beyond the target) |
| H6 | > | 0.3875 | NOT SUPPORTED (raw change at the target sites, matched controls change as much) |
| H7 | > | 0.66 | NOT SUPPORTED |
| H8 | < | 0.00525 | SUPPORTED |
| H9 | < | 0.0064 | NOT SUPPORTED (Delta passes; placebo control fails: 1 of 10 placebo means at or beyond the target) |
| H10 | < | 0.2 | NOT SUPPORTED |

`tests/test_stat_checks.py` expresses the study in this format and shows that:

- the pre-registration record is read-only, named by its UTC time, and tamper-evident (an
  edited rule or document breaks its hash);
- one statistical check per prediction passes `check_claim` against the record, and the rules
  reproduce all ten recorded verdicts, so every rule fires, in both directions;
- BH and Holm recomputed from the ten p values match the recorded values within 0.1%;
- `check_claim` names each inconsistency: interval, verdict, timing, rules, dataset, seal,
  control outcome, unknown prediction, changed input;
- under `sharur reverify`, a replay within tolerance passes, a drifted mean fails, a q that
  crosses 0.05 changes H4's verdict and fails, and an edited input, a new seal, an edited
  record and a missing record are flagged.
