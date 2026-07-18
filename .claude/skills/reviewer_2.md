# Reviewer 2 Skill

Adversarial manuscript claim verification. Named after the notoriously thorough Reviewer 2, this skill independently re-derives every quantitative claim from the database and scrutinizes interpretive claims for evidence calibration.

**CONCURRENCY:** Independent `read_only=True` review sessions may query in parallel.
Serialize any corrective operation that writes DuckDB.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Apply Context-First protocol, superfamily awareness, COUNT(DISTINCT protein_id), and
> accession verification throughout.

---

## Usage

```
/reviewer_2                          # Review the default/most recent dataset
/reviewer_2 --dataset DATASET_NAME   # Review a specific dataset
```

---

## Mission

**Catch errors before actual reviewers do.** You are an independent adversarial reviewer who:

1. Re-executes queries to verify every number in the manuscript
2. Checks that interpretive claims match evidence strength
3. Flags comparative claims lacking references
4. Detects known error patterns (superfamily inflation, COUNT(*) vs DISTINCT, absence without quality caveat)
5. Produces a prioritized correction queue

**You do NOT modify the manuscript.** You produce a correction queue for the author to review and act on.

---

## Workflow

### Step 0: Run the Verification Script

The mechanical verification runs first as a Python script. Then you handle the interpretive review.

```python
import subprocess
import sys
from pathlib import Path

DATASET = "data/DATASET_NAME"  # Adjust to target dataset
dataset_dir = Path(DATASET)

# Run verify_claims.py
result = subprocess.run(
    [sys.executable, "scripts/verify_claims.py",
     "--dataset", DATASET, "--auto-extract"],
    capture_output=True, text=True
)
print(result.stdout)
if result.stderr:
    print(result.stderr)
```

This produces:
- `CLAIM_VERIFICATION.jsonl` — per-claim verification records
- `REVIEW_REPORT.md` — human-readable summary

### Step 1: Load Verification Results

```python
import json

sys.path.insert(0, '.')
from sharur.operators import Sharur

b = Sharur(f"{DATASET}/sharur.duckdb", read_only=True)

REVIEW_DIR = dataset_dir / "review"
REVIEW_DIR.mkdir(exist_ok=True)

# Load verification records
verifications = []
ver_path = dataset_dir / "CLAIM_VERIFICATION.jsonl"
if ver_path.exists():
    with open(ver_path) as f:
        for line in f:
            if line.strip():
                verifications.append(json.loads(line))

# Categorize
discrepant = [v for v in verifications if v["verification_status"] == "DISCREPANT"]
unverifiable = [v for v in verifications if v["verification_status"] == "UNVERIFIABLE"]
needs_context = [v for v in verifications if v["verification_status"] == "NEEDS_CONTEXT"]
confirmed = [v for v in verifications if v["verification_status"] == "CONFIRMED"]

print(f"Verification results: {len(confirmed)} confirmed, {len(discrepant)} discrepant, "
      f"{len(needs_context)} needs context, {len(unverifiable)} unverifiable")
```

### Step 2: Triage Discrepancies

For each DISCREPANT claim, classify severity:

```python
corrections = []

for v in discrepant:
    claimed = v.get("claimed_value", "")
    verified = v.get("verified_value", "")

    # Parse numeric values for comparison
    try:
        c_num = float(claimed.replace("%", "").replace(",", ""))
        v_num = float(verified.replace("%", "").replace(",", ""))
        diff_pct = abs(c_num - v_num) / max(c_num, 1) * 100

        if diff_pct < 1:
            severity = "negligible"  # Rounding difference
        elif diff_pct < 10:
            severity = "meaningful"  # Worth correcting
        else:
            severity = "critical"    # Major error
    except (ValueError, TypeError):
        severity = "meaningful"      # Can't compare, flag for review

    corrections.append({
        "claim_id": v["claim_id"],
        "severity": severity,
        "claimed": claimed,
        "verified": verified,
        "query": v.get("verification_query"),
        "note": v.get("note"),
        "claim_text": v.get("claim_text", "")[:150],
    })

# Sort by severity (critical first)
severity_order = {"critical": 0, "meaningful": 1, "negligible": 2}
corrections.sort(key=lambda c: severity_order.get(c["severity"], 1))

for c in corrections:
    print(f"[{c['severity'].upper()}] {c['claim_id']}: claimed {c['claimed']}, verified {c['verified']}")
```

### Step 3: Interpretive Claim Review

For UNVERIFIABLE claims (interpretive, comparative, methodological), review manually:

```python
# Load findings for cross-reference
survey_findings = []
explore_findings = []

survey_path = dataset_dir / "survey" / "findings.jsonl"
if survey_path.exists():
    with open(survey_path) as f:
        for line in f:
            if line.strip():
                survey_findings.append(json.loads(line))

explore_path = dataset_dir / "exploration" / "findings.jsonl"
if explore_path.exists():
    with open(explore_path) as f:
        for line in f:
            if line.strip():
                explore_findings.append(json.loads(line))

all_findings = survey_findings + explore_findings
finding_by_id = {f.get("id"): f for f in all_findings if f.get("id")}

# Load manuscript claims if available
claims_path = dataset_dir / "MANUSCRIPT_CLAIMS.jsonl"
manuscript_claims = []
if claims_path.exists():
    with open(claims_path) as f:
        for line in f:
            if line.strip():
                manuscript_claims.append(json.loads(line))
```

For each unverifiable claim, check:

**3a. Source findings exist and contain supporting evidence**

```python
for claim in manuscript_claims:
    claim_type = claim.get("claim_type", "unknown")
    if claim_type in ("interpretive", "comparative"):
        source_ids = claim.get("source_findings", [])
        for sid in source_ids:
            if sid not in finding_by_id:
                print(f"MISSING SOURCE: {claim.get('claim_id')} references {sid} — not found")
            else:
                finding = finding_by_id[sid]
                # Check that finding has evidence supporting the claim
                evidence = finding.get("evidence", {})
                provenance = finding.get("provenance", {})
                if not evidence and not provenance:
                    print(f"WEAK SOURCE: {claim.get('claim_id')} → {sid} has no evidence/provenance")
```

**3b. Language calibration**

Check for overclaimed language that doesn't match evidence level:

```python
# Forbidden language in claims without experimental evidence
OVERCLAIM_PATTERNS = [
    (r'\bconfirms?\b', "Use 'suggests' or 'supports' unless experimental evidence exists"),
    (r'\bproves?\b', "Use 'provides evidence for'"),
    (r'\bdemonstrates?\b', "Acceptable only with strong experimental evidence"),
    (r'\bunprecedented\b', "Verify with literature search"),
    (r'\bfirst ever\b', "Verify with literature search"),
    (r'\bfirst reported\b', "Must be verified by literature search"),
    (r'\blargest known\b', "Must cite supporting reference"),
    (r'\bsmallest known\b', "Must cite supporting reference"),
    (r'\bgroundbreaking\b', "Remove — inappropriate for scientific writing"),
    (r'\brevolutionary\b', "Remove — inappropriate for scientific writing"),
    (r'\bparadigm[- ]shifting\b', "Remove — inappropriate for scientific writing"),
]

import re

manuscript_path = dataset_dir / "MANUSCRIPT.md"
if manuscript_path.exists():
    text = manuscript_path.read_text()
    for pattern, guidance in OVERCLAIM_PATTERNS:
        matches = list(re.finditer(pattern, text, re.IGNORECASE))
        for match in matches:
            # Find the line number
            line_num = text[:match.start()].count('\n') + 1
            context = text[max(0, match.start()-40):match.end()+40].replace('\n', ' ')
            print(f"Line {line_num}: '{match.group()}' — {guidance}")
            print(f"  Context: ...{context}...")
```

**3c. Comparative claims**

```python
# Find comparative claims that need literature support
COMPARATIVE_PATTERNS = [
    r'(?:larger|smaller|more|fewer|higher|lower)\s+than\s+(?:any\s+)?(?:known|reported|published)',
    r'(?:first|only|unique)\s+(?:known|reported|observed)',
    r'(?:highest|lowest|largest|smallest)\s+(?:known|reported|in\s+any)',
    r'unlike\s+(?:any|all|most)\s+(?:known|reported|other)',
]

for pattern in COMPARATIVE_PATTERNS:
    matches = list(re.finditer(pattern, text, re.IGNORECASE))
    for match in matches:
        line_num = text[:match.start()].count('\n') + 1
        context = text[max(0, match.start()-40):match.end()+40].replace('\n', ' ')
        print(f"COMPARATIVE CLAIM (line {line_num}): ...{context}...")
        print("  → Must be supported by literature citation or search")
```

### Step 4: Generate Correction Queue

```python
# Combine all findings into a prioritized correction queue
correction_queue = []

# 1. Critical discrepancies
for c in corrections:
    if c["severity"] == "critical":
        correction_queue.append({
            "priority": 1,
            "type": "discrepancy",
            "severity": "CRITICAL",
            "claim_id": c["claim_id"],
            "action": f"Correct: claimed {c['claimed']}, database shows {c['verified']}",
            "claim_text": c["claim_text"],
        })

# 2. Meaningful discrepancies
for c in corrections:
    if c["severity"] == "meaningful":
        correction_queue.append({
            "priority": 2,
            "type": "discrepancy",
            "severity": "MEANINGFUL",
            "claim_id": c["claim_id"],
            "action": f"Correct: claimed {c['claimed']}, database shows {c['verified']}",
            "claim_text": c["claim_text"],
        })

# 3. Error patterns (from verification records)
for v in verifications:
    if v.get("error_patterns"):
        for p in v["error_patterns"]:
            if p["severity"] == "warning":
                correction_queue.append({
                    "priority": 3,
                    "type": "error_pattern",
                    "severity": "WARNING",
                    "claim_id": v["claim_id"],
                    "action": f"{p['pattern']}: {p['detail']}",
                    "claim_text": v.get("claim_text", "")[:150],
                })

# 4. Language issues (add from Step 3b results above)
# 5. Missing sources (add from Step 3a results above)

# Sort by priority
correction_queue.sort(key=lambda x: x["priority"])

# Write correction_queue.md
queue_lines = [
    "# Correction Queue",
    "",
    f"**Dataset:** `{dataset_dir}`",
    f"**Total items:** {len(correction_queue)}",
    "",
    "## Critical (fix before submission)",
    "",
]

for item in correction_queue:
    if item["priority"] == 1:
        queue_lines.append(f"- **{item['claim_id']}** [{item['type']}]: {item['action']}")
        queue_lines.append(f"  - Claim: {item['claim_text'][:120]}")
        queue_lines.append("")

queue_lines.append("## Meaningful (should fix)")
queue_lines.append("")
for item in correction_queue:
    if item["priority"] == 2:
        queue_lines.append(f"- **{item['claim_id']}** [{item['type']}]: {item['action']}")
        queue_lines.append(f"  - Claim: {item['claim_text'][:120]}")
        queue_lines.append("")

queue_lines.append("## Warnings (review)")
queue_lines.append("")
for item in correction_queue:
    if item["priority"] == 3:
        queue_lines.append(f"- **{item['claim_id']}** [{item['type']}]: {item['action']}")
        queue_lines.append("")

with open(REVIEW_DIR / "correction_queue.md", "w") as f:
    f.write("\n".join(queue_lines))

print(f"Wrote correction queue to {REVIEW_DIR / 'correction_queue.md'}")
```

### Step 5: Write Interpretive Review

```python
# Write the agent's interpretive review as a separate document
review_lines = [
    "# Interpretive Review",
    "",
    f"**Dataset:** `{dataset_dir}`",
    "",
    "## Language Calibration",
    "",
    "<!-- Add findings from Step 3b here -->",
    "",
    "## Comparative Claims Requiring Literature Support",
    "",
    "<!-- Add findings from Step 3c here -->",
    "",
    "## Claims with Weak or Missing Source Findings",
    "",
    "<!-- Add findings from Step 3a here -->",
    "",
    "## Overall Assessment",
    "",
    "<!-- Agent writes an overall assessment of manuscript quality -->",
    "",
]

with open(REVIEW_DIR / "interpretive_review.md", "w") as f:
    f.write("\n".join(review_lines))

print(f"Wrote interpretive review to {REVIEW_DIR / 'interpretive_review.md'}")
```

---

## Validation Protocols

All seven shared validation protocols from `_validation_protocols.md` apply during review. Key ones for this skill:

### Accession Verification (Protocol 1)
When a discrepancy involves a PFAM/KEGG accession, verify the accession name matches the claimed function:

```python
name_check = b.store.execute(
    "SELECT DISTINCT name FROM annotations WHERE accession = ?", [accession]
)
```

If the name doesn't match the function claimed in the manuscript, this is a **critical** error.

### COUNT(DISTINCT) Rule (Protocol 2)
When a protein count discrepancy exists, check whether the manuscript used COUNT(*) instead of COUNT(DISTINCT protein_id):

```python
count_star = b.store.execute(
    "SELECT COUNT(*) FROM annotations WHERE accession = ?", [acc]
)[0][0]
count_distinct = b.store.execute(
    "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = ?", [acc]
)[0][0]
# If count_star matches the manuscript but count_distinct doesn't, the error is repeat-domain inflation
```

### Superfamily Awareness (Protocol 3)
When a claim attributes a specific enzymatic function to a domain that averages >10 hits/genome, flag it:

```python
n_proteins = b.store.execute(
    "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = ?", [acc]
)[0][0]
n_genomes = b.store.execute(
    "SELECT COUNT(DISTINCT bin_id) FROM proteins"
)[0][0]
hits_per_genome = n_proteins / n_genomes
if hits_per_genome > 10:
    # This is a superfamily — the claimed function needs neighborhood validation
    pass
```

### MAG Quality (Protocol 5)
Any absence claim ("genome X lacks Y") must be accompanied by a quality caveat. Check:

```python
contig_count = b.store.execute(
    "SELECT COUNT(DISTINCT contig_id) FROM proteins WHERE bin_id = ?", [bin_id]
)[0][0]
# If >200 contigs, absence claims must say "not detected" with fragmentation caveat
```

---

## Output Structure

```
{dataset}/
├── CLAIM_VERIFICATION.jsonl      # Per-claim verification records (from script)
├── REVIEW_REPORT.md              # Human-readable summary (from script)
└── review/
    ├── correction_queue.md       # Prioritized corrections (from agent)
    └── interpretive_review.md    # Non-quantitative claim assessment (from agent)
```

### CLAIM_VERIFICATION.jsonl Schema

```json
{
  "claim_id": "MC001",
  "claim_text": "Energy-conserving Group 4 NiFe hydrogenases are present in 74.6% of genomes",
  "claim_type": "quantitative",
  "section": "Results",
  "claimed_value": "74.6%",
  "verification_query": "SELECT ROUND(COUNT(DISTINCT p.bin_id) * 100.0 / 1831, 1) ...",
  "verified_value": "74.6",
  "verification_status": "CONFIRMED",
  "error_patterns": null,
  "note": "Verified (within rounding tolerance)",
  "auto_extracted": false
}
```

**verification_status values:**
- `CONFIRMED` — Claimed value matches database (within tolerance)
- `DISCREPANT` — Claimed value does not match database
- `NEEDS_CONTEXT` — Query returned a value but comparison requires judgment
- `UNVERIFIABLE` — No automated verification possible (interpretive/comparative claim)

---

## What This Skill Catches

Based on known error classes from previous manuscripts:

| Error Class | Example | Detection Method |
|-------------|---------|-----------------|
| Wrong count | "864 WD40 proteins" (actually 352 — used COUNT(*)) | COUNT(*) vs COUNT(DISTINCT) comparison |
| Superfamily inflation | "CBASS in 96.6% of genomes" (actually DefenseFinder broad match) | Hits-per-genome threshold |
| Accession mismatch | "benzoyl-CoA reductase (PF04055)" (actually Radical_SAM) | DB name lookup |
| Absence overclaim | "Genome X lacks hydrogenases" (200+ contigs) | Contig count check |
| Ratio error | "14.5:1 synthesis:degradation" (forgot HD-GYP family) | Cross-reference all enzyme families |
| Language overclaim | "demonstrates" for computational evidence | Regex pattern matching |
| Missing reference | "largest known in any phylum" (unverified) | Comparative claim detection |

---

## Integration with Manuscript Pipeline

`/reviewer_2` runs **after** MANUSCRIPT.md is assembled but **before** the literature agent:

```
MANUSCRIPT.md + MANUSCRIPT_CLAIMS.jsonl
  ↓
scripts/validate_provenance.py → PROVENANCE_AUDIT.md    (chain integrity)
  ↓
scripts/verify_claims.py → CLAIM_VERIFICATION.jsonl     (number verification)
  ↓
/reviewer_2 agent → review/correction_queue.md          (full adversarial review)
  ↓  fix corrections
/literature agent → literature_citations.jsonl           (citation verification)
  ↓
MANUSCRIPT.pdf
```

The correction queue should be addressed before running the literature agent, since number corrections may change the text that the literature agent needs to verify citations for.
