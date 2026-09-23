#!/usr/bin/env python3
"""Adversarial claim verification: re-derive manuscript numbers from the database.

Complements validate_provenance.py (which checks chain integrity) by actually
re-executing queries to verify that claimed numbers match the database.

Usage:
    python scripts/verify_claims.py --dataset data/DATASET/
    python scripts/verify_claims.py --dataset data/DATASET/ --auto-extract

Outputs:
    CLAIM_VERIFICATION.jsonl  — per-claim verification record
    REVIEW_REPORT.md          — human-readable summary
"""

import argparse
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_jsonl(path: Path) -> list[dict]:
    """Load a JSONL file, returning empty list if missing."""
    if not path.exists():
        return []
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    entries.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return entries


def connect_db(dataset_dir: Path):
    """Connect to the dataset's DuckDB and return the connection."""
    import duckdb
    db_path = dataset_dir / "sharur.duckdb"
    if not db_path.exists():
        print(f"Error: database not found at {db_path}", file=sys.stderr)
        return None
    return duckdb.connect(str(db_path), read_only=True)


def get_dataset_totals(conn) -> dict:
    """Get fundamental dataset statistics for verification."""
    totals = {}
    try:
        totals["n_proteins"] = conn.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
    except Exception:
        totals["n_proteins"] = None
    try:
        totals["n_genomes"] = conn.execute(
            "SELECT COUNT(DISTINCT bin_id) FROM proteins"
        ).fetchone()[0]
    except Exception:
        totals["n_genomes"] = None
    return totals


# ---------------------------------------------------------------------------
# Claim extraction from MANUSCRIPT.md
# ---------------------------------------------------------------------------

# Patterns that capture numeric claims in scientific prose
_NUM_PATTERNS = [
    # "1,831 MAGs" / "2,921,111 proteins" / "45 genomes"
    re.compile(
        r'([\d,]+)\s+'
        r'(MAGs?|genomes?|proteins?|contigs?|genes?|enzymes?|families|systems?|'
        r'hydrogenases?|operons?|loci|arrays?|spacers?|islands?|species|'
        r'clusters?|domains?|blocks?|hits?|members?|copies|instances?)',
        re.IGNORECASE,
    ),
    # "74.6%" / "34.2% of"
    re.compile(r'([\d.]+)\s*%'),
    # ">1,000 aa" / ">5,000 amino acids"
    re.compile(r'[><≥≤]([\d,]+)\s*(aa\b|amino acids?|bp\b|kb\b|Mb\b)', re.IGNORECASE),
    # "N/M genomes" ratio pattern
    re.compile(r'(\d+)\s*/\s*(\d+)\s+(genomes?|MAGs?)'),
]


def extract_claims_from_manuscript(manuscript_path: Path) -> list[dict]:
    """Auto-extract quantitative claims from MANUSCRIPT.md lines."""
    if not manuscript_path.exists():
        return []

    claims = []
    claim_counter = 0

    with open(manuscript_path) as f:
        for line_num, line in enumerate(f, 1):
            text = line.strip()
            # Skip headings, figure captions, references, empty, metadata
            if not text:
                continue
            if text.startswith("#"):
                continue
            if text.startswith("!"):
                continue
            if text.startswith("**Figure"):
                continue
            if text.startswith("---"):
                continue
            if text.startswith("[CITE:"):
                continue

            # Check each pattern
            has_number = False
            for pat in _NUM_PATTERNS:
                if pat.search(text):
                    has_number = True
                    break

            if has_number:
                claim_counter += 1
                claims.append({
                    "claim_id": f"AUTO-{claim_counter:03d}",
                    "line_number": line_num,
                    "claim_text": text[:200],
                    "claim_type": "quantitative",
                    "source_findings": [],
                    "auto_extracted": True,
                })

    return claims


# ---------------------------------------------------------------------------
# Query generation from claim text
# ---------------------------------------------------------------------------

def parse_claimed_value(claim_text: str) -> tuple[str | None, float | None]:
    """Extract the numeric value and its type from a claim.

    Returns (value_string, numeric_value) or (None, None).
    """
    # Percentage
    m = re.search(r'([\d.]+)\s*%', claim_text)
    if m:
        return f"{m.group(1)}%", float(m.group(1))

    # Comma-separated large numbers with context
    m = re.search(r'([\d,]+)\s+(MAGs?|genomes?|proteins?|contigs?)', claim_text, re.IGNORECASE)
    if m:
        val = m.group(1).replace(",", "")
        return m.group(0), float(val)

    # N/M ratio
    m = re.search(r'(\d+)\s*/\s*(\d+)', claim_text)
    if m:
        return f"{m.group(1)}/{m.group(2)}", float(m.group(1))

    # Bare large number
    m = re.search(r'([\d,]{3,})', claim_text)
    if m:
        val = m.group(1).replace(",", "")
        return m.group(1), float(val)

    return None, None


def _term_to_predicate(text_lower: str) -> str | None:
    """Map a biological term in claim text to a predicate name for querying.

    Returns predicate string or None if no mapping found.
    """
    TERM_MAP = {
        "nife hydrogenase": "nife_hydrogenase",
        "nife group 4": "nife_group4",
        "nife group 3": "nife_group3",
        "nife group 2": "nife_group2",
        "nife group 1": "nife_group1",
        "fefe hydrogenase": "fefe_hydrogenase",
        "hydrogenase": "hydrogenase",
        "crispr": "crispr_associated",
        "crispr-cas": "crispr_associated",
        "glycosyltransferase": "glycosyltransferase_cazy",
        "transposase": "transposase",
        "adhesin": "adhesin",
        "flagell": "flagellum",
        "wood-ljungdahl": "wood_ljungdahl",
        "wood ljungdahl": "wood_ljungdahl",
        "calvin cycle": "calvin_cycle",
        "methanogenesis": "methanogenesis",
        "cytochrome oxidase": "aerobic_respiration",
        "rnf complex": "electron_transport",
        "restriction-modification": "restriction_modification",
        "toxin-antitoxin": "toxin_antitoxin",
        "defense system": "defense_system",
        "type iv pil": "type_iv_pilus",
        "s-layer": "s_layer",
        "s_layer": "s_layer",
    }
    for term, pred in TERM_MAP.items():
        if term in text_lower:
            return pred
    return None


def _term_to_annotation_name(text_lower: str) -> str | None:
    """Map a biological term to an annotation name pattern for LIKE queries."""
    NAME_MAP = {
        "glycosyltransferase": "Glyco_trans%",
        "pilz": "PilZ",
        "dockerin": "Dockerin%",
        "cohesin": "Cohesin%",
        "wd40": "WD40",
        "tpr": "TPR_%",
        "radical sam": "Radical_SAM",
    }
    for term, name_pat in NAME_MAP.items():
        if term in text_lower:
            return name_pat
    return None


def generate_verification_query(claim: dict, totals: dict) -> str | None:
    """Attempt to generate a SQL verification query for a claim.

    This uses heuristic pattern matching on claim text. Returns None if
    no query can be generated (claim requires interpretive review).
    """
    text = claim.get("claim_text", "")
    text_lower = text.lower()

    n_genomes = totals.get("n_genomes")
    n_proteins = totals.get("n_proteins")

    # --- Extract PFAM/KEGG accession if present ---
    accession = None
    m = re.search(r'(PF\d{5}|K\d{5})', text)
    if m:
        accession = m.group(1)

    # --- Extract biological term mappings ---
    predicate = _term_to_predicate(text_lower)
    annot_name = _term_to_annotation_name(text_lower)

    # === SIZE-BASED PATTERNS ===

    # "N proteins exceeding X aa" / "N exceeding X aa"
    m = re.search(
        r'([\d,]+)\s+(?:proteins?\s+)?exceeding\s+([\d,]+)\s*(?:aa|amino)',
        text_lower,
    )
    if m:
        threshold = m.group(2).replace(",", "")
        return f"SELECT COUNT(*) FROM proteins WHERE sequence_length > {threshold}"

    # "proteins exceed X aa" / "larger than X aa" / "> X aa"
    m = re.search(
        r'(?:exceed|larger than|greater than|>)\s*([\d,]+)\s*(?:aa|amino)',
        text_lower,
    )
    if m:
        threshold = m.group(1).replace(",", "")
        return f"SELECT COUNT(*) FROM proteins WHERE sequence_length > {threshold}"

    # "largest protein" / "at N aa"
    if "largest protein" in text_lower:
        return "SELECT MAX(sequence_length) FROM proteins"

    # === PERCENTAGE-BASED PATTERNS ===

    # "X% of genomes" — with predicate or accession
    if re.search(r'[\d.]+\s*%\s*(?:of\s+)?(?:genomes?|mags?)', text_lower):
        if accession and n_genomes:
            return (
                f"SELECT ROUND(COUNT(DISTINCT p.bin_id) * 100.0 / {n_genomes}, 1) "
                f"FROM annotations a "
                f"JOIN proteins p ON a.protein_id = p.protein_id "
                f"WHERE a.accession = '{accession}'"
            )
        if predicate and n_genomes:
            return (
                f"SELECT ROUND(COUNT(DISTINCT p.bin_id) * 100.0 / {n_genomes}, 1) "
                f"FROM proteins p "
                f"JOIN protein_predicates pp ON p.protein_id = pp.protein_id "
                f"WHERE '{predicate}' = ANY(pp.predicates)"
            )

    # "X% of genomes" — present in N% (another phrasing)
    if re.search(r'present in.*[\d.]+\s*%', text_lower) or \
       re.search(r'[\d.]+\s*%.*present', text_lower):
        if predicate and n_genomes:
            return (
                f"SELECT ROUND(COUNT(DISTINCT p.bin_id) * 100.0 / {n_genomes}, 1) "
                f"FROM proteins p "
                f"JOIN protein_predicates pp ON p.protein_id = pp.protein_id "
                f"WHERE '{predicate}' = ANY(pp.predicates)"
            )

    # "X% of the proteome" — fraction of total proteins
    if re.search(r'[\d.]+\s*%\s*of\s+the\s+proteome', text_lower):
        if predicate and n_proteins:
            return (
                f"SELECT ROUND(COUNT(DISTINCT p.protein_id) * 100.0 / {n_proteins}, 1) "
                f"FROM proteins p "
                f"JOIN protein_predicates pp ON p.protein_id = pp.protein_id "
                f"WHERE '{predicate}' = ANY(pp.predicates)"
            )
        if annot_name and n_proteins:
            return (
                f"SELECT ROUND(COUNT(DISTINCT a.protein_id) * 100.0 / {n_proteins}, 1) "
                f"FROM annotations a "
                f"WHERE a.name LIKE '{annot_name}'"
            )

    # === TOTAL COUNT PATTERNS ===

    # "N predicted proteins" / "N proteins in total"
    if re.search(r'[\d,]+\s+predicted\s+proteins?', text_lower) or \
       re.search(r'[\d,]+\s+proteins?\s+(were|across|in total)', text_lower) or \
       re.search(r'containing\s+[\d,]+\s+proteins?', text_lower):
        return "SELECT COUNT(*) FROM proteins"

    # "N MAGs" / "N genomes" with context
    if re.search(r'[\d,]+\s+(MAGs?|genomes?)\s+(containing|comprising|across|in\s)', text_lower):
        return "SELECT COUNT(DISTINCT bin_id) FROM proteins"

    # === ACCESSION-BASED PATTERNS ===

    if accession:
        # "in N genomes" / "in N/M genomes"
        if re.search(r'in\s+\d+\s*(?:/\s*\d+\s*)?genomes?', text_lower):
            return (
                f"SELECT COUNT(DISTINCT p.bin_id) "
                f"FROM annotations a "
                f"JOIN proteins p ON a.protein_id = p.protein_id "
                f"WHERE a.accession = '{accession}'"
            )
        # Generic accession count
        return (
            f"SELECT COUNT(DISTINCT protein_id) FROM annotations "
            f"WHERE accession = '{accession}'"
        )

    # === PREDICATE-BASED PATTERNS ===

    # "present in N genomes" / "in N/M genomes"
    if predicate and re.search(r'(?:present\s+)?in\s+[\d,]+\s*(?:/\s*[\d,]+\s*)?genomes?', text_lower):
        return (
            f"SELECT COUNT(DISTINCT p.bin_id) "
            f"FROM proteins p "
            f"JOIN protein_predicates pp ON p.protein_id = pp.protein_id "
            f"WHERE '{predicate}' = ANY(pp.predicates)"
        )

    return None


# ---------------------------------------------------------------------------
# Error pattern detectors
# ---------------------------------------------------------------------------

def detect_error_patterns(claim: dict, conn, totals: dict) -> list[dict]:
    """Run heuristic error pattern detectors on a claim.

    Returns list of {pattern_name, detail, severity} dicts.
    """
    patterns_found = []
    text = claim.get("claim_text", "")

    # 1. Superfamily inflation: accession with >10 hits/genome claimed as specific function
    accessions = re.findall(r'(PF\d{5}|K\d{5})', text)
    n_genomes = totals.get("n_genomes", 1) or 1
    for acc in accessions:
        try:
            n_proteins = conn.execute(
                "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = ?",
                [acc]
            ).fetchone()[0]
            hits_per_genome = n_proteins / n_genomes
            if hits_per_genome > 10:
                patterns_found.append({
                    "pattern": "superfamily_inflation",
                    "detail": f"{acc} averages {hits_per_genome:.1f} hits/genome — likely superfamily-level annotation",
                    "severity": "warning",
                })
        except Exception:
            pass

    # 2. COUNT(*) vs COUNT(DISTINCT) discrepancy hint
    for acc in accessions:
        try:
            count_star = conn.execute(
                "SELECT COUNT(*) FROM annotations WHERE accession = ?", [acc]
            ).fetchone()[0]
            count_distinct = conn.execute(
                "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = ?", [acc]
            ).fetchone()[0]
            if count_distinct > 0 and count_star / count_distinct > 2.0:
                patterns_found.append({
                    "pattern": "repeat_domain_inflation",
                    "detail": (
                        f"{acc}: COUNT(*)={count_star} vs COUNT(DISTINCT)={count_distinct} "
                        f"(ratio {count_star/count_distinct:.1f}x) — "
                        f"if claim uses raw count, it's inflated by repeat domains"
                    ),
                    "severity": "warning",
                })
        except Exception:
            pass

    # 3. Accession name mismatch — check if claimed function matches DB name
    for acc in accessions:
        try:
            db_name = conn.execute(
                "SELECT DISTINCT name FROM annotations WHERE accession = ? LIMIT 1", [acc]
            ).fetchone()
            if db_name:
                db_name = db_name[0]
                # Simple check: if the claim text mentions a specific function
                # that doesn't appear in the DB name at all
                # (We can't do full NLP here; just flag for agent review)
                patterns_found.append({
                    "pattern": "accession_name_check",
                    "detail": f"{acc} = '{db_name}' in database — verify this matches the claimed function",
                    "severity": "info",
                })
        except Exception:
            pass

    # 4. MAG absence without quality caveat
    absence_words = ["absent", "lacks", "missing", "devoid", "without"]
    if any(w in text.lower() for w in absence_words):
        if not any(q in text.lower() for q in ["not detected", "contigs", "fragmented", "assembly"]):
            patterns_found.append({
                "pattern": "absence_without_caveat",
                "detail": "Absence claim without MAG quality caveat — should use 'not detected' and note fragmentation",
                "severity": "warning",
            })

    return patterns_found


# ---------------------------------------------------------------------------
# Claim verification
# ---------------------------------------------------------------------------

def verify_claim(claim: dict, conn, totals: dict) -> dict:
    """Verify a single claim against the database.

    Returns a verification record.
    """
    claim_text = claim.get("claim_text", "")
    claimed_str, claimed_num = parse_claimed_value(claim_text)
    query = generate_verification_query(claim, totals)
    error_patterns = detect_error_patterns(claim, conn, totals)

    record = {
        "claim_id": claim.get("claim_id", "?"),
        "claim_text": claim_text[:200],
        "claim_type": claim.get("claim_type", "unknown"),
        "section": claim.get("section", ""),
        "claimed_value": claimed_str,
        "verification_query": query,
        "verified_value": None,
        "verification_status": "UNVERIFIABLE",
        "error_patterns": error_patterns if error_patterns else None,
        "note": None,
        "auto_extracted": claim.get("auto_extracted", False),
    }

    if not query:
        record["note"] = "No verification query generated — requires interpretive review"
        return record

    try:
        result = conn.execute(query).fetchone()
        if result is not None:
            verified_val = result[0]
            record["verified_value"] = str(verified_val)

            # Compare
            if claimed_num is not None and verified_val is not None:
                try:
                    verified_num = float(verified_val)
                    diff = abs(verified_num - claimed_num)
                    # Percentage claims: compare absolute difference
                    if "%" in (claimed_str or ""):
                        if diff < 0.2:
                            record["verification_status"] = "CONFIRMED"
                            record["note"] = "Verified (within rounding tolerance)"
                        elif diff < 1.0:
                            record["verification_status"] = "CONFIRMED"
                            record["note"] = f"Verified (difference: {diff:.1f} percentage points)"
                        else:
                            record["verification_status"] = "DISCREPANT"
                            record["note"] = f"Discrepancy: claimed {claimed_str}, verified {verified_val}"
                    else:
                        # Count claims: exact or within 1%
                        if claimed_num == 0:
                            if verified_num == 0:
                                record["verification_status"] = "CONFIRMED"
                            else:
                                record["verification_status"] = "DISCREPANT"
                                record["note"] = f"Claimed 0, found {verified_val}"
                        elif diff / max(claimed_num, 1) < 0.01:
                            record["verification_status"] = "CONFIRMED"
                            record["note"] = "Exact match"
                        elif diff / max(claimed_num, 1) < 0.05:
                            record["verification_status"] = "CONFIRMED"
                            record["note"] = f"Within 5% (claimed {claimed_str}, verified {verified_val})"
                        else:
                            record["verification_status"] = "DISCREPANT"
                            record["note"] = f"Discrepancy: claimed {claimed_str}, verified {verified_val}"
                except (ValueError, TypeError):
                    record["verification_status"] = "NEEDS_CONTEXT"
                    record["note"] = f"Query returned {verified_val}, manual comparison needed"
            else:
                record["verification_status"] = "NEEDS_CONTEXT"
                record["note"] = f"Query returned {verified_val}, could not parse claimed value"
    except Exception as e:
        record["verification_status"] = "UNVERIFIABLE"
        record["note"] = f"Query execution failed: {e}"

    return record


# ---------------------------------------------------------------------------
# Cross-reference: find uncovered manuscript lines
# ---------------------------------------------------------------------------

def find_uncovered_statements(manuscript_path: Path, claims: list[dict]) -> list[dict]:
    """Find lines in manuscript with numbers but no claim entry covering them."""
    if not manuscript_path.exists():
        return []

    claim_texts = {c.get("claim_text", "")[:60] for c in claims}
    uncovered = []

    with open(manuscript_path) as f:
        for line_num, line in enumerate(f, 1):
            text = line.strip()
            if not text or text.startswith("#") or text.startswith("!") or text.startswith("---"):
                continue
            if text.startswith("**Figure"):
                continue

            # Check if line has a quantitative assertion
            has_number = False
            for pat in _NUM_PATTERNS:
                if pat.search(text):
                    has_number = True
                    break

            if not has_number:
                continue

            # Check if any claim covers this line
            covered = any(
                ct and (ct[:40] in text or text[:40] in ct)
                for ct in claim_texts
            )
            if not covered:
                uncovered.append({
                    "line_number": line_num,
                    "text": text[:200],
                })

    return uncovered


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_review_report(
    verifications: list[dict],
    uncovered: list[dict],
    totals: dict,
    dataset_dir: Path,
    had_claims_file: bool,
) -> str:
    """Generate REVIEW_REPORT.md from verification results."""

    confirmed = [v for v in verifications if v["verification_status"] == "CONFIRMED"]
    discrepant = [v for v in verifications if v["verification_status"] == "DISCREPANT"]
    unverifiable = [v for v in verifications if v["verification_status"] == "UNVERIFIABLE"]
    needs_context = [v for v in verifications if v["verification_status"] == "NEEDS_CONTEXT"]

    # Collect all error patterns
    all_patterns = []
    for v in verifications:
        if v.get("error_patterns"):
            for p in v["error_patterns"]:
                p["claim_id"] = v["claim_id"]
                all_patterns.append(p)
    warnings = [p for p in all_patterns if p["severity"] == "warning"]

    lines = [
        "# Claim Verification Report",
        "",
        f"**Dataset:** `{dataset_dir}`",
        f"**Claims source:** {'MANUSCRIPT_CLAIMS.jsonl' if had_claims_file else 'Auto-extracted from MANUSCRIPT.md'}",
        f"**Database totals:** {totals.get('n_genomes', '?')} genomes, "
        f"{totals.get('n_proteins', '?'):,} proteins" if totals.get('n_proteins') else "",
        "",
        "## Summary",
        "",
        f"| Status | Count |",
        f"|--------|-------|",
        f"| CONFIRMED | {len(confirmed)} |",
        f"| DISCREPANT | {len(discrepant)} |",
        f"| NEEDS_CONTEXT | {len(needs_context)} |",
        f"| UNVERIFIABLE | {len(unverifiable)} |",
        f"| **Total claims** | **{len(verifications)}** |",
        "",
    ]

    if verifications:
        confirmed_pct = len(confirmed) / len(verifications) * 100
        lines.append(
            f"**Verification rate:** {len(confirmed)}/{len(verifications)} "
            f"({confirmed_pct:.0f}%) claims confirmed by database query."
        )
        lines.append("")

    # Discrepancies table (the critical section)
    if discrepant:
        lines.append("## Discrepancies")
        lines.append("")
        lines.append("*These claims do not match database values. Review and correct.*")
        lines.append("")
        lines.append("| Claim ID | Claimed | Verified | Note |")
        lines.append("|----------|---------|----------|------|")
        for v in discrepant:
            cid = v["claim_id"]
            claimed = v.get("claimed_value", "?")
            verified = v.get("verified_value", "?")
            note = (v.get("note") or "").replace("|", "/")
            lines.append(f"| {cid} | {claimed} | {verified} | {note} |")
        lines.append("")

        # Discrepancy details
        lines.append("### Discrepancy Details")
        lines.append("")
        for v in discrepant:
            lines.append(f"**{v['claim_id']}**: {v['claim_text'][:120]}")
            lines.append(f"- Claimed: {v.get('claimed_value', '?')}")
            lines.append(f"- Verified: {v.get('verified_value', '?')}")
            if v.get("verification_query"):
                lines.append(f"- Query: `{v['verification_query']}`")
            lines.append("")

    # Error patterns
    if warnings:
        lines.append("## Error Patterns Detected")
        lines.append("")
        lines.append("| Claim | Pattern | Detail |")
        lines.append("|-------|---------|--------|")
        for p in warnings:
            detail = p["detail"].replace("|", "/")
            lines.append(f"| {p['claim_id']} | {p['pattern']} | {detail} |")
        lines.append("")

    # Needs context
    if needs_context:
        lines.append("## Claims Needing Context")
        lines.append("")
        lines.append("*Query returned a value but could not be automatically compared to the claim.*")
        lines.append("")
        for v in needs_context:
            lines.append(f"- **{v['claim_id']}**: {v['claim_text'][:100]}")
            if v.get("note"):
                lines.append(f"  - {v['note']}")
        lines.append("")

    # Unverifiable claims
    if unverifiable:
        lines.append("## Unverifiable Claims")
        lines.append("")
        lines.append(
            f"*{len(unverifiable)} claims could not be verified by automated query. "
            f"These require interpretive review by the `/reviewer_2` agent.*"
        )
        lines.append("")
        # Show first 20
        for v in unverifiable[:20]:
            ctype = v.get("claim_type", "?")
            lines.append(f"- **{v['claim_id']}** ({ctype}): {v['claim_text'][:100]}")
        if len(unverifiable) > 20:
            lines.append(f"- *... and {len(unverifiable) - 20} more*")
        lines.append("")

    # Uncovered manuscript statements
    if uncovered:
        lines.append("## Uncovered Quantitative Statements")
        lines.append("")
        lines.append(
            f"*{len(uncovered)} lines in MANUSCRIPT.md contain numbers but have no "
            f"corresponding claim entry. Consider adding claim entries for traceability.*"
        )
        lines.append("")
        lines.append("| Line | Text |")
        lines.append("|------|------|")
        for item in uncovered[:30]:
            text = item["text"].replace("|", "\\|")
            lines.append(f"| {item['line_number']} | {text[:120]} |")
        if len(uncovered) > 30:
            lines.append(f"| ... | *{len(uncovered) - 30} more* |")
        lines.append("")

    # All confirmed claims (for completeness)
    if confirmed:
        lines.append("## Confirmed Claims")
        lines.append("")
        lines.append(f"*{len(confirmed)} claims verified successfully.*")
        lines.append("")
        for v in confirmed:
            lines.append(f"- **{v['claim_id']}**: {v['claim_text'][:80]} — {v.get('note', 'OK')}")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def verify(
    dataset_dir: Path,
    auto_extract: bool = False,
    claims_file: str | None = None,
    manuscript_file: str | None = None,
) -> tuple[list[dict], str]:
    """Run claim verification on a dataset. Returns (verifications, report)."""

    # Connect to database
    conn = connect_db(dataset_dir)
    if conn is None:
        return [], "# Error\n\nDatabase not found."

    totals = get_dataset_totals(conn)

    # Resolve file paths — support non-standard naming (e.g. BATHYARCHAEIA_MANUSCRIPT.md)
    if claims_file:
        claims_path = dataset_dir / claims_file
    else:
        claims_path = dataset_dir / "MANUSCRIPT_CLAIMS.jsonl"
        # Auto-detect prefixed claims files if standard name missing
        if not claims_path.exists():
            candidates = list(dataset_dir.glob("*_MANUSCRIPT_CLAIMS.jsonl"))
            if candidates:
                claims_path = candidates[0]

    if manuscript_file:
        manuscript_path = dataset_dir / manuscript_file
    else:
        manuscript_path = dataset_dir / "MANUSCRIPT.md"
        if not manuscript_path.exists():
            candidates = list(dataset_dir.glob("*_MANUSCRIPT.md"))
            if candidates:
                manuscript_path = candidates[0]

    had_claims_file = claims_path.exists() and load_jsonl(claims_path)

    if had_claims_file:
        claims = load_jsonl(claims_path)
        # Normalize field names (same as validate_provenance.py)
        for claim in claims:
            if "id" in claim and "claim_id" not in claim:
                claim["claim_id"] = claim["id"]
            if "claim" in claim and "claim_text" not in claim:
                claim["claim_text"] = claim["claim"]
            if "type" in claim and "claim_type" not in claim:
                claim["claim_type"] = claim["type"]
        print(f"Loaded {len(claims)} claims from {claims_path.name}")
    elif auto_extract or not claims_path.exists():
        claims = extract_claims_from_manuscript(manuscript_path)
        print(f"Auto-extracted {len(claims)} quantitative claims from MANUSCRIPT.md")
    else:
        claims = []
        print("No MANUSCRIPT_CLAIMS.jsonl found and --auto-extract not set.")

    if not claims:
        conn.close()
        return [], "# Claim Verification Report\n\nNo claims to verify."

    # Verify each claim
    verifications = []
    for claim in claims:
        record = verify_claim(claim, conn, totals)
        verifications.append(record)

    # Find uncovered statements
    uncovered = find_uncovered_statements(manuscript_path, claims)

    # Generate report
    report = generate_review_report(
        verifications, uncovered, totals, dataset_dir, bool(had_claims_file)
    )

    conn.close()
    return verifications, report


def main():
    parser = argparse.ArgumentParser(
        description="Verify manuscript claims against the database."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Path to the dataset directory",
    )
    parser.add_argument(
        "--auto-extract",
        action="store_true",
        help="Auto-extract claims from MANUSCRIPT.md if no MANUSCRIPT_CLAIMS.jsonl exists",
    )
    parser.add_argument(
        "--verification-output",
        default="CLAIM_VERIFICATION.jsonl",
        help="Output JSONL filename (default: CLAIM_VERIFICATION.jsonl)",
    )
    parser.add_argument(
        "--report-output",
        default="REVIEW_REPORT.md",
        help="Output report filename (default: REVIEW_REPORT.md)",
    )
    parser.add_argument(
        "--claims-file",
        default=None,
        help="Claims JSONL filename (default: auto-detect MANUSCRIPT_CLAIMS.jsonl or *_MANUSCRIPT_CLAIMS.jsonl)",
    )
    parser.add_argument(
        "--manuscript-file",
        default=None,
        help="Manuscript filename (default: auto-detect MANUSCRIPT.md or *_MANUSCRIPT.md)",
    )
    args = parser.parse_args()

    dataset_dir = Path(args.dataset)
    if not dataset_dir.exists():
        print(f"Error: dataset directory not found: {dataset_dir}", file=sys.stderr)
        return 1

    verifications, report = verify(
        dataset_dir,
        auto_extract=args.auto_extract,
        claims_file=args.claims_file,
        manuscript_file=args.manuscript_file,
    )

    # Write outputs
    verification_path = dataset_dir / args.verification_output
    with open(verification_path, "w") as f:
        for v in verifications:
            f.write(json.dumps(v) + "\n")
    print(f"Wrote {len(verifications)} verification records to {verification_path}")

    report_path = dataset_dir / args.report_output
    report_path.write_text(report)
    print(f"Wrote review report to {report_path}")

    # Print summary
    confirmed = sum(1 for v in verifications if v["verification_status"] == "CONFIRMED")
    discrepant = sum(1 for v in verifications if v["verification_status"] == "DISCREPANT")
    unverifiable = sum(1 for v in verifications if v["verification_status"] == "UNVERIFIABLE")
    needs_context = sum(1 for v in verifications if v["verification_status"] == "NEEDS_CONTEXT")

    print(f"\nResults: {confirmed} confirmed, {discrepant} discrepant, "
          f"{needs_context} needs context, {unverifiable} unverifiable")

    if discrepant > 0:
        print(f"\n!! {discrepant} DISCREPANCIES found — review REVIEW_REPORT.md")
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
