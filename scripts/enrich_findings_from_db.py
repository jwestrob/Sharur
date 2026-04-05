#!/usr/bin/env python3
"""Bootstrap verification blocks and representative protein IDs from DuckDB."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb


@dataclass(frozen=True)
class VerificationSpec:
    claim: str
    query: str
    expected_ref: Any


@dataclass(frozen=True)
class FindingRecipe:
    verification_specs: tuple[VerificationSpec, ...]
    protein_ids_query: str | None = None


def _load_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    if not path.exists():
        return records
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def _write_jsonl(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, default=str) + "\n")


def _scalar(con: duckdb.DuckDBPyConnection, query: str) -> Any:
    row = con.execute(query).fetchone()
    return row[0] if row else None


def _string_list(
    con: duckdb.DuckDBPyConnection,
    query: str,
) -> list[str]:
    return [str(row[0]) for row in con.execute(query).fetchall() if row and row[0]]


def _merge_verification(
    existing: list[dict[str, Any]] | None,
    new_entries: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    merged: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()

    for entry in existing or []:
        claim = str(entry.get("claim") or "")
        query = str(entry.get("query") or "")
        key = (claim, query)
        if key not in seen:
            seen.add(key)
            merged.append(entry)

    for entry in new_entries:
        key = (entry["claim"], entry["query"])
        if key not in seen:
            seen.add(key)
            merged.append(entry)

    return merged


def _merge_string_list(
    existing: list[str] | None,
    new_values: list[str],
) -> list[str]:
    merged: list[str] = []
    seen: set[str] = set()

    for value in (existing or []) + new_values:
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            merged.append(text)
    return merged


def _contigs_for_protein_ids(
    con: duckdb.DuckDBPyConnection,
    protein_ids: list[str],
) -> list[str]:
    if not protein_ids:
        return []
    placeholders = ", ".join("?" for _ in protein_ids)
    rows = con.execute(
        f"""
        SELECT DISTINCT contig_id
        FROM proteins
        WHERE protein_id IN ({placeholders})
        ORDER BY contig_id
        """,
        protein_ids,
    ).fetchall()
    return [str(row[0]) for row in rows if row and row[0]]


def enrich_record_with_recipe(
    con: duckdb.DuckDBPyConnection,
    record: dict[str, Any],
    recipe: FindingRecipe,
) -> tuple[dict[str, Any], dict[str, Any]]:
    updated = dict(record)
    added_verification: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []

    for spec in recipe.verification_specs:
        actual = _scalar(con, spec.query)
        if actual != spec.expected_ref:
            mismatches.append(
                {
                    "claim": spec.claim,
                    "expected_ref": spec.expected_ref,
                    "actual": actual,
                }
            )
            continue
        added_verification.append(
            {
                "claim": spec.claim,
                "query": spec.query,
                "expected": actual,
            }
        )

    if added_verification:
        updated["verification"] = _merge_verification(
            updated.get("verification"),
            added_verification,
        )

    added_protein_ids: list[str] = []
    if recipe.protein_ids_query:
        added_protein_ids = _string_list(con, recipe.protein_ids_query)
        if added_protein_ids:
            updated["protein_ids"] = _merge_string_list(
                updated.get("protein_ids"),
                added_protein_ids,
            )
            updated["contigs"] = _merge_string_list(
                updated.get("contigs"),
                _contigs_for_protein_ids(con, added_protein_ids),
            )

    changed = (
        updated.get("verification") != record.get("verification")
        or updated.get("protein_ids") != record.get("protein_ids")
        or updated.get("contigs") != record.get("contigs")
    )

    return updated, {
        "changed": changed,
        "added_verification": len(added_verification),
        "added_protein_ids": len(added_protein_ids),
        "mismatches": mismatches,
    }


def omni_production_recipes() -> dict[str, FindingRecipe]:
    return {
        "survey-001": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Group 4 NiFe hydrogenases are present in 2,990 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT pp.protein_id)
                    FROM protein_predicates pp,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'nife_group4'
                    """,
                    expected_ref=2990,
                ),
                VerificationSpec(
                    claim="Group 4 NiFe hydrogenases occur in 1,366 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'nife_group4'
                    """,
                    expected_ref=1366,
                ),
                VerificationSpec(
                    claim="Group 4 NiFe hydrogenases occur in 74.6% of bins.",
                    query="""
                    SELECT ROUND(
                        100.0 * COUNT(DISTINCT p.bin_id) / (SELECT COUNT(*) FROM bins),
                        1
                    )
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'nife_group4'
                    """,
                    expected_ref=74.6,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM protein_predicates pp
            JOIN proteins p ON pp.protein_id = p.protein_id,
                 UNNEST(pp.predicates) AS t(pred)
            WHERE pred = 'nife_group4'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-002": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Glycosyltransferase-associated PFAM hits occur in 68,575 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND (a.name LIKE 'Glycos_transf%' OR a.name LIKE 'Glyco_trans%')
                    """,
                    expected_ref=68575,
                ),
                VerificationSpec(
                    claim="Glycosyltransferase-associated PFAM hits are present in all 1,831 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND (a.name LIKE 'Glycos_transf%' OR a.name LIKE 'Glyco_trans%')
                    """,
                    expected_ref=1831,
                ),
                VerificationSpec(
                    claim="Glycosyltransferase-associated PFAM hits account for 2.3% of proteins in the dataset.",
                    query="""
                    SELECT ROUND(
                        100.0 * COUNT(DISTINCT a.protein_id) / (SELECT COUNT(*) FROM proteins),
                        1
                    )
                    FROM annotations a
                    WHERE a.source = 'pfam'
                      AND (a.name LIKE 'Glycos_transf%' OR a.name LIKE 'Glyco_trans%')
                    """,
                    expected_ref=2.3,
                ),
                VerificationSpec(
                    claim="Omnitrophota average 37.5 glycosyltransferase-associated PFAM hits per genome.",
                    query="""
                    SELECT ROUND(
                        COUNT(DISTINCT a.protein_id) * 1.0 / (SELECT COUNT(*) FROM bins),
                        1
                    )
                    FROM annotations a
                    WHERE a.source = 'pfam'
                      AND (a.name LIKE 'Glycos_transf%' OR a.name LIKE 'Glyco_trans%')
                    """,
                    expected_ref=37.5,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'pfam'
              AND (a.name LIKE 'Glycos_transf%' OR a.name LIKE 'Glyco_trans%')
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-003": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="The largest protein in the dataset is 85,804 aa long.",
                    query="SELECT MAX(sequence_length) FROM proteins",
                    expected_ref=85804,
                ),
                VerificationSpec(
                    claim="There are 2,034 proteins longer than 5,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT protein_id)
                    FROM proteins
                    WHERE sequence_length > 5000
                    """,
                    expected_ref=2034,
                ),
                VerificationSpec(
                    claim="Proteins longer than 5,000 aa occur in 792 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT bin_id)
                    FROM proteins
                    WHERE sequence_length > 5000
                    """,
                    expected_ref=792,
                ),
            ),
            protein_ids_query="""
            SELECT protein_id
            FROM proteins
            ORDER BY sequence_length DESC, protein_id
            LIMIT 5
            """,
        ),
        "survey-004": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="PF00990 GGDEF domains occur in 11,773 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF00990'
                    """,
                    expected_ref=11773,
                ),
                VerificationSpec(
                    claim="PF07238 PilZ domains occur in 16,324 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF07238'
                    """,
                    expected_ref=16324,
                ),
                VerificationSpec(
                    claim="PF07238 PilZ domains occur in 1,801 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF07238'
                    """,
                    expected_ref=1801,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'pfam'
              AND a.accession = 'PF07238'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-005": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="At least five Rnf subunits co-occur in 493 genomes.",
                    query="""
                    WITH rnf AS (
                      SELECT p.bin_id, COUNT(DISTINCT a.accession) AS n_rnf
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN (
                        'K03613','K03614','K03615','K03616','K03617','K03618'
                      )
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM rnf WHERE n_rnf >= 5
                    """,
                    expected_ref=493,
                ),
            ),
            protein_ids_query="""
            WITH complete AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.accession IN (
                'K03613','K03614','K03615','K03616','K03617','K03618'
              )
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) >= 5
              ORDER BY p.bin_id
              LIMIT 1
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM complete)
              AND a.accession IN (
                'K03613','K03614','K03615','K03616','K03617','K03618'
              )
            ORDER BY p.gene_index, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-006": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="PF04055 Radical SAM occurs in 46,646 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF04055'
                    """,
                    expected_ref=46646,
                ),
                VerificationSpec(
                    claim="PF04055 Radical SAM occurs in all 1,831 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF04055'
                    """,
                    expected_ref=1831,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'pfam'
              AND a.accession = 'PF04055'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-007": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Flagellin K02406 occurs in 168 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K02406'
                    """,
                    expected_ref=168,
                ),
                VerificationSpec(
                    claim="Flagellin K02406 is detected in 66 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K02406'
                    """,
                    expected_ref=66,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.accession = 'K02406'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-009": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="There are 194 near-identical genome pairs under the duplicate-pair heuristic.",
                    query="""
                    SELECT COUNT(*) FROM (
                      SELECT b1.bin_id, b2.bin_id
                      FROM bins b1
                      JOIN bins b2
                        ON b1.bin_id < b2.bin_id
                       AND ABS(b1.total_length - b2.total_length) < 1000
                       AND b1.n_contigs = b2.n_contigs
                    ) dup
                    """,
                    expected_ref=194,
                ),
            ),
        ),
        "survey-010": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Wood-Ljungdahl core markers (CODH + ACS + MTHFR) co-occur in 174 genomes.",
                    query="""
                    WITH wl AS (
                      SELECT p.bin_id,
                             BOOL_OR(a.accession = 'K00198') AS has_codh,
                             BOOL_OR(a.accession = 'K14138') AS has_acs,
                             BOOL_OR(a.accession = 'K00297') AS has_mthfr
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN ('K00198','K14138','K00297')
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM wl
                    WHERE has_codh AND has_acs AND has_mthfr
                    """,
                    expected_ref=174,
                ),
                VerificationSpec(
                    claim="Wood-Ljungdahl-complete genomes make up 9.5% of bins.",
                    query="""
                    WITH wl AS (
                      SELECT p.bin_id,
                             BOOL_OR(a.accession = 'K00198') AS has_codh,
                             BOOL_OR(a.accession = 'K14138') AS has_acs,
                             BOOL_OR(a.accession = 'K00297') AS has_mthfr
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN ('K00198','K14138','K00297')
                      GROUP BY p.bin_id
                    )
                    SELECT ROUND(
                      100.0 * COUNT(*) / (SELECT COUNT(*) FROM bins),
                      1
                    )
                    FROM wl
                    WHERE has_codh AND has_acs AND has_mthfr
                    """,
                    expected_ref=9.5,
                ),
                VerificationSpec(
                    claim="K00297 methylene-THF reductase is present in 1,285 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K00297'
                    """,
                    expected_ref=1285,
                ),
            ),
            protein_ids_query="""
            WITH complete AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.accession IN ('K00198','K14138','K00297')
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) = 3
              ORDER BY p.bin_id
              LIMIT 1
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM complete)
              AND a.accession IN ('K00198','K14138','K00297')
            ORDER BY p.gene_index, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-011": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Rubisco predicates occur in 55 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'rubisco'
                    """,
                    expected_ref=55,
                ),
                VerificationSpec(
                    claim="Rubisco predicates occur in 58 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT pp.protein_id)
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'rubisco'
                    """,
                    expected_ref=58,
                ),
                VerificationSpec(
                    claim="No proteins carry the prk predicate.",
                    query="""
                    SELECT COUNT(DISTINCT pp.protein_id)
                    FROM protein_predicates pp,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'prk'
                    """,
                    expected_ref=0,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM protein_predicates pp
            JOIN proteins p ON pp.protein_id = p.protein_id,
                 UNNEST(pp.predicates) AS t(pred)
            WHERE pred = 'rubisco'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-012": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="DUF1015 PF06245 occurs in 1,629 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF06245'
                    """,
                    expected_ref=1629,
                ),
                VerificationSpec(
                    claim="DUF1015 PF06245 occurs in 1,557 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF06245'
                    """,
                    expected_ref=1557,
                ),
                VerificationSpec(
                    claim="DUF1015 PF06245 occurs in 85.0% of genomes.",
                    query="""
                    SELECT ROUND(
                      100.0 * COUNT(DISTINCT p.bin_id) / (SELECT COUNT(*) FROM bins),
                      1
                    )
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF06245'
                    """,
                    expected_ref=85.0,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.accession = 'PF06245'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-013": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="There are 13,422 unannotated proteins longer than 1,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT p.protein_id)
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'unannotated'
                      AND p.sequence_length > 1000
                    """,
                    expected_ref=13422,
                ),
                VerificationSpec(
                    claim="The longest unannotated protein is 39,880 aa.",
                    query="""
                    SELECT MAX(p.sequence_length)
                    FROM protein_predicates pp
                    JOIN proteins p ON pp.protein_id = p.protein_id,
                         UNNEST(pp.predicates) AS t(pred)
                    WHERE pred = 'unannotated'
                    """,
                    expected_ref=39880,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM protein_predicates pp
            JOIN proteins p ON pp.protein_id = p.protein_id,
                 UNNEST(pp.predicates) AS t(pred)
            WHERE pred = 'unannotated'
              AND p.sequence_length > 1000
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-014": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="At least five F-type ATPase subunits co-occur in 842 genomes.",
                    query="""
                    WITH atp AS (
                      SELECT p.bin_id,
                             COUNT(DISTINCT a.accession) FILTER (
                               WHERE a.accession IN (
                                 'K02108','K02109','K02110','K02111',
                                 'K02112','K02113','K02114','K02115'
                               )
                             ) AS n_f
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN (
                        'K02108','K02109','K02110','K02111',
                        'K02112','K02113','K02114','K02115'
                      )
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM atp WHERE n_f >= 5
                    """,
                    expected_ref=842,
                ),
                VerificationSpec(
                    claim="At least four V/A-type ATPase subunits co-occur in 191 genomes.",
                    query="""
                    WITH atp AS (
                      SELECT p.bin_id,
                             COUNT(DISTINCT a.accession) FILTER (
                               WHERE a.accession IN (
                                 'K02117','K02118','K02119','K02120',
                                 'K02121','K02122','K02123','K02124'
                               )
                             ) AS n_v
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN (
                        'K02117','K02118','K02119','K02120',
                        'K02121','K02122','K02123','K02124'
                      )
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM atp WHERE n_v >= 4
                    """,
                    expected_ref=191,
                ),
            ),
            protein_ids_query="""
            WITH complete AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.accession IN (
                'K02108','K02109','K02110','K02111',
                'K02112','K02113','K02114','K02115'
              )
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) >= 5
              ORDER BY p.bin_id
              LIMIT 1
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM complete)
              AND a.accession IN (
                'K02108','K02109','K02110','K02111',
                'K02112','K02113','K02114','K02115'
              )
            ORDER BY p.gene_index, p.protein_id
            LIMIT 5
            """,
        ),
        "survey-015": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="NifH, NifD, and NifK co-occur in 11 genomes.",
                    query="""
                    WITH nif AS (
                      SELECT p.bin_id,
                             BOOL_OR(a.accession = 'K02588') AS has_h,
                             BOOL_OR(a.accession = 'K02586') AS has_d,
                             BOOL_OR(a.accession = 'K02591') AS has_k
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN ('K02588','K02586','K02591')
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM nif
                    WHERE has_h AND has_d AND has_k
                    """,
                    expected_ref=11,
                ),
            ),
            protein_ids_query="""
            WITH complete AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.accession IN ('K02588','K02586','K02591')
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) = 3
              ORDER BY p.bin_id
              LIMIT 1
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM complete)
              AND a.accession IN ('K02588','K02586','K02591')
            ORDER BY p.gene_index, p.protein_id
            LIMIT 5
            """,
        ),
        "E001": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Three PacBio bins are single-contig genomes.",
                    query="""
                    SELECT COUNT(*) FROM (
                      SELECT bin_id
                      FROM proteins
                      WHERE bin_id LIKE '%PACBIO%'
                      GROUP BY bin_id
                      HAVING COUNT(DISTINCT contig_id) = 1
                    ) pacbio
                    """,
                    expected_ref=3,
                ),
                VerificationSpec(
                    claim="The largest protein in the PacBio genomes is 78,142 aa.",
                    query="""
                    SELECT MAX(sequence_length)
                    FROM proteins
                    WHERE bin_id LIKE '%PACBIO%'
                    """,
                    expected_ref=78142,
                ),
            ),
            protein_ids_query="""
            SELECT protein_id
            FROM proteins
            WHERE bin_id LIKE '%PACBIO%'
            ORDER BY sequence_length DESC, protein_id
            LIMIT 5
            """,
        ),
        "E003": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="GGDEF domains occur in 52 proteins longer than 5,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF00990'
                      AND p.sequence_length > 5000
                    """,
                    expected_ref=52,
                ),
                VerificationSpec(
                    claim="GGDEF domains in proteins longer than 5,000 aa occur in 50 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF00990'
                      AND p.sequence_length > 5000
                    """,
                    expected_ref=50,
                ),
                VerificationSpec(
                    claim="Fifteen GGDEF proteins are longer than 20,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF00990'
                      AND p.sequence_length > 20000
                    """,
                    expected_ref=15,
                ),
                VerificationSpec(
                    claim="Peptidase_C39-family domains occur in 113 proteins longer than 5,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession IN ('PF03412','PF13529')
                      AND p.sequence_length > 5000
                    """,
                    expected_ref=113,
                ),
                VerificationSpec(
                    claim="Peptidase_C39-family domains in proteins longer than 5,000 aa occur in 112 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.source = 'pfam'
                      AND a.accession IN ('PF03412','PF13529')
                      AND p.sequence_length > 5000
                    """,
                    expected_ref=112,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'pfam'
              AND a.accession = 'PF00990'
              AND p.sequence_length > 5000
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "E004": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="TIM+PGK glycolytic fusion proteins occur in 82 proteins.",
                    query="""
                    SELECT COUNT(*)
                    FROM (
                      SELECT a.protein_id
                      FROM annotations a
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    ) tim_pgk
                    """,
                    expected_ref=82,
                ),
                VerificationSpec(
                    claim="TIM+PGK glycolytic fusion proteins occur in 81 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT bin_id)
                    FROM (
                      SELECT a.protein_id,
                             MAX(p.bin_id) AS bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    ) tim_pgk
                    """,
                    expected_ref=81,
                ),
                VerificationSpec(
                    claim="The longest TIM+PGK glycolytic fusion protein is 15,058 aa.",
                    query="""
                    SELECT MAX(sequence_length)
                    FROM proteins
                    WHERE protein_id IN (
                      SELECT a.protein_id
                      FROM annotations a
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    )
                    """,
                    expected_ref=15058,
                ),
                VerificationSpec(
                    claim="Twenty-nine TIM+PGK glycolytic fusion proteins also carry Gp_dh_C domains.",
                    query="""
                    WITH tim_pgk AS (
                      SELECT a.protein_id
                      FROM annotations a
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN tim_pgk t ON a.protein_id = t.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'Gp_dh_C'
                    """,
                    expected_ref=29,
                ),
                VerificationSpec(
                    claim="Twenty-four TIM+PGK glycolytic fusion proteins also carry GPI domains.",
                    query="""
                    WITH tim_pgk AS (
                      SELECT a.protein_id
                      FROM annotations a
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN tim_pgk t ON a.protein_id = t.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'GPI'
                    """,
                    expected_ref=24,
                ),
                VerificationSpec(
                    claim="Thirteen TIM+PGK glycolytic fusion proteins also carry F_bP_aldolase domains.",
                    query="""
                    WITH tim_pgk AS (
                      SELECT a.protein_id
                      FROM annotations a
                      WHERE a.source = 'pfam'
                        AND a.name IN ('TIM', 'PGK')
                      GROUP BY a.protein_id
                      HAVING COUNT(DISTINCT a.name) = 2
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN tim_pgk t ON a.protein_id = t.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'F_bP_aldolase'
                    """,
                    expected_ref=13,
                ),
            ),
            protein_ids_query="""
            WITH glycolytic_fusions AS (
              SELECT a.protein_id
              FROM annotations a
              WHERE a.source = 'pfam'
                AND a.name IN ('TIM', 'PGK')
              GROUP BY a.protein_id
              HAVING COUNT(DISTINCT a.name) = 2
            )
            SELECT p.protein_id
            FROM proteins p
            JOIN glycolytic_fusions gf ON p.protein_id = gf.protein_id
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "E007": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Complete V/A-type ATPases (>=4 KEGG subunits) occur in 191 genomes.",
                    query="""
                    SELECT COUNT(*)
                    FROM (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    ) va_bins
                    """,
                    expected_ref=191,
                ),
                VerificationSpec(
                    claim="Complete F-type ATPases (>=5 KEGG subunits) occur in 842 genomes.",
                    query="""
                    SELECT COUNT(*)
                    FROM (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ) f_bins
                    """,
                    expected_ref=842,
                ),
                VerificationSpec(
                    claim="Twenty-six genomes encode both F-type and V/A-type ATPases.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    )
                    SELECT COUNT(*) FROM f_bins
                    WHERE bin_id IN (SELECT bin_id FROM va_bins)
                    """,
                    expected_ref=26,
                ),
                VerificationSpec(
                    claim="Eight hundred twenty-four genomes encode neither complete F-type nor complete V/A-type ATPases.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    )
                    SELECT COUNT(*)
                    FROM bins
                    WHERE bin_id NOT IN (SELECT bin_id FROM f_bins)
                      AND bin_id NOT IN (SELECT bin_id FROM va_bins)
                    """,
                    expected_ref=824,
                ),
                VerificationSpec(
                    claim="V/A-type-only genomes average 1,594 proteins.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    ),
                    stats AS (
                      SELECT b.bin_id,
                             COUNT(p.protein_id) AS protein_count
                      FROM bins b
                      LEFT JOIN proteins p ON b.bin_id = p.bin_id
                      GROUP BY b.bin_id
                    )
                    SELECT ROUND(AVG(protein_count), 0)
                    FROM stats
                    WHERE bin_id IN (SELECT bin_id FROM va_bins)
                      AND bin_id NOT IN (SELECT bin_id FROM f_bins)
                    """,
                    expected_ref=1594.0,
                ),
                VerificationSpec(
                    claim="F-type-only genomes average 1,770 proteins.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    ),
                    stats AS (
                      SELECT b.bin_id,
                             COUNT(p.protein_id) AS protein_count
                      FROM bins b
                      LEFT JOIN proteins p ON b.bin_id = p.bin_id
                      GROUP BY b.bin_id
                    )
                    SELECT ROUND(AVG(protein_count), 0)
                    FROM stats
                    WHERE bin_id IN (SELECT bin_id FROM f_bins)
                      AND bin_id NOT IN (SELECT bin_id FROM va_bins)
                    """,
                    expected_ref=1770.0,
                ),
                VerificationSpec(
                    claim="V/A-type-only genomes average 197 contigs.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    )
                    SELECT ROUND(AVG(n_contigs), 0)
                    FROM bins
                    WHERE bin_id IN (SELECT bin_id FROM va_bins)
                      AND bin_id NOT IN (SELECT bin_id FROM f_bins)
                    """,
                    expected_ref=197.0,
                ),
                VerificationSpec(
                    claim="F-type-only genomes average 180 contigs.",
                    query="""
                    WITH f_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02111','K02112','K02113','K02114',
                          'K02115','K02108','K02109','K02110'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 5
                    ),
                    va_bins AS (
                      SELECT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession IN (
                          'K02117','K02118','K02119','K02120',
                          'K02121','K02122','K02123','K02124'
                        )
                      GROUP BY p.bin_id
                      HAVING COUNT(DISTINCT a.accession) >= 4
                    )
                    SELECT ROUND(AVG(n_contigs), 0)
                    FROM bins
                    WHERE bin_id IN (SELECT bin_id FROM f_bins)
                      AND bin_id NOT IN (SELECT bin_id FROM va_bins)
                    """,
                    expected_ref=180.0,
                ),
            ),
            protein_ids_query="""
            WITH va_bins AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.source = 'kegg'
                AND a.accession IN (
                  'K02117','K02118','K02119','K02120',
                  'K02121','K02122','K02123','K02124'
                )
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) >= 4
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM va_bins)
              AND a.accession IN (
                'K02117','K02118','K02119','K02120',
                'K02121','K02122','K02123','K02124'
              )
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "E009": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="PF00990 GGDEF domains occur in 11,773 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    WHERE a.source = 'pfam'
                      AND a.accession = 'PF00990'
                    """,
                    expected_ref=11773,
                ),
                VerificationSpec(
                    claim="3,492 GGDEF proteins carry GAF domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'GAF'
                    """,
                    expected_ref=3492,
                ),
                VerificationSpec(
                    claim="2,560 GGDEF proteins carry Response_reg domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'Response_reg'
                    """,
                    expected_ref=2560,
                ),
                VerificationSpec(
                    claim="936 GGDEF proteins carry HAMP domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'HAMP'
                    """,
                    expected_ref=936,
                ),
                VerificationSpec(
                    claim="858 GGDEF proteins carry HD domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'HD'
                    """,
                    expected_ref=858,
                ),
                VerificationSpec(
                    claim="435 GGDEF proteins carry EAL domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'EAL'
                    """,
                    expected_ref=435,
                ),
                VerificationSpec(
                    claim="515 GGDEF proteins carry CHASE2 domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'CHASE2'
                    """,
                    expected_ref=515,
                ),
                VerificationSpec(
                    claim="392 GGDEF proteins carry HisKA domains.",
                    query="""
                    WITH ggdef AS (
                      SELECT DISTINCT protein_id
                      FROM annotations
                      WHERE source = 'pfam'
                        AND accession = 'PF00990'
                    )
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN ggdef g ON a.protein_id = g.protein_id
                    WHERE a.source = 'pfam'
                      AND a.name = 'HisKA'
                    """,
                    expected_ref=392,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'pfam'
              AND a.accession = 'PF00990'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "E015": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="There are 2,034 proteins longer than 5,000 aa.",
                    query="""
                    SELECT COUNT(DISTINCT protein_id)
                    FROM proteins
                    WHERE sequence_length > 5000
                    """,
                    expected_ref=2034,
                ),
                VerificationSpec(
                    claim="Proteins longer than 5,000 aa occur in 792 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT bin_id)
                    FROM proteins
                    WHERE sequence_length > 5000
                    """,
                    expected_ref=792,
                ),
                VerificationSpec(
                    claim="There are 6,547 proteins with at least 10 distinct PFAM domains.",
                    query="""
                    WITH pfam_counts AS (
                      SELECT a.protein_id,
                             COUNT(DISTINCT a.accession) AS n_pfam,
                             MAX(p.bin_id) AS bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'pfam'
                      GROUP BY a.protein_id
                    )
                    SELECT COUNT(*) FROM pfam_counts WHERE n_pfam >= 10
                    """,
                    expected_ref=6547,
                ),
                VerificationSpec(
                    claim="Proteins with at least 10 distinct PFAM domains occur in 1,672 genomes.",
                    query="""
                    WITH pfam_counts AS (
                      SELECT a.protein_id,
                             COUNT(DISTINCT a.accession) AS n_pfam,
                             MAX(p.bin_id) AS bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'pfam'
                      GROUP BY a.protein_id
                    )
                    SELECT COUNT(DISTINCT bin_id) FROM pfam_counts WHERE n_pfam >= 10
                    """,
                    expected_ref=1672,
                ),
            ),
            protein_ids_query="""
            SELECT protein_id
            FROM proteins
            WHERE sequence_length > 5000
            ORDER BY sequence_length DESC, protein_id
            LIMIT 5
            """,
        ),
        "E027": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="There are 533 candidate multi-heme cytochromes with >=10 CxxCH motifs and length >=200 aa.",
                    query="""
                    SELECT COUNT(*)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    """,
                    expected_ref=533,
                ),
                VerificationSpec(
                    claim="Candidate multi-heme cytochromes occur in 241 raw genomes before duplicate collapse.",
                    query="""
                    SELECT COUNT(DISTINCT bin_id)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    """,
                    expected_ref=241,
                ),
                VerificationSpec(
                    claim="Fifteen candidate multi-heme cytochromes have at least 30 CxxCH motifs.",
                    query="""
                    SELECT COUNT(*)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 30
                    """,
                    expected_ref=15,
                ),
                VerificationSpec(
                    claim="Forty-three candidate multi-heme cytochromes have at least 20 CxxCH motifs.",
                    query="""
                    SELECT COUNT(*)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 20
                    """,
                    expected_ref=43,
                ),
                VerificationSpec(
                    claim="One hundred thirty candidate multi-heme cytochromes have at least 15 CxxCH motifs.",
                    query="""
                    SELECT COUNT(*)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 15
                    """,
                    expected_ref=130,
                ),
                VerificationSpec(
                    claim="Candidate multi-heme cytochromes range from 229 aa to 2,149 aa with median length 627 aa.",
                    query="""
                    SELECT
                      CAST(MIN(sequence_length) AS VARCHAR) || '-' ||
                      CAST(MAX(sequence_length) AS VARCHAR) || ' median ' ||
                      CAST(CAST(MEDIAN(sequence_length) AS INTEGER) AS VARCHAR)
                    FROM proteins
                    WHERE sequence_length >= 200
                      AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    """,
                    expected_ref="229-2149 median 627",
                ),
                VerificationSpec(
                    claim="MtrC-MtrF domains occur in 63 multi-heme genomes, corresponding to a 34.6x enrichment over non-multi-heme genomes.",
                    query="""
                    WITH mh AS (
                      SELECT DISTINCT bin_id
                      FROM proteins
                      WHERE sequence_length >= 200
                        AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    ),
                    non_mh AS (
                      SELECT DISTINCT bin_id FROM proteins
                      EXCEPT
                      SELECT bin_id FROM mh
                    ),
                    dom AS (
                      SELECT DISTINCT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'pfam'
                        AND a.name LIKE 'Mtrc-MtrF%'
                    )
                    SELECT
                      CAST((SELECT COUNT(*) FROM mh WHERE bin_id IN (SELECT bin_id FROM dom)) AS VARCHAR) ||
                      ' multi-heme genomes; fold=' ||
                      CAST(
                        ROUND(
                          (100.0 * (SELECT COUNT(*) FROM mh WHERE bin_id IN (SELECT bin_id FROM dom)) / (SELECT COUNT(*) FROM mh)) /
                          NULLIF(
                            (100.0 * (SELECT COUNT(*) FROM non_mh WHERE bin_id IN (SELECT bin_id FROM dom)) / (SELECT COUNT(*) FROM non_mh)),
                            0
                          ),
                          1
                        ) AS VARCHAR
                      )
                    """,
                    expected_ref="63 multi-heme genomes; fold=34.6",
                ),
                VerificationSpec(
                    claim="GSu_C4xC__C2xCH domains are exclusive to multi-heme genomes (14 vs 0).",
                    query="""
                    WITH mh AS (
                      SELECT DISTINCT bin_id
                      FROM proteins
                      WHERE sequence_length >= 200
                        AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    ),
                    non_mh AS (
                      SELECT DISTINCT bin_id FROM proteins
                      EXCEPT
                      SELECT bin_id FROM mh
                    ),
                    gsu AS (
                      SELECT DISTINCT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'pfam'
                        AND a.name = 'GSu_C4xC__C2xCH'
                    )
                    SELECT
                      CAST((SELECT COUNT(*) FROM mh WHERE bin_id IN (SELECT bin_id FROM gsu)) AS VARCHAR) ||
                      ' vs ' ||
                      CAST((SELECT COUNT(*) FROM non_mh WHERE bin_id IN (SELECT bin_id FROM gsu)) AS VARCHAR)
                    """,
                    expected_ref="14 vs 0",
                ),
                VerificationSpec(
                    claim="coxI occurs in 60 multi-heme genomes, corresponding to 24.9% prevalence.",
                    query="""
                    WITH mh AS (
                      SELECT DISTINCT bin_id
                      FROM proteins
                      WHERE sequence_length >= 200
                        AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
                    ),
                    cox AS (
                      SELECT DISTINCT p.bin_id
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.source = 'kegg'
                        AND a.accession = 'K02274'
                    )
                    SELECT
                      CAST((SELECT COUNT(*) FROM mh WHERE bin_id IN (SELECT bin_id FROM cox)) AS VARCHAR) ||
                      ' genomes; pct=' ||
                      CAST(
                        ROUND(
                          100.0 * (SELECT COUNT(*) FROM mh WHERE bin_id IN (SELECT bin_id FROM cox)) /
                          (SELECT COUNT(*) FROM mh),
                          1
                        ) AS VARCHAR
                      )
                    """,
                    expected_ref="60 genomes; pct=24.9",
                ),
            ),
            protein_ids_query="""
            SELECT protein_id
            FROM proteins
            WHERE sequence_length >= 200
              AND array_length(regexp_extract_all(sequence, 'C..CH')) >= 10
            ORDER BY array_length(regexp_extract_all(sequence, 'C..CH')) DESC,
                     sequence_length DESC,
                     protein_id
            LIMIT 5
            """,
        ),
        "D008": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Cas1 and Cas2 co-occur in 156 genomes.",
                    query="""
                    WITH crispr AS (
                      SELECT p.bin_id,
                             BOOL_OR(a.accession = 'K15342') AS has_cas1,
                             BOOL_OR(a.accession = 'K09951') AS has_cas2
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN ('K15342','K09951')
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM crispr WHERE has_cas1 AND has_cas2
                    """,
                    expected_ref=156,
                ),
                VerificationSpec(
                    claim="Cas3 markers occur in 125 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K07012'
                    """,
                    expected_ref=125,
                ),
                VerificationSpec(
                    claim="Cas9 markers occur in 13 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K09952'
                    """,
                    expected_ref=13,
                ),
            ),
            protein_ids_query="""
            WITH complete AS (
              SELECT p.bin_id
              FROM annotations a
              JOIN proteins p ON a.protein_id = p.protein_id
              WHERE a.accession IN ('K15342','K09951')
              GROUP BY p.bin_id
              HAVING COUNT(DISTINCT a.accession) = 2
              ORDER BY p.bin_id
              LIMIT 1
            )
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id IN (SELECT bin_id FROM complete)
              AND a.accession IN ('K15342','K09951','K07012')
            ORDER BY p.gene_index, p.protein_id
            LIMIT 5
            """,
        ),
        "D031": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Rubrerythrin PF02915 occurs in 2,852 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF02915'
                    """,
                    expected_ref=2852,
                ),
                VerificationSpec(
                    claim="Rubrerythrin PF02915 occurs in 1,233 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF02915'
                    """,
                    expected_ref=1233,
                ),
                VerificationSpec(
                    claim="Catalase K03781 occurs in only 20 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'K03781'
                    """,
                    expected_ref=20,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.accession = 'PF02915'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
        "D032": FindingRecipe(
            verification_specs=(
                VerificationSpec(
                    claim="Membrane H_PPase PF03030 occurs in 1,870 proteins.",
                    query="""
                    SELECT COUNT(DISTINCT a.protein_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF03030'
                    """,
                    expected_ref=1870,
                ),
                VerificationSpec(
                    claim="Membrane H_PPase occurs in 1,548 genomes.",
                    query="""
                    SELECT COUNT(DISTINCT p.bin_id)
                    FROM annotations a
                    JOIN proteins p ON a.protein_id = p.protein_id
                    WHERE a.accession = 'PF03030'
                    """,
                    expected_ref=1548,
                ),
                VerificationSpec(
                    claim="Membrane H_PPase occurs without soluble pyrophosphatase in 1,473 genomes.",
                    query="""
                    WITH flags AS (
                      SELECT p.bin_id,
                             BOOL_OR(a.accession = 'PF03030') AS has_hppase,
                             BOOL_OR(a.accession IN ('PF00719','K01507')) AS has_soluble
                      FROM annotations a
                      JOIN proteins p ON a.protein_id = p.protein_id
                      WHERE a.accession IN ('PF03030','PF00719','K01507')
                      GROUP BY p.bin_id
                    )
                    SELECT COUNT(*) FROM flags
                    WHERE has_hppase AND NOT has_soluble
                    """,
                    expected_ref=1473,
                ),
            ),
            protein_ids_query="""
            SELECT DISTINCT p.protein_id
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.accession = 'PF03030'
            ORDER BY p.sequence_length DESC, p.protein_id
            LIMIT 5
            """,
        ),
    }


def enrich_findings_dataset(
    dataset_dir: str | Path,
    *,
    apply_changes: bool = False,
) -> dict[str, Any]:
    dataset_dir = Path(dataset_dir)
    db_path = dataset_dir / "sharur.duckdb"
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")

    recipes = omni_production_recipes()
    summaries: list[dict[str, Any]] = []
    total_changed = 0
    total_verified = 0
    total_proteins = 0
    missing_recipe: list[str] = []
    mismatches: dict[str, list[dict[str, Any]]] = {}

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        for phase in ("survey", "exploration"):
            path = dataset_dir / phase / "findings.jsonl"
            records = _load_jsonl(path)
            updated_records: list[dict[str, Any]] = []
            changed = 0
            verified = 0
            proteins = 0

            for record in records:
                recipe = recipes.get(str(record.get("id")))
                if recipe is None:
                    updated_records.append(record)
                    missing_recipe.append(str(record.get("id")))
                    continue

                updated, summary = enrich_record_with_recipe(con, record, recipe)
                updated_records.append(updated)
                if summary["changed"]:
                    changed += 1
                verified += summary["added_verification"]
                proteins += summary["added_protein_ids"]
                if summary["mismatches"]:
                    mismatches[str(record.get("id"))] = summary["mismatches"]

            if apply_changes:
                _write_jsonl(path, updated_records)

            summaries.append(
                {
                    "phase": phase,
                    "path": str(path),
                    "records": len(records),
                    "changed": changed,
                    "verification_entries_added": verified,
                    "protein_ids_added": proteins,
                }
            )
            total_changed += changed
            total_verified += verified
            total_proteins += proteins
    finally:
        con.close()

    return {
        "dataset_dir": str(dataset_dir),
        "apply_changes": apply_changes,
        "recipes_available": len(recipes),
        "records_without_recipe": sorted(set(missing_recipe)),
        "mismatches": mismatches,
        "phases": summaries,
        "changed": total_changed,
        "verification_entries_added": total_verified,
        "protein_ids_added": total_proteins,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Bootstrap verification blocks and representative protein IDs from DuckDB "
            "for findings that have safe, reproducible query recipes."
        )
    )
    parser.add_argument("--dataset", required=True, help="Path to dataset directory.")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write updated findings.jsonl files in place.",
    )
    args = parser.parse_args()

    summary = enrich_findings_dataset(args.dataset, apply_changes=args.apply)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
