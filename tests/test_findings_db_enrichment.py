from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import duckdb

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script(relative_path: str, module_name: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


ENRICH_MODULE = _load_script(
    "scripts/enrich_findings_from_db.py",
    "test_enrich_findings_from_db_script",
)


def test_enrich_record_with_recipe_adds_verification_and_protein_ids(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "mini.duckdb"
    con = duckdb.connect(str(db_path))
    con.execute(
        """
        CREATE TABLE bins (
            bin_id VARCHAR,
            completeness FLOAT,
            contamination FLOAT,
            taxonomy VARCHAR,
            n_contigs INTEGER,
            total_length INTEGER
        )
        """
    )
    con.execute(
        """
        CREATE TABLE proteins (
            protein_id VARCHAR,
            contig_id VARCHAR,
            bin_id VARCHAR,
            start INTEGER,
            end_coord INTEGER,
            strand VARCHAR,
            gene_index INTEGER,
            sequence VARCHAR,
            sequence_length INTEGER,
            gc_content FLOAT
        )
        """
    )
    con.execute(
        """
        CREATE TABLE protein_predicates (
            protein_id VARCHAR,
            predicates VARCHAR[],
            updated_at TIMESTAMP
        )
        """
    )
    con.execute(
        "INSERT INTO bins VALUES ('bin_a', 90, 0, NULL, 1, 1000), ('bin_b', 90, 0, NULL, 1, 900)"
    )
    con.execute(
        """
        INSERT INTO proteins VALUES
        ('prot_1', 'contig_1', 'bin_a', 1, 300, '+', 1, NULL, 100, 0.5),
        ('prot_2', 'contig_2', 'bin_b', 1, 450, '+', 1, NULL, 150, 0.5)
        """
    )
    con.execute(
        """
        INSERT INTO protein_predicates VALUES
        ('prot_1', ['example_pred'], NOW()),
        ('prot_2', ['example_pred'], NOW())
        """
    )

    recipe = ENRICH_MODULE.FindingRecipe(
        verification_specs=(
            ENRICH_MODULE.VerificationSpec(
                claim="example_pred occurs in 2 proteins",
                query="""
                SELECT COUNT(DISTINCT pp.protein_id)
                FROM protein_predicates pp,
                     UNNEST(pp.predicates) AS t(pred)
                WHERE pred = 'example_pred'
                """,
                expected_ref=2,
            ),
        ),
        protein_ids_query="""
        SELECT protein_id
        FROM proteins
        ORDER BY sequence_length DESC, protein_id
        LIMIT 2
        """,
    )

    record = {
        "id": "test-001",
        "title": "Test finding",
        "category": "general",
        "description": "Test finding body.",
        "evidence": {"source_format": "test"},
    }

    updated, summary = ENRICH_MODULE.enrich_record_with_recipe(con, record, recipe)

    assert summary["changed"] is True
    assert summary["added_verification"] == 1
    assert summary["added_protein_ids"] == 2
    assert updated["verification"][0]["expected"] == 2
    assert updated["protein_ids"] == ["prot_2", "prot_1"]
    assert updated["contigs"] == ["contig_1", "contig_2"]

    con.close()


def test_enrich_record_with_recipe_reports_mismatch(tmp_path: Path) -> None:
    db_path = tmp_path / "mini.duckdb"
    con = duckdb.connect(str(db_path))
    con.execute("CREATE TABLE proteins (protein_id VARCHAR, contig_id VARCHAR, bin_id VARCHAR, start INTEGER, end_coord INTEGER, strand VARCHAR, gene_index INTEGER, sequence VARCHAR, sequence_length INTEGER, gc_content FLOAT)")
    con.execute("INSERT INTO proteins VALUES ('prot_1', 'contig_1', 'bin_a', 1, 300, '+', 1, NULL, 100, 0.5)")

    recipe = ENRICH_MODULE.FindingRecipe(
        verification_specs=(
            ENRICH_MODULE.VerificationSpec(
                claim="wrong count",
                query="SELECT COUNT(*) FROM proteins",
                expected_ref=2,
            ),
        ),
    )

    updated, summary = ENRICH_MODULE.enrich_record_with_recipe(
        con,
        {"id": "test-002", "evidence": {"source_format": "test"}},
        recipe,
    )

    assert summary["changed"] is False
    assert summary["added_verification"] == 0
    assert summary["mismatches"][0]["actual"] == 1
    assert "verification" not in updated

    con.close()


def test_omni_production_recipes_include_second_tranche_ids() -> None:
    recipes = ENRICH_MODULE.omni_production_recipes()

    assert {"E004", "E007", "E009", "E027"}.issubset(recipes)
    assert len(recipes["E004"].verification_specs) >= 3
    assert len(recipes["E007"].verification_specs) >= 4
    assert len(recipes["E009"].verification_specs) >= 4
    assert len(recipes["E027"].verification_specs) >= 4
