"""Dataset health checks, the `sharur health` command and the /health page."""

import json

import pytest
from fastapi.testclient import TestClient
from typer.testing import CliRunner

from sharur.assemblies import gene_call_deficits
from sharur.health import run_checks
from sharur.storage.duckdb_store import DuckDBStore


@pytest.fixture
def dataset(tmp_path):
    """Four genomes, each with one problem:
    g1: 4 proteins on a 50 kb assembly (gene calls missing), the only MinCED report;
    g2: no completeness, numeric strands, an unscanned 2 kb assembly;
    g3: no proteins;
    g4: no annotation, one protein its own contig, one on a contig shared with g1."""
    d = tmp_path / "ds"
    d.mkdir()
    path = d / "sharur.duckdb"
    s = DuckDBStore(str(path))
    c = s.conn
    for b, comp in (("g1", 92.0), ("g2", None), ("g3", 80.0), ("g4", 88.0)):
        c.execute("INSERT INTO bins (bin_id, completeness, contamination) VALUES (?, ?, 1.0)", [b, comp])
    for contig, b in (("g1_c1", "g1"), ("g2_c1", "g2"), ("shared", "g1"), ("g4_p1", "g4"), ("empty_c", "g3")):
        c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 5000)", [contig, b])
    proteins = [
        ("g1_1", "g1_c1", "g1", 1, 300, "+"), ("g1_2", "g1_c1", "g1", 400, 900, "-"), ("g1_3", "g1_c1", "g1", 1000, 1500, "+"),
        ("g2_1", "g2_c1", "g2", 1, 300, "1"), ("g2_2", "g2_c1", "g2", 400, 900, "-1"),
        ("g4_p1", "g4_p1", "g4", 0, 300, "+"),
        ("g4_s1", "shared", "g4", 0, 300, "+"), ("g1_s", "shared", "g1", 0, 300, "+"),
    ]
    c.executemany("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) "
                  "VALUES (?, ?, ?, ?, ?, ?)", proteins)
    ann = [(i, pid, "pfam", "PF00001", "A") for i, pid in enumerate(["g1_1", "g1_2", "g2_1", "g2_2"], 1)]
    c.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name) VALUES (?, ?, ?, ?, ?)",
                  ann)
    s.close()
    fna = d / "genomes_fna"
    fna.mkdir()
    (fna / "g1.fna").write_text(">c\n" + "A" * 50_000 + "\n")
    (fna / "g2.fna").write_text(">c\n" + "A" * 2_000 + "\n")
    crispr = d / "stage05c_crispr"
    crispr.mkdir()
    (crispr / "g1_crispr.txt").write_text("")
    seal = {"generated_at": "2026-01-01T00:00:00", "seal_strength": "structural",
            "identity": {"database": {"schema_version": 9, "tables": [{"table": "proteins", "rows": 2}]}}}
    (d / "dataset.seal.json").write_text(json.dumps(seal))
    return path


def test_checks_find_each_problem(dataset):
    report = run_checks(dataset)
    by = {f.key: f for f in report.findings}
    assert by["schema"].status == "ok"
    assert by["seal"].status == "warn" and any("proteins: 2 → 8" in e["id"] for e in by["seal"].examples)
    assert by["completeness"].status == "warn" and [e["id"] for e in by["completeness"].examples] == ["g2"]
    assert by["empty_genomes"].status == "fail" and by["empty_genomes"].examples[0]["id"] == "g3"
    gene = by["gene_calls"]
    assert gene.status == "warn" and gene.examples[0]["id"] == "g1" and "4 proteins on" in gene.examples[0]["note"]
    assert by["strand"].status == "warn" and "'1' 1" in by["strand"].summary
    pos = by["positionless"]
    assert pos.status == "warn" and "2 proteins share 1 contig ID across genomes" in pos.summary
    assert "1 protein in 1 genome is its own contig" in pos.summary and "stacked" not in pos.summary
    cov = by["annotation_coverage"]
    assert cov.status == "warn" and "no annotation at all" in cov.summary and cov.table["rows"][0][0] == "pfam"
    assert by["crispr_scan"].status == "warn" and by["crispr_scan"].examples[0]["id"] == "g2"
    assert report.worst == "fail"
    json.dumps(report.to_dict())
    assert "## Genes" in report.to_markdown()


def test_annotation_coverage_flags_a_dense_source_missing_from_a_large_genome(tmp_path):
    """defensefinder hits 10 of 200 proteins in g1 and g2, so g3's 200 proteins predict 10 hits: g3 was never
    searched. hyddb (one hit) and caller rows (defensefinder_system) predict too few to read anything into."""
    path = tmp_path / "sharur.duckdb"
    s = DuckDBStore(str(path))
    c = s.conn
    for b in ("g1", "g2", "g3"):
        c.execute("INSERT INTO bins (bin_id) VALUES (?)", [b])
        c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 300000)", [f"{b}_c", b])
        c.executemany("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) VALUES (?, ?, ?, ?, ?, '+')",
                      [(f"{b}_{i}", f"{b}_c", b, 1000 * i + 1, 1000 * i + 900) for i in range(200)])
    rows = [(f"{b}_{i}", "pfam") for b in ("g1", "g2", "g3") for i in range(100)]
    rows += [(f"{b}_{i}", "defensefinder") for b in ("g1", "g2") for i in range(10)]
    rows += [("g1_0", "hyddb")] + [(f"g1_{i}", "defensefinder_system") for i in range(10)]
    c.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name) VALUES (?, ?, ?, 'X', 'X')",
                  [(n, pid, src) for n, (pid, src) in enumerate(rows, 1)])
    s.close()

    cov = {f.key: f for f in run_checks(path).findings}["annotation_coverage"]

    assert cov.status == "warn" and "1 genome lacks defensefinder hits that 2 others have" in cov.summary
    assert [e["note"] for e in cov.examples] == ["no defensefinder"]


def test_gene_call_deficits_prefer_assemblies():
    proteins = {f"g{i}": 1000 for i in range(6)} | {"small": 20}
    completeness = dict.fromkeys(proteins, 90.0)
    assert gene_call_deficits(proteins, completeness) == {"small"}
    assert gene_call_deficits(proteins, completeness, {"small": 30.0}) == set()


def test_health_command(dataset):
    from sharur.cli import app

    runner = CliRunner()
    out = runner.invoke(app, ["health", "--db", str(dataset), "--format", "json"])
    assert out.exit_code == 0 and json.loads(out.output)["counts"]["fail"] == 1
    assert runner.invoke(app, ["health", "--db", str(dataset), "--strict"]).exit_code == 1


def test_health_page_and_rail_pill(dataset):
    from sharur.browser import create_app

    client = TestClient(create_app(dataset, background=False))
    page = client.get("/health")
    assert page.status_code == 200 and "Gene calls per genome" in page.text and 'href="/genome/g1"' in page.text
    assert "health-pill hpill-fail" in client.get("/").text
    assert client.get("/api/health").json()["worst"] == "fail"
