"""Per-sample coverage import, genome and feature abundance, coverage outliers, stage 08."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from sharur.abundance import (
    coverage_outliers,
    feature_abundance,
    feature_abundance_markdown,
    genome_abundance,
    genome_abundance_markdown,
    import_coverage,
    parse_coverm,
    samples,
)
from sharur.storage.duckdb_store import DuckDBStore


REPO_ROOT = Path(__file__).resolve().parents[1]

# bin A: four 10 kb contigs, a3 sits at 4x the others; bin B: one 20 kb contig
CONTIGS = {"a1": ("A", 10000), "a2": ("A", 10000), "a3": ("A", 10000), "a4": ("A", 10000), "b1": ("B", 20000)}


def _core(store: DuckDBStore) -> DuckDBStore:
    values = ", ".join(f"('{c}', '{b}', {n})" for c, (b, n) in CONTIGS.items())
    store.conn.execute(f"""
        INSERT INTO bins (bin_id) VALUES ('A'), ('B');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES {values};
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES ('a1_1', 'a1', 'A', 1, 900, '+', 0, 300), ('b1_1', 'b1', 'B', 1, 900, '+', 0, 300);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue)
        VALUES (1, 'b1_1', 'kofam', 'K00370', 'narG', 1e-50);
        INSERT INTO protein_predicates (protein_id, predicates) VALUES
            ('a1_1', ['hydrogenase']), ('b1_1', ['nitrate_reductase']);
    """)
    return store


def _coverm(path: Path, depths: dict[str, dict[str, float]]) -> Path:
    names = list(depths)
    header = ["Contig"] + [f"{s} {m}" for s in names for m in ("Mean", "Covered Fraction", "Read Count")]
    lines = ["\t".join(header)]
    for contig, (_, length) in CONTIGS.items():
        row = [contig]
        for s in names:
            depth = depths[s][contig]
            row += [str(depth), "0.95", str(int(depth * length / 150))]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")
    return path


DEPTHS = {
    "s1": {"a1": 10, "a2": 10, "a3": 40, "a4": 10, "b1": 5},
    "s2": {"a1": 2, "a2": 2, "a3": 8, "a4": 2, "b1": 20},
}


@pytest.fixture
def loaded(tmp_path):
    store = _core(DuckDBStore())
    sidecar = tmp_path / "abundance.duckdb"
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\thabitat\ns1\tsediment\ns2\twater\n")
    result = import_coverage(store, sidecar, [_coverm(tmp_path / "cov.tsv", DEPTHS)], sample_metadata=meta)
    return store, sidecar, result


def test_coverm_parsing_handles_multiple_samples(tmp_path):
    rows = parse_coverm(_coverm(tmp_path / "cov.tsv", DEPTHS))
    assert len(rows) == 10
    row = next(r for r in rows if r["sample_id"] == "s1" and r["contig_id"] == "a3")
    assert (row["mean_depth"], row["covered_fraction"], row["read_count"]) == (40.0, 0.95, 2666)


def test_import_records_samples_and_metadata(loaded):
    _, sidecar, result = loaded
    assert result == {"sidecar": str(sidecar), "samples": ["s1", "s2"], "rows": 10, "unknown_contigs": 0}
    assert samples(sidecar)[0] == {"sample_id": "s1", "metadata": {"habitat": "sediment"}, "contigs": 5}


def test_genome_abundance_weights_by_length_and_reads(loaded):
    _, sidecar, _ = loaded
    rows = {(r["sample_id"], r["bin_id"]): r for r in genome_abundance(sidecar)}
    a1 = rows[("s1", "A")]
    assert a1["mean_depth"] == pytest.approx(17.5)            # (10+10+40+10)/4 with equal lengths
    reads_a, reads_b = sum(int(DEPTHS["s1"][c] * 10000 / 150) for c in ("a1", "a2", "a3", "a4")), int(5 * 20000 / 150)
    assert a1["relative_abundance"] == pytest.approx(reads_a / (reads_a + reads_b))
    assert rows[("s2", "B")]["relative_abundance"] > rows[("s2", "A")]["relative_abundance"]
    assert "of mapped reads" in genome_abundance_markdown(list(rows.values()))


def test_feature_abundance_by_predicate_and_annotation(loaded):
    store, sidecar, _ = loaded
    by_predicate = feature_abundance(store, sidecar, predicate="nitrate_reductase")
    by_annotation = feature_abundance(store, sidecar, annotation="K00370")
    assert by_predicate["carrier_genomes"] == by_annotation["carrier_genomes"] == 1
    s2 = next(s for s in by_predicate["samples"] if s["sample_id"] == "s2")
    assert s2["share_in_carrier_genomes"] > 0.5 and s2["top_carriers"][0]["bin_id"] == "B"
    assert "1 genomes carry it" in feature_abundance_markdown(by_predicate)
    with pytest.raises(ValueError):
        feature_abundance(store, sidecar)


def test_coverage_outliers_flag_the_high_contig(loaded):
    _, sidecar, _ = loaded
    result = coverage_outliers(sidecar, "A")
    assert [(o["contig_id"], o["direction"]) for o in result["outliers"]] == [("a3", "high")]
    assert result["outliers"][0]["median_log2_ratio"] == pytest.approx(2.0)
    assert result["samples_used"] == 2


def test_unknown_contigs_are_rejected_unless_allowed(tmp_path):
    store = _core(DuckDBStore())
    table = tmp_path / "long.tsv"
    table.write_text("sample_id\tcontig_id\tmean_depth\ns1\ta1\t3\ns1\tzz\t4\n")
    with pytest.raises(ValueError, match="absent from the dataset"):
        import_coverage(store, tmp_path / "ab.duckdb", [table], fmt="long")
    result = import_coverage(store, tmp_path / "ab.duckdb", [table], fmt="long", allow_unknown_contigs=True)
    assert (result["rows"], result["unknown_contigs"]) == (1, 1)


def test_stage08_runs_coverm_per_sample_and_imports(tmp_path):
    data = tmp_path / "data"
    asm = data / "stage00_prepared" / "genomes" / "A.fna"
    asm.parent.mkdir(parents=True)
    asm.write_text(">a1\nACGT\n")
    (data / "stage00_prepared" / "processing_manifest.json").write_text(
        json.dumps({"genomes": [{"genome_id": "A", "filename": "A.fna", "output_path": str(asm)}]}))
    db = data / "sharur.duckdb"
    _core(DuckDBStore(db)).close()
    reads = tmp_path / "reads.tsv"
    r1 = tmp_path / "s1_R1.fq"
    r1.write_text("@r\nA\n+\nI\n")
    reads.write_text(f"sample_id\tread1\thabitat\ns1\t{r1}\tsoil\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "coverm"
    fake.write_text(
        "#!/bin/sh\n"
        "while [ $# -gt 0 ]; do if [ \"$1\" = -o ]; then out=$2; fi; shift; done\n"
        "printf 'Contig\\treference.fna/s1_R1.fq Mean\\treference.fna/s1_R1.fq Length\\n"
        "a1\\t7.5\\t10000\\n' > \"$out\"\n")
    fake.chmod(0o755)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}
    subprocess.run([sys.executable, str(REPO_ROOT / "src/ingest/08_coverage.py"), "--data-dir", str(data),
                    "--db", str(db), "--reads", str(reads), "--threads", "2"], check=True, env=env)
    (row,) = genome_abundance(data / "abundance.duckdb")
    assert (row["sample_id"], row["bin_id"], row["mean_depth"]) == ("s1", "A", 7.5)
    assert samples(data / "abundance.duckdb")[0]["metadata"] == {"habitat": "soil"}
    manifest = json.loads((data / "stage08_coverage" / "processing_manifest.json").read_text())
    assert manifest["threads"] == 2 and manifest["samples"] == ["s1"]
