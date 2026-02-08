"""Run core tools against a DuckDB built from dummy_dataset."""

from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

from bennu.core.session import ExplorationSession
from bennu.tools.find_proteins import FindProteinsTool, FindProteinsParams
from bennu.tools.get_context import GetContextTool, GetContextParams
from bennu.tools.manage_sets import ManageSetsTool, ManageSetsParams
from bennu.tools.export import ExportTool, ExportParams


REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"
kb_module = SourceFileLoader("kb_build_full", str(KB_PATH)).load_module()
KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder
PipelineOutputs = kb_module.PipelineOutputs


def _build_db(tmp_path: Path) -> tuple[Path, list[str], int]:
    """Build DuckDB from dummy_dataset, returning path, protein ids, annotated count."""
    dataset_dir = REPO_ROOT / "dummy_dataset"
    if not dataset_dir.exists():
        pytest.skip("dummy_dataset not available")

    # Reuse integration builder to construct minimal stage outputs
    from tests.test_ingest_integration_dataset import _parse_contigs  # type: ignore

    data_dir = tmp_path / "data"
    stage02 = data_dir / "stage02_dfast_qc"
    stage03 = data_dir / "stage03_prodigal"
    stage04 = data_dir / "stage04_astra"
    stage05a = data_dir / "stage05a_gecco"
    stage05b = data_dir / "stage05b_dbcan"
    stage05c = data_dir / "stage05c_crispr"

    genomes_manifest = {"genomes": []}
    pfam_rows = []
    gecco_clusters = []
    protein_ids: list[str] = []

    for fasta_path in sorted(dataset_dir.glob("*.fna")):
        raw_bin_id = fasta_path.name.replace(".contigs.fna", "")
        bin_id = raw_bin_id.replace("_", "")
        contigs = list(_parse_contigs(fasta_path))
        genomes_manifest["genomes"].append(
            {
                "genome_id": bin_id,
                "taxonomy": {"name": "unknown", "completeness": 90.0, "contamination": 1.0},
                "n_contigs": len(contigs),
                "total_length": sum(len(seq) for _, seq in contigs),
            }
        )

        bin_dir = stage03 / "genomes" / bin_id
        bin_dir.mkdir(parents=True, exist_ok=True)
        faa_path = bin_dir / f"{bin_id}.faa"
        with open(faa_path, "w") as faa:
            for idx, (contig_id, seq) in enumerate(contigs):
                pid = f"{bin_id}_ctg{idx:05d}_00001"
                protein_ids.append(pid)
                faa.write(f">{pid} # 1 # {min(len(seq), 300)} # 1 # ID={pid};partial=00\n")
                faa.write("M" * 30 + "\n")
                if idx == 0:
                    pfam_rows.append(
                        {
                            "sequence_id": pid,
                            "hmm_name": "PF00001",
                            "e_value": 1e-5,
                            "bit_score": 42.0,
                            "ali_from": 1,
                            "ali_to": 30,
                            "name": "MockDomain",
                            "description": "Synthetic annotation",
                        }
                    )
                    gecco_clusters.append(
                        {
                            "cluster_id": f"{bin_id}_cluster1",
                            "contig": pid.rsplit("_", 1)[0],
                            "start": 1,
                            "end": 100,
                            "bgc_type": "nrps",
                            "protein_list": [pid],
                        }
                    )

    stage02.mkdir(parents=True, exist_ok=True)
    (stage02 / "processing_manifest.json").write_text(__import__("json").dumps(genomes_manifest))

    if pfam_rows:
        header = [
            "sequence_id",
            "hmm_name",
            "e_value",
            "bit_score",
            "ali_from",
            "ali_to",
            "name",
            "description",
        ]
        lines = ["\t".join(header)]
        for row in pfam_rows:
            lines.append("\t".join(str(row[h]) for h in header))
        stage04.mkdir(parents=True, exist_ok=True)
        (stage04 / "synthetic_hits_df.tsv").write_text("\n".join(lines) + "\n")

    stage05a.mkdir(parents=True, exist_ok=True)
    (stage05a / "combined_bgc_data.json").write_text(__import__("json").dumps({"clusters": gecco_clusters}))
    stage05b.mkdir(parents=True, exist_ok=True)
    (stage05b / "synthetic_cazyme_results.json").write_text(__import__("json").dumps({"annotations": []}))
    stage05c.mkdir(parents=True, exist_ok=True)
    (stage05c / "synthetic_crispr_arrays.json").write_text(__import__("json").dumps({"arrays": []}))

    outputs = PipelineOutputs(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=stage02,
        stage03_dir=stage03,
        stage04_dir=stage04,
        stage05a_dir=stage05a,
        stage05b_dir=stage05b,
        stage05c_dir=stage05c,
        stage06_dir=data_dir / "stage06_embeddings",
    )

    db_path = data_dir / "bennu.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    builder.build()
    return db_path, protein_ids, len(pfam_rows)


def test_tool_invocations_on_dataset(tmp_path):
    db_path, protein_ids, annotated_count = _build_db(tmp_path)
    session = ExplorationSession(db_path=db_path)

    # find_proteins by PFAM
    fp_params = FindProteinsParams(domains=["PF00001"], limit=10)
    fp_result = FindProteinsTool().execute(fp_params, session)
    assert fp_result.success
    assert fp_result.count >= annotated_count  # annotated per-bin proteins retrieved

    # push focus from result and get_context
    if fp_result.data:
        session.push_focus("protein", fp_result.data[0].protein_id)
    ctx_params = GetContextParams(window_genes=1)
    ctx_result = GetContextTool().execute(ctx_params, session)
    assert ctx_result.success
    assert ctx_result.count >= 1

    # manage_sets create and export
    ms_params = ManageSetsParams(action="create", set_name="hit_set", protein_ids=[p.protein_id for p in fp_result.data[:2]])
    ms_result = ManageSetsTool().execute(ms_params, session)
    assert ms_result.success
    assert ms_result.data.name == "hit_set"

    ex_params = ExportParams(format="tsv", source="set", source_id="hit_set", include_sequences=False, include_annotations=False)
    ex_result = ExportTool().execute(ex_params, session)
    assert ex_result.success
    assert "protein_id" in ex_result.data.split("\n")[0]
