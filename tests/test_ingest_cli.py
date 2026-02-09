from pathlib import Path

from scripts import ingest


def test_ingest_dry_run_plans_all_stages(tmp_path):
    plan = ingest.run(
        input_dir=Path("dummy_dataset"),
        data_dir=tmp_path,
        output=tmp_path / "sharur.duckdb",
        mode="tools",
        force=False,
        skip_quast=False,
        skip_dfast=False,
        skip_prodigal=False,
        skip_astra=False,
        skip_gecco=False,
        skip_dbcan=False,
        skip_crispr=True,  # no runner bundled
        skip_embeddings=False,
        dry_run=True,
    )
    assert plan is not None
    stages_seen = " ".join(" ".join(cmd) for cmd in plan)
    expected = [
        "00_prepare_inputs.py",
        "01_run_quast.py",
        "02_dfast_qc.py",
        "03_prodigal.py",
        "04_astra_scan.py",
        "gecco_bgc.py",
        "dbcan_cazyme.py",
        "06_esm2_embeddings.py",
    ]
    for stage in expected:
        assert stage in stages_seen
    # Order sanity: first is stage 00, last is stage 06
    assert plan[0][0].endswith("00_prepare_inputs.py")
    assert plan[-1][0].endswith("06_esm2_embeddings.py")
    # Ensure dbCAN input path points to all_protein_symlinks
    dbcan_cmd = [cmd for cmd in plan if "dbcan_cazyme.py" in cmd[0]][0]
    assert "all_protein_symlinks" in " ".join(dbcan_cmd)
