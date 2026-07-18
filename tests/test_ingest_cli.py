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
        "07_build_knowledge_base.py",
        "06_esm2_embeddings.py",
        "vector_index_runner.py",
    ]
    for stage in expected:
        assert stage in stages_seen
    # Order sanity: build the core DB before optional embeddings.
    assert plan[0][0].endswith("00_prepare_inputs.py")
    assert plan[-2][0].endswith("06_esm2_embeddings.py")
    assert plan[-1][0].endswith("vector_index_runner.py")
    stage07_index = next(
        i for i, command in enumerate(plan) if command[0].endswith("07_build_knowledge_base.py")
    )
    stage06_index = next(
        i for i, command in enumerate(plan) if command[0].endswith("06_esm2_embeddings.py")
    )
    assert stage07_index < stage06_index
    assert plan[-2][-1] == str(tmp_path / "embeddings")
    # Ensure dbCAN input path points to all_protein_symlinks
    dbcan_cmd = next(cmd for cmd in plan if "dbcan_cazyme.py" in cmd[0])
    assert "all_protein_symlinks" in " ".join(dbcan_cmd)


def test_dry_run_does_not_create_dataset_directory(tmp_path):
    data_dir = tmp_path / "planned_dataset"

    ingest.run(
        input_dir=Path("dummy_dataset"),
        data_dir=data_dir,
        output=data_dir / "sharur.duckdb",
        mode="tools",
        skip_quast=True,
        skip_dfast=True,
        skip_prodigal=False,
        skip_astra=False,
        skip_gecco=True,
        skip_dbcan=True,
        skip_crispr=True,
        skip_embeddings=True,
        dry_run=True,
    )

    assert not data_dir.exists()


def test_default_plan_skips_optional_and_deprecated_stages(tmp_path):
    plan = ingest.run(
        input_dir=Path("dummy_dataset"),
        data_dir=tmp_path / "dataset",
        output=tmp_path / "dataset" / "sharur.duckdb",
        mode="tools",
        dry_run=True,
    )

    assert plan is not None
    stages_seen = " ".join(" ".join(command) for command in plan)
    assert "01_run_quast.py" not in stages_seen
    assert "02_dfast_qc.py" not in stages_seen
    assert "gecco_bgc.py" not in stages_seen
    assert "dbcan_cazyme.py" not in stages_seen
    assert "03_prodigal.py" in stages_seen
    assert "04_astra_scan.py" in stages_seen
    assert "minced_crispr.py" in stages_seen
    assert "07_build_knowledge_base.py" in stages_seen
    assert "06_esm2_embeddings.py" in stages_seen


def test_cazyme_consensus_opt_in_reaches_stage07(tmp_path):
    plan = ingest.run(
        input_dir=Path("dummy_dataset"),
        data_dir=tmp_path / "dataset",
        output=tmp_path / "dataset" / "sharur.duckdb",
        mode="tools",
        enable_cazymes=True,
        dry_run=True,
    )

    builder = next(
        command for command in plan or [] if command[0].endswith("07_build_knowledge_base.py")
    )
    assert "--enable-cazymes" in builder
