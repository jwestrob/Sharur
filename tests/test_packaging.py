"""Guards for what the wheel must carry beyond Python modules."""

import re
import tomllib
from importlib.machinery import SourceFileLoader
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PYPROJECT = tomllib.loads((REPO / "pyproject.toml").read_text())
FORCE = PYPROJECT["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]


def test_every_stage_the_ingest_cli_runs_is_packaged():
    cli = (REPO / "sharur/ingest_cli.py").read_text()
    stages = set(re.findall(r'stage_dir / "([0-9a-z_]+\.py)"', cli))
    assert {"07_build_knowledge_base.py", "08_coverage.py"} <= stages
    packaged = {Path(target).name for target in FORCE.values() if target.startswith("sharur/ingest/stages/")}
    assert stages <= packaged, stages - packaged
    for source, target in FORCE.items():
        assert (REPO / source).exists(), source
        assert target.startswith("sharur/"), target


def test_stage_helpers_imported_at_runtime_are_packaged():
    # stage 07 imports classify_cazymes for --enable-cazymes
    assert FORCE.get("scripts/classify_cazymes.py") == "sharur/ingest/stages/classify_cazymes.py"


def test_provider_attribution_stays_attached_to_the_wheel():
    for notice in ("DATA_LICENSES.md", "CITATIONS.md"):
        assert FORCE.get(notice) == f"sharur/{notice}"
        assert (REPO / notice).is_file()


def test_package_data_lives_inside_the_package_and_is_not_excluded():
    wheel = PYPROJECT["tool"]["hatch"]["build"]["targets"]["wheel"]
    assert wheel["packages"] == ["sharur"]
    assert not wheel.get("exclude") and not wheel.get("only-include")
    for directory, pattern in [("sharur/browser/templates", "*.html"), ("sharur/browser/static", "*"),
                               ("sharur/predicates/mappings/data", "*.tsv")]:
        assert list((REPO / directory).glob(pattern)), directory
    assert FORCE.get("config/predicates_v2") == "sharur/predicates_v2/config"


def test_stage07_finds_reference_maps_beside_the_dataset(tmp_path, monkeypatch):
    module = SourceFileLoader("kb_refs", str(REPO / "src/ingest/07_build_knowledge_base.py")).load_module()
    monkeypatch.delenv("SHARUR_REFERENCE_DIR", raising=False)
    monkeypatch.chdir(tmp_path)
    reference = tmp_path / "data" / "reference"
    reference.mkdir(parents=True)
    (reference / "ko_list").write_text("knum\tdefinition\n")
    dirs = module.reference_dirs(tmp_path / "data" / "my_dataset")
    assert dirs[0] == reference
    assert module.find_reference(dirs, "ko_list") == reference / "ko_list"
    # a file missing beside the dataset falls through to later locations, never to the dataset's own reference dir
    found = module.find_reference(dirs, "pfam_id_desc.tsv")
    assert found is None or found.parent != reference
    override = tmp_path / "elsewhere"
    override.mkdir()
    monkeypatch.setenv("SHARUR_REFERENCE_DIR", str(override))
    assert module.reference_dirs(tmp_path / "data" / "my_dataset")[0] == override
