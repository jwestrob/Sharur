"""Guards for the files a new user installs from: environment.yml and the example genome set."""

import csv
import subprocess
from pathlib import Path

import pytest
import yaml
from packaging.requirements import Requirement

from sharur.ingest.input_integrity import generate_genome_id


REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "complete_genomes"


def _pip_lines() -> list[str]:
    env = yaml.safe_load((REPO / "environment.yml").read_text())
    for dependency in env["dependencies"]:
        if isinstance(dependency, dict) and "pip" in dependency:
            return dependency["pip"]
    raise AssertionError("environment.yml has no pip section")


def test_environment_pip_lines_are_valid_requirement_file_lines():
    # conda writes these lines into a pip requirements file; a requirement carrying its own
    # option (for example "--no-deps pkg") makes pip abort the whole pip section.
    for line in _pip_lines():
        if line.startswith("-e "):
            continue
        assert not line.startswith("-"), line
        Requirement(line)


def test_environment_installs_sharur_and_its_stage_tools():
    assert _pip_lines()[0] == '-e ".[all,dev]"'
    env = yaml.safe_load((REPO / "environment.yml").read_text())
    conda = {d.split("=")[0].split(" ")[0] for d in env["dependencies"] if isinstance(d, str)}
    # Stage 03, 04 (hmmer for profile tools), 05c, 09
    assert {"prodigal", "hmmer", "minced", "mmseqs2"} <= conda
    assert any(Requirement(line).name == "aksha" for line in _pip_lines() if not line.startswith("-"))


def _example_rows() -> list[dict[str, str]]:
    with (EXAMPLE / "genomes.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_example_genome_table_is_consistent():
    rows = _example_rows()
    assert 12 <= len(rows) <= 16
    ids = [row["genome_id"] for row in rows]
    assert len(set(ids)) == len(ids)
    assert len({row["assembly"] for row in rows}) == len(rows)
    for row in rows:
        # Stage 00 derives the genome id from the file name download.sh writes.
        assert generate_genome_id(Path(f"{row['genome_id']}.fna")) == row["genome_id"]
        assert row["assembly"].startswith("GCA_")
        assert row["assembly_level"] == "Complete Genome"
        replicons = row["replicon_detail"].split(";")
        assert len(replicons) == int(row["replicons"])
        assert sum(int(r.split(":")[3]) for r in replicons) == int(row["total_length"])
    assert {row["domain"] for row in rows} == {"Archaea", "Bacteria"}
    assert sum(int(row["plasmids"]) > 0 for row in rows) >= 2
    assert sum(int(row["total_length"]) for row in rows) <= 50_000_000


@pytest.mark.skipif(not Path("/bin/bash").exists(), reason="needs bash")
def test_example_download_script_parses():
    subprocess.run(["/bin/bash", "-n", str(EXAMPLE / "download.sh")], check=True)
