"""Guards for the example genome set a new user downloads."""

import csv
import subprocess
from pathlib import Path

import pytest
from sharur.ingest.input_integrity import generate_genome_id


REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "complete_genomes"


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
