#!/usr/bin/env python3
"""Download the reference inputs the predicate maps are built from (maintainer step).

Writes a dated snapshot directory with a ``manifest.json`` (URL, retrieval
time, sha256, size, upstream release) for:

- Pfam-A.clans.tsv.gz and Pfam.version.gz (EBI FTP, CC0)
- pfam2go (Gene Ontology external2go) and go-basic.obo (CC BY 4.0)
- enzyme.dat and enzclass.txt (Expasy ENZYME, CC BY 4.0)
- uniprot_sprot.dat.gz and reldate.txt (UniProtKB/Swiss-Prot, CC BY 4.0)
- ``kegg_inputs/``: the KEGG REST files ``sharur setup-kegg`` builds from
  (academic use only; subject to KEGG's terms)
- with ``--kegg-links``: KEGG gene -> KO links for every organism Swiss-Prot
  references (about 2,200 requests, resumable) and the joined
  ``swissprot_kegg_ko.tsv``

Then run ``scripts/rebuild_predicate_maps.py --snapshot <dir>``.

Usage:
    python scripts/fetch_reference_snapshots.py [--dest data/reference/snapshots/2026-09-24] [--kegg-links]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.predicates.mappings.kegg_build import REQUEST_INTERVAL, REST, fetch_inputs


SOURCES = {
    "Pfam-A.clans.tsv.gz": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz",
    "Pfam.version.gz": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz",
    "pfam2go": "https://current.geneontology.org/ontology/external2go/pfam2go",
    "go-basic.obo": "https://purl.obolibrary.org/obo/go/go-basic.obo",
    "enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
    "enzclass.txt": "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
    "uniprot_sprot.dat.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
    "reldate.txt": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as response, open(tmp, "wb") as out:
        shutil.copyfileobj(response, out, 1 << 20)
    tmp.replace(dest)


def releases(dest: Path) -> dict[str, str]:
    found = {}
    if (dest / "reldate.txt").exists():
        first = (dest / "reldate.txt").read_text().splitlines()
        found["swissprot"] = next((line.split("Release")[1].split()[0] for line in first if "Swiss-Prot" in line), "")
    if (dest / "Pfam.version.gz").exists():
        text = gzip.decompress((dest / "Pfam.version.gz").read_bytes()).decode()
        found["pfam"] = next((line.split(":", 1)[1].strip() for line in text.splitlines()
                              if line.lower().startswith("pfam release")), "")
    return found


def fetch_kegg_links(dest: Path, log=print) -> None:
    organisms = set()
    with gzip.open(dest / "uniprot_sprot.dat.gz", "rt", encoding="latin-1") as handle:
        for line in handle:
            if line.startswith("DR   KEGG; "):
                organisms.add(line[11:].split(";")[0].split(":")[0])
    links = dest / "kegg_org_ko"
    links.mkdir(exist_ok=True)
    todo = sorted(o for o in organisms if not (links / f"{o}.tsv").exists())
    log(f"KEGG gene -> KO links: {len(organisms)} organisms, {len(todo)} to fetch")
    for i, org in enumerate(todo, 1):
        try:
            with urllib.request.urlopen(f"{REST}/link/ko/{org}", timeout=300) as response:
                (links / f"{org}.tsv").write_bytes(response.read())
        except OSError:
            (links / f"{org}.tsv").write_bytes(b"")  # KEGG returns nothing for some organisms
        time.sleep(REQUEST_INTERVAL)
        if i % 200 == 0:
            log(f"  {i}/{len(todo)}")
    subprocess.run([sys.executable, str(Path(__file__).with_name("link_swissprot_kegg.py")),
                    "--swissprot", str(dest / "uniprot_sprot.dat.gz"), "--org-ko-dir", str(links),
                    "--out", str(dest / "swissprot_kegg_ko.tsv"), "--retrieved", date.today().isoformat()],
                   check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dest", type=Path, default=Path("data/reference/snapshots") / date.today().isoformat())
    parser.add_argument("--kegg-links", action="store_true",
                        help="Also fetch KEGG gene -> KO links for Swiss-Prot organisms (slow; academic use)")
    parser.add_argument("--skip-kegg", action="store_true", help="Skip KEGG REST files")
    args = parser.parse_args()

    dest = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    manifest: dict = {"created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "files": {}}
    for name, url in SOURCES.items():
        path = dest / name
        if not path.exists():
            print(f"Fetching {name}")
            download(url, path)
        manifest["files"][name] = {"url": url, "sha256": sha256(path), "bytes": path.stat().st_size,
                                   "retrieved": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
                                   .isoformat(timespec="seconds")}
    manifest["releases"] = releases(dest)
    if not args.skip_kegg:
        print("KEGG REST files are for academic use and subject to KEGG's terms (https://www.kegg.jp/kegg/legal.html).")
        if not (dest / "kegg_inputs" / "ko_list.tsv").exists():
            fetch_inputs(dest / "kegg_inputs")
        manifest["kegg_inputs"] = "kegg_inputs/"
        if args.kegg_links:
            fetch_kegg_links(dest)
            manifest["files"]["swissprot_kegg_ko.tsv"] = {"sha256": sha256(dest / "swissprot_kegg_ko.tsv")}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Snapshot in {dest}: {json.dumps(manifest['releases'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
