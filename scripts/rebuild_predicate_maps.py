#!/usr/bin/env python3
"""Rebuild every predicate map from one reference snapshot, in dependency order.

0. SAM-dependent EC table from ENZYME reactions (``build_ec_sam_table.py``)
1. HydDB x KOfam snapshot (``build_kegg_hyddb_snapshot.py``; needs hmmsearch and
   KOfam profiles, skipped otherwise)
2. KO reviewed-protein consensus (``build_kegg_swissprot_consensus.py``; needs
   ``swissprot_kegg_ko.tsv`` from ``fetch_reference_snapshots.py --kegg-links``)
3. Local KEGG map (``sharur setup-kegg --inputs``)
4. Pfam map (``build_pfam_predicate_map.py``; uses the local KEGG map for the
   reviewed proteins' KO predicates)
5. CAZy map (``build_cazy_predicate_map.py``)

Dropped-proposal reports land in the snapshot directory. Afterwards run the
integrity recounts: ``make recount SNAPSHOT=<dir>``.

Usage:
    python scripts/rebuild_predicate_maps.py --snapshot data/reference/snapshots/2026-09-24 \\
        [--kofam-dir ~/.config/Astra/KOFAM]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parent


def run(*cmd: str) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--kofam-dir", type=Path, default=Path.home() / ".config/Astra/KOFAM")
    parser.add_argument("--kofam-ko-list", type=Path, default=Path("data/reference/ko_list"))
    parser.add_argument("--hmm", type=Path, default=Path.home() / ".config/Astra/PFAM/Pfam-A.hmm",
                        help="Installed Pfam-A.hmm (families reported by annotation)")
    args = parser.parse_args()

    snap = args.snapshot
    manifest = json.loads((snap / "manifest.json").read_text()) if (snap / "manifest.json").exists() else {}
    release = manifest.get("releases", {}).get("swissprot") or "unknown"
    py = sys.executable
    swissprot, obo = snap / "uniprot_sprot.dat.gz", snap / "go-basic.obo"
    links = snap / "swissprot_kegg_ko.tsv"
    kofam_list = args.kofam_ko_list if args.kofam_ko_list.exists() else snap / "kegg_inputs/kofam_ko_list"

    run(py, str(SCRIPTS / "build_ec_sam_table.py"), "--enzyme-dat", str(snap / "enzyme.dat"))

    if shutil.which("hmmsearch") and args.kofam_dir.is_dir():
        run(py, str(SCRIPTS / "build_kegg_hyddb_snapshot.py"), "--ko-list", str(snap / "kegg_inputs/ko_list.tsv"),
            "--kofam-dir", str(args.kofam_dir), "--kofam-ko-list", str(kofam_list))
    else:
        print("Skipping HydDB x KOfam snapshot (needs hmmsearch and KOfam profiles); keeping the shipped one.")

    if links.exists():
        run(py, str(SCRIPTS / "build_kegg_swissprot_consensus.py"), "--swissprot", str(swissprot),
            "--swissprot-release", release, "--swissprot-kegg", str(links), "--go-obo", str(obo))
    else:
        print("Skipping KO consensus (no swissprot_kegg_ko.tsv; fetch with --kegg-links); keeping the shipped table.")

    run("sharur", "setup-kegg", "--inputs", str(snap / "kegg_inputs"), "--kofam-ko-list", str(kofam_list),
        "--report", str(snap / "kegg_dropped_pairs.tsv"))

    pfam = [py, str(SCRIPTS / "build_pfam_predicate_map.py"), "--pfam-hmm", str(args.hmm),
            "--pfam-clans", str(snap / "Pfam-A.clans.tsv.gz"), "--pfam2go", str(snap / "pfam2go"),
            "--go-obo", str(obo), "--enzyme-dat", str(snap / "enzyme.dat"), "--swissprot", str(swissprot),
            "--swissprot-release", release, "--report", str(snap / "pfam_dropped_pairs.tsv")]
    if links.exists():
        pfam += ["--swissprot-kegg", str(links)]
    run(*pfam)

    run(py, str(SCRIPTS / "build_cazy_predicate_map.py"), "--swissprot", str(swissprot),
        "--swissprot-release", release, "--go-obo", str(obo), "--report", str(snap / "cazy_dropped_pairs.tsv"))
    print(f"\nMaps rebuilt from {snap}. Verify: make recount SNAPSHOT={snap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
