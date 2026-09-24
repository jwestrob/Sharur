#!/usr/bin/env python3
"""Link Swiss-Prot's ``DR KEGG`` genes to KEGG orthologs.

Reads the KEGG gene IDs referenced by a UniProtKB/Swiss-Prot flat file and the
KEGG REST ``link/ko/<organism>`` files for those organisms, and writes
``gene<TAB>KO[,KO]`` rows for the referenced genes. The output feeds
``--swissprot-kegg`` in ``scripts/build_pfam_predicate_map.py`` and
``scripts/build_kegg_predicate_map.py``.

Fetch the organism files first, one request at a time, e.g.:
    for org in $(organisms); do curl -sf https://rest.kegg.jp/link/ko/$org -o org_ko/$org.tsv; done

Usage:
    python scripts/link_swissprot_kegg.py --swissprot uniprot_sprot.dat.gz \\
        --org-ko-dir org_ko/ --out swissprot_kegg_ko.tsv
"""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path


def referenced_genes(swissprot: Path) -> set[str]:
    genes = set()
    with gzip.open(swissprot, "rt", encoding="latin-1") as handle:
        for line in handle:
            if line.startswith("DR   KEGG; "):
                genes.add(line[11:].split(";")[0])
    return genes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--swissprot", type=Path, required=True)
    parser.add_argument("--org-ko-dir", type=Path, required=True, help="Directory of <organism>.tsv link/ko files")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--retrieved", default="", help="Retrieval date recorded in the header")
    args = parser.parse_args()

    genes = referenced_genes(args.swissprot)
    organisms = sorted({g.split(":")[0] for g in genes})
    missing = [o for o in organisms if not (args.org_ko_dir / f"{o}.tsv").exists()]
    if missing:
        raise SystemExit(f"{len(missing)} organisms lack link files, e.g. {missing[:5]}")

    links: dict[str, set[str]] = defaultdict(set)
    for org in organisms:
        with open(args.org_ko_dir / f"{org}.tsv") as handle:
            for line in handle:
                gene, _, ko = line.rstrip("\n").partition("\t")
                if gene in genes:
                    links[gene].add(ko.removeprefix("ko:"))

    with open(args.out, "w") as out:
        out.write(f"# KEGG REST link/ko for {len(organisms)} organisms referenced by {args.swissprot.name}"
                  f"{'; retrieved ' + args.retrieved if args.retrieved else ''}\n")
        for gene in sorted(links):
            out.write(f"{gene}\t{','.join(sorted(links[gene]))}\n")
    print(f"{len(genes):,} referenced genes; {len(links):,} linked to KOs across {len(organisms):,} organisms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
