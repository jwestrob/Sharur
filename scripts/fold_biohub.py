#!/usr/bin/env python3
"""Fold proteins with ESMFold2 through the Biohub batch API and write PDB files.

No GPU and no third-party packages: the script only talks HTTP (stdlib urllib), so it
runs on the Mac or the biotite login node without loading either.

Usage:
    python scripts/fold_biohub.py OUT_DIR input.faa [more.faa ...] [--pae] [--model NAME]

Auth: BIOHUB_API_KEY, or ESM_API_KEY (the same key works). On Jacob's Mac ESM_API_KEY is
set in the shell environment; on biotite it is not set (as of 2026-11-03).

Output, per sequence (named from the first header token, sanitized):
    OUT_DIR/<name>.pdb    atom37 coordinates, pLDDT x 100 in the B-factor column
    OUT_DIR/<name>.json   model, pTM, per-residue pLDDT, sequence (and PAE with --pae)
    OUT_DIR/summary.tsv   name, length, mean pLDDT, pTM, full header

API facts (docs.biohub.ai/api/protein/models/batch_overview, checked 2026-11-03):
    - base https://biohub.ai/api/v1/batch, endpoints submit/list/status/cancel, Bearer auth
    - model esmfold2-fast-2026-05; up to 4,096 tokens per sequence; 2,000 items per submit
    - quota 5,000 fold items per user per day; 429 responses carry Retry-After
    - each finished task returns a signed URL to a JSON with `coordinates` (atom37, null for
      absent atoms), `plddt` (0-1 scale) and `ptm`; there is no PDB string
    - sequences may contain only amino-acid letters (plus X, B, Z, U, O, ':' and '|'); a
      trailing '*' from gene callers is rejected with 422, so this script strips it
    - results are deleted 31 days after submission
"""
import argparse
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.request

BASE = "https://biohub.ai/api/v1/batch"
ATOM37 = ["N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD", "CD1", "CD2", "ND1", "ND2",
          "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2",
          "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT"]
AA3 = dict(A="ALA", R="ARG", N="ASN", D="ASP", C="CYS", Q="GLN", E="GLU", G="GLY", H="HIS", I="ILE", L="LEU",
           K="LYS", M="MET", F="PHE", P="PRO", S="SER", T="THR", W="TRP", Y="TYR", V="VAL")
MAX_ITEMS = 2000


def headers():
    key = os.environ.get("BIOHUB_API_KEY") or os.environ.get("ESM_API_KEY")
    if not key:
        sys.exit("set BIOHUB_API_KEY or ESM_API_KEY")
    return {"Authorization": f"Bearer {key}", "Content-Type": "application/json", "User-Agent": "sharur-fold/0.1"}


def post(path, body, H):
    for _ in range(8):
        req = urllib.request.Request(BASE + path, data=json.dumps(body).encode(), headers=H, method="POST")
        try:
            return json.loads(urllib.request.urlopen(req, timeout=120).read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(int(e.headers.get("Retry-After", "30")))
                continue
            sys.exit(f"HTTP {e.code} on {path}: {e.read().decode()[:300]}")
    sys.exit("rate-limited repeatedly")


def read_fasta(paths):
    recs, h, s = [], None, []
    for p in paths:
        for line in open(p):
            line = line.strip()
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(s)))
                h, s = line[1:], []
            elif line:
                s.append(line)
        if h is not None:
            recs.append((h, "".join(s)))
            h, s = None, []
    return recs


def to_pdb(seq, coords, plddt):
    out, n = [], 1
    for i, (aa, res) in enumerate(zip(seq, coords)):
        for name, xyz in zip(ATOM37, res):
            if not xyz or any(v is None or (isinstance(v, float) and math.isnan(v)) for v in xyz):
                continue
            nm = f" {name:<3}" if len(name) < 4 else name
            out.append(f"ATOM  {n:5d} {nm} {AA3.get(aa, 'UNK')} A{i + 1:4d}    "
                       f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{plddt[i] * 100:6.2f}           {name[0]}")
            n += 1
    return "\n".join(out + ["TER", "END"]) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out_dir")
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--model", default="esmfold2-fast-2026-05")
    ap.add_argument("--pae", action="store_true", help="request and save the PAE matrix")
    ap.add_argument("--poll", type=int, default=15, help="seconds between status checks")
    a = ap.parse_args()
    H = headers()
    os.makedirs(a.out_dir, exist_ok=True)

    recs = [(h, s.rstrip("*")) for h, s in read_fasta(a.fasta)]
    bad = [h for h, s in recs if "*" in s]
    if bad:
        sys.exit(f"internal stop codon in {len(bad)} sequences, e.g. {bad[0][:60]}")
    names, seen = [], {}
    for h, _ in recs:
        nm = re.sub(r"[^A-Za-z0-9_.-]", "_", h.split()[0])[:100] or "seq"
        seen[nm] = seen.get(nm, 0) + 1
        names.append(nm if seen[nm] == 1 else f"{nm}_{seen[nm]}")

    tasks = []  # (index, task_id)
    for start in range(0, len(recs), MAX_ITEMS):
        chunk = recs[start:start + MAX_ITEMS]
        body = {"endpoint": "fold",
                "payload": [{"model": a.model, "sequence": s, "include_pae": a.pae} for _, s in chunk]}
        batch = post("/submit", body, H)["data"]["batch_id"]
        print(f"batch {batch}: {len(chunk)} items", flush=True)
        tasks += [(start + i, f"{batch}-{i:05d}") for i in range(len(chunk))]

    pending, rows, failed = dict(tasks), [], []
    while pending:
        time.sleep(a.poll)
        for i, tid in list(pending.items()):
            d = post("/status", {"task_id": tid}, H)["data"]
            st = d.get("status")
            if st in ("failed", "cancelled"):
                failed.append((names[i], d.get("error", st)))
                print(f"FAILED {names[i]}: {d.get('error', st)}", flush=True)
                del pending[i]
                continue
            if st != "done":
                continue
            j = json.loads(urllib.request.urlopen(d["response"], timeout=120).read())
            h, s = recs[i]
            open(f"{a.out_dir}/{names[i]}.pdb", "w").write(to_pdb(s, j["coordinates"], j["plddt"]))
            keep = {k: j.get(k) for k in ("model", "ptm", "plddt", "sequence", "pae")}
            json.dump(keep, open(f"{a.out_dir}/{names[i]}.json", "w"))
            mp = 100 * sum(j["plddt"]) / len(j["plddt"])
            rows.append((names[i], len(s), round(mp, 1), round(j["ptm"], 3), h))
            print(f"done {names[i]}  L={len(s)}  pLDDT={mp:.1f}  pTM={j['ptm']:.2f}", flush=True)
            del pending[i]

    with open(f"{a.out_dir}/summary.tsv", "w") as f:
        f.write("name\tlength\tmean_plddt\tptm\theader\n")
        for r in sorted(rows):
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} structures, {len(failed)} failed -> {a.out_dir}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
