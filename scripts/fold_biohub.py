#!/usr/bin/env python3
"""Fold proteins with ESMFold2 through the Biohub batch API, optionally with an MSA and ligands.

No GPU and no third-party packages: the script only talks HTTP (stdlib urllib), so it
runs on the Mac or the biotite login node without loading either.

Usage:
    # plain fold (fast model, endpoint "fold")
    python scripts/fold_biohub.py OUT_DIR input.faa [more.faa ...]

    # protein + ligands and/or MSA (full model, endpoint "fold_all_atom")
    python scripts/fold_biohub.py OUT_DIR one_protein.faa --msa family.a3m --ligand HEC:auto
    python scripts/fold_biohub.py OUT_DIR many.faa --msa-dir msas/ --ligand HEC:auto --ligand CA:2
    python scripts/fold_biohub.py OUT_DIR p.faa --smiles "CC(=O)O:1"

Options:
    --msa FILE          one MSA for a single input sequence (aligned FASTA or A3M)
    --msa-dir DIR       per-sequence MSAs named <name>.a3m / .fa / .fasta / .aln (name = sanitized
                        first header token, as used for the output files)
    --ligand CCD:N      add N copies of a CCD ligand (repeatable). HEC:auto adds one c-type heme
                        per CxxCH motif in each sequence.
    --smiles SMI:N      add N copies of a SMILES ligand (repeatable)
    --mod POS:CCD       residue modification, 1-based position, e.g. 349:PYL for pyrrolysine on a K
                        placeholder (repeatable). Writing O/U in the sequence gives an UNK residue with
                        backbone + CB only; the modification builds the real side chain.
    --model NAME        default esmfold2-fast-2026-05 for plain folds, esmfold2-2026-05 for all-atom
    --pae               request and save PAE
    --all-atom          use the all-atom endpoint and full model even for a plain fold (for like-for-like
                        comparisons with oligomer or ligand folds)
    --resume BATCH_ID   collect an earlier submission of the same inputs (after an interruption) without
                        resubmitting; models already in OUT_DIR are skipped
    --copies N          fold N copies of each protein as one homo-oligomer (chains A, B, C, ...); the
                        complex must stay within 4,096 tokens. Uses the all-atom endpoint.

MSA handling: lowercase A3M insertions are dropped; the row whose ungapped sequence equals the
query (else the first row) becomes the anchor; columns where the anchor has a gap are removed, so
every row has the query's length. Build MSAs on biotite (mafft on a protein family, or
colabfold_search); this script only reads them.

Auth: BIOHUB_API_KEY, or ESM_API_KEY (the same key works). On Jacob's Mac ESM_API_KEY is set in
the shell environment; on biotite it is not set (as of 2026-11-03).

Output per sequence: <name>.pdb (protein ATOM records plus ligand HETATM records on chains B, C, ...;
pLDDT x 100 in the B-factor column), <name>.json (model, pTM, ipTM, pLDDT, PAE if requested), and
summary.tsv. With HEC ligands, summary.tsv also reports how many hemes sit on a CxxCH motif
(motif His axial to Fe within 2.6 A and both Cys SG within 2.7 A of the heme vinyl carbons).

API facts (docs.biohub.ai, checked 2026-11-03): base https://biohub.ai/api/v1/batch; endpoints
fold / fold_all_atom; up to 4,096 tokens per item (ligand atoms count); 2,000 items per submission;
5,000 items per user per day shared across fold endpoints; a trailing '*' is rejected (stripped
here); results expire after 31 days. Covalent bonds are supported by the API but not exposed here;
hemes placed without them landed on their CxxCH motifs in 32/34 test cases.
"""
import argparse
import json
import math
import os
import re
import string
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
AA1 = {v: k for k, v in AA3.items()}
MAX_ITEMS = 2000
CHAIN_IDS = string.ascii_uppercase[1:] + string.ascii_lowercase + string.digits


def headers():
    key = os.environ.get("BIOHUB_API_KEY") or os.environ.get("ESM_API_KEY")
    if not key:
        sys.exit("set BIOHUB_API_KEY or ESM_API_KEY")
    return {"Authorization": f"Bearer {key}", "Content-Type": "application/json", "User-Agent": "sharur-fold/0.2"}


TRANSIENT_HTTP = (500, 502, 503, 504)


def post(path, body, H, attempts=10):
    """POST with retries: 429 waits for Retry-After; 5xx and dropped connections back off (15 s, doubling, max 5 min)."""
    for attempt in range(attempts):
        req = urllib.request.Request(BASE + path, data=json.dumps(body).encode(), headers=H, method="POST")
        wait = min(300, 15 * 2 ** attempt)
        try:
            return json.loads(urllib.request.urlopen(req, timeout=120).read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(int(e.headers.get("Retry-After", "30")))
                continue
            if e.code in TRANSIENT_HTTP and attempt < attempts - 1:
                print(f"HTTP {e.code} on {path}; retrying in {wait} s", file=sys.stderr, flush=True)
                time.sleep(wait)
                continue
            sys.exit(f"HTTP {e.code} on {path}: {e.read().decode()[:400]}")
        except (urllib.error.URLError, ConnectionError, TimeoutError) as e:
            if attempt < attempts - 1:
                print(f"{type(e).__name__} on {path}; retrying in {wait} s", file=sys.stderr, flush=True)
                time.sleep(wait)
                continue
            sys.exit(f"{type(e).__name__} on {path}: {e}")
    sys.exit(f"{path}: gave up after {attempts} attempts")


def read_fasta(paths):
    recs, h, s = [], None, []
    for p in paths:
        for line in open(p):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(s)))
                h, s = line[1:].strip(), []
            elif line.strip():
                s.append(line.strip())
        if h is not None:
            recs.append((h, "".join(s)))
            h, s = None, []
    return recs


def load_msa(path, query):
    rows = [re.sub(r"[a-z.]", "", s) for _, s in read_fasta([path])]  # drop A3M insertions
    rows = [r.upper().replace(".", "-") for r in rows if r]
    if not rows:
        sys.exit(f"empty MSA: {path}")
    anchor = next((i for i, r in enumerate(rows) if r.replace("-", "") == query), 0)
    cols = [k for k, a in enumerate(rows[anchor]) if a != "-"]
    out = ["".join(r[k] if k < len(r) else "-" for k in cols) for r in rows]
    out.insert(0, out.pop(anchor))
    if out[0] != query:
        print(f"  note: MSA anchor row differs from the query in {sum(a != b for a, b in zip(out[0], query))} positions "
              f"({path}); using the query as row 1", flush=True)
        out[0] = query
    if len({len(r) for r in out}) != 1 or len(out[0]) != len(query):
        sys.exit(f"MSA {path} does not map onto the query (length {len(query)})")
    return out


def parse_counts(specs, kind):
    out = []
    for s in specs or []:
        val, _, n = s.rpartition(":")
        if not val:
            sys.exit(f"--{kind} expects VALUE:N, got {s}")
        out.append((val, n))
    return out


def to_pdb_atom37(seq, coords, plddt):
    out, n = [], 1
    for i, (aa, res) in enumerate(zip(seq, coords)):
        for name, xyz in zip(ATOM37, res):
            if not xyz or any(v is None or (isinstance(v, float) and math.isnan(v)) for v in xyz):
                continue
            nm = f" {name:<3}" if len(name) < 4 else name
            out.append(f"ATOM  {n:5d} {nm} {AA3.get(aa, 'UNK')} A{i + 1:4d}    "
                       f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{plddt[i] * 100:6.2f}           {name[0]}")
            n += 1
    return "\n".join(out + ["TER", "END"]) + "\n", None


def to_pdb_complex(c, lig_names, protein_chains=("A",)):
    """Write the fold_all_atom `complex` object; return (pdb_text, atoms) for scoring."""
    chains = [c["metadata"]["chain_lookup"][str(x)] for x in c["chain_id"]]
    lines, atoms, n, resno = [], [], 1, {}
    for t, (a0, a1) in enumerate(c["token_to_atoms"]):
        ch = chains[t]
        resno[ch] = resno.get(ch, 0) + 1
        tok = c["sequence"][t]
        if ch in protein_chains:
            rn = tok if len(tok) == 3 else AA3.get(tok, "UNK")
        else:
            rn = lig_names.get(ch, "LIG")[:3]
        for k in range(a0, a1):
            x, y, z = c["atom_positions"][k]
            an, el = c["atom_names"][k], c["atom_elements"][k]
            nm = f" {an:<3}" if len(an) < 4 else an[:4]
            rec = "HETATM" if c["atom_hetero"][k] else "ATOM  "
            lines.append(f"{rec}{n:5d} {nm} {rn:>3} {ch}{resno[ch]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00"
                         f"{c['plddt'][t] * 100:6.2f}          {el:>2}")
            atoms.append((ch, resno[ch], rn, an, (x, y, z), c["plddt"][t]))
            n += 1
    return "\n".join(lines + ["END"]) + "\n", atoms


def heme_check(atoms, query):
    motif = {m.start() + 5: (m.start() + 1, m.start() + 4) for m in re.finditer(r"(?=C..CH)", query)}
    prot = [a for a in atoms if a[0] == "A"]
    ne2 = [(a[1], a[4]) for a in prot if a[2] == "HIS" and a[3] == "NE2"]
    sg = {a[1]: a[4] for a in prot if a[2] == "CYS" and a[3] == "SG"}
    n_hec, n_ok = 0, 0
    for ch in sorted({a[0] for a in atoms if a[0] != "A"}):
        het = {a[3]: a[4] for a in atoms if a[0] == ch}
        if "FE" not in het or not ({"CAB", "CAC"} & het.keys()):
            continue
        n_hec += 1
        ax = [r for r, p in ne2 if math.dist(het["FE"], p) < 2.6]
        mh = [r for r in ax if r in motif]
        if mh:
            vin = [het[k] for k in ("CAB", "CAC") if k in het]
            d = [min(math.dist(sg[c], v) for v in vin) for c in motif[mh[0]] if c in sg]
            n_ok += len(d) == 2 and max(d) < 2.7
    return n_hec, n_ok


def main():
    ap = argparse.ArgumentParser(description="Fold with ESMFold2 via the Biohub batch API (see module docstring).")
    ap.add_argument("out_dir")
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--model", default=None)
    ap.add_argument("--msa", help="MSA for a single input sequence (aligned FASTA or A3M)")
    ap.add_argument("--msa-dir", help="directory of per-sequence MSAs named <name>.a3m/.fa/.fasta/.aln")
    ap.add_argument("--ligand", action="append", help="CCD:N, repeatable; HEC:auto = one heme per CxxCH")
    ap.add_argument("--smiles", action="append", help="SMILES:N, repeatable")
    ap.add_argument("--mod", action="append", help="POS:CCD residue modification (1-based), repeatable")
    ap.add_argument("--pae", action="store_true", help="request and save the PAE matrix")
    ap.add_argument("--copies", type=int, default=1, help="fold N copies of each protein as a homo-oligomer")
    ap.add_argument("--all-atom", action="store_true", help="force the all-atom endpoint and full model")
    ap.add_argument("--poll", type=int, default=15, help="seconds between status checks")
    ap.add_argument("--resume", action="append", metavar="BATCH_ID",
                    help="collect an earlier submission of the same inputs instead of submitting again "
                         "(one batch id per 2,000-item chunk, in order); models already in OUT_DIR are skipped")
    a = ap.parse_args()
    H = headers()
    os.makedirs(a.out_dir, exist_ok=True)

    recs = [(h, s.rstrip("*").upper()) for h, s in read_fasta(a.fasta)]
    bad = [h for h, s in recs if "*" in s]
    if bad:
        sys.exit(f"internal stop codon in {len(bad)} sequences, e.g. {bad[0][:60]}")
    if a.msa and len(recs) != 1:
        sys.exit("--msa takes one input sequence; use --msa-dir for several")
    names, seen = [], {}
    for h, _ in recs:
        nm = re.sub(r"[^A-Za-z0-9_.-]", "_", h.split()[0])[:100] or "seq"
        seen[nm] = seen.get(nm, 0) + 1
        names.append(nm if seen[nm] == 1 else f"{nm}_{seen[nm]}")

    ligs, smis = parse_counts(a.ligand, "ligand"), parse_counts(a.smiles, "smiles")
    mods = [{"position": int(p) - 1, "ccd": c} for p, c in (m.split(":", 1) for m in (a.mod or []))]
    all_atom = bool(ligs or smis or a.msa or a.msa_dir or mods or a.copies > 1 or a.all_atom)
    prot_ids = list(string.ascii_uppercase[:max(1, a.copies)])
    lig_ids = [x for x in CHAIN_IDS if x not in prot_ids]
    model = a.model or ("esmfold2-2026-05" if all_atom else "esmfold2-fast-2026-05")
    endpoint = "fold_all_atom" if all_atom else "fold"

    payload, lig_maps = [], []
    for nm, (h, s) in zip(names, recs):
        if not all_atom:
            payload.append({"model": model, "sequence": s, "include_pae": a.pae})
            lig_maps.append({})
            continue
        prot = {"type": "protein", "id": "A", "sequence": s}
        if mods:
            prot["modifications"] = mods
        mpath = a.msa
        if a.msa_dir:
            mpath = next((os.path.join(a.msa_dir, nm + ext) for ext in (".a3m", ".fa", ".fasta", ".aln")
                          if os.path.exists(os.path.join(a.msa_dir, nm + ext))), None)
            if not mpath:
                print(f"  note: no MSA for {nm} in {a.msa_dir}; folding single-sequence", flush=True)
        if mpath:
            prot["msa"] = {"sequences": load_msa(mpath, s)}
        if a.copies * len(s) > 4096:
            sys.exit(f"{nm}: {a.copies} x {len(s)} residues exceeds the API's 4,096 tokens per complex")
        ents = [prot] + [{**prot, "id": cid} for cid in prot_ids[1:]]
        lmap, k = {}, 0
        for code, n in ligs:
            cnt = len(re.findall(r"(?=C..CH)", s)) if n == "auto" else int(n)
            for _ in range(cnt):
                cid = lig_ids[k]; k += 1
                ents.append({"type": "ligand", "id": cid, "ccd": [code]}); lmap[cid] = code
        for smi, n in smis:
            for _ in range(int(n)):
                cid = lig_ids[k]; k += 1
                ents.append({"type": "ligand", "id": cid, "smiles": smi}); lmap[cid] = "LIG"
        if len(ents) > 36:
            sys.exit(f"{nm}: {len(ents)} chains; the API allows at most 36 per complex (1 protein + 35 ligands). "
                     f"Lower the ligand count for this sequence.")
        payload.append({"model": model, "all_atom_input": {"sequences": ents}, "include_pae": a.pae,
                        "include_pair_chains_iptm": True})
        lig_maps.append(lmap)
        desc = ", ".join(f"{v}x{sum(1 for x in lmap.values() if x == v)}" for v in dict.fromkeys(lmap.values()))
        print(f"  {nm}: L={len(s)}{f' x {a.copies} copies' if a.copies > 1 else ''}{' + MSA ' + str(len(prot['msa']['sequences'])) + ' rows' if 'msa' in prot else ''}"
              f"{' + ' + desc if desc else ''}", flush=True)

    tasks = []
    resume = [b for arg in (a.resume or []) for b in arg.split(",") if b]
    chunks = list(range(0, len(payload), MAX_ITEMS))
    if resume and len(resume) != len(chunks):
        sys.exit(f"--resume needs {len(chunks)} batch id(s) for {len(payload)} inputs, got {len(resume)}")
    for k, start in enumerate(chunks):
        chunk = payload[start:start + MAX_ITEMS]
        if resume:
            batch = resume[k]
            print(f"batch {batch}: resuming {len(chunk)} items", flush=True)
        else:
            batch = post("/submit", {"endpoint": endpoint, "payload": chunk}, H)["data"]["batch_id"]
            print(f"batch {batch}: {len(chunk)} items, endpoint {endpoint}, model {model}", flush=True)
        tasks += [(start + i, f"{batch}-{i:05d}") for i in range(len(chunk))
                  if not (resume and os.path.exists(f"{a.out_dir}/{names[start + i]}.pdb"))]
    if resume:
        print(f"{len(payload) - len(tasks)} already in {a.out_dir}; collecting {len(tasks)}", flush=True)

    pending, rows, failed, retries = dict(tasks), [], [], {}
    while pending:
        time.sleep(a.poll)
        for i, tid in list(pending.items()):
            d = post("/status", {"task_id": tid}, H)["data"]
            st = d.get("status")
            if st in ("failed", "cancelled"):
                err = str(d.get("error", st))
                # Batch items are sometimes rejected against the 768-token on-demand limit; a resubmission succeeds.
                if "maximum allowed sequence length (768)" in err and retries.get(i, 0) < 3:
                    retries[i] = retries.get(i, 0) + 1
                    b2 = post("/submit", {"endpoint": endpoint, "payload": [payload[i]]}, H)["data"]["batch_id"]
                    print(f"retry {retries[i]} {names[i]} (transient 768-token rejection) -> {b2}", flush=True)
                    pending[i] = f"{b2}-00000"
                    continue
                failed.append((names[i], err))
                print(f"FAILED {names[i]}: {err}", flush=True)
                del pending[i]
                continue
            if st != "done":
                continue
            j = json.loads(urllib.request.urlopen(d["response"], timeout=120).read())
            h, s = recs[i]
            if all_atom:
                pdb, atoms = to_pdb_complex(j["complex"], lig_maps[i], tuple(prot_ids))
                ca = [x[5] for x in atoms if x[0] == "A" and x[3] == "CA"]
                mp = 100 * sum(ca) / max(1, len(ca))
                n_hec, n_ok = heme_check(atoms, s) if "HEC" in lig_maps[i].values() else (0, 0)
            else:
                pdb, _ = to_pdb_atom37(s, j["coordinates"], j["plddt"])
                mp = 100 * sum(j["plddt"]) / len(j["plddt"])
                n_hec = n_ok = 0
            open(f"{a.out_dir}/{names[i]}.pdb", "w").write(pdb)
            keep = {k: j.get(k) for k in ("model", "ptm", "interface_ptm", "plddt", "pae", "pair_chains_iptm")}
            json.dump(keep, open(f"{a.out_dir}/{names[i]}.json", "w"))
            iptm = j.get("interface_ptm")
            rows.append((names[i], len(s), round(mp, 1), round(j["ptm"], 3), "" if iptm is None else round(iptm, 3),
                         f"{n_ok}/{n_hec}" if n_hec else "", h))
            print(f"done {names[i]}  L={len(s)}  pLDDT={mp:.1f}  pTM={j['ptm']:.2f}"
                  f"{f'  ipTM={iptm:.2f}' if iptm is not None else ''}{f'  hemes on CxxCH {n_ok}/{n_hec}' if n_hec else ''}",
                  flush=True)
            del pending[i]

    summary = f"{a.out_dir}/summary.tsv"
    if resume:
        # Rows for models collected before the interruption, from their saved JSON.
        fresh = {r[0] for r in rows}
        for i, (h, s) in enumerate(recs):
            path = f"{a.out_dir}/{names[i]}.json"
            if names[i] in fresh or not os.path.exists(path):
                continue
            j = json.load(open(path))
            pl = j.get("plddt") or []
            flat = [x for row in pl for x in (row if isinstance(row, list) else [row])]
            mp = 100 * sum(flat) / max(1, len(flat))
            iptm = j.get("interface_ptm")
            rows.append((names[i], len(s), round(mp, 1), round(j["ptm"], 3), "" if iptm is None else round(iptm, 3), "", h))
    with open(summary, "w") as f:
        f.write("name\tlength\tprotein_plddt\tptm\tiptm\themes_on_cxxch\theader\n")
        for r in sorted(rows):
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} structures, {len(failed)} failed -> {a.out_dir}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
