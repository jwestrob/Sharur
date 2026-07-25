#!/usr/bin/env python
"""Build per-bin dossiers for the SR-VP archaeal Stage 2 swarm.

Bulk-query strategy: one query per metric across all bins, then assemble
dossiers in-memory. Much faster than per-bin querying (14k+ small queries).

Outputs:
  data/srvp_archaea/swarm/_data/bin_dossiers.jsonl    (all 967 bins)
  data/srvp_archaea/swarm/_data/chunks/chunk_NN.jsonl (20 chunks)
  data/srvp_archaea/swarm/_context/shared_context.md
"""
from __future__ import annotations
import duckdb
import json
import math
import re
import sys
from pathlib import Path
from collections import defaultdict

DATASET = Path("data/srvp_archaea")
DB = DATASET / "sharur.duckdb"
OUT_DATA = DATASET / "swarm" / "_data"
OUT_CHUNKS = OUT_DATA / "chunks"
OUT_CTX = DATASET / "swarm" / "_context"
N_READERS = 20

ARCHAEAL_MARKERS = {
    "K03044": "rpoA1", "K13800": "rpoB", "K03046": "rpoC1",
    "K02945": "rpsA", "K02982": "rpsC", "K02935": "rpsL",
    "K02863": "rpl1", "K03553": "recA", "K07466": "radA",
    "K02356": "efp", "K02117": "atpA", "K07151": "stt3-OST",
}
MCR_SUBUNIT_PFAMS = ["MCR_alpha", "MCR_beta_N", "MCR_gamma", "MCR_anchor"]
ASGARD_ESPS = ('Actin', 'Profilin', 'Gelsolin', 'ADF', 'Tubulin',
               'ESCRT-III', 'Snf7', 'Vps20_like', 'Bro1', 'Vps4_C',
               'Ras', 'RhoGAP', 'RhoGEF', 'ArfGap', 'Sec7', 'Roadblock',
               'Ubiquitin', 'Rad60-SLD', 'UBA_e1_C', 'UQ_con', 'Ubc', 'OTU')
DPANN_MARKERS = ("K02982", "K02935", "K02356", "K03470")


def collect_findings_index(survey_dir: Path, exploration_dir: Path) -> dict:
    idx = defaultdict(list)
    for d in (survey_dir, exploration_dir):
        if not d.exists():
            continue
        for p in sorted(d.glob("*_findings.jsonl")):
            for line in p.read_text().splitlines():
                if not line.strip():
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                fid = rec.get("id") or rec.get("finding_id")
                if not fid:
                    continue
                headline = (rec.get("title") or rec.get("description", ""))[:120]
                bin_ids = []
                for k in ("bin_ids", "affected_bins", "bins"):
                    v = rec.get(k)
                    if isinstance(v, list):
                        bin_ids += [str(b) for b in v]
                desc = rec.get("description", "") + " " + str(rec.get("evidence", ""))
                for m in re.finditer(r"SR_VP_[A-Za-z0-9_]+", desc):
                    bin_ids.append(m.group(0))
                for bid in set(bin_ids):
                    idx[bid].append({"finding_id": fid, "note": headline.strip()})
    return idx


def main():
    if not DB.exists():
        sys.exit(f"DB not found: {DB}")
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_CHUNKS.mkdir(parents=True, exist_ok=True)
    OUT_CTX.mkdir(parents=True, exist_ok=True)

    print(f"[1/12] connect {DB} read-only", flush=True)
    con = duckdb.connect(str(DB), read_only=True)

    print("[2/12] index Stage 1 findings → bin_id caveats", flush=True)
    findings_idx = collect_findings_index(DATASET / "survey", DATASET / "exploration")
    print(f"       caveats for {len(findings_idx)} bin_ids "
          f"(total {sum(len(v) for v in findings_idx.values())} caveat instances)", flush=True)

    print("[3/12] fetch bin list", flush=True)
    bins = con.execute("""
        SELECT bin_id, COALESCE(taxonomy,'unclassified') as taxonomy,
               total_length, n_contigs
        FROM bins ORDER BY taxonomy NULLS LAST, bin_id
    """).fetchall()
    n_bins = len(bins)
    print(f"       {n_bins} bins", flush=True)

    # Bulk: protein counts per bin
    print("[4/12] bulk: protein counts", flush=True)
    pc = dict(con.execute("SELECT bin_id, COUNT(*) FROM proteins GROUP BY bin_id").fetchall())

    # Bulk: top PFAM / KEGG per bin (use a window function to limit to top 10 per bin)
    print("[5/12] bulk: top-10 PFAM per bin", flush=True)
    pfam_rows = con.execute("""
        WITH ranked AS (
          SELECT p.bin_id, a.name, COUNT(DISTINCT a.protein_id) n,
                 ROW_NUMBER() OVER (PARTITION BY p.bin_id ORDER BY COUNT(DISTINCT a.protein_id) DESC) rk
          FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
          WHERE a.source = 'pfam'
          GROUP BY p.bin_id, a.name
        )
        SELECT bin_id, name, n FROM ranked WHERE rk <= 10
        ORDER BY bin_id, rk
    """).fetchall()
    top_pfam = defaultdict(list)
    for bid, name, n in pfam_rows:
        top_pfam[bid].append({"name": name, "n": int(n)})

    print("[6/12] bulk: top-10 KEGG per bin", flush=True)
    kegg_rows = con.execute("""
        WITH ranked AS (
          SELECT p.bin_id, a.name, COUNT(DISTINCT a.protein_id) n,
                 ROW_NUMBER() OVER (PARTITION BY p.bin_id ORDER BY COUNT(DISTINCT a.protein_id) DESC) rk
          FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
          WHERE a.source = 'kofam'
          GROUP BY p.bin_id, a.name
        )
        SELECT bin_id, name, n FROM ranked WHERE rk <= 10
        ORDER BY bin_id, rk
    """).fetchall()
    top_kegg = defaultdict(list)
    for bid, name, n in kegg_rows:
        top_kegg[bid].append({"name": name, "n": int(n)})

    # Bulk: archaeal marker presence per bin
    print("[7/12] bulk: archaeal marker presence", flush=True)
    marker_kos = list(ARCHAEAL_MARKERS.keys())
    placeholders = ",".join("?" * len(marker_kos))
    marker_rows = con.execute(f"""
        SELECT p.bin_id, a.name
        FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = 'kofam' AND a.name IN ({placeholders})
        GROUP BY p.bin_id, a.name
    """, marker_kos).fetchall()
    marker_present = defaultdict(set)
    for bid, ko in marker_rows:
        marker_present[bid].add(ko)

    # Bulk: defense systems
    print("[8/12] bulk: defense systems", flush=True)
    ds_rows = con.execute("""
        SELECT genome_id, system_type, COUNT(*) FROM defense_systems
        GROUP BY genome_id, system_type
    """).fetchall()
    defense_by_bin = defaultdict(lambda: {"total": 0, "types": set()})
    for gid, stype, n in ds_rows:
        defense_by_bin[gid]["total"] += int(n)
        defense_by_bin[gid]["types"].add(stype)

    # Bulk: HydDB subgroups
    print("[9/12] bulk: HydDB subgroups", flush=True)
    hyddb_rows = con.execute("""
        SELECT p.bin_id, sa.atom_id, COUNT(DISTINCT sa.protein_id)
        FROM semantic_atoms sa JOIN proteins p ON sa.protein_id = p.protein_id
        WHERE sa.source_db = 'hyddb_subgroup'
        GROUP BY p.bin_id, sa.atom_id
    """).fetchall()
    hyddb_by_bin = defaultdict(dict)
    for bid, atom, n in hyddb_rows:
        hyddb_by_bin[bid][atom] = int(n)

    # Bulk: CRISPR arrays
    print("[10/12] bulk: CRISPR arrays", flush=True)
    crispr_rows = con.execute("""
        SELECT c.bin_id, COUNT(*) FROM loci l JOIN contigs c ON l.contig_id=c.contig_id
        WHERE l.locus_type = 'crispr' GROUP BY c.bin_id
    """).fetchall()
    crispr_by_bin = {bid: int(n) for bid, n in crispr_rows}

    # Bulk: archaea-aware PFAM facets in one query
    print("[11/12] bulk: archaea-aware PFAM facets (MCR, AMO, F420, Molybdopterin, HAO, ESPs)", flush=True)
    arch_rows = con.execute("""
        SELECT p.bin_id, a.name, COUNT(DISTINCT a.protein_id)
        FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = 'pfam' AND (
              a.name IN ('MCR_alpha','MCR_beta_N','MCR_gamma','MCR_anchor',
                         'Archaeal_AmoA','F420_oxidored','LLM',
                         'Actin','Profilin','Gelsolin','ADF','Tubulin',
                         'ESCRT-III','Snf7','Vps20_like','Bro1','Vps4_C',
                         'Ras','RhoGAP','RhoGEF','ArfGap','Sec7','Roadblock',
                         'Ubiquitin','Rad60-SLD','UBA_e1_C','UQ_con','Ubc','OTU')
              OR a.name LIKE 'Molybdopterin%'
              OR a.name LIKE 'Hydroxylam%'
        )
        GROUP BY p.bin_id, a.name
    """).fetchall()
    arch_pfam = defaultdict(dict)
    for bid, name, n in arch_rows:
        arch_pfam[bid][name] = int(n)

    # Bulk: HAO via KEGG K10535
    hao_kegg_rows = con.execute("""
        SELECT p.bin_id, COUNT(DISTINCT a.protein_id)
        FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = 'kofam' AND a.name = 'K10535'
        GROUP BY p.bin_id
    """).fetchall()
    hao_kegg_by_bin = {bid: int(n) for bid, n in hao_kegg_rows}

    # Bulk: DPANN markers
    dpann_rows = con.execute(f"""
        SELECT p.bin_id, COUNT(DISTINCT a.name)
        FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = 'kofam' AND a.name IN ({','.join(['?']*len(DPANN_MARKERS))})
        GROUP BY p.bin_id
    """, list(DPANN_MARKERS)).fetchall()
    dpann_by_bin = {bid: int(n) for bid, n in dpann_rows}

    # Assemble + write
    print("[12/12] assemble + write dossiers", flush=True)
    dossier_path = OUT_DATA / "bin_dossiers.jsonl"
    with dossier_path.open("w") as f:
        for i, (bid, tax, tlen, nc) in enumerate(bins):
            tlen = int(tlen or 0)
            nc = int(nc or 0)
            if nc == 1:
                closure = "single_contig"
            elif nc <= 10:
                closure = "lightly_fragmented"
            elif nc <= 50:
                closure = "fragmented"
            else:
                closure = "heavily_fragmented"

            ap = arch_pfam.get(bid, {})
            mcr_counts = {k: ap[k] for k in MCR_SUBUNIT_PFAMS if k in ap}
            amo = ap.get("Archaeal_AmoA", 0)
            f420 = ap.get("F420_oxidored", 0) + ap.get("LLM", 0)
            moly = sum(v for k, v in ap.items() if k.startswith("Molybdopterin"))
            hao_pfam = sum(v for k, v in ap.items() if k.startswith("Hydroxylam"))
            hao_total = hao_pfam + hao_kegg_by_bin.get(bid, 0)
            esp_count = sum(ap.get(k, 0) for k in ASGARD_ESPS)

            ds = defense_by_bin.get(bid, {"total": 0, "types": set()})

            caveats = findings_idx.get(bid, [])[:8]

            dossier = {
                "bin_id": bid,
                "taxonomy": tax or "unclassified",
                "total_length_bp": tlen,
                "n_contigs": nc,
                "n_proteins": int(pc.get(bid, 0)),
                "closure": closure,
                "top_pfam": top_pfam.get(bid, []),
                "top_kegg": top_kegg.get(bid, []),
                "marker_score": len(marker_present.get(bid, set())),
                "marker_kos_present": sorted(marker_present.get(bid, set())),
                "defense_systems": {
                    "total": int(ds["total"]),
                    "n_distinct_types": len(ds["types"]),
                    "types": sorted(ds["types"]),
                },
                "hyddb_subgroups": hyddb_by_bin.get(bid, {}),
                "n_crispr_arrays": crispr_by_bin.get(bid, 0),
                "mcr_subunit_counts": mcr_counts,
                "amo_pfam_count": amo,
                "hao_count": hao_total,
                "f420_dependent_enzyme_count": f420,
                "molybdopterin_alpha_count": moly,
                "asgard_esp_signature_count": esp_count,
                "dpann_marker_count": dpann_by_bin.get(bid, 0),
                "caveat_notes": caveats,
                "n_findings_cited": len(caveats),
            }
            f.write(json.dumps(dossier) + "\n")
    print(f"       wrote {n_bins} dossiers to {dossier_path}", flush=True)

    # Build chunks (group bins phylum-balanced for better reader workload)
    bins_per_chunk = math.ceil(n_bins / N_READERS)
    print(f"[13/13] build {N_READERS} chunks (~{bins_per_chunk} bins each)", flush=True)
    with dossier_path.open() as src:
        lines = src.readlines()
    for r in range(N_READERS):
        chunk_path = OUT_CHUNKS / f"chunk_{r+1:02d}.jsonl"
        with chunk_path.open("w") as out:
            out.writelines(lines[r*bins_per_chunk:(r+1)*bins_per_chunk])
        n_in_chunk = sum(1 for _ in chunk_path.open())
        print(f"       chunk_{r+1:02d}: {n_in_chunk} bins", flush=True)

    # Shared context
    ctx_path = OUT_CTX / "shared_context.md"
    ctx_path.write_text(f"""# SR-VP Archaeal Swarm — Shared Context

You are one of {N_READERS} parallel readers. Each handles ~{bins_per_chunk} bins from the Serpens Ridge Vernal Pool archaeal metagenome and produces a per-organism digest. NO topic specialization — every reader uses the same prompt. Read all bins in your assigned dossier carefully; don't skim.

## Dataset facts

- **{n_bins} archaeal MAGs** from SR-VP, mixed Illumina + PacBio HiFi (Abawaca2, Maxbin2, MetaBAT2, hifiasm-meta).
- **GTDB-Tk r220 ar53** (Luis, 2024-07-10): Thermoproteota 599 (62%), Thermoplasmatota 139, Halobacteriota 105 (mostly Methanoperedens / ANME-2c), Aenigmatarchaeota 41, Nanoarchaeota 23, Micrarchaeota 19, **Asgardarchaeota 14**, Hadarchaeota 12, B1Sed10-29 7, EX4484-52 4, Methanobacteriota_B 2, SpSt-1190 1, p__ 1.
- **Bin sizes**: 0.3–5.4 Mb. Asgard largest (median 2.80 Mb), Nanoarchaeota smallest (0.63 Mb).
- **Total proteins**: 1,988,970. **V2 atoms**: 12.5M.
- **Defense** (MacSyFinder-validated): 1,902 across 666 bins.
- **CRISPR** (MinCED): 410 arrays across 234 bins.

## Stage 1 already analyzed

Reports + JSONLs in `data/srvp_archaea/survey/` and `data/srvp_archaea/exploration/`. Narrative synthesis: `data/srvp_archaea/reports/STORY.md`. Read STORY.md before starting.

Headline biology you should know:
- **Methanoperedens (ANME-2c)** dominates the mcrA signal — 71 bins, reverse methanogenesis, defense champions
- **Helarchaeales (Asgard) bin with 3 independent mcrBGA operons** — candidate short-chain alkane oxidizer
- **44 Nitrosotalea AOA bins** with conserved 121 aa unannotated AMO factor in amoA-amoC gap
- **Pan-Asgard gelsolin** + actin+gelsolin+Snf7 trifecta + Njordarchaeia sampylation operon
- **Hadarchaeota outlier** — 12 metabolically self-sufficient bins; 2 candidate sulfate reducers
- **DPANN mega-adhesins** (up to 11,391 aa NosD/Beta_helix); **13,831 aa fully-dark Aenigmatarchaeota surface glycoprotein** with 11,595 aa paralog

## Annotation traps to avoid

- **Cas12f / TnpB / IS-OrfB**: PFAM `Cas12f*` is TnpB by default — DO NOT report as Cas. See `_validation_protocols.md` §8.
- **Molybdopterin α-subunit ≠ NXR/NAR/FDH function call** — fold superfamily; substrate needs phylogeny.
- **HydDB → Complex I (NuoH-N) cross-classification** in aerobic AOA bins — check ±8 gene window for NuoH/L/M.
- **Cytochrom_C_asm (CcmF/H/I)** is cyt c maturation machinery, not EET.
- **K10944** matches both archaeal amoA and bacterial pmoA — use PFAM `Archaeal_AmoA` to disambiguate.

## Dossier fields

| Field | Meaning |
|---|---|
| `bin_id` | bin name |
| `taxonomy` | GTDB r220 ar53 lineage |
| `total_length_bp` / `n_contigs` / `n_proteins` | basic stats |
| `closure` | single_contig / lightly_fragmented / fragmented / heavily_fragmented |
| `marker_score` | archaeal SCG count present (0-12) |
| `marker_kos_present` | list of present marker KO IDs |
| `top_pfam` / `top_kegg` | top-10 by distinct protein count |
| `defense_systems` | MacSyFinder-validated totals + types |
| `hyddb_subgroups` | NiFe groups 1-4, FeFe A-C |
| `n_crispr_arrays` | MinCED arrays |
| `mcr_subunit_counts` | dict of MCR_alpha/MCR_beta_N/MCR_gamma/MCR_anchor → count |
| `amo_pfam_count` | `Archaeal_AmoA` proteins |
| `hao_count` | PFAM Hydroxylam% + KEGG K10535 |
| `f420_dependent_enzyme_count` | PFAM F420_oxidored + LLM |
| `molybdopterin_alpha_count` | PFAM Molybdopterin% — FOLD count, not function |
| `asgard_esp_signature_count` | ESP-family PFAM proteins |
| `dpann_marker_count` | small DPANN-core marker count |
| `caveat_notes` | Stage 1 finding IDs + headlines involving this bin |
| `n_findings_cited` | length of caveat_notes |

The `caveat_notes` field is your guard against propagating Stage 1 omissions. When a finding flagged a KO/annotation gap for a bin, it appears here. Address those caveats directly in your digest.

## Your job

For each bin, write a short paragraph (5-12 sentences) describing what's interesting about that organism. Combine dossier fields with knowledge of archaeal biology. If the bin has Stage 1 caveat notes, address them. Flag anomalies not yet in the Stage 1 catalog.

Output:
- `data/srvp_archaea/swarm/level1_digests/reader_NN.md` — per-bin digests
- `data/srvp_archaea/swarm/level1_digests/reader_NN_findings.jsonl` — anything new (anomaly, candidate, gap)

Follow `docs/findings_spec.md` for new findings — verification SQL, falsification for novelty≥2, no forbidden language.
""")
    print(f"       wrote shared context at {ctx_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
