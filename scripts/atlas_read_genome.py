#!/usr/bin/env python3
"""Atlas genome reader — exhaustive single-genome inventory for the Atlas skill."""
import json
import os
import sys

# Known superfamily traps: annotations that describe protein folds, not specific functions
SUPERFAMILY_ACCESSIONS = {
    "PF00069", "PF07714",  # Pkinase / PK_Tyr_Ser-Thr
    "PF00515",  # TPR
    "PF00400",  # WD40
    "PF04055",  # Radical_SAM
    "PF00346",  # Complex1_49kDa (MBL-fold)
    "PF13432", "PF14559", "PF13181",  # TPR variants
    "PF00023",  # Ankyrin
    "PF07719",  # TPR_2
}

SUPERFAMILY_DEFENSE = {
    "Mokosh_type_I__MkoA_B",
    "BREX__pglW",
}


def read_genome(genome, short_name, db_path="data/susan_genomes/sharur.duckdb"):
    from sharur.operators import Sharur
    b = Sharur(db_path)

    out_dir = "data/susan_genomes/atlas"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Basic stats
    total_proteins = b.store.execute(
        f"SELECT COUNT(*) FROM proteins WHERE bin_id = '{genome}'"
    )[0][0]
    n_contigs = b.store.execute(
        f"SELECT COUNT(DISTINCT contig_id) FROM proteins WHERE bin_id = '{genome}'"
    )[0][0]
    avg_len = b.store.execute(
        f"SELECT ROUND(AVG(sequence_length), 0) FROM proteins WHERE bin_id = '{genome}'"
    )[0][0]
    max_len = b.store.execute(
        f"SELECT MAX(sequence_length) FROM proteins WHERE bin_id = '{genome}'"
    )[0][0]

    # 2. Annotation source census
    census = b.store.execute(f"""
        SELECT a.source, COUNT(DISTINCT a.protein_id) as n_proteins
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = '{genome}'
        GROUP BY a.source
        ORDER BY n_proteins DESC
    """)

    # 3. Top 30 annotations
    top_annots = b.store.execute(f"""
        SELECT a.source, a.accession, a.name,
               COUNT(DISTINCT a.protein_id) as n
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = '{genome}'
        GROUP BY a.source, a.accession, a.name
        ORDER BY n DESC
        LIMIT 30
    """)

    # 4. Context-first checks for high-hit annotations
    context_checks = []
    flags = []

    for t in top_annots:
        source, accession, name, count = t
        if count < 10:
            break  # sorted desc, so all remaining are <10 too

        # Check if this is a known superfamily
        is_superfamily = accession in SUPERFAMILY_ACCESSIONS or accession in SUPERFAMILY_DEFENSE

        if is_superfamily or count > total_proteins * 0.02:  # >2% of proteome
            co_annots = b.store.execute(f"""
                SELECT a2.source, a2.accession, a2.name,
                       COUNT(DISTINCT a1.protein_id) as n
                FROM annotations a1
                JOIN annotations a2 ON a1.protein_id = a2.protein_id
                JOIN proteins p ON a1.protein_id = p.protein_id
                WHERE p.bin_id = '{genome}'
                  AND a1.accession = '{accession}'
                  AND a2.accession != '{accession}'
                GROUP BY a2.source, a2.accession, a2.name
                ORDER BY n DESC LIMIT 10
            """)

            co_list = [{"source": c[0], "accession": c[1], "name": c[2], "count": c[3]} for c in co_annots]

            # Determine verdict
            if accession in SUPERFAMILY_DEFENSE:
                # Check if co-annotated with Pkinase
                kinase_co = any(c[1] in ("PF00069", "PF07714") for c in co_annots)
                if kinase_co:
                    verdict = "superfamily_fp"
                    flags.append(f"{accession} ({count} hits) = Ser/Thr kinase superfamily, not defense")
                else:
                    verdict = "ambiguous_annotation"
            elif accession in SUPERFAMILY_ACCESSIONS:
                verdict = "superfamily_fp"
                flags.append(f"{accession} {name} ({count} hits) = fold/scaffold, not specific function")
            elif count > total_proteins * 0.05:  # >5% seems high
                verdict = "ambiguous_annotation"
                flags.append(f"{accession} {name} ({count} hits, {round(count/total_proteins*100,1)}%) — unusually prevalent, check specificity")
            else:
                verdict = "validated"

            context_checks.append({
                "accession": accession,
                "name": name,
                "source": source,
                "count": count,
                "co_annotations": co_list[:5],
                "verdict": verdict,
            })

    # 5. Unannotated proteins
    unannotated = b.store.execute(f"""
        SELECT p.protein_id, p.sequence_length
        FROM proteins p
        LEFT JOIN annotations a ON p.protein_id = a.protein_id
        WHERE p.bin_id = '{genome}'
          AND a.protein_id IS NULL
        ORDER BY p.sequence_length DESC
    """)

    # 6. Giant proteins (>1000 aa)
    giants = b.store.execute(f"""
        SELECT p.protein_id, p.sequence_length,
               COUNT(DISTINCT a.accession) as n_annotations
        FROM proteins p
        LEFT JOIN annotations a ON p.protein_id = a.protein_id
        WHERE p.bin_id = '{genome}' AND p.sequence_length > 1000
        GROUP BY p.protein_id, p.sequence_length
        ORDER BY p.sequence_length DESC
    """)
    giant_unannotated = [g for g in giants if g[2] == 0]

    if giant_unannotated:
        flags.append(f"{len(giant_unannotated)} giant proteins (>1000 aa) with zero annotations")

    # 7. Defense systems (actual, after superfamily filtering)
    defense_raw = b.store.execute(f"""
        SELECT a.accession, a.name, COUNT(DISTINCT a.protein_id) as n
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = '{genome}' AND a.source = 'defensefinder'
        GROUP BY a.accession, a.name
        ORDER BY n DESC
    """)
    defense_filtered = []
    for d in defense_raw:
        if d[0] not in SUPERFAMILY_DEFENSE:
            defense_filtered.append({"accession": d[0], "name": d[1], "count": d[2]})
        else:
            # Check kinase co-annotation
            kinase_check = b.store.execute(f"""
                SELECT COUNT(DISTINCT a1.protein_id)
                FROM annotations a1
                JOIN annotations a2 ON a1.protein_id = a2.protein_id
                JOIN proteins p ON a1.protein_id = p.protein_id
                WHERE p.bin_id = '{genome}'
                  AND a1.accession = '{d[0]}'
                  AND a2.accession IN ('PF00069', 'PF07714')
            """)[0][0]
            if kinase_check < d[2] * 0.5:  # Less than half are kinases
                defense_filtered.append({"accession": d[0], "name": d[1], "count": d[2], "note": "partially validated"})

    # 8. HydDB hits
    hyddb_hits = b.store.execute(f"""
        SELECT a.accession, a.name, COUNT(DISTINCT a.protein_id) as n
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = '{genome}' AND a.source = 'hyddb'
        GROUP BY a.accession, a.name
        ORDER BY n DESC
    """)

    # 9. Sample neighborhoods
    neighborhoods_sampled = 0
    neighborhood_notes = []

    # Pick targets: largest unannotated, a giant with annotations, an interesting annotation
    nbr_targets = []

    # Largest unannotated
    if unannotated:
        nbr_targets.append(("largest_unannotated", unannotated[0][0], unannotated[0][1]))

    # Giant with most annotations
    annotated_giants = [g for g in giants if g[2] > 0]
    if annotated_giants:
        nbr_targets.append(("most_annotated_giant", annotated_giants[0][0], annotated_giants[0][1]))

    # HydDB hit if any
    if hyddb_hits:
        hyddb_pid = b.store.execute(f"""
            SELECT a.protein_id FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id = '{genome}' AND a.source = 'hyddb'
            LIMIT 1
        """)[0][0]
        nbr_targets.append(("hyddb_hit", hyddb_pid, None))

    # DefenseFinder hit (non-superfamily)
    if defense_filtered:
        def_pid = b.store.execute(f"""
            SELECT a.protein_id FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE p.bin_id = '{genome}' AND a.source = 'defensefinder'
              AND a.accession NOT IN ('Mokosh_type_I__MkoA_B', 'BREX__pglW')
            LIMIT 1
        """)
        if def_pid:
            nbr_targets.append(("defense_hit", def_pid[0][0], None))

    # Cluster of consecutive unannotated proteins
    consec_unannotated = b.store.execute(f"""
        WITH unanno AS (
            SELECT p.protein_id, p.contig_id, p.gene_index
            FROM proteins p
            LEFT JOIN annotations a ON p.protein_id = a.protein_id
            WHERE p.bin_id = '{genome}' AND a.protein_id IS NULL
        ),
        gaps AS (
            SELECT u1.protein_id, u1.contig_id, u1.gene_index,
                   u1.gene_index - ROW_NUMBER() OVER (PARTITION BY u1.contig_id ORDER BY u1.gene_index) as grp
            FROM unanno u1
        ),
        clusters AS (
            SELECT contig_id, grp, MIN(gene_index) as start_idx, MAX(gene_index) as end_idx,
                   COUNT(*) as cluster_size
            FROM gaps
            GROUP BY contig_id, grp
            HAVING COUNT(*) >= 3
        )
        SELECT c.contig_id, c.start_idx, c.end_idx, c.cluster_size
        FROM clusters c
        ORDER BY c.cluster_size DESC
        LIMIT 3
    """)

    if consec_unannotated:
        # Get a protein from the middle of the largest cluster
        cl = consec_unannotated[0]
        mid_idx = (cl[1] + cl[2]) // 2
        mid_pid = b.store.execute(f"""
            SELECT protein_id FROM proteins
            WHERE contig_id = '{cl[0]}' AND gene_index = {mid_idx}
        """)
        if mid_pid:
            nbr_targets.append(("unannotated_cluster", mid_pid[0][0], None))

    for label, pid, size_info in nbr_targets[:6]:
        try:
            nbr = b.get_neighborhood(pid, window=8, all_annotations=True)
            neighborhoods_sampled += 1
            # Summarize neighborhood
            n_genes = len(nbr) if isinstance(nbr, list) else "unknown"
            neighborhood_notes.append({
                "label": label,
                "protein_id": pid,
                "n_neighbors": n_genes,
            })
        except Exception as e:
            neighborhood_notes.append({
                "label": label,
                "protein_id": pid,
                "error": str(e)[:100],
            })

    # 10. Predicate summary (protein_predicates uses array column)
    predicates = b.store.execute(f"""
        SELECT pred, COUNT(DISTINCT pp.protein_id) as n
        FROM protein_predicates pp
        JOIN proteins p ON pp.protein_id = p.protein_id,
        UNNEST(pp.predicates) AS t(pred)
        WHERE p.bin_id = '{genome}'
          AND pred NOT IN ('small', 'medium', 'large', 'giant', 'confident_hit',
                           'multi_source', 'pfam_annotated', 'kegg_annotated', 'unannotated')
          AND pred NOT LIKE 'pfam:%' AND pred NOT LIKE 'kegg:%'
        GROUP BY pred
        ORDER BY n DESC
        LIMIT 30
    """)

    # 11. Build summary sentence
    phylum = genome.split("_")[-2] if "vamb" in genome or "metabat" in genome or "concoct" in genome or "maxbin" in genome else "Unknown"
    # Extract phylum from bin_id more carefully
    for p in ["Myxococcota", "Planctomycetota", "Proteobacteria", "Bacteroidota", "Verrucomicrobiota", "Omnitrophota"]:
        if p in genome:
            phylum = p
            break
    else:
        phylum = "Unknown"

    unanno_pct = round(len(unannotated) / total_proteins * 100, 1)
    summary = (
        f"{phylum} MAG with {total_proteins} proteins across {n_contigs} contigs; "
        f"{unanno_pct}% unannotated; {len(giants)} giant proteins (>1000 aa), "
        f"{len(giant_unannotated)} unannotated giants; "
        f"{len(defense_filtered)} defense system types detected (after superfamily filtering); "
        f"{len(hyddb_hits)} HydDB classifications."
    )

    # 12. Write inventory
    inventory = {
        "genome": genome,
        "short_name": short_name,
        "phylum": phylum,
        "n_proteins": total_proteins,
        "n_contigs": n_contigs,
        "avg_protein_length": avg_len,
        "max_protein_length": max_len,
        "annotation_coverage": {row[0]: row[1] for row in census},
        "annotation_pcts": {row[0]: round(row[1] / total_proteins * 100, 1) for row in census},
        "unannotated_count": len(unannotated),
        "unannotated_pct": unanno_pct,
        "giant_proteins": len(giants),
        "giant_unannotated": [g[0] for g in giant_unannotated],
        "giant_unannotated_count": len(giant_unannotated),
        "top_annotations": [
            {"source": t[0], "accession": t[1], "name": t[2], "count": t[3]}
            for t in top_annots[:15]
        ],
        "predicates_top15": [
            {"predicate": pr[0], "count": pr[1]} for pr in predicates[:15]
        ],
        "defense_systems": defense_filtered,
        "hyddb_hits": [{"accession": h[0], "name": h[1], "count": h[2]} for h in hyddb_hits],
        "neighborhoods_sampled": neighborhoods_sampled,
        "neighborhood_notes": neighborhood_notes,
        "unannotated_clusters": [
            {"contig": c[0], "start": c[1], "end": c[2], "size": c[3]}
            for c in (consec_unannotated or [])
        ],
        "flags": flags,
        "context_first_checks": context_checks,
        "summary": summary,
    }

    inv_path = os.path.join(out_dir, f"inventory_{short_name}.jsonl")
    with open(inv_path, "w") as f:
        f.write(json.dumps(inventory) + "\n")
    print(f"Wrote {inv_path}")

    # 13. Log notable findings
    findings = []

    # Giant unannotated proteins >2000 aa
    big_unanno = [g for g in giant_unannotated if g[1] > 2000]
    if big_unanno:
        findings.append({
            "id": f"atlas-{short_name}-001",
            "title": f"{len(big_unanno)} unannotated giant proteins >2000 aa in {phylum} MAG {short_name}",
            "category": "unknown_proteins",
            "description": (
                f"Genome {genome} contains {len(big_unanno)} proteins larger than 2000 aa "
                f"with zero HMM annotations from any source (PFAM, KEGG, DefenseFinder, HydDB). "
                f"Sizes: {', '.join(str(g[1]) + ' aa' for g in big_unanno[:5])}. "
                f"These may represent novel protein families or require E-value based re-annotation."
            ),
            "evidence": {
                "genome": genome,
                "protein_ids": [g[0] for g in big_unanno[:10]],
                "sizes": [g[1] for g in big_unanno[:10]],
            },
            "n_genomes": 1,
            "provenance": {
                "query": f"SELECT protein_id, sequence_length FROM proteins LEFT JOIN annotations ... WHERE bin_id='{genome}' AND sequence_length > 2000 AND annotations IS NULL",
                "raw_result": f"{len(big_unanno)} proteins found",
                "interpretation": "Large proteins with no HMM hits warrant structural prediction or E-value based re-search",
            },
            "figures": [],
            "related_findings": [],
            "phase": "atlas",
        })

    # Unusual defense content
    real_defense_count = sum(d["count"] for d in defense_filtered)
    if real_defense_count > 20:
        findings.append({
            "id": f"atlas-{short_name}-002",
            "title": f"Dense defense system repertoire in {phylum} MAG {short_name} ({real_defense_count} proteins)",
            "category": "defense_systems",
            "description": (
                f"After filtering superfamily false positives (Mokosh/BREX kinases), "
                f"genome {genome} retains {real_defense_count} defense-associated proteins "
                f"across {len(defense_filtered)} system types: "
                f"{', '.join(d['accession'] + ' (' + str(d['count']) + ')' for d in defense_filtered[:5])}."
            ),
            "evidence": {
                "genome": genome,
                "defense_systems": defense_filtered,
            },
            "n_genomes": 1,
            "provenance": {
                "query": "DefenseFinder annotations with kinase superfamily filtering",
                "raw_result": f"{real_defense_count} proteins in {len(defense_filtered)} systems",
                "interpretation": "Dense defense repertoire after superfamily correction",
            },
            "figures": [],
            "related_findings": [],
            "phase": "atlas",
        })

    # Very high unannotated fraction
    if unanno_pct > 50:
        findings.append({
            "id": f"atlas-{short_name}-003",
            "title": f"Majority of proteome unannotated in {phylum} MAG {short_name} ({unanno_pct}%)",
            "category": "annotation_gaps",
            "description": (
                f"Over half ({unanno_pct}%) of the {total_proteins} proteins in {genome} "
                f"have no HMM annotations from any source. This suggests either a highly divergent "
                f"organism with novel protein families, or assembly/gene-calling artifacts. "
                f"The genome has {n_contigs} contigs."
            ),
            "evidence": {
                "genome": genome,
                "unannotated_count": len(unannotated),
                "total_proteins": total_proteins,
            },
            "n_genomes": 1,
            "provenance": {
                "query": f"LEFT JOIN annotations ... WHERE bin_id='{genome}' AND annotations IS NULL",
                "raw_result": f"{len(unannotated)}/{total_proteins} unannotated",
                "interpretation": "High unannotated fraction may indicate divergent lineage or annotation gaps",
            },
            "figures": [],
            "related_findings": [],
            "phase": "atlas",
        })

    # HydDB hits
    if hyddb_hits:
        findings.append({
            "id": f"atlas-{short_name}-004",
            "title": f"Hydrogenase(s) detected in {phylum} MAG {short_name}: {', '.join(h[1] for h in hyddb_hits)}",
            "category": "energy_metabolism",
            "description": (
                f"HydDB classified {sum(h[2] for h in hyddb_hits)} proteins in {genome}: "
                f"{'; '.join(h[0] + ' ' + h[1] + ' (' + str(h[2]) + ' hits)' for h in hyddb_hits)}. "
                f"These require neighborhood validation to distinguish true hydrogenases from Complex I / NuoD false positives."
            ),
            "evidence": {
                "genome": genome,
                "hyddb_hits": [{"accession": h[0], "name": h[1], "count": h[2]} for h in hyddb_hits],
            },
            "n_genomes": 1,
            "provenance": {
                "query": f"SELECT accession, name, COUNT(DISTINCT protein_id) FROM annotations WHERE source='hyddb' AND bin_id='{genome}'",
                "raw_result": f"{len(hyddb_hits)} HydDB classifications",
                "interpretation": "HydDB hits need neighborhood curation to confirm",
            },
            "figures": [],
            "related_findings": [],
            "phase": "atlas",
        })

    if findings:
        findings_path = os.path.join(out_dir, f"findings_{short_name}.jsonl")
        with open(findings_path, "w") as f:
            for finding in findings:
                f.write(json.dumps(finding) + "\n")
        print(f"Wrote {len(findings)} findings to {findings_path}")
    else:
        print(f"No notable findings for {short_name}")

    return inventory


if __name__ == "__main__":
    genome = sys.argv[1]
    short_name = sys.argv[2]
    read_genome(genome, short_name)
