# Maintenance after the 2026-10-23 rule fixes

Run from the main checkout after merging. Serialize every database write in a maintenance
window, then reseal.

## Hydrogenase neighborhood rule

The neighborhood verdicts are computed at read time (`/hydrogenases` in the browser, Step 3 of
`.claude/skills/hydrogenase.md`). Stored data stays valid as it is:

- `hydrogenase_classifications` (`ko_support` from the protein's own KOfam hits,
  `curation_status` from its own Pfam domains), `scripts/classify_hydrogenases.py`, the
  `hyddb_subgroup` labels and the V2 hydrogenase atoms all read the protein itself.
- Removing K07464 (Cas4, EC 3.1.12.1, which the old EC pattern read as EC 1.12) from
  `kegg_hyddb_snapshot.tsv` leaves the KEGG map unchanged: the KO captured no HydDB reference.

Follow-through:

1. Restart any running `sharur browse` session; it caches hydrogenase calls per process.
2. Re-measure any finding or report that cites neighborhood-supported Group 4 counts.
   `IMPACT_hydrogenase.tsv` gives the SR-VP changes; for another dataset:

   ```bash
   python scripts/compare_hydrogenase_neighborhood.py --old-rev 38b0093 --db data/DATASET/sharur.duckdb
   ```

## cas_nuclease and phage_terminase (GO:0004519 anchors `endonuclease` alone)

Scope, measured read-only on 2026-10-23:

- KEGG map: rejudging a scratch copy of the local build changes one KO. K05986 drops
  `cas_nuclease` and `phage_terminase` and keeps `endonuclease`, `enzyme`, `hydrolase` and
  `nuclease` (map sha256 05c167889aaeb14d -> 98e78993f0b4a5bd). K05986 was the only KO whose
  Swiss-Prot consensus carried either predicate.
- Pfam map: unchanged. No Pfam pair for either predicate rests on GO:0004519; the three
  `phage_terminase` pairs (PF03354, PF04466, PF17289) rest on the family name.
- Datasets: `srvp_bacteria_pb` holds 109 K05986 proteins in 98 genomes, each carrying
  `cas_domain`, `cas_nuclease`, `crispr_associated`, `phage_terminase` and `phage_related`;
  none sits in `crispr_cas_systems`, `defense_systems` or `loci`. `srvp_archaea`, `dpann`,
  `dpann_reingest_2026-10-14` and `omni_production` carry no K05986 annotation.

Commands:

1. Rejudge the local KEGG map from its evidence snapshot (writes `data/reference/kegg/`):

   ```bash
   sharur setup-kegg --rejudge
   python -m pytest tests/test_kegg_map_integrity.py tests/test_kegg_rules.py -q
   ```

   Expect `map sha256 05c167889aaeb14d -> 98e78993f0b4a5bd`. `test_kegg_map_integrity.py`
   fails on K05986 until this step runs.

2. Refresh predicates for the 109 affected proteins:

   ```bash
   python - <<'EOF'
   from sharur.operators import Sharur
   b = Sharur("data/srvp_bacteria_pb/sharur.duckdb")
   ids = [r[0] for r in b.store.execute(
       "SELECT DISTINCT protein_id FROM annotations WHERE source = 'kegg' AND accession = 'K05986'")]
   assert len(ids) == 109, len(ids)
   b.generate_v2(protein_ids=ids, update_legacy_predicates=True)
   b.store.close()
   EOF
   ```

   Verify (expect 0; a scratch copy of the 109 proteins gave 0 and kept `endonuclease` and
   `nuclease` on all 109):

   ```sql
   SELECT COUNT(DISTINCT pp.protein_id) FROM protein_predicates pp
   JOIN annotations a USING (protein_id)
   WHERE a.source = 'kegg' AND a.accession = 'K05986'
     AND (list_contains(pp.predicates, 'cas_domain') OR list_contains(pp.predicates, 'cas_nuclease')
          OR list_contains(pp.predicates, 'phage_terminase'));
   ```

3. Reseal:

   ```bash
   sharur seal --db data/srvp_bacteria_pb/sharur.duckdb
   sharur verify-seal data/srvp_bacteria_pb/dataset.seal.json
   ```

Predicate stamps: both commits edit KEGG rule files (`kegg_hyddb_snapshot.tsv`,
`kegg_swissprot_consensus.tsv`) and step 1 changes the KEGG map, so `sharur describe` reports
every dataset's predicates as stale against the installed maps. The content difference is the
K05986 change above. A full regeneration clears the stamp for one dataset:

```bash
sharur compute-predicates --db data/DATASET/sharur.duckdb --chunk-size 100000 --workers N
```
