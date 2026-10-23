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
