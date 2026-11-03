# Staged writes to a canonical database

Every write to a dataset's `sharur.duckdb` (a loader, a caller rerun, a predicate regeneration)
can run against a staged clone first. `sharur stage run` shows what the write changes and which
recorded findings it moves before the canonical file changes.

```bash
sharur stage run --db data/DATASET/sharur.duckdb --label cazy-2026-11 \
    --findings data/DATASET --findings data/JOINT \
    --reports 'data/JOINT/reports/*.md' -- \
    python scripts/classify_cazymes.py load --db {staged} ...
```

1. **Stage.** The database is checkpointed and cloned beside itself. On APFS (macOS) and
   reflink filesystems the clone takes milliseconds and no space until pages diverge.
   `<staged>.stage.json` records the source's file state.
2. **Run.** The command after `--` runs with `{staged}` replaced by the staged path, also
   exported as `SHARUR_STAGED_DB`.
3. **Diff.** `sharur diff-db` compares every table: schema, row counts, an order-free content
   hash, rows present on only one side, and added/removed/changed keys for keyed tables. For a
   table with a `source` column it names the sources of the differing rows.
4. **Reverify.** With `--findings`, `sharur reverify --compare` replays the findings that read
   what changed (changed table names, or the changed `source` values for `annotations`) against
   the live and staged databases, and lists every recorded number that moves, with the report
   lines that cite the finding.
5. **Publish.** `sharur stage publish STAGED` replaces the canonical file when it is unchanged
   since staging. It holds DuckDB's exclusive lock through the rename, keeps the prior file as
   `sharur.duckdb.pre-LABEL-TIMESTAMP` (a hard link), and re-seals with the existing seal
   policy. `--publish` on `stage run` publishes directly when no recorded verification result
   changes. `sharur stage discard STAGED` drops a staged copy.

Reports land in `--report-dir` (default: beside the staged file): `LABEL.diff.md/json` and
`LABEL.reverify.md/tsv`.

`sharur diff-db BEFORE AFTER` also works on its own, for instance between a pre-publish backup
and the live database.
