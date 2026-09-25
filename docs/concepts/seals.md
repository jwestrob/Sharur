# Seals and provenance

## Dataset seals

A seal is a receipt for a dataset. After ingest finishes, run:

```bash
sharur seal --db data/my_dataset/sharur.duckdb
```

This writes `data/my_dataset/dataset.seal.json`, which records:

- the input genome files;
- every database table, its columns and its row count;
- the ingest stages that ran and their manifests;
- the annotation sources and curated callers present;
- fingerprints of key files (embeddings, indexes, findings);
- the Sharur version, git commit and package versions that produced it.

All of this is summarized in one fingerprint, the **dataset ID**.

Later, check that the dataset still matches its receipt:

```bash
sharur verify-seal data/my_dataset/dataset.seal.json
```

`verify-seal` rebuilds the receipt from what is on disk and compares the two. A match means you are analyzing exactly the dataset the seal describes. A mismatch names what changed: a table's row count, a replaced input file, a new annotation source.

### When to use it

- **Resuming work** after a break, to confirm nothing changed underneath you.
- **Sharing or copying** a dataset to another machine or a collaborator.
- **Reporting results**: cite the dataset ID so each number points at one exact dataset.
- **Multi-agent campaigns**: `sharur-query` serves only sealed databases, so every agent reads the same data.

By default a seal fingerprints large files by sampling their contents, which takes seconds even for very large databases. `sharur seal --full` hashes every byte of every file, which is slower and stronger. After changing a dataset on purpose (adding annotations, recomputing predicates), write a new seal with `--force`.

## Predicate provenance

Each time predicates are computed, Sharur records which predicate maps, vocabulary and configuration produced them, with their hashes and the git commit. `sharur describe` and `sharur preflight` report whether a database's predicates match the maps currently installed, so an upgrade that changes the maps is visible before you compare results. See [How predicates are built](../predicate_construction.md#provenance-and-verification).

## Capability preflight

```bash
sharur preflight --db data/my_dataset/sharur.duckdb --format json
```

`preflight` reports what a dataset and the current machine can do: which annotation sources, callers, embeddings and external tools are available. Agents read it before planning an analysis.
