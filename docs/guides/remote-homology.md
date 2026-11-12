# Remote homology

Pfam at gathering thresholds, KOfam and the system callers leave many proteins without an
annotation. Many of them still resemble a characterized family when a profile of the protein's
family is compared with profiles of known families. Sharur records that similarity in two
optional layers:

| Layer | Unit | Search | Stored in |
|---|---|---|---|
| CDD per protein | every unique protein sequence | RPS-BLAST against CDD, hits at E ≤ 0.01 | the best model of each region in `protein_resemblance` (`source = 'cdd'`); every hit in the stage's staging file |
| Family tier | one profile per protein family | MAFFT alignment of up to 250 unique members (representative first), optional HHblits enrichment against UniRef30, HHsearch against Pfam and PDB profile databases | `family_annotations` |

Both layers are **observed similarity**. A row says that a protein, or the family it belongs
to, resembles a profile; it names no system, pathway or function. The SQL surface labels these
rows `evidence = 'resembles'`. Both layers live in their own tables, so `annotations` holds the
curated per-domain hits alone: predicates, `unannotated` ("no hit in the curated sources") and
every query that counts annotation rows read the same rows with or without the layers.

## The "resembles" rule

Every family-tier hit with probability ≥ 50% is stored, so the threshold can change later. The
default rule, applied by the `family_resemblance` view and the SQL macros below:

- **resembles a characterized family**: a Pfam or PDB hit with probability ≥ 90% over ≥ 30
  aligned columns whose name and description carry no unknown-function marker;
- **resembles a family of unknown function**: the same, for hits to DUF or UPF families,
  "uncharacterized", "hypothetical" or "unknown function" proteins and structures
  (`unknown_function = true`, `hit_class = 'unknown_function'`).

UniRef30 hits (`source = 'hhblits_uniref'`, the first 20 at probability ≥ 50%) are kept as
context: they show whether homologs exist outside the dataset, and their names are automatic.

## Configuration

Database paths are site-specific and live in a JSON configuration, never in code. Every
database carries a name and a version, and every stored row carries both, so a refresh with
newer databases is a rerun:

```json
{
  "hhsuite_bin": "/path/to/hhsuite/bin",
  "mafft": "mafft",
  "rpsblast": "rpsblast",
  "databases": {
    "pfam":   {"prefix": "/path/to/pfam_hh/pfam", "name": "Pfam", "version": "39.0", "kind": "pfam"},
    "pdb":    {"prefix": "/path/to/pdb100/pdb100", "name": "PDB100", "version": "YYYY-MM-DD", "kind": "pdb",
               "titles": "/path/to/pdb_titles.tsv"},
    "uniref": {"prefix": "/path/to/UniRef30_YYYY_MM", "name": "UniRef30", "version": "YYYY_MM", "kind": "uniref"},
    "cdd":    {"prefix": "/path/to/cdd/Cdd", "name": "CDD", "version": "3.xx"}
  },
  "family_targets": ["pfam", "pdb"],
  "msa_kinds": ["family_only"],
  "search": {"hhblits_iterations": 2, "hhblits_evalue": 0.001, "max_hits": 100,
             "min_probability": 50, "uniref_max_hits": 20, "cdd_evalue": 0.01, "pdb_tool": "hhsearch"}
}
```

`config/remote_homology.example.json` holds the same template. `sharur remote-homology check
--config CONFIG` reports missing tools or databases and the tool versions.

- **`msa_kinds`.** `family_only` searches the family alignment itself (the cheap mode).
  `enriched` first runs HHblits against UniRef30 (`-n 2 -e 0.001`) and searches the enriched
  alignment; it finds more distant relatives and costs several times more. List both to keep
  both sets of hits (`family_annotations.msa_kind`).
- **`search.pdb_tool`.** `hhsearch` (default) compares the family profile with every PDB
  profile; `hhblits` runs one HHblits iteration with its column-state prefilter, which needs the
  database's `_a3m` and `_cs219` files and is far faster. Rows are labelled `hhsearch_pdb` or
  `hhblits_pdb`, and both count for the rule. Measure the switch on a pilot before adopting it:
  compare the families meeting the rule and the time per family under both settings.
- **PDB titles.** HH-suite PDB databases often name templates by chain only (`1abc_A`). A
  two-column `pdb code<TAB>title` file (`titles`) supplies entry titles, which the
  unknown-function flag reads. Build it once from a local PDB mirror (the `_struct.title` of
  each mmCIF) or from a database's header file.
- **Pfam as an HH-suite database.** HH-suite 3.1 `hhsearch` cannot read Pfam's HMMER3 models.
  Build the database from the seed alignments of the same Pfam release: one alignment per
  family (`.` replaced by `-`, upper case), `hhmake -i fam.fa -o fam.hhm -name "<ACC> <ID>; <DE>"
  -M 50` (most seed alignments mark no insert columns, so `-M a2m` fails), then `ffindex_build
  -s -f list PREFIX_hhm.ffdata PREFIX_hhm.ffindex`. That suffices for HH-suite 3.1; HH-suite
  3.2 and later (bioconda installs 3.3) need `PREFIX_a3m` and `PREFIX_cs219` as well
  (`cstranslate -I a3m` with HH-suite's `cs219.lib` and `context_data.lib`), and `sharur
  remote-homology check` names any that are missing. Hit names then come from the first record of
  each A3M, so give it the `"<ACC> <ID>; <DE>"` name. Self-test the result: a seed sequence of a
  family should rank its own family first.

## Family tier

The family tier needs a protein-family clustering ([Protein families](protein-families.md)):
a Sharur build, or an external clustering registered with `sharur families import`. Several
clusterings can sit side by side; every family-tier row carries its `clustering_run`.

**Scope.** By default, families found in three or more genomes whose members carry no Pfam
hit (`family_stats.n_genomes >= 3 AND family_stats.pfam_dark`). `--min-genomes` and
`--include-pfam-annotated` change the rule; `--family-id` (repeatable) selects named families.

```bash
sharur remote-homology export --db data/DS/sharur.duckdb \
  --families data/DS/stage09_families/protein_families.duckdb --run RUN --out-dir work/
sharur remote-homology search families --work work/ --config rh.json          # this host
sharur remote-homology load-families --work work/ --config rh.json --staging work/remote_homology.duckdb
```

`export` reads the database and the family tables only. Each selected family gets one FASTA:
the representative's sequence first, then up to 249 further unique member sequences in MD5
order of their sequence (a fixed, content-defined sample). Identical member sequences are written
once.

Per family the runner writes `out/<bucket>/<family>/`: the `.hhr` files (gzipped once
parsed), `hits.tsv` and `done.json` with each step's time. The query length of every result is
checked against the representative's length. A family with `done.json` for the current
configuration is skipped, so a rerun continues where a job stopped; a configuration change
(new database version, other MSA kinds) reruns what it affects. `load-families` refuses a
partial run unless `--allow-partial` (pending and failed families then get
`family_annotation_status` rows).

**On a Slurm cluster.** `runner.py` needs only Python 3, HH-suite, MAFFT and BLAST+, so the
searches run on any node that sees the work directory and the databases:

```bash
scripts/submit_remote_homology_shards.sh families -n 40 -m 10 -w work/ -c rh.json \
  -s "-A ACCOUNT -p PARTITION --qos QOS -c 24 --exclude NODES"
```

Each array task searches every 40th family (largest families first) with all of its CPUs:
`$SLURM_CPUS_ON_NODE // 4` workers of four threads. Resubmitting the same command continues an
interrupted run. Measure the cost first with a pilot: `-- --limit 100` runs the first 100
families of each shard (or `--only FILE` a list of family keys), and `done.json` records each
step's time. HHblits against UniRef30 and HHsearch against a full PDB profile database take most
of the time; Pfam searches are fast.

## CDD per protein

```bash
sharur remote-homology export-cdd --db data/DS/sharur.duckdb --out-dir work/
sharur remote-homology search cdd --work work/ --config rh.json        # or: submit_remote_homology_shards.sh cdd ...
sharur remote-homology compact-cdd --work work/                         # beside the raw RPS-BLAST tables
sharur remote-homology load-cdd --work work/ --config rh.json --staging work/remote_homology.duckdb
```

Each unique sequence is searched once (named by its MD5). `compact-cdd` reduces the raw
RPS-BLAST tables to two Parquet files (hits at E ≤ 0.01 per sequence, sorted by MD5, and the
CDD models hit), so on a cluster only `work/cdd/compact/` needs to travel.

CDD bundles models from several collections (NCBI-curated `cd` models, Pfam, SMART, COG, PRK,
TIGRFAMs and others), so one domain region usually hits several models. Two tiers keep this
manageable:

- **The staging file holds every hit** at E ≤ 0.01, per unique sequence (`cdd_hits`,
  `cdd_models`), with its query and model ranges.
- **Activation writes the best model of each region** to `protein_resemblance`. Per sequence, hits are
  taken by bit score (ties: lower E-value, then lower model id); a hit is kept when, against
  every hit kept before it, the two query ranges share at most `--max-overlap` (default 0.5) of
  the shorter range. This is the overlap measure `sharur architecture` uses for Pfam domains.
  The best-scoring hit always stays, so every protein with a CDD hit keeps at least one row.
  `load-cdd --max-overlap` sets the rule (`cdd_hits.best_in_region`); `activate --cdd-rows all`
  writes every hit instead.

The full table stays queryable by attaching the staging file read-only:

```sql
ATTACH 'stage10_cdd/remote_homology_cdd.duckdb' AS cdd (READ_ONLY);
SELECT p.protein_id, m.accession, m.name, h.evalue, h.score, h.start_aa, h.end_aa, h.best_in_region
FROM cdd.cdd_hits h
JOIN cdd.cdd_models m USING (pssm_id)
JOIN proteins p ON md5(rtrim(p.sequence, '*')) = h.sequence_md5
WHERE p.protein_id = 'PROTEIN_ID'
ORDER BY h.score DESC;
```

## During ingest

`sharur-ingest` runs the layers as optional stages:

```bash
sharur-ingest ... --with-families --with-cdd --with-family-homology --remote-homology-config rh.json
```

Stage 10 (`stage10_cdd/remote_homology_cdd.duckdb`) needs only Stage 07; Stage 11
(`stage11_family_homology/remote_homology_families.duckdb`) searches Stage 09's families
(`--family-homology-run RUN` picks a run other than the default family run). The dataset
database stays unchanged until activation.

## Activating

```bash
sharur remote-homology activate --db data/DS/sharur.duckdb \
  --staging work/remote_homology.duckdb \
  --families-staging data/DS/stage09_families/protein_families.duckdb --publish
```

This follows the staged-write protocol ([Staged writes](staged-writes.md)), as `sharur families
activate` does: the database is cloned, the rows are loaded into the clone, the changed tables
are diffed, and `--publish` swaps the clone in, keeps the previous file as
`sharur.duckdb.pre-remote-homology-*` and reseals the dataset. Without `--publish` the staged
copy waits for `sharur stage publish` or `sharur stage discard`.

- `--staging` is repeatable (for example Stage 10 and Stage 11 files).
- `--families-staging` (repeatable) loads clustering runs before the family rows, file by file in
  the order given; a later file's run replaces a same-named one. By default only the searched runs
  are taken (`--family-runs searched`); `--family-runs all` takes every run of every file, so one
  staged write can bring in a dataset's own family build, imported clusterings and the
  remote-homology rows together:

  ```bash
  sharur remote-homology activate --db data/DS/sharur.duckdb \
    --families-staging data/DS/stage09_families/protein_families.duckdb \
    --families-staging work/imported_clusterings.duckdb --family-runs all \
    --staging work/remote_homology.duckdb --label families-and-remote-homology --publish
  ```

  A database that already holds the searched runs needs no family file. Every searched run must
  exist with the clustering input digest the search recorded, and every searched family must
  exist in it.
- Family rows of the same clustering run and source are replaced; other runs stay.
- CDD rows replace the `protein_resemblance` rows of source `cdd`: the best model of each
  region by default, every hit with `--cdd-rows all`. Before that, the database's set of unique
  protein sequences must hash to the digest recorded when they were exported for the search; a
  database with added, removed or changed proteins needs a new CDD run. Every row carries the
  search and the database name and version; `protein_resemblance_runs` records the E-value
  cutoff, the region rule, the rows stored and the protein-level count of every hit, and the
  staging file that holds the full table.
- Activation (re)creates the views `family_resemblance` and `family_annotations_carried`
  (below).

A staging file that was loaded more than once keeps free blocks; `sharur remote-homology
compact-staging --staging FILE` rewrites it whole (every table checked by row count and content
fingerprint, then one atomic rename) before it is copied or activated.

`sharur remote-homology status --db ...` and `sharur preflight` (the `remote_homology`
capability) report the searches, how many families meet the rule, and CDD coverage: the rows
stored in `protein_resemblance` and the protein-level count of every hit in the staging file.

### Databases activated by earlier releases

Earlier releases wrote the CDD rows to `annotations` (`source = 'cdd'`) and the database name
and version to `annotation_source_versions`. `status` reports such rows
(`annotation_rows_to_move`), `preflight` marks the capability `stale`, and the SQL surface
leaves them out until they move:

```bash
sharur remote-homology migrate --db data/DS/sharur.duckdb            # staged: diff, then publish
sharur stage publish data/DS/sharur.duckdb.resemblance-table.staging
```

`migrate` follows the staged-write protocol (label `resemblance-table`; `--publish` publishes
directly). In one transaction on the clone it copies the rows to `protein_resemblance` with the
recorded database name, version and search, checks the copy against the originals by row count
and content fingerprint, moves the `annotation_source_versions` entry to
`protein_resemblance_runs`, deletes the rows from `annotations` and creates the family-tier
views. `annotations` is only deleted from, so its other rows stay exactly as they were. A
database without such rows gets the views and nothing else.

## Tables

| Table / view | One row per | Columns |
|---|---|---|
| `family_annotations` | family, source, MSA kind and hit | `family_id`, `clustering_run`, `search_run`, `source` (`hhsearch_pfam`, `hhsearch_pdb` or `hhblits_pdb`, `hhblits_uniref`), `db_name`, `db_version`, `msa_kind` (`enriched`, `family_only`), `hit_rank`, `hit_accession`, `hit_name`, `hit_description`, `probability` (%), `evalue`, `score`, `aligned_cols`, `identities`, `q_start`, `q_end`, `q_length` (representative length), `t_start`, `t_end`, `t_length`, `unknown_function`, `hit_label` |
| `family_annotation_status` | searched family | `status` (`ok`, `error`, `pending`), member, genome and sequence counts, representative, `n_hits`, `verified` (query length matched), `seconds`, `steps` |
| `family_annotation_runs` | search | `search_run`, `clustering_run`, selection rule, `msa_kinds`, `targets`, `databases` (names, versions, paths), parameters, tool versions, family counts, export and configuration digests |
| `family_resemblance` (view) | hit meeting the rule | the hit columns plus `hit_class` (`characterized`, `unknown_function`) and `evidence` = `resembles` |
| `protein_resemblance` | best CDD model per region and protein (default) | `protein_id`, `source` (`cdd`), `search_run`, `db_name`, `db_version`, `accession` (model, e.g. `cd00009`), `name`, `description`, `evalue`, `score` (bits), `start_aa`, `end_aa` |
| `protein_resemblance_runs` | source | `search_run`, `loaded_at`, `db_name`, `db_version`, `db_prefix`, `selection_rule`, `rows_stored`, `protein_rows_all` (every hit, expanded to proteins), `rows_loaded`, `proteins`, `staging_path`, `parameters` |
| `family_annotations_carried` (view) | family hit and family of another clustering run with the same representative | `family_id` and `clustering_run` of the target family, `source_run`, `source_family_id`, `source_members`, `target_members`, the hit columns, `unknown_function`, `evidence` = `resembles` |
| `cdd_hits`, `cdd_models`, `cdd_runs` (staging file) | every hit per unique sequence; model; search | `hit_id`, `sequence_md5`, `pssm_id`, E-value, bits, ranges, `best_in_region`; accession, name, description, length; database name and version, E-value cutoff, region rule and `max_overlap`, counts, export digest |

Query coordinates (`q_start`, `q_end`) refer to the family representative's sequence, which
leads every family alignment.

## Querying

Each SQL-surface macro below needs the remote-homology tables ([SQL surface](sql-surface.md)):

```sql
-- what one protein resembles: hits of its families in every clustering run, and its own protein_resemblance rows
FROM resembles_of('PROTEIN_ID');                               -- min_probability := 90, min_cols := 30
-- proteins resembling a target, in the feature_proteins shape (evidence = 'resembles')
FROM resembles_proteins('PF00004');
-- genome x target grid in the feature_presence shape
FROM resembles_presence('PF00004', min_completeness := 90);
-- CDD rows through the general feature specs: name the source
FROM feature_proteins('cdd:cd00009');
```

An unqualified spec (`PF00004`, `AAA`) matches `annotations` only; CDD rows appear when the
spec names `cdd:`. `protein_summary`, `neighbors` and `slot` list CDD hits in their own
`resembles` column, apart from `other` and `n_annotations`. `protein_resemblance_rows()` is
the table itself, and an empty table of the same shape on a database without the layer, so
these macros work on every database.

### Hits under another clustering run

Family-tier hits belong to the clustering run that was searched. A later run that only adds
members to existing families (for example one that moves fragments out of very small families
into larger ones) keeps each searched family's representative, so its hits carry over through the
shared representative:

```sql
-- hits of every searched run, on the families of TARGET_RUN that share the searched family's representative
FROM family_annotations_via('TARGET_RUN');
-- the stored view: the same rows for every target run, readable by any DuckDB client
FROM family_annotations_carried WHERE clustering_run = 'TARGET_RUN';
```

`source_members` and `target_members` show how much a family grew. Check the premise before
relying on the join: every member of a searched family should sit in its target counterpart.

```sql
WITH src AS (
  SELECT f.protein_id, f.family_id, s.rep_protein_id
  FROM protein_families f JOIN family_stats s USING (clustering_run, family_id)
  WHERE f.clustering_run = 'SOURCE_RUN'
    AND f.family_id IN (SELECT family_id FROM family_annotation_status WHERE status = 'ok'))
SELECT COUNT(DISTINCT src.family_id) FILTER (WHERE t.family_id IS NULL) AS no_counterpart,
       COUNT(DISTINCT src.family_id) FILTER (WHERE m.family_id <> t.family_id) AS members_elsewhere
FROM src
LEFT JOIN family_stats t ON t.rep_protein_id = src.rep_protein_id AND t.clustering_run = 'TARGET_RUN'
JOIN protein_families m ON m.protein_id = src.protein_id AND m.clustering_run = 'TARGET_RUN';
```

Families of the target run that enter a scope only through added members have no searched
counterpart; search them under the target run when they matter.

Families of unknown function stay distinguishable:

```sql
SELECT clustering_run, hit_class, COUNT(DISTINCT family_id) AS families
FROM family_resemblance GROUP BY ALL ORDER BY ALL;
```

## Reading the rows

A resemblance is evidence of homology to a family, at the strength the probability and the
aligned length give. It supports wording such as "the family resembles Pfam family X (HHsearch
probability P% over N aligned columns)" and a hypothesis tagged UNVERIFIED. A named function,
system or pathway still needs a curated caller's call; domains compatible with a name remain a
hypothesis. A family profile summarizes its members: in a permissive clustering a few members
can carry an extra domain, so a hit covers the region given by `q_start`–`q_end` of the
representative, and a family whose members differ in architecture can mix them.
