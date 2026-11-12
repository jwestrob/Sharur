# Protein families

Pfam and KOfam describe the proteins they recognize. Protein families give every protein a
handle, recognized or dark: Sharur clusters all proteins of a dataset with
[MMseqs2](https://github.com/soedinglab/MMseqs2) into families and, inside each family,
tighter subfamilies. The layer is optional. It is computed into a staging DuckDB and copied into
the dataset database with one command.

| Level | Criterion (every pair a member and its representative) | Default |
|---|---|---|
| Family | ≥ 30% identity over the alignment, alignment covers ≥ 80% of both sequences | `--family-min-seq-id 0.3` |
| Subfamily | the same at a tighter identity | `--subfamily-min-seq-id 0.5` |

Both levels run `mmseqs cluster` with `-c 0.8 --cov-mode 0 -s 6 --max-seqs 50 --cluster-mode 0
--cluster-reassign 1`. Clustering is hierarchical: subfamilies are clustered from all proteins, and
families from the subfamily representatives, so every subfamily lies inside exactly one family.

## Building families

**During ingest.** `--with-families` adds Stage 09 after the database build:

```bash
sharur-ingest --input-dir genomes/ --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb --with-families
```

**For an existing dataset.**

```bash
sharur families run --db data/my_dataset/sharur.duckdb --out-dir data/my_dataset/stage09_families
```

Either way the dataset database is only read. `stage09_families/` receives the staging database
`protein_families.duckdb`, the clustering tables (`clustering/level0.tsv.gz` for subfamilies,
`level1.tsv.gz` for families of subfamily representatives, `clustering.json` with every MMseqs2
command, version and timing) and `export.json` (protein count and sequence digest). MMseqs2 and
DuckDB use every allocated CPU (`$SLURM_CPUS_ON_NODE` under Slurm).

**On a cluster node.** `sharur/protein_families/mmseqs.py` needs only Python and MMseqs2, so the
clustering runs on any host that holds the exported FASTA:

```bash
sharur families export --db data/my_dataset/sharur.duckdb --out-dir work/       # proteins.faa, export.json
# on the cluster, in a job with N CPUs:
python3 mmseqs.py cluster --fasta proteins.faa --out-dir clustering --levels 0.5,0.3 --threads "$SLURM_CPUS_ON_NODE"
# back beside the database (copy export.json into clustering/):
sharur families build --db data/my_dataset/sharur.duckdb --clustering clustering \
  --staging work/protein_families.duckdb
```

On Slurm give the job as many CPUs as the cluster allows per job, and exclude nodes whose CPUs
miss the instruction set the MMseqs2 binary was built for (`SHARUR_SLURM_EXCLUDE` for the ingest
bundle). Clustering 3.6 million proteins at both levels took about 30 minutes on 64 CPUs; the
four-genome `dummy_dataset` (10,102 proteins) takes about 30 seconds on 4 CPUs.

`build` checks that the database's protein sequences still match the exported digest and refuses
otherwise.

## Activating

```bash
sharur families activate --db data/my_dataset/sharur.duckdb \
  --staging data/my_dataset/stage09_families/protein_families.duckdb --publish
```

This follows the staged-write protocol ([Staged writes](staged-writes.md)): the database is
cloned, the runs of the staging file are copied into the clone, the family tables are diffed, and
`--publish` swaps the clone in, keeps the previous file as `sharur.duckdb.pre-protein-families-*`
and reseals the dataset. Without `--publish` the staged copy stays for review
(`sharur stage publish` or `sharur stage discard`). Before copying, every staged protein must
exist in `proteins` with the same `bin_id`, and the sequences of each run's proteins must match
the digest recorded at clustering. Runs already in the database under other names are kept; a run
of the same name is replaced. `--run NAME` (repeatable) copies selected runs.

`sharur families status --db ...` lists the runs and the proteins that lack a family;
`sharur preflight` reports the layer as the `protein_families` capability.

## Tables

Several clusterings of the same proteins live side by side, each a **clustering run**. A Sharur
build records two runs: `<build_id>-family` and `<build_id>-subfamily`, the second naming the first
as its `parent_run`. The optional truncated-gene rescue adds two derived runs on top,
`<build_id>-family-rescued` and `<build_id>-subfamily-rescued` (see below).

| Table / view | One row per | Columns |
|---|---|---|
| `protein_families` | protein and run | `protein_id`, `family_id`, `is_representative`, `clustering_run`, `bin_id`, `assignment` (`cluster`, `import`, for later additions `search`, `new_subfamily`, `new_family`, and `rescue` in rescued runs), `batch_id` |
| `family_stats` | family and run | `family_id`, `clustering_run`, `n_members`, `n_genomes`, `rep_protein_id`, `rep_length`, `parent_family_id` (subfamily runs), `n_subfamilies` (family runs of a build), `n_orders` (GTDB orders in `bins.taxonomy`), `median_length`, `pfam_fraction` (members with any Pfam hit), `dominant_pfam_architecture` and `dominant_pfam_architecture_fraction` (share of all members), `ko_fraction`, `dominant_ko` and `dominant_ko_fraction`, `pfam_dark` (no member has a Pfam hit), `n_rescued` (members placed by the rescue) |
| `protein_family_runs` | run | `clustering_run`, `level`, `parent_run`, `build_id`, `origin` (`sharur`, `import` or `rescue`), `is_default`, `tool`, `tool_version`, `parameters`, `input_proteins`, `input_digest`, `families`, `source_db`, `details` |
| `protein_family_assignments` | assignment batch | `batch_id`, `build_id`, `family_run`, `subfamily_run`, parameters, protein count and digest |
| `protein_family_rescues` | rescued protein | `rescue_id`, `family_run`, `subfamily_run`, `base_family_run`, `protein_id`, `bin_id`, `partial`, `original_family_id`, `original_subfamily_id`, `family_id`, `subfamily_id`, `joined_subfamily`, `hit_protein_id`, `identity`, `query_coverage`, `target_coverage`, `bits`, `evalue`, `query_length`, `target_length` |
| `genome_protein_families` (view) | genome and family of the default family run | `genome_id`, `feature_id`, `n_proteins`, `n_rescued` |
| `genome_protein_subfamilies` (view) | genome and subfamily of the default subfamily run | `genome_id`, `feature_id`, `n_proteins`, `n_rescued` |
| `genome_protein_families_clustered`, `genome_protein_subfamilies_clustered` (views) | the same over the clustering partition (rescued members left out) | `genome_id`, `feature_id`, `n_proteins` |
| `genome_family_presence` (view) | run, genome and family | `clustering_run`, `genome_id`, `feature_id`, `n_proteins`, `n_rescued` |
| `protein_family_levels` (view) | protein in the default build | `protein_id`, `bin_id`, `family_id`, `subfamily_id`, `family_run`, `subfamily_run`, `rescued` |

Dominant architectures follow the rules of `sharur architecture` (Pfam hits, overlaps resolved,
N- to C-terminal order). KOs come from the `kofam` and `kegg` sources; a protein carrying two KOs
counts toward both. The tables carry their primary keys only: DuckDB scans `protein_families`
for one protein or family in a fraction of a second, and secondary ART indexes cost about 110 bytes
per row each on disk (files from before 2026-11-12 lose theirs at the next activation or rescue).

### Identifiers

- **Build:** `pfr-` plus ten hex digits of a SHA-256 over both levels' criteria, the MMseqs2
  version, the input sequence digest and the resulting partition. The same proteins clustered with
  the same parameters into the same partition give the same build id.
- **Family / subfamily:** `fam:<tag>:<representative protein id>` and
  `sub:<tag>:<representative protein id>`, `<tag>` being the build id's ten hex digits. The part
  after the second colon is the representative's protein id.
- **Imported run:** the name given at import; its family ids are `fam:<tag>:<representative>`
  with `<tag>` hashed from the name, parameters, sequence digest and partition.
- **Sequence digest:** SHA-256 over `protein_id<TAB>sequence<LF>` lines in `protein_id` order, with
  trailing `*` removed.

## Querying

Per-genome presence has the shape of any per-genome feature table, so family presence can sit
beside Pfam or KO presence in a screen:

```sql
-- families, Pfam families and KOs as (genome_id, feature_id)
SELECT genome_id, feature_id FROM genome_protein_families
UNION ALL
SELECT DISTINCT p.bin_id, a.accession FROM annotations a JOIN proteins p USING (protein_id)
WHERE a.source = 'pfam'
UNION ALL
SELECT DISTINCT p.bin_id, a.accession FROM annotations a JOIN proteins p USING (protein_id)
WHERE a.source IN ('kofam', 'kegg');

-- recurrent Pfam-dark families
SELECT family_id, n_members, n_genomes, n_orders, rep_protein_id, median_length
FROM family_stats JOIN protein_family_runs r USING (clustering_run)
WHERE r.level = 'family' AND r.is_default AND pfam_dark AND n_genomes >= 5
ORDER BY n_genomes DESC;

-- the family and subfamily of one protein, and its relatives in other genomes
SELECT * FROM protein_family_levels WHERE protein_id = 'PROTEIN_ID';
SELECT l.protein_id, l.bin_id FROM protein_family_levels l
WHERE l.family_id = (SELECT family_id FROM protein_family_levels WHERE protein_id = 'PROTEIN_ID');
```

The [SQL surface](sql-surface.md) carries three family macros on every Sharur connection:
`family_proteins(families)` lists members in the `feature_proteins` columns,
`family_presence(families)` gives the genome × family grid of `feature_presence` (completeness,
contamination and GTDB ranks included), and `families_of(ids)` gives each protein's family in every
clustering run with the family's statistics:

```sql
FROM families_of('PROTEIN_ID');
FROM family_presence('fam:0123456789:PROTEIN_ID', min_completeness := 90);
```

A family is a sequence-similarity cluster. `dominant_pfam_architecture` and `dominant_ko` report
what annotations its members carry; they describe the members. A function name for a family needs
the same curated evidence as for a single protein.

## Several runs and imported clusterings

A clustering made elsewhere (for example a permissive one used as tokens by another tool) can be
registered as a run directly from its two-column `mmseqs createtsv` output:

```bash
sharur families import --db data/my_dataset/sharur.duckdb --tsv clu.tsv --name permissive-e3c50s75 \
  --level permissive --tool-version 18-8cc5c \
  --parameters "mmseqs cluster -s 7.5 -e 1e-3 --min-seq-id 0 -c 0.5 --cov-mode 0 --cluster-mode 0" \
  --staging data/my_dataset/stage09_families/protein_families.duckdb
sharur families activate --db data/my_dataset/sharur.duckdb \
  --staging data/my_dataset/stage09_families/protein_families.duckdb --run permissive-e3c50s75 --publish
```

Every member must be a protein of the database and appear once, and each representative must be a
member of its own cluster. The import records the sequence digest of the clustered proteins, so
activation checks the run like any other. `--default` makes an imported run the default of its
`--level`, which moves the per-level views to it; otherwise the views keep following the Sharur
build and `genome_family_presence` lists every run.

## Adding genomes

Proteins that arrive after clustering get their rows by assignment against the build's representatives:

```bash
sharur families assign --db data/my_dataset/sharur.duckdb --publish
```

1. New proteins are searched against the subfamily representatives with the subfamily criteria;
   a hit joins that subfamily (`search`).
2. The rest are clustered with the subfamily criteria; each new subfamily's representative is
   searched against the family representatives with the family criteria, and a hit joins that
   family (`new_subfamily`).
3. The remaining new subfamily representatives are clustered with the family criteria into new
   families (`new_family`).

Existing families keep their members and ids; new ones are named after their own representatives
in the build's namespace. Each batch is recorded in `protein_family_assignments`. On
`dummy_dataset`, clustering three genomes and assigning the fourth (826 proteins: 18 `search`, 137
`new_subfamily`, 671 `new_family`) agreed with a full re-cluster of all four genomes at an adjusted
Rand index of 0.986 for families and 0.994 for subfamilies; 7,618 of 7,674 families were identical.

**When to re-cluster.** Run `sharur families run` and activate the new build when

- more than 10% of a build's proteins were placed by assignment (`status` and `preflight` flag it:
  assignment never merges two existing families that new genomes connect, and representatives stay
  those of the original proteins);
- proteins were removed or their sequences changed (activation and `assign` refuse: the digests
  or representatives no longer match);
- the parameters or the MMseqs2 version change, or ids must be comparable with a fresh clustering
  of the whole dataset.

## Choosing the thresholds

The family level uses 30% identity over 80% of both sequences, a common setting for family-level
MMseqs2 clustering. The subfamily level should split what a family merges (paralogs and
functionally distinct groups) while keeping each orthologous group whole. Its default was chosen
by clustering a dataset of 3.6 million archaeal proteins at candidate identities 0.4 to 0.8 and
scoring each partition two ways, with the selection rule fixed before the results were read:

1. **KO agreement** on proteins that carry exactly one KO (and KOs carried by at least two such
   proteins): pairwise precision (pairs in one subfamily that share their KO), pairwise recall
   (pairs sharing a KO that share a subfamily), V-measure, the share of multi-member subfamilies
   holding more than one KO, and the share of KOs split across subfamilies.
2. **Curated clades.** Membership of phylogenetically defined clades of one enzyme family (clades
   with strong support, against their paralogs and nearest outgroup branches): how many subfamilies
   each clade spans, the share of each clade in its largest subfamily, and outsiders that share a
   subfamily with clade members.

Rule: the lowest identity whose subfamilies are KO-pure (pairwise precision ≥ 0.995 and at most 1%
of multi-member subfamilies mixing KOs) and keep the curated clades apart from each other and from
their outgroup. Lower identities keep more of each orthologous group together, so the lowest
identity that meets the purity and separation conditions gives the most complete subfamilies. 0.5
was the lowest candidate that met the rule. Repeat the exercise with `--subfamily-min-seq-id` for a
dataset whose diversity differs much from that of the test set.

Hierarchical families depend a little on the subfamily threshold, because families are clustered
from the subfamily representatives. Two independent clusterings at 0.3 and 0.5 nest poorly (a
direct 0.5 cluster often spans several direct 0.3 families), so the hierarchy is what keeps every
subfamily inside one family.

## Truncated gene calls

Coverage of both sequences keeps fragments out of families: a gene call truncated by a contig end
covers too little of its full-length relatives and often ends up a singleton. On `dummy_dataset`,
84% of truncated calls are family singletons against 61% of complete calls (four genomes from four
phyla, so most families there are singletons anyway). Measure it for a dataset with Prodigal's
`partial` flags:

```sql
SELECT CASE WHEN p.partial = '00' THEN 'complete' ELSE 'truncated' END AS call,
       COUNT(*) AS proteins, AVG((s.n_members = 1)::INT) AS singleton_share
FROM protein_family_levels l JOIN proteins p USING (protein_id)
JOIN family_stats s ON s.family_id = l.family_id AND s.clustering_run = l.family_run
WHERE p.partial IS NOT NULL GROUP BY 1;
```

### The rescue layer (optional, off by default)

The rescue places such fragments into existing families as a separate, labeled layer. Gene calls
are never modified, and the clustering partition and build id stay as they were.

- **Candidates:** proteins whose Prodigal partial flag marks a contig-edge truncation (`10`, `01`,
  `11`) and whose family in the build has one member.
- **Targets:** every protein of a family with at least `--min-family-size` members (default 2).
- **Search:** MMseqs2 with the build's alignment settings; the best hit by bit score is accepted when
  the alignment covers at least `--min-query-coverage` of the candidate (default 0.8; coverage of
  the candidate only, `--cov-mode 2`) at identity of at least `--min-seq-id` (default: the family
  floor, 0.3).
- **Result:** the candidate joins the hit's family, and the hit's subfamily when identity also
  reaches the subfamily floor. Derived runs `<build_id>-family-rescued` and
  `<build_id>-subfamily-rescued` hold the build's membership with these proteins moved
  (`assignment = 'rescue'`); `protein_family_rescues` records each one with its original family, the
  hit, identity and coverages. `family_stats` has rows for both the base and the rescued runs.

```bash
sharur-ingest ... --with-families --family-rescue            # Stage 09 adds the layer
sharur families run ... --rescue                             # the same for an existing dataset
sharur families rescue --db data/my_dataset/sharur.duckdb \
  --staging data/my_dataset/stage09_families/protein_families.duckdb             # add to a staging file
sharur families rescue --db data/my_dataset/sharur.duckdb --include --publish   # add to a dataset, as default
sharur families rescue ... --drop                                               # remove the layer
```

On a staging file the layer is written into a copy-on-write clone, the build's base runs are
checked unchanged by fingerprint, the clone replaces the file with one rename, and the result is
rewritten without free blocks. For a large dataset the search can run elsewhere:
`--export-inputs DIR` writes `candidates.faa` and `targets.faa`,
`python3 mmseqs.py search --query candidates.faa --target targets.faa --cov-mode 2 --top 5 ...`
runs on a cluster node, and `--hits hits.tsv --tool-version VERSION` applies the result. `sharur
families assign` recomputes an existing layer with its stored parameters, so newly assigned
fragments are covered too.

**Including or excluding rescued members.** By default an added layer is excluded: it is stored,
and the defaults stay on the base runs. With `--include` (`--rescue-include` for Stage 09 and
`sharur families run`) the rescued runs become the defaults: `genome_protein_families`,
`genome_protein_subfamilies` and `protein_family_levels` then include rescued members, flagged by
`n_rescued` and `rescued`. `genome_protein_families_clustered` and
`genome_protein_subfamilies_clustered` always give the clustering partition, and
`genome_family_presence` and `family_stats` carry the base and the rescued runs side by side. The
SQL macros read the rescued runs wherever a layer exists unless given `include_rescued := false`,
and mark rescued members with `evidence = 'rescued'`.

**Why it is off by default.** The rescue was validated against acceptance criteria written before
the runs, on a dataset of 3.6 million archaeal proteins:

1. **Simulated truncation:** complete proteins from multi-member families, stratified by family
   size and length, cut as a contig end would cut them (an N- or C-terminal piece, kept fraction
   uniform between 0.15 and 0.95), with the intact original removed from the targets. Criterion:
   at most 1% of fragments join a family other than their own.
2. **Negative control:** the same for complete proteins that are family singletons. Criterion: at
   most 2% join any family.
3. **Real truncated genes:** among rescued genes carrying a Pfam hit, every domain of the gene in the
   dominant Pfam architecture of the family it joined. Criterion: at least 95%.

The rule met none of the three. Most simulated fragments returned to their own family, but joins to
another family ran at several percent in every kept-fraction and length bin, and fragments of
genuine singletons joined a family often. Coverage of the candidate alone cannot separate a gene cut
by a contig end from a protein that shares a domain or a stretch with another family. The layer
therefore stays opt-in: use it where more complete families matter more than a few percent of
misplaced fragments, and read `rescue` rows as observed homology to a family, a weaker statement
than membership from clustering.
