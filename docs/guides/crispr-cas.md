# CRISPR-Cas subtypes

`sharur cas-type` calls CRISPR-Cas systems by subtype with the profiles and scoring scheme of [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper) (Russel et al. 2020, *The CRISPR Journal* 3:462). It runs on the dataset's own gene calls and MinCED arrays, so every call points at stored protein and locus IDs.

## Setup

The Cas profiles, subtype scoring table, effector cutoffs and repeat model come as one Aksha database:

```bash
aksha initialize --hmms CCTyper
```

Sharur finds it through Aksha's install record. `--cctyper-db DIR` or `SHARUR_CCTYPER_DB` points at another copy. Repeat typing also needs `xgboost`.

## Running

```bash
sharur cas-type --db data/my_dataset/sharur.duckdb --dry-run   # report only
sharur cas-type --db data/my_dataset/sharur.duckdb             # write the calls
```

Each genome is searched separately on every CPU (`--workers`), so E-values match a per-assembly CRISPRCasTyper run. `--genome` limits the run to named genomes. Writing replaces earlier calls from this caller; re-seal the dataset afterwards.

## What it does

1. **Profiles.** About 700 Cas profiles run against every protein. Per protein, the best-scoring profile is kept if it passes its thresholds: effector-specific E-value and coverage cutoffs for Cas9, Cas12 and Cas13 families, and E-value < 0.01 with 30% sequence and profile coverage for the rest.
2. **Operons.** Hits with at most three other genes between them form one operon.
3. **Subtypes.** Each operon is scored against the subtype table, using each gene's best profile.
   - Operons of three or more genes take the top-scoring subtype.
   - A score of 5 or below needs a signature effector.
   - Tied subtypes give `Ambiguous`.
   - Six or more genes from two subtypes with type-unique genes give `Hybrid(...)`.
   - One- and two-gene operons need a single-effector signature gene.
   - Interference and adaptation completeness list the share of each subtype's core gene sets found.
4. **Arrays.** Each MinCED array's consensus repeat gets a subtype prediction from CRISPRCasTyper's repeat classifier. It is named when the probability is ≥ 0.75. An array is trusted when its repeats are conserved (identity > 70%) and its spacers are diverse (identity < 55%) and even in length (SEM < 3.5), or when the repeat prediction reaches 0.9.
5. **Linking.** Arrays within 10 kb of an operon join it. The joint call reconciles the operon subtype with the nearest array's repeat subtype.

## Tables

| Table | Contents |
|---|---|
| `crispr_cas_systems` | One row per operon: `status`, joint `prediction`, operon call `prediction_cas`, `best_type`, `best_score`, completeness, member proteins and profiles, linked array loci and distances |
| `crispr_array_types` | Per array locus: repeat `subtype`, `probability`, `prediction`, identity statistics, `trusted`, `near_cas` |
| `system_proteins` | Operon members with `system_source = 'cctyper'`, profile and score |

`status` separates confident calls from candidates:

| status | Meaning |
|---|---|
| `crispr_cas` | Typed operon with a linked array and an agreeing joint call |
| `cas` | Typed operon, no array within 10 kb |
| `crispr_cas_putative` | Linked array, joint call `Unknown` or `(Putative)` |
| `cas_putative` | Operon call `False` (too little Cas evidence) or `Ambiguous` |

Named subtype claims rest on `crispr_cas` and `cas` rows. Putative rows are candidates: report their genes and profiles.

## Agreement with CRISPRCasTyper

On 30 DPANN genomes (20 with MinCED arrays, 10 with Cas-domain hits only), run through CRISPRCasTyper 1.9.0 with the same proteins and gene coordinates, the port reproduces all 527 operons field for field. That covers span, prediction, best subtype, score, completeness and strand.

For the 32 arrays found by both, repeat subtypes and probabilities are identical. Trust status matches for 31: CRISPRCasTyper estimates repeat identity from a random sample of ten repeats, and the port uses an evenly spaced one.

Arrays come from MinCED, so genomes MinCED did not scan have no arrays to link. The port also skips CRISPRCasTyper's search for known repeats near operons that lack an array.
