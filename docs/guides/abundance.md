# Abundance across samples

With read-mapping coverage, Sharur answers which genomes, and which functions, are abundant in which samples. Coverage is optional and lives in a sidecar, `abundance.duckdb`, beside the dataset database; the core database and its seal stay unchanged.

## Getting coverage in

**During ingest.** Give `sharur-ingest` a reads table and it adds Stage 08, which maps each sample to the dataset's assemblies with [CoverM](https://github.com/wwood/CoverM) using every allocated CPU:

```bash
sharur-ingest --input-dir genomes/ --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb --reads reads.tsv
```

`reads.tsv` is tab-separated with a header: `sample_id`, `read1`, and `read2` for paired reads or `interleaved` (`true`) for interleaved files. Any further columns (habitat, depth, date) are stored as sample metadata.

**From existing coverage tables.** Import CoverM `contig` output, or a long table, directly:

```bash
coverm contig --reference all_contigs.fna -1 s1_R1.fq.gz -2 s1_R2.fq.gz \
  -m mean covered_fraction count length -o s1.tsv
sharur import-coverage s1.tsv s2.tsv --db data/my_dataset/sharur.duckdb --samples samples.tsv
sharur import-coverage coverage_long.tsv --format long --db data/my_dataset/sharur.duckdb
```

The long format has columns `sample_id`, `contig_id` and any of `mean_depth`, `covered_fraction`, `read_count`, `contig_length`. Contigs must exist in the dataset; `--allow-unknown-contigs` keeps the matching rows when the mapping reference held more. Re-importing a sample replaces its rows.

## Questions you can ask

```bash
sharur abundance --db data/my_dataset/sharur.duckdb                    # genomes per sample
sharur abundance --db data/my_dataset/sharur.duckdb --predicate nife_group1
sharur abundance --db data/my_dataset/sharur.duckdb --annotation K00370
sharur coverage-outliers GENOME_ID --db data/my_dataset/sharur.duckdb
```

```python
b.abundance(samples=["s1"])
b.feature_abundance(predicate="nife_group1")
b.coverage_outliers("GENOME_ID")
```

| Quantity | Definition |
|---|---|
| Genome mean depth | Contig mean depths weighted by contig length |
| Covered fraction | Contig covered fractions weighted by contig length |
| Relative abundance | The genome's share of reads mapped to dataset contigs in that sample (depth × length when read counts are absent) |
| Feature abundance | Summed relative abundance of genomes that carry the predicate or annotation, with the top carriers |
| Coverage outlier | A contig whose depth differs from its genome's median contig depth by at least 2× (`--min-log2 1`) in most samples where the genome is covered |

Relative abundance is a share of *mapped* reads, so it describes the community represented by the dataset's genomes. Coverage outliers point at contigs worth a second look: binning errors, multicopy elements, or strain-level variation.
