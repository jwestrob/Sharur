# Expression from metatranscriptomes

With metatranscriptome reads mapped to a dataset's contigs, Sharur answers which genomes are transcriptionally active in which samples, which genes a genome transcribes most, and how much of a sample's transcription goes to a feature. Counts are optional and live in a sidecar, `expression.duckdb`, beside the dataset database; the core database and its seal stay unchanged.

## Getting counts in

Map each RNA library to all dataset contigs at once (competitive mapping), keep primary alignments that pass your identity and aligned-fraction filters, and count fragments per CDS with featureCounts. Build the annotation from the dataset itself so gene IDs are protein IDs:

```python
import duckdb
con = duckdb.connect("data/my_dataset/sharur.duckdb", read_only=True)
con.execute("""COPY (SELECT protein_id AS GeneID, contig_id AS Chr, start AS Start, end_coord AS "End",
                            strand AS Strand FROM proteins ORDER BY contig_id, start)
               TO 'genes.saf' (DELIMITER '\t', HEADER)""")
```

```bash
featureCounts -p --countReadPairs -s 1 -F SAF -a genes.saf --largestOverlap -o lib.s1.txt lib.bam
featureCounts -p --countReadPairs -s 2 -F SAF -a genes.saf --largestOverlap -o lib.s2.txt lib.bam
```

Import either a long table or featureCounts tables:

```bash
# long: sample_id, protein_id, count_s1, count_s2 (or sense[, antisense]); TSV, CSV or Parquet
sharur import-expression gene_counts.parquet --db data/my_dataset/sharur.duckdb --samples rna_samples.tsv
# featureCounts tables of one orientation
sharur import-expression lib.s2.txt --format featurecounts --orientation sense --db data/my_dataset/sharur.duckdb
```

With `count_s1`/`count_s2`, the orientation that holds more fragments in a library is taken as sense, unless the samples table gives `sense_orientation` (`s1` or `s2`). The sample record keeps the strand specificity (sense / (sense + antisense)); a value near 0.5 means an unstranded library. Proteins must exist in the dataset; `--allow-unknown-proteins` keeps the matching rows. Re-importing a sample replaces its rows.

## Questions you can ask

```bash
sharur expression --db data/my_dataset/sharur.duckdb                    # genome transcript share per sample
sharur expression --db data/my_dataset/sharur.duckdb --bin GENOME_ID    # a genome's most expressed genes
sharur expression --db data/my_dataset/sharur.duckdb --annotation K00399
sharur expression --db data/my_dataset/sharur.duckdb --predicate nife_group1
```

```python
b.expression(samples=["lib1"])
b.genome_genes("GENOME_ID", top=20)
b.feature_expression(annotation="K00399")
```

## Reading the numbers

- **Transcript share** is a genome's share of the sample's sense fragments on dataset genes. It mixes abundance and per-cell activity; compare it with the genome's DNA share from `sharur abundance` for an activity ratio.
- **Genome TPM and percentile** normalize within one genome (fragments per kb of CDS over all of its CDS, genes without reads included), so they rank a genome's own genes and are comparable across samples with different depths. A percentile near 1.0 means one of the genome's most transcribed genes.
- **Near-identical genomes** split reads from shared sequence, so their counts are lower bounds; aggregate them at the population level before comparing.
- rRNA reads fall outside CDS and are left out by construction; the per-library fraction of fragments assigned to CDS reflects rRNA depletion and the share of the community in the dataset.
