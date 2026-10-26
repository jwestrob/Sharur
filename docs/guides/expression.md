# Expression from metatranscriptomes

With metatranscriptome reads mapped to a dataset's contigs, Sharur answers which genomes are transcriptionally active in which samples, which genes a genome transcribes most, and how much of a sample's transcription goes to a feature. Counts are optional and live in a sidecar, `expression.duckdb`, beside the dataset database; the core database and its seal stay unchanged.

## Getting counts in

Map each RNA library to all dataset contigs at once (competitive mapping), keep primary alignments that pass your identity and aligned-fraction filters, remove rRNA fragments, and count fragments per CDS with featureCounts. Build the annotation from the dataset itself so gene IDs are protein IDs:

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

**Remove rRNA and count structural RNAs alongside the CDS.** Run SortMeRNA (for example 7.0 with `smr_v4.3_default_db`, `--no-best --num_alignments 1`) on the reads in the filtered BAM, and drop every fragment with an rRNA hit on either mate. Each pair aligns independently, so filtering the existing BAM by read name matches a remap of the cleaned reads. Then append the structural-RNA loci to the SAF under one meta-feature ID per family (for example `RNA__tmRNA`), so `--largestOverlap` assigns reads that straddle a CDS edge to the RNA. Genes that sit beside an abundant RNA need this. For example, `smpB` often overlaps tmRNA by a few tens of bp, below any useful overlap-flag threshold, and otherwise collects the tmRNA reads. Drop the `RNA__` rows before import. They are the per-family structural-RNA totals for QC. When reading featureCounts tables with DuckDB, declare `Chr`, `Start`, `End` and `Strand` as VARCHAR, since meta-feature rows join their coordinates with `;`.

Import either a long table or featureCounts tables:

```bash
# long: sample_id, protein_id, count_s1, count_s2 (or sense[, antisense]); TSV, CSV or Parquet
sharur import-expression gene_counts.parquet --db data/my_dataset/sharur.duckdb --samples rna_samples.tsv
# featureCounts tables of one orientation
sharur import-expression lib.s2.txt --format featurecounts --orientation sense --db data/my_dataset/sharur.duckdb
```

**Flag CDS inside structural RNAs.** Gene callers place ORFs inside rRNA, tRNA, tmRNA, RNase P, SRP and 6S RNA genes, and reads from those abundant RNAs then count as "expression" of the ORFs, on either strand. Annotate the RNAs (for example Infernal `cmsearch --cut_ga` with those Rfam families), intersect them with the CDS, and pass the overlap table at import:

```bash
sharur import-expression gene_counts.parquet --db data/my_dataset/sharur.duckdb \
  --gene-flags cds_rna_overlap.tsv --min-overlap-bp 60
```

The table needs `protein_id` and `overlap_bp` (`rna_name` is kept as detail). Flagged genes leave every query's numerators, denominators and rankings; `--include-flagged` (or `include_flagged=True`) keeps them.

With `count_s1`/`count_s2`, sense is the orientation favored by most genes with at least 10 fragments (fragment totals decide when fewer than 20 genes qualify), unless the samples table gives `sense_orientation` (`s1` or `s2`). The per-gene vote matters for metatranscriptomes: a handful of ORFs inside abundant non-coding RNAs can collect more one-strand fragments than every mRNA in the library together. The sample record keeps the strand specificity, the fraction of those genes with more sense than antisense fragments; about 0.5 means an unstranded library. Proteins must exist in the dataset; `--allow-unknown-proteins` keeps the matching rows. Re-importing a sample replaces its rows.

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
- After rRNA removal, the per-library fraction of fragments assigned to CDS reflects how much of the transcribing community the dataset covers, along with the share of other non-coding RNAs (tmRNA, RNase P and SRP RNA stay abundant).
- **Activity per cell** compares a genome's share of ribosomal-protein transcripts with its share of genome copies (DNA reads / genome length) in the paired DNA sample. Ribosomal-protein genes are conserved, so at population level recount with a tighter identity filter (for example NM <= 1% of the reference span) to keep reads from uncatalogued relatives out. Weight the ribosomal-protein families equally: compute each family's length-normalized share separately and average across families. Summing fragments over the whole panel lets the head genes of ribosomal-protein operons dominate (uL3 can hold a quarter of archaeal ribosomal-protein reads), and that shifts archaeal against bacterial activity.
