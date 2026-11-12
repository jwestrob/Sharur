# Example dataset: 16 complete public genomes

A small, closed-genome test set for a first ingest and for checking an installation against
known biology. All 16 are GenBank assemblies at level "Complete Genome"; every one of their 27
replicons is circular in its NCBI record. Together they hold 39.4 Mb; Prodigal calls 38,554 proteins on them.

| Genome | Assembly | Size (Mb) | Replicons | Why it is here |
|---|---|---:|---:|---|
| Nanoarchaeum equitans Kin4-M | GCA_000008085.1 | 0.49 | 1 | DPANN episymbiont, reduced genome |
| Ignicoccus hospitalis KIN4/I | GCA_000017945.1 | 1.30 | 1 | Host of N. equitans |
| Nanosynbacter lyticus TM7x | GCA_000803625.1 | 0.71 | 1 | CPR (Saccharibacteria) episymbiont |
| Schaalia odontolytica XH001 | GCA_005696695.1 | 2.36 | 1 | Host of TM7x |
| Methanocaldococcus jannaschii DSM 2661 | GCA_000091665.1 | 1.74 | 3 | Two extrachromosomal elements; CRISPR arrays; hydrogenases |
| Pyrococcus furiosus DSM 3638 | GCA_008245085.1 | 1.89 | 1 | CRISPR-Cas types I and III; NiFe hydrogenases |
| Haloferax volcanii DS2 | GCA_000025685.1 | 4.01 | 5 | Four plasmids; CRISPR-Cas type I-B on pHV4 |
| Sulfolobus acidocaldarius DSM 639 | GCA_000012285.1 | 2.23 | 1 | CRISPR arrays |
| Escherichia coli K-12 MG1655 | GCA_000005845.2 | 4.64 | 1 | CRISPR-Cas I-E; four NiFe hydrogenases; type I restriction-modification |
| Salmonella Typhimurium LT2 | GCA_000006945.2 | 4.95 | 2 | Virulence plasmid pSLT; type III secretion |
| Streptococcus pyogenes SF370 | GCA_000006785.2 | 1.85 | 1 | Source of SpCas9 (CRISPR-Cas II-A) |
| Desulfovibrio vulgaris Hildenborough | GCA_000195755.1 | 3.77 | 2 | Plasmid pDV1; FeFe and NiFe hydrogenases |
| Bacillus subtilis 168 | GCA_000009045.1 | 4.22 | 1 | Nonribosomal peptide and polyketide gene clusters |
| Buchnera aphidicola APS | GCA_000009605.1 | 0.66 | 3 | Reduced endosymbiont; plasmids pLeu and pTrp |
| Mycoplasmoides genitalium G37 | GCA_000027325.1 | 0.58 | 1 | Genetic code 4 (UGA = Trp) |
| Vibrio cholerae N16961 | GCA_000006745.1 | 4.03 | 2 | Two chromosomes; type VI secretion |

`genomes.tsv` lists each assembly with its replicons (`accession:role:name:length:topology`),
taken from the NCBI Datasets and nuccore records.

## Download

```bash
bash examples/complete_genomes/download.sh genomes/
```

The script needs `curl`, `gzip` and `md5sum`. It fetches each assembly's `*_genomic.fna.gz`
from the NCBI FTP site one at a time, checks it against NCBI's `md5checksums.txt` and writes
`genomes/<genome_id>.fna` (39 MB in total). Stage 00 derives each genome ID from that file
name.

## Ingest

```bash
sharur-ingest --input-dir genomes --data-dir data/complete_genomes \
  --output data/complete_genomes/sharur.duckdb --skip-embeddings
```

Add `--profile slurm --slurm-cpus N` on a cluster ([QUICKSTART.md](../../QUICKSTART.md)).
`--skip-embeddings` leaves out Stage 06 for machines without a GPU.

## Expected answers

These hold for a correct installation and make quick checks of a finished database:

- Every genome carries the universal ribosomal proteins (for example Pfam `Ribosomal_L2_C`,
  `Ribosomal_S3_C`, `Ribosomal_L5`, `Ribosomal_L14`).
- MinCED finds CRISPR arrays in N. equitans (two arrays), H. volcanii (one on the main
  chromosome, two on pHV4), E. coli K-12, S. pyogenes SF370, P. furiosus, M. jannaschii and
  S. acidocaldarius.
- `sharur cas-type` calls subtype II-A in S. pyogenes SF370, I-E in E. coli K-12 and I-B on
  pHV4 in H. volcanii.
- The hydrogenase classifier assigns NiFe group 1 and group 4 enzymes in E. coli, group 3 and
  group 4 in P. furiosus, and FeFe enzymes in D. vulgaris.
- Stage 03 calls genes with genetic code 11, so M. genitalium (code 4) shows more and shorter
  proteins than its GenBank annotation.
