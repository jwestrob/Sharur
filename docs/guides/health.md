# Dataset health

`sharur health` checks a dataset for the problems that quietly distort analyses. Each check reports exact counts, a few example genomes, proteins or contigs, and the command that fixes it. The whole report takes about a second on a few million proteins and only reads the database.

```bash
sharur health --db data/my_dataset/sharur.duckdb                 # markdown report
sharur health --db data/my_dataset/sharur.duckdb --format json   # for agents and scripts
sharur health --db data/my_dataset/sharur.duckdb --strict        # exit 1 when a check fails
```

The browser shows the same report at `/health`. A pill in the left rail counts the checks that need care and links there.

## Statuses

| Status | Meaning |
|---|---|
| **fail** | Fix before analysing |
| **warn** | Results that depend on it need care |
| **info** | Worth knowing |
| **ok** | All clear |

## Checks

| Area | Check | What it looks at |
|---|---|---|
| Database | Schema version | The schema against the version this Sharur expects |
| Database | Dataset seal | Table row counts and schema against `dataset.seal.json`; `sharur verify-seal` gives the full comparison |
| Genomes | Completeness estimates | Genomes with completeness, which the presence/absence matrix uses to weigh absences |
| Genomes | Contamination | Genomes above 10% contamination |
| Genomes | Genomes without proteins | Genomes with no gene calls at all |
| Genes | Gene calls per genome | Genomes with under 0.5 proteins per kb of assembly, where most gene calls are missing. Without an assembly file, the cutoff is half the proteins its completeness implies |
| Genes | Strand encoding | `+`/`-` and `1`/`-1` mixed in one dataset |
| Genes | Genes without a genomic position | Proteins that are their own contig, contig IDs shared across genomes, and proteins stacked at coordinate 0; neighborhoods and operons skip them |
| Genes | Proteins sharing coordinates | Several accessions mapped to one locus |
| Genes | Contig lengths | Contigs without genes, and lengths taken from the last gene, which `sharur backfill-contig-context` replaces with assembly lengths |
| Annotations | Annotation coverage | Hits per source. Each source's typical share of a genome's proteins predicts how many hits a genome of a given size should carry; a genome expecting more than seven hits from a source and holding none was most likely never searched with it (the pattern an interrupted search leaves). Sources that hit about one gene per genome, such as hydrogenases, stay below that expectation, and rows written by system callers record calls, so a genome without them holds a result |
| Annotations | Function labels vs installed maps | Whether stored labels came from the maps installed now |
| Callers | CRISPR array scan | Genomes with a MinCED report among those with assemblies |
| Callers | Curated system callers | Defense, secretion and CRISPR-Cas calls, given the hits they need |
| References | Local KEGG build | Module completeness and KO names, from `sharur setup-kegg` |

Assemblies are found where the browser finds them: the stage 00 manifest, then `genomes_fna/`, `genomes_fna_new/`, `source/` or `assemblies/` beside the database, named after their genome.
