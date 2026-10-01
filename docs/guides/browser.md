# Browsing a dataset

`sharur browse` serves a read-only website over one dataset, for looking around without a terminal and for sharing what you find.

```bash
sharur browse --db data/my_dataset/sharur.duckdb          # http://localhost:8800/
```

Browse without knowing an identifier. The left rail offers five ways in:

| Section | What you find |
|---|---|
| **Tree of life** | A treemap of the dataset's clades; each clade page shows its subclades, genome size and assembly quality, the functions carried more often than in the rest of the dataset, its most common systems, and a sortable genome table |
| **Functions** | The labels that vary most between genomes, every label by category, and a page per label with its prevalence across the tree and example proteins |
| **Pathways** | KEGG modules with the share of genomes carrying each; a page per module with completeness by clade and how often each step is found |
| **Systems** | Defense and secretion system types and mobile-element regions, with their distribution across the tree |
| **Families** | Every Pfam domain and VOG family in the dataset. Pfam pages show the functional labels the family maps to with their evidence, its architectures and partner domains; VOG pages show the VOGdb consensus description and functional category, whether carrier genes sit in prophage regions or islands, and the Pfam domains on the same proteins. Both show protein lengths and prevalence across the tree |
| **Discover** | The dataset's extremes and oddities, each with a full list: giant proteins drawn to scale (marked when the gene runs off its contig) and the clades richest in them; Pfam domains fused on one protein only within one order, family or genus while common apart elsewhere; the longest stretches of unannotated genes (whole unannotated contigs flagged); the largest unannotated proteins; the longest domain repeats; curated system types seen in three genomes or fewer; and domain-pattern search. *Surprise me* opens a random entry |

Genome pages draw the contig landscape (click a contig to open it in the contig viewer: genes on both strands colored by function, with prophages, islands and systems marked above, and controls to pan and zoom), compare the genome's function profile with the dataset average, and list its pathways, systems and largest proteins. Protein pages draw the domain architecture to scale and the gene neighborhood colored by function category, with each label linking to its evidence.

VOG descriptions and categories come from VOGdb's `vog.annotations.tsv`, found in `data/reference/vogdb/`, `~/.sharur/vogdb/`, the installed VOGdb HMM directory (`~/.config/Astra/VOGdb`), or `$SHARUR_VOG_ANNOTATIONS`; without it VOG families appear by identifier.

**Compare** (`/compare`, or the *Compare with* box on any genome or clade page) sets two genomes or clades side by side: genome statistics, then the functions, pathways, systems and Pfam domains more common on each side, each linking to the proteins behind it. The **function heatmap** (`/heatmap`) shows chosen labels, or presets such as hydrogenases and terminal oxidases, across every clade at a chosen rank; each cell opens the proteins in that clade.

The search box suggests taxa, functions, pathways, systems, genomes and proteins as you type (press `/` to focus it). Protein pages show the amino-acid sequence with copy and FASTA download buttons; the CLI, Python API and agent-facing cards stay sequence-free.

## Presence/absence matrix

`/matrix` draws genomes as rows and features as columns: KEGG orthologs, Pfam families, KEGG modules (cells hold completeness), curated systems, or function labels. Open it from any clade page, from **Compare** (*genome-by-genome matrix*), or from a pathway page (*Steps across genomes*), which lays out the module's KOs step by step.

- **One clade.** *Most variable* shows the features carried by about half the genomes, clustered so co-inherited sets sit together. *Absences beyond incompleteness* shows the most common features whose absences exceed what genome completeness explains.
- **Two clades.** The largest differences in prevalence, each with a two-sided Fisher's exact p and a Benjamini-Hochberg q across every feature either clade carries. Related genomes share features by descent, so these values rank features and overstate independent evidence.
- **Rows** sort by taxonomy, by clustering (Jaccard distance on the shown features), or by completeness. A colour strip marks the first rank that splits the clade, and a bar shows each genome's completeness.
- **Cells** open the genome's hits for that feature. **Download TSV** saves the matrix with completeness and lineage.

**Completeness and absences.** Genome completeness comes from DFAST_QC at ingest, or from a CheckM, CheckM2 or GTDB table:

```bash
sharur import-quality checkm2/quality_report.tsv --db data/my_dataset/sharur.duckdb
sharur import-quality ar53_metadata.tsv.gz --db data/my_dataset/sharur.duckdb   # GTDB; CheckM2 columns preferred
```

Absences in genomes below 70% complete are hatched. For KOs, Pfam families and function labels in at least half of a clade, the matrix asks whether incompleteness explains the absences. If every genome carried the feature and reported it with probability equal to its completeness, the carrier count would follow a Poisson-binomial distribution. The table reports the expected and observed carriers, with q from the distribution's lower tail.

Absences beyond incompleteness have three sources: gene loss, gene calls missing from the dataset, and annotation thresholds that miss divergent homologs. Genomes with under 0.5 proteins per kb of assembly have most of their gene calls missing. They are marked in red and left out of the test. Without an assembly file, the cutoff is half the proteins their completeness implies.


## Similar proteins, structures and findings

- **Similar proteins** on protein pages come from the dataset's persistent embedding index. Build it once with `sharur build-vector-index --db data/my_dataset/sharur.duckdb`; the panel loads after the page and shows the nearest neighbors with their genome, architecture and best hit. For ESM-2 mean-pooled embeddings, cosine values crowd near 1, so the ranking carries more information than the score. `$SHARUR_BROWSE_EMBEDDINGS` points the browser at an embedding file elsewhere.
- **Structures** appear for proteins with a model in `structures/`: either a file named after the protein ID (characters outside `A-Z a-z 0-9 _ . -` written as `_`), or a JSON result file whose records name both `protein_id` and `pdb_path`, as the structure and Foldseek workflows write. The viewer colors by pLDDT, lists the Foldseek matches recorded with the model, and labels models shorter than the protein as partial.
- **Agent findings** from `findings.jsonl` (in the dataset directory or one level below) appear on the proteins, genomes, clades and system types they reference, and under **Findings** in the left rail. A finding's page reruns its verification queries against the database: a single read-only SELECT each, with file-reading functions refused and a 10-second timeout. Python and shell checks are listed and skipped.

## Loci side by side

Every system type has a **View every call as its locus** page: each call drawn with a few flanking genes, aligned on the system's core component and turned so it points right, with system genes outlined. Domain, VOG and function pages have **Stack neighborhoods**, the same view for every carrier of the family, sampled one genome at a time. Genes take their family's color in every row, so conserved neighborhoods read as vertical bands and the odd ones out stand apart. Filter by subtype or clade and widen the flank as needed.

## Curation

Flag any protein, genome, system call, family, pathway or clade as **verified**, **suspicious**, **interesting** or **follow up**, and add free-text notes. Set your name once (left rail, "Curating as"); every flag and note records its author and time, and deletions keep the history. Notes live in `browser_notes.sqlite` beside the dataset (or `--notes PATH`), apart from the evidence-backed labels. **Flags & notes** lists them with filters and exports them as TSV or JSONL.

**Triage** walks a list one item at a time: every call of a system type, every carrier of a family, or everything carrying a flag. Each item shows its locus and the curation panel; keys `1`–`4` toggle flags, `n` writes a note, `j`/`k` move on and back, `o` opens the item.

## Working faster

- **Collection:** "＋ Collect" (or `c`) on protein and genome pages builds a set kept in this browser; the Collection page downloads protein FASTA and copies a ready prompt for an agent.
- **Previews:** hover over a protein or genome link for its architecture, lineage, labels and flags.
- **Tables:** sort by any column, filter, download the rows shown as TSV, and move with `j`/`k` and Enter.
- **Genes:** `[` and `]` step to the previous and next gene on a protein page.
- **Figures:** hover over any track, neighborhood, contig view, stack or histogram for SVG and PNG downloads, drawn in the light theme for papers and slides.
- **Recently viewed** pages sit in the left rail; `?` lists every shortcut.

## Sharing

The server binds to this machine by default. To let collaborators on your network open it:

```bash
sharur browse --db data/my_dataset/sharur.duckdb --host 0.0.0.0 --share
```

This prints a link containing a random access token. The first visit exchanges the token for a browser cookie; after that, links such as `/protein/<id>` work as they are, so a page can be pasted into a chat for anyone who has opened the token link. Binding beyond this machine always requires the token.

The database opens read-only. Dataset-wide summaries (taxonomy, function prevalence, systems) load at startup, and pathway completeness and discovery lists fill in moments later; a few seconds for about 2,000 genomes. Requests share one connection, which suits a few people browsing. Agents and large campaigns use `sharur-query`.
