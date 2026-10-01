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
| **Discover** | Giant proteins with their domain architecture, the largest unannotated proteins, the longest domain repeats, and domain-pattern search |

Genome pages draw the contig landscape (click a contig to open it in the contig viewer: genes on both strands colored by function, with prophages, islands and systems marked above, and controls to pan and zoom), compare the genome's function profile with the dataset average, and list its pathways, systems and largest proteins. Protein pages draw the domain architecture to scale and the gene neighborhood colored by function category, with each label linking to its evidence.

VOG descriptions and categories come from VOGdb's `vog.annotations.tsv`, found in `data/reference/vogdb/`, `~/.sharur/vogdb/`, the Astra VOGdb directory, or `$SHARUR_VOG_ANNOTATIONS`; without it VOG families appear by identifier.

The search box suggests taxa, functions, pathways, systems, genomes and proteins as you type (press `/` to focus it). Protein pages show the amino-acid sequence with copy and FASTA download buttons; the CLI, Python API and agent-facing cards stay sequence-free.

## Similar proteins, structures and findings

- **Similar proteins** on protein pages come from the dataset's persistent embedding index. Build it once with `sharur build-vector-index --db data/my_dataset/sharur.duckdb`; the panel loads after the page and shows the nearest neighbors with their genome, architecture and best hit. For ESM-2 mean-pooled embeddings, cosine values crowd near 1, so the ranking carries more information than the score. `$SHARUR_BROWSE_EMBEDDINGS` points the browser at an embedding file elsewhere.
- **Structures** appear for proteins with a model in `structures/`: either a file named after the protein ID (characters outside `A-Z a-z 0-9 _ . -` written as `_`), or a JSON result file whose records name both `protein_id` and `pdb_path`, as the structure and Foldseek workflows write. The viewer colors by pLDDT, lists the Foldseek matches recorded with the model, and labels models shorter than the protein as partial.
- **Agent findings** from `findings.jsonl` (in the dataset directory or one level below) appear on the proteins, genomes, clades and system types they reference, and under **Findings** in the left rail. A finding's page reruns its verification queries against the database: a single read-only SELECT each, with file-reading functions refused and a 10-second timeout. Python and shell checks are listed and skipped.

## Sharing

The server binds to this machine by default. To let collaborators on your network open it:

```bash
sharur browse --db data/my_dataset/sharur.duckdb --host 0.0.0.0 --share
```

This prints a link containing a random access token. The first visit exchanges the token for a browser cookie; after that, links such as `/protein/<id>` work as they are, so a page can be pasted into a chat for anyone who has opened the token link. Binding beyond this machine always requires the token.

The database opens read-only. Dataset-wide summaries (taxonomy, function prevalence, systems) load at startup, and pathway completeness and discovery lists fill in moments later; a few seconds for about 2,000 genomes. Requests share one connection, which suits a few people browsing. Agents and large campaigns use `sharur-query`.
