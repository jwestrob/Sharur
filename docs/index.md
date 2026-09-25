<p align="center">
  <img src="assets/logo.png" alt="Sharur logo" width="360">
</p>

# Sharur

Sharur turns a folder of metagenome-assembled genomes into a database that you, and AI agents working with you, can question directly: which genomes carry a Group 3 NiFe hydrogenase, what sits next to this giant unannotated protein, how complete is the Wood–Ljungdahl pathway in this bin, and why Sharur thinks a given protein binds SAM.

It is built for datasets of thousands of genomes and millions of proteins, and for analysis by agents with limited context: every answer is bounded, typed, and free of raw sequences, and every functional label traces back to a recorded reason.

## What you get

- **One database per dataset.** Gene calls, Pfam, KOfam, HydDB, VOGdb, CAZy, DefenseFinder and TXSScan annotations, CRISPR arrays, biosynthetic gene clusters and protein embeddings, loaded into DuckDB by one resumable command.
- **Readable functional labels.** Annotations become *predicates* such as `nife_group3`, `sam_binding` or `crispr_associated`. Each one carries its evidence, and `sharur why` shows the chain for any protein. See [Predicates](predicate_construction.md).
- **Questions as commands.** Search by predicate, walk genomic neighborhoods, summarize a protein on one screen, score KEGG module completeness, find similar proteins by embedding.
- **A receipt for every dataset.** Seals record exactly what a dataset contains, so results can be tied to the data that produced them. See [Seals and provenance](concepts/seals.md).
- **Agent-ready.** Claude Code and Codex read the project's instructions and skills directly; any agent that runs shell commands or Python can use the same interface. See [Working with agents](concepts/agents.md).

## Where to start

| If you want to… | Read |
|---|---|
| Install Sharur and its tools | [Installation](getting-started/installation.md) |
| Load your genomes | [Ingest a dataset](getting-started/quickstart.md) |
| Look around a finished dataset | [First look](getting-started/first-look.md) |
| Understand a functional label | [How predicates are built](predicate_construction.md) |
| Look up a command | [Command-line reference](reference/cli.md) |
