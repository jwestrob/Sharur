# Browsing a dataset

`sharur browse` serves a read-only website over one dataset, for looking around without a terminal and for sharing what you find.

```bash
sharur browse --db data/my_dataset/sharur.duckdb          # http://localhost:8800/
```

| Page | Shows |
|---|---|
| `/` | Proteins, genomes, annotation sources, curated callers, predicate map status |
| `/protein/<id>` | Location and contig position, domain architecture drawn to scale, a clickable gene neighborhood, predicates, annotations, KEGG module steps nearby |
| `/protein/<id>/why/<predicate>` | The evidence behind one predicate |
| `/genome/<id>` | Genome metadata, largest proteins, KEGG modules at least half complete, abundance per sample |
| `/predicate/<id>` | Proteins carrying a predicate |
| `/architecture?pattern=…` | [Domain-architecture search](architecture-search.md) results |

The search box accepts a protein ID, a genome ID, a predicate, or a domain pattern. Pages contain no sequences.

## Sharing

The server binds to this machine by default. To let collaborators on your network open it:

```bash
sharur browse --db data/my_dataset/sharur.duckdb --host 0.0.0.0 --share
```

This prints a link containing a random access token. The first visit exchanges the token for a browser cookie; after that, links such as `/protein/<id>` work as they are, so a page can be pasted into a chat for anyone who has opened the token link. Binding beyond this machine always requires the token.

The database opens read-only and requests run one at a time on a single connection, which suits a few people browsing. Agents and large campaigns use `sharur-query`.
