# Domain-architecture search

Search proteins by the order of their domains: repeat expansions, fusions, N- or C-terminal domains, and combinations that a single-domain query cannot express.

```bash
sharur architecture "Big_* {5,} . ( VWA | VWA_2 | VWA_3 )" --db data/my_dataset/sharur.duckdb
```

```python
b.search_architecture("^ SP ( Cadherin | Big_2 )+", limit=20)
b.architecture(protein_id)          # one protein's domains, N- to C-terminal
```

Results list each matching protein with its genome, length, a compact architecture (`SP - Big_2 x12 - VWA`) and the matched domains with their protein coordinates. `total` counts every match; `--limit` controls how many are shown. Protein cards show the Pfam architecture too.

## Pattern language

| Token | Matches |
|---|---|
| `Big_2`, `PF00092` | One domain, by name or accession (case-sensitive) |
| `Big_*`, `TPR_?` | A glob over domain names |
| `.` | Any one domain |
| `?` `*` `+` `{m}` `{m,}` `{m,n}` | Repeats of the preceding token or group |
| `( A \| B )` | Either alternative; groups can be repeated |
| `^`, `$` | The N-terminal or C-terminal end of the architecture |

Separate tokens with spaces. `+` can follow a name directly (`Big_2+`); put a space before `?` and `*`, because attached to a name they are glob characters (`Big_* +` is one or more domains whose names start with `Big_`). Without anchors a pattern matches anywhere in the architecture.

Examples:

| Pattern | Finds |
|---|---|
| `TPR_* {10,}` | Ten or more consecutive TPR-family domains |
| `^ . * VWA_2 $` | Proteins ending in a VWA_2 domain |
| `^ . {40,} $` | Proteins with at least 40 resolved domains |
| `Big_* {3,} . ? ( VWA \| VWA_2 \| VWA_3 )` | Three or more Big domains, then a VWA domain, optionally one domain apart |

## How architectures are built

A protein's architecture is its domain hits with protein coordinates from the chosen sources (`--source`, Pfam by default; repeatable). Overlapping hits are resolved by E-value: hits are accepted best first, and a hit that overlaps an accepted one by more than half of the shorter domain is dropped (`--max-overlap`). Accepted domains are ordered by start position.

Resolution keeps one call per region, which is what makes repeat counts meaningful, and it means a weaker family hit under a stronger one is absent from the architecture. Use `--max-overlap 1` to keep every hit when you need all family calls in a region. Sources without coordinates (CAZy families, system callers) are left out.

Repeat counts reflect resolved Pfam hits. Divergent repeats that fall below the profile thresholds are absent, so a count is a lower bound on the repeats a protein carries.
