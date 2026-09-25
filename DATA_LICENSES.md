# Data Licenses and Attribution

Sharur's code is MIT-licensed (`LICENSE`). Some shipped data files are derived
from third-party resources; this file lists them, their licenses, and the
resources Sharur uses without redistributing.

## Shipped derived data

| File(s) under `sharur/predicates/mappings/data/` | Derived from | License |
|---|---|---|
| `pfam_predicates.tsv`, `pfam_evidence_snapshot.tsv` | Pfam family names and descriptions (EMBL-EBI) | CC0 1.0 |
| (GO evidence in the above) | Gene Ontology and InterPro `pfam2go` | CC BY 4.0 |
| (ENZYME evidence in the above), `ec_sam_dependent.tsv`, EC class entries in `kegg_map.py` | Expasy ENZYME (`enzyme.dat`, `enzclass.txt`), SIB Swiss Institute of Bioinformatics | CC BY 4.0 |
| `swissprot:` evidence counts in the Pfam/CAZy maps, `kegg_swissprot_consensus.tsv` | UniProtKB/Swiss-Prot, The UniProt Consortium | CC BY 4.0 |
| `kegg_hyddb_snapshot.tsv`, `sharur/hydrogenase/subgroups.py` | HydDB reference labels and Table 1 interpretations (Søndergaard, Pedersen & Greening 2016, *Sci Rep* 6:34212) | CC BY 4.0 |
| `cazy_predicates.tsv` | CAZy family identifiers and class definitions; Swiss-Prot consensus (above) | Identifiers; see above |
| `kegg_brite_predicates.tsv`, `kegg_module_predicates.tsv`, `kegg_predicate_proposals.tsv` | Sharur's own rules, keyed by KEGG identifiers and BRITE category names | MIT (Sharur) |

## Built on the user's machine (never redistributed)

- **KEGG** (Kanehisa Laboratories): `sharur setup-kegg` fetches KEGG data
  through KEGG REST, which is for academic use, or builds from a licensed KEGG
  copy with `--inputs`. The resulting `kegg_predicates.tsv`,
  `kegg_evidence_snapshot.tsv` and `kegg_modules.tsv` stay local.
  Terms: https://www.kegg.jp/kegg/legal.html
- **KOfam** profiles and `ko_list`, used for annotation (via Astra) and by the
  HydDB x KOfam snapshot builder.

## Downloaded by the user for annotation

Pfam-A HMMs, KOfam, HydDB references, dbCAN/CAZy databases, VOGdb,
DefenseFinder and TXSScan models, and other reference databases are fetched
by the user from their providers under the providers' terms; see
`CITATIONS.md` for citations.
