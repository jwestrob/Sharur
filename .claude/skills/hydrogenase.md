# Hydrogenase Skill

Inventory, curate, and interpret [NiFe]-, [FeFe]-, and [Fe]-hydrogenase calls in a
Sharur dataset, keeping observed evidence, reference assignments, and physiological
interpretation as separate layers.

**CONCURRENCY:** Independent `read_only=True` sessions may query in parallel.
Serialize subgroup refreshes and every other DuckDB write.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

> **Literature dispatch:** For ambiguous annotations, unknown Foldseek hits, or
> comparative claims ("first known", "largest"), dispatch a literature agent.
> Read `.claude/skills/literature.md` for protocols.

---

## What Sharur's subgroup call is

Proteins with a HydDB HMM hit (`source = 'hyddb'`, class `NiFe` / `FeFe` / `Fe_only`)
receive the subgroup of their best DIAMOND match among the installed HydDB reference
sequences. This is a **Sharur nearest-reference assignment**. The published HydDB
classifier (Søndergaard et al. 2016) additionally votes over k=4 neighbors, screens
homologous non-hydrogenase families, and uses the downstream gene to separate [FeFe]
Group A subtypes. Describe Sharur calls as "assigned to [NiFe] Group 1h by nearest
HydDB reference" and treat them as provisional.

Each HydDB-hit protein has one row in `hydrogenase_classifications`:

| Column | Meaning |
|---|---|
| `outcome` | `assigned`, `class_conflict`, `no_reference_hit`, `missing_sequence`, `unparsed_reference_label` |
| `discovery_classes`, `discovery_hit_count`, `discovery_best_score` | HMM evidence (raw `hyddb` rows stay in `annotations`) |
| `reference_accession`, `reference_label`, `reference_class`, `reference_subgroup` | Nearest reference and its exact label |
| `pident`, `evalue`, `bitscore` | Best alignment (no coverage is recorded yet) |
| `interpretation_status`, `reference_role` | Table 1 status (`characterized` / `putative` / `unresolved` / `unverified`) and role text |
| `has_nifese_hases`, `has_fe_hyd`, `has_complex1`, `has_hmd` | Pfam domain observations on the same protein |
| `curation_status`, `curation_reason` | `domain_check_cleared` or `needs_curation`, with the reason |
| `ko_support`, `ko_support_detail` | How the protein's KOfam hits relate to the assigned subgroup (below), with each associated KO's HydDB label counts |
| `reference_release`, `reference_sha256`, `classifier_version` | Provenance |

Derived labels appear as `hyddb_subgroup` annotations and V2 atoms with relation
`supports`:

- **Structural:** `nife_group1`–`nife_group4`, `fefe_groupA`–`fefe_groupC`, `fe_only_hydrogenase`.
- **Functional:** emitted only for subgroups Table 1 characterizes (for example
  `uptake_hydrogenase` for 1h and 2a, `h2_sensor` for 2b, `nad_coupled` for 3d).
  Putative, unresolved, and unverified subgroups carry structure only. [FeFe] Group A
  subtypes carry structure only because A1–A4 require gene organization.
- **Review flags:** `hyddb_needs_curation` (catalytic-domain check not met),
  `hyddb_class_conflict` (HMM class and reference class disagree; no subgroup labels
  are emitted), `hyddb_ko_conflict` (every associated KOfam hit captures HydDB
  references outside the assigned group), and `hydrogenase_complex1_review` (HydDB
  NiFe hit with Complex I-superfamily domains and no NiFeSe_Hases).
- **Supporting flag:** `hyddb_ko_supported` (a KOfam hit captures HydDB references of
  the assigned subgroup at ≥80%). It is a quality flag and adds no functional claim.

### KOfam support (KO → HydDB subgroup associations)

`sharur/predicates/mappings/data/kegg_hyddb_snapshot.tsv` records, for every KO KEGG
names as a hydrogenase, the HydDB labels of the reference hydrogenases that score at or
above that KO's KOfam threshold (`scripts/build_kegg_hyddb_snapshot.py`). For example
K15830 (hycE) captures `[NiFe]_Group_4a` 47/47; K14090 (echE) captures 4e 130, 4c 27
and 4g 13 of 170. `sharur.hydrogenase.ko_association` exposes these associations.

They describe which references a KOfam profile captures, so they support or question a
nearest-reference call and never change it. `ko_support` grades each assignment:
`subgroup` (≥80% of an associated KO's references carry the assigned label),
`compatible` (some do), `group` (none do; some share the group), `conflict` (none
share the group or type, for every associated KO), `none` (no associated KO with ≥5
references). [FeFe] KOs carry no associations: HydDB [FeFe] references are
catalytic-domain segments that score below KOfam's full-length thresholds.

Report KOfam support as corroboration beside the reference label ("assigned to [NiFe]
Group 4a by nearest reference; its hycE KOfam hit captures Group 4a references 47/47").

The exact subgroup label is always in `reference_label` and the direct-access
predicate `hyddb_subgroup:<Group_xx>`. Report that label alongside any interpretation.

## Why domain evidence needs care

[NiFe]-hydrogenase large subunits and respiratory Complex I subunit D share the
Complex1_49kDa superfamily; HydDB's own NiFe catalytic-domain check uses that
superfamily. The Pfam profile NiFeSe_Hases is typical of Groups 1–3 and can also occur
on Group 4 assignments. Complex1_49kDa/30kDa domains without NiFeSe_Hases are
compatible with either a Group 4 hydrogenase or a respiratory Complex I subunit.
Domain observations therefore set a review state; they neither validate a
physiological role nor exclude a hydrogenase. Operon context resolves the case.

Pfam hits may store either the accession (`PF00374`) or the profile name
(`NiFeSe_Hases`) in `accession`; match both columns.

---

## Usage

```
/hydrogenase                     # Full inventory with curation
/hydrogenase --group4-only       # Focus on Group 4 curation
/hydrogenase --genome GENOME_ID  # Single genome analysis
```

---

## Prompt

You are curating hydrogenase calls in a metagenomic dataset. Every call that carries a
review flag needs neighborhood evidence before it appears as a hydrogenase in a finding.
Follow the steps in order.

### Step 1: Inventory

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

print(b.store.execute("""
    SELECT outcome, reference_class, curation_status, COUNT(*) AS proteins
    FROM hydrogenase_classifications
    GROUP BY ALL ORDER BY ALL
"""))

print(b.store.execute("""
    SELECT reference_label, interpretation_status, curation_status, COUNT(*) AS proteins
    FROM hydrogenase_classifications WHERE outcome = 'assigned'
    GROUP BY ALL ORDER BY reference_label
"""))
```

If `hydrogenase_classifications` is absent, the database predates the reconciled
classifier: run `python scripts/classify_hydrogenases.py --db ...` (dry run) and ask
before publishing.

### Step 2: Partition by domain evidence

```python
rows = b.store.execute("""
    SELECT protein_id, reference_class, reference_label, curation_status, curation_reason,
           has_nifese_hases, has_fe_hyd, has_complex1
    FROM hydrogenase_classifications
    WHERE outcome = 'assigned'
""")
cleared = [r for r in rows if r[3] == "domain_check_cleared"]
review = [r for r in rows if r[3] == "needs_curation"]
complex1_only = [r for r in review if r[7] and not r[5]]
print(f"catalytic domain observed: {len(cleared)}; needs curation: {len(review)} "
      f"(Complex I-superfamily only: {len(complex1_only)})")

conflicts = b.store.execute("""
    SELECT protein_id, discovery_classes, reference_label, curation_reason
    FROM hydrogenase_classifications WHERE outcome = 'class_conflict'
""")
print(f"class conflicts: {len(conflicts)}")
```

A cleared domain check means the catalytic domain for the assigned type is present.
The subgroup remains a nearest-reference assignment, and the physiological role remains
an interpretation.

### Step 3: Neighborhood evidence for every review-flagged call

Classify each review-flagged protein as **supported**, **Complex I context**, or
**ambiguous** from its ±6-gene neighborhood.

```python
COMPLEX_I_KOS = {f"K{n:05d}" for n in range(330, 344)}          # nuoA-N
HYC_KOS = {f"K{n:05d}" for n in range(15827, 15834)}            # hycB-G, hycA (formate hydrogenlyase)
HYF_KOS = {f"K{n:05d}" for n in range(12136, 12146)}            # hyfA-J (hydrogenase-4)
ECH_KOS = {f"K{n:05d}" for n in range(14086, 14092)}            # echA-F
MATURATION_KOS = {"K04651", "K04652", "K04653", "K04654", "K04655", "K04656", "K03605"}  # hypA-F, hyaD/hybD
COMPLEX_I_PFAM_NAMES = {"Complex1_30kDa", "Complex1_49kDa", "Oxidored_q6", "Oxidored_q5_N", "Oxidored_q4"}
HYDROGENASE_PFAM_NAMES = {"NiFeSe_Hases", "Fe_hyd_lg_C", "Fe_hyd_SSU"}
NAME_PREFIXES = {"ech": "Ech", "eha": "Eha", "ehb": "Ehb", "hyc": "Hyc", "hyf": "Hyf", "coo": "Coo"}


def neighborhood_evidence(protein_id, window=6):
    nbr = b.get_neighborhood(protein_id, window=window, all_annotations=True)
    complex_i, hydrogenase = set(), set()
    for gene in nbr:
        if gene.get("protein_id") == protein_id:
            continue
        for ann in gene.get("annotations", []):
            acc, name = ann.get("accession") or "", ann.get("name") or ""
            if acc in COMPLEX_I_KOS or {acc, name} & COMPLEX_I_PFAM_NAMES:
                complex_i.add(f"{acc} {name}")
            if acc in HYC_KOS | HYF_KOS | ECH_KOS | MATURATION_KOS or {acc, name} & HYDROGENASE_PFAM_NAMES:
                hydrogenase.add(f"{acc} {name}")
            for prefix, complex_name in NAME_PREFIXES.items():
                if name.lower().startswith(prefix):
                    hydrogenase.add(f"{complex_name}: {name}")
    if hydrogenase and not complex_i:
        return "supported", complex_i, hydrogenase
    if complex_i and not hydrogenase:
        return "complex_i_context", complex_i, hydrogenase
    return "ambiguous", complex_i, hydrogenase


verdicts = {pid: neighborhood_evidence(pid) for pid, *_ in review}
for label in ("supported", "complex_i_context", "ambiguous"):
    print(label, sum(v[0] == label for v in verdicts.values()))
```

Report the counts per verdict for this dataset. Rates differ widely between datasets
and lineages; compute them, and cite them only for the dataset measured.

### Step 4: Interpret by subgroup

Use the Table 1 interpretation stored in `reference_role` and `interpretation_status`.
For Group 4 complexes, name a complex (Ech, Hyc, Hyf, Eha, Ehb, Coo) only when its
subunit genes appear in the neighborhood, and pair the name with the HydDB subgroup
label from `reference_label`. When neighborhood genes and the reference label point to
different subgroups, report both and mark the case unresolved.

### Step 5: Metabolic context

For supported calls, record genome-level context (methanogenesis, Wood-Ljungdahl,
fermentation, respiratory chains) from predicates in the same bin. Directionality of
Group 3 and Group 4 enzymes depends on that context.

### Step 6: Figures

> **Read `.claude/skills/visualize.md` before generating figures.** Use
> `plot_locus_multisource.py` for report figures and `b.visualize_neighborhood()` for
> quick checks. Include representative supported and Complex I-context loci.

### Step 7: Inventory finding

```python
summary = {
    "assigned": len(rows),
    "domain_check_cleared": len(cleared),
    "needs_curation": len(review),
    "neighborhood": {k: sum(v[0] == k for v in verdicts.values())
                     for k in ("supported", "complex_i_context", "ambiguous")},
    "class_conflicts": len(conflicts),
}
```

---

## Subgroup reference (Søndergaard et al. 2016, Table 1)

| Subgroup | Table 1 name | Sharur functional predicates |
|---|---|---|
| NiFe 1a, 1b, 1c, 1d, 1e, 1g, 1h, 1j, 1k | Respiratory H2-uptake | `uptake_hydrogenase` |
| NiFe 1f | Oxygen-protecting; role unresolved | — |
| NiFe 1i | Coriobacteria-type (putative) | — |
| NiFe 1l | Present in the installed reference; interpretation unverified | — |
| NiFe 2a | Cyanobacteria-type; respiratory uptake | `uptake_hydrogenase` |
| NiFe 2b | Histidine kinase-linked; H2 sensing | `h2_sensor` |
| NiFe 2c, 2d, 2e | Putative or unresolved roles | — |
| NiFe 3a | F420-coupled | `bidirectional_hydrogenase`, `f420_reducing` |
| NiFe 3b | NADP-coupled | `bidirectional_hydrogenase`, `nadp_coupled` |
| NiFe 3c | Heterodisulfide reductase-linked; bifurcating | `bidirectional_hydrogenase`, `heterodisulfide_reductase_linked`, `bifurcating_hydrogenase` |
| NiFe 3d | NAD-coupled | `bidirectional_hydrogenase`, `nad_coupled` |
| NiFe 4a | Formate hydrogenlyase | `h2_evolving`, `formate_coupled` |
| NiFe 4b | Formate-respiring, Mrp-linked | `h2_evolving`, `energy_conserving_hydrogenase`, `formate_coupled` |
| NiFe 4c | Carbon monoxide-respiring | `h2_evolving`, `energy_conserving_hydrogenase`, `co_coupled` |
| NiFe 4d | Ferredoxin-coupled, Mrp-linked | `h2_evolving`, `energy_conserving_hydrogenase`, `ferredoxin_coupled` |
| NiFe 4e | Ferredoxin-coupled, Ech-type | `h2_evolving`, `energy_conserving_hydrogenase`, `ferredoxin_coupled`, `ech_hydrogenase` |
| NiFe 4f | Formate-coupled (putative) | — |
| NiFe 4g | Ferredoxin-coupled (putative) | — |
| NiFe 4h | Ferredoxin-coupled, Eha-type | `h2_evolving`, `energy_conserving_hydrogenase`, `ferredoxin_coupled` |
| NiFe 4i | Ferredoxin-coupled, Ehb-type | `h2_evolving`, `energy_conserving_hydrogenase`, `ferredoxin_coupled` |
| FeFe A1–A4 | Prototypical / glutamate synthase-linked (putative) / bifurcating / formate dehydrogenase-linked; assigned from the downstream gene | — |
| FeFe B | Colonic-type (putative) | — |
| FeFe C1–C3 | Putative sensory/regulatory | — |
| Fe | Methenyl-H4MPT dehydrogenase | `methanogen_hydrogenase` |

A nearest-reference A3 match is a provisional Group A assignment. A bifurcation claim
for an [FeFe] enzyme needs the downstream NuoF-domain gene, as in HydDB.

## Key Pfam profiles

| Profile | Accession | Reading |
|---|---|---|
| NiFeSe_Hases | PF00374 | NiFe catalytic domain; clears the NiFe domain check |
| Fe_hyd_lg_C / Fe_hyd_SSU | PF02906 / PF02256 | FeFe catalytic domains; clear the FeFe domain check |
| HMD | PF03201 | [Fe]-hydrogenase domain |
| Complex1_49kDa / Complex1_30kDa | PF00346 / PF00329 | Shared by Group 4 hydrogenases and Complex I; review state |
| Proton_antipo_M | PF00361 | Shared by Complex I, Ech, Mrp antiporters |
| Ni_hydr_CYTB | PF01292 | Cytochrome b; needs a catalytic subunit for context |

---

## Language

- Say: "N proteins were assigned to [NiFe] Group 1h by nearest HydDB reference; M carry
  the NiFeSe_Hases catalytic domain."
- Say: "Group 1h enzymes characterized to date provide electrons for aerobic
  respiration; H2 uptake is a hypothesis for these proteins."
- Say: "no hydrogenase was detected in this MAG (N contigs)" for absences.
- Report raw call counts with curated counts, name the curation method, and keep
  subgroup labels beside any interpretation.

### Claim escalation

| Claim | Evidence required |
|---|---|
| "HydDB HMM hit" | Raw `hyddb` row |
| "Assigned to [NiFe] Group X by nearest reference" | `hydrogenase_classifications.outcome = 'assigned'` |
| "Catalytic domain observed" | `curation_status = 'domain_check_cleared'` |
| "KOfam hit corroborates the subgroup" | `ko_support = 'subgroup'` |
| "Supported Group 4 hydrogenase" | Step 3 neighborhood verdict `supported` |
| Named complex (Ech, Hyc, ...) | Subunit genes of that complex in the neighborhood |
| Physiological role in this organism | Above, plus genome-level metabolic context |

---

## Output specification

Each analysis appends at least one inventory finding:

```jsonl
{
  "id": "survey-NNN",
  "title": "Hydrogenase inventory: A assigned (B catalytic-domain observed); C review-flagged, of which D supported by neighborhood context",
  "category": "energy_metabolism",
  "evidence": "HydDB HMM discovery; Sharur nearest-reference subgroup (HydDB <release>); Pfam catalytic-domain check; neighborhood curation for review-flagged calls",
  "verification": [
    {"claim": "A assigned", "query": "SELECT COUNT(*) FROM hydrogenase_classifications WHERE outcome = 'assigned'", "expected": "A"}
  ],
  "phase": "survey"
}
```

Include: raw and curated counts, the curation method, subgroup-level breakdown with exact
labels, and metabolic context for supported calls.
