"""
Co-location engine for MacSyFinder-style gene system validation.

Replaces the MacSyFinder subprocess calls in validate_defense_systems.py and
validate_secretion_systems.py with an in-process implementation that reads
HMM hits directly from DuckDB and performs co-localization using the same
XML model definitions.

Why: MacSyFinder takes 20-30 minutes on large datasets (3,500+ genomes) even
with --previous-run. This engine does the same work in <10 seconds by reading
hits from the annotations table and protein positions from the proteins table,
both already indexed in DuckDB.

The algorithm is:
  1. Parse XML model definitions (DefenseFinder or TXSScan)
  2. Load HMM hits from DuckDB (annotations + proteins JOIN)
  3. For each model, filter hits to relevant profiles, cluster on each contig
     within inter_gene_max_space, validate clusters against quorum rules
  4. Handle loner genes and multi_loci merging (TXSScan only)
  5. Resolve overlapping candidates with an exact, deterministic conflict solver

Output: (systems_df, genes_df) DataFrames matching the schema produced by
validate_defense_systems.py and validate_secretion_systems.py.
"""

from __future__ import annotations

import functools
import logging
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import duckdb
import pandas as pd

from sharur.storage.migrations import get_current_version, run_migrations


logger = logging.getLogger(__name__)

_AMBIGUOUS_HMM_NAME_PREFIX = "__sharur_ambiguous_hmm_name__:"


# ------------------------------------------------------------------ #
# Data classes
# ------------------------------------------------------------------ #


@dataclass
class GeneSpec:
    """One gene slot in a system model definition."""

    name: str  # HMM profile name (= annotations.accession)
    presence: str  # mandatory | accessory | neutral | forbidden
    loner: bool = False
    multi_system: bool = False
    multi_model: bool = False
    inter_gene_max_space: int | None = None
    exchangeables: list[str] = field(default_factory=list)


@dataclass
class SystemModel:
    """A system definition parsed from a MacSyFinder XML model file."""

    name: str  # e.g. "Gabija", "T3SS"
    family: str  # e.g. "defense-finder-models", "TXSScan"
    inter_gene_max_space: int
    min_mandatory_genes_required: int
    min_genes_required: int
    max_nb_genes: int
    multi_loci: bool
    genes: list[GeneSpec]

    # Derived lookup tables — built by _build_lookups() after parsing
    n_mandatory_genes: int = 0  # total mandatory gene slots in model
    n_accessory_genes: int = 0
    profile_to_gene_slot: dict[str, GeneSpec] = field(default_factory=dict)
    all_profiles: set[str] = field(default_factory=set)
    mandatory_profiles: set[str] = field(default_factory=set)
    accessory_profiles: set[str] = field(default_factory=set)
    neutral_profiles: set[str] = field(default_factory=set)
    forbidden_profiles: set[str] = field(default_factory=set)
    loner_profiles: set[str] = field(default_factory=set)
    multi_system_profiles: set[str] = field(default_factory=set)

    def _build_lookups(self) -> None:
        """Populate derived lookup tables from gene specs."""
        self.profile_to_gene_slot.clear()
        self.all_profiles.clear()
        self.mandatory_profiles.clear()
        self.accessory_profiles.clear()
        self.neutral_profiles.clear()
        self.forbidden_profiles.clear()
        self.loner_profiles.clear()
        self.multi_system_profiles.clear()
        self.n_mandatory_genes = sum(1 for g in self.genes if g.presence == "mandatory")
        self.n_accessory_genes = sum(1 for g in self.genes if g.presence == "accessory")
        for gene in self.genes:
            # The parent profile itself
            all_names = [gene.name, *gene.exchangeables]
            for profile_name in all_names:
                # Map every profile (parent + exchangeables) to the parent gene slot
                self.profile_to_gene_slot[profile_name] = gene
                self.all_profiles.add(profile_name)

                if gene.presence == "mandatory":
                    self.mandatory_profiles.add(profile_name)
                elif gene.presence == "accessory":
                    self.accessory_profiles.add(profile_name)
                elif gene.presence == "neutral":
                    self.neutral_profiles.add(profile_name)
                elif gene.presence == "forbidden":
                    self.forbidden_profiles.add(profile_name)

                if gene.loner:
                    self.loner_profiles.add(profile_name)
                if gene.multi_system:
                    self.multi_system_profiles.add(profile_name)


@dataclass
class SystemHit:
    """A validated system instance — one cluster of genes matching a model."""

    model_name: str
    model_family: str
    genome_id: str
    contig_ids: list[str]  # exactly one replicon/contig
    protein_ids: list[str]
    gene_indices: list[int]
    profile_names: list[str]
    scores: list[float]
    gene_presences: list[str]  # per-gene: mandatory/accessory/neutral/loner
    mandatory_satisfied: int
    total_satisfied: int  # mandatory + accessory (not neutral)
    score: float  # sum of HMM scores
    primary_mandatory_hits: int = 0  # Bug 2 fix: mandatory slots hit by their PRIMARY profile
    n_mandatory_total: int = 0  # Bug 2 fix: total mandatory slots in the model
    max_genes_total: int = 0
    model_score: float = 0.0

    @property
    def protein_set(self) -> set[str]:
        return set(self.protein_ids)

    @property
    def multi_system_proteins(self) -> set[str]:
        """Proteins that are flagged multi_system and should not block conflicts."""
        # Populated during validation; stored as a set alongside the hit
        return getattr(self, "_multi_system_proteins", set())

    @multi_system_proteins.setter
    def multi_system_proteins(self, val: set[str]) -> None:
        self._multi_system_proteins = val

    @property
    def multi_model_proteins(self) -> set[str]:
        """Proteins explicitly allowed to participate in different models."""
        return getattr(self, "_multi_model_proteins", set())

    @multi_model_proteins.setter
    def multi_model_proteins(self, val: set[str]) -> None:
        self._multi_model_proteins = val

    @property
    def true_loner_proteins(self) -> set[str]:
        """Out-of-cluster loners, shareable between same-model candidates."""
        return getattr(self, "_true_loner_proteins", set())

    @true_loner_proteins.setter
    def true_loner_proteins(self, val: set[str]) -> None:
        self._true_loner_proteins = val

    @property
    def rescued_multi_system_proteins(self) -> set[str]:
        """Multi-system hits copied into an otherwise incomplete candidate."""
        return getattr(self, "_rescued_multi_system_proteins", set())

    @rescued_multi_system_proteins.setter
    def rescued_multi_system_proteins(self, val: set[str]) -> None:
        self._rescued_multi_system_proteins = val

    @property
    def contig_id(self) -> str:
        """Return the single contig backing this structured call."""
        if len(self.contig_ids) != 1:
            raise ValueError(
                f"Structured system {self.model_name!r} spans {len(self.contig_ids)} contigs"
            )
        return self.contig_ids[0]


# ------------------------------------------------------------------ #
# XML model parser
# ------------------------------------------------------------------ #


def _parse_bool(val: str | None) -> bool:
    """Parse XML boolean attribute (1, True, true → True; everything else → False)."""
    if val is None:
        return False
    return val.strip().lower() in ("1", "true")


def _parse_model_xml(xml_path: Path, family: str) -> SystemModel | None:
    """Parse a single MacSyFinder XML model file into a SystemModel."""
    try:
        tree = ET.parse(xml_path)
    except ET.ParseError as e:
        logger.warning(f"Failed to parse XML {xml_path}: {e}")
        return None

    root = tree.getroot()
    if root.tag != "model":
        logger.warning(f"Unexpected root tag '{root.tag}' in {xml_path}")
        return None

    # Model-level attributes
    model_inter_gene_max_space = int(root.get("inter_gene_max_space", "5"))
    min_mandatory = int(root.get("min_mandatory_genes_required", "1"))
    min_genes = int(root.get("min_genes_required", "1"))
    max_nb_genes_text = root.get("max_nb_genes")
    multi_loci = _parse_bool(root.get("multi_loci"))

    # Model name = XML filename stem
    model_name = xml_path.stem

    genes: list[GeneSpec] = []
    for gene_elem in root.findall("gene"):
        name = gene_elem.get("name", "")
        presence = gene_elem.get("presence", "mandatory")
        loner = _parse_bool(gene_elem.get("loner"))
        multi_system = _parse_bool(gene_elem.get("multi_system"))
        multi_model = _parse_bool(gene_elem.get("multi_model"))
        gene_max_space = gene_elem.get("inter_gene_max_space")
        gene_inter_gene_max_space = (
            int(gene_max_space) if gene_max_space is not None else None
        )

        # Parse exchangeables
        exchangeables: list[str] = []
        exch_elem = gene_elem.find("exchangeables")
        if exch_elem is not None:
            for exch_gene in exch_elem.findall("gene"):
                exch_name = exch_gene.get("name", "")
                if exch_name:
                    exchangeables.append(exch_name)

        genes.append(
            GeneSpec(
                name=name,
                presence=presence,
                loner=loner,
                multi_system=multi_system,
                multi_model=multi_model,
                inter_gene_max_space=gene_inter_gene_max_space,
                exchangeables=exchangeables,
            )
        )

    model = SystemModel(
        name=model_name,
        family=family,
        inter_gene_max_space=model_inter_gene_max_space,
        min_mandatory_genes_required=min_mandatory,
        min_genes_required=min_genes,
        max_nb_genes=(
            int(max_nb_genes_text)
            if max_nb_genes_text is not None
            else sum(gene.presence in {"mandatory", "accessory"} for gene in genes)
        ),
        multi_loci=multi_loci,
        genes=genes,
    )
    model._build_lookups()
    return model


def parse_models(models_dir: Path) -> dict[str, SystemModel]:
    """Parse all XML model files under a MacSyFinder models directory.

    Handles both DefenseFinder and TXSScan directory structures:
      - DefenseFinder: models_dir/defense-finder-models/definitions/{category}/{family}/{model}.xml
      - TXSScan: models_dir/TXSScan/definitions/{domain}/{membrane}/{model}.xml

    Args:
        models_dir: Root models directory (e.g. ~/.macsyfinder/models).

    Returns:
        Dict mapping model name to SystemModel. Names are unique per model file.
        If two models have the same filename, later ones overwrite earlier ones
        (this shouldn't happen in practice).
    """
    models: dict[str, SystemModel] = {}

    for family_dir in sorted(models_dir.iterdir()):
        if not family_dir.is_dir():
            continue
        family_name = family_dir.name  # "defense-finder-models" or "TXSScan"

        definitions_dir = family_dir / "definitions"
        if not definitions_dir.is_dir():
            continue

        # Recursively find all XML files under definitions/
        for xml_path in sorted(definitions_dir.rglob("*.xml")):
            model = _parse_model_xml(xml_path, family=family_name)
            if model is not None:
                models[model.name] = model

    return models


def parse_models_for_family(models_dir: Path, model_family: str) -> dict[str, SystemModel]:
    """Parse XML models for a specific family (e.g. 'defense-finder-models', 'TXSScan').

    Args:
        models_dir: Root models directory (e.g. ~/.macsyfinder/models).
        model_family: Directory name under models_dir (e.g. 'defense-finder-models').

    Returns:
        Dict mapping model name to SystemModel.
    """
    family_dir = models_dir / model_family
    if not family_dir.is_dir():
        logger.error(f"Model family directory not found: {family_dir}")
        return {}

    definitions_dir = family_dir / "definitions"
    if not definitions_dir.is_dir():
        logger.error(f"No definitions/ subdirectory in {family_dir}")
        return {}

    models: dict[str, SystemModel] = {}
    for xml_path in sorted(definitions_dir.rglob("*.xml")):
        model = _parse_model_xml(xml_path, family=model_family)
        if model is not None:
            models[model.name] = model

    return models


# ------------------------------------------------------------------ #
# Hit loader
# ------------------------------------------------------------------ #


@dataclass
class HitRecord:
    """A single HMM hit with protein position context."""

    protein_id: str
    accession: str
    score: float
    contig_id: str
    bin_id: str
    gene_index: int
    evalue: float = float("inf")
    annotation_id: int = 0

    @property
    def replicon_key(self) -> tuple[str, str]:
        """Database-stable replicon identity.

        ``contig_id`` is not guaranteed to be globally unique across bins
        (reference accessions such as ``WP`` can recur), so every spatial
        operation must include the owning genome.
        """
        return self.bin_id, self.contig_id


def load_hits(conn: duckdb.DuckDBPyConnection, source: str) -> list[HitRecord]:
    """Load HMM hits from the DuckDB annotations + proteins tables.

    Args:
        conn: Open DuckDB connection (read-only is fine).
        source: Annotation source filter (e.g. 'defensefinder', 'txsscan').

    Returns:
        List of HitRecord, sorted by contig_id then gene_index.
    """
    rows = conn.execute(
        """
        SELECT a.protein_id, a.accession, a.score,
               p.contig_id, p.bin_id, p.gene_index,
               a.evalue, a.annotation_id
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = ?
        ORDER BY p.bin_id, p.contig_id, p.gene_index
    """,
        [source],
    ).fetchall()

    return [
        HitRecord(
            protein_id=r[0],
            accession=r[1],
            score=float(r[2]) if r[2] is not None else 0.0,
            contig_id=r[3],
            bin_id=r[4],
            gene_index=int(r[5]) if r[5] is not None else 0,
            evalue=float(r[6]) if r[6] is not None else float("inf"),
            annotation_id=int(r[7]),
        )
        for r in rows
    ]


def _best_hit_sort_key(hit: HitRecord) -> tuple[float, float, str, int, str]:
    """Deterministic ordering equivalent to MacSyFinder's score-first filter."""
    return (
        -hit.score,
        hit.evalue,
        hit.accession,
        hit.annotation_id,
        hit.protein_id,
    )


def _select_best_hits_per_position(hits: list[HitRecord]) -> list[HitRecord]:
    """Keep one best profile hit per protein position.

    MacSyFinder applies ``get_best_hits(..., key="score")`` before filtering
    hits into individual models. That global competition is important: related
    secretion-system profiles must not let one protein satisfy several model
    components or competing system subtypes.

    Sharur uses ``(bin_id, contig_id, gene_index)`` as the
    replicon-position key because contig labels are not globally unique. Ties
    are resolved by e-value, accession, annotation ID, and protein ID so the
    result is stable across processes and Python hash seeds.
    """
    best_by_position: dict[tuple[str, str, int], HitRecord] = {}
    for hit in hits:
        key = (*hit.replicon_key, hit.gene_index)
        current = best_by_position.get(key)
        if current is None or _best_hit_sort_key(hit) < _best_hit_sort_key(current):
            best_by_position[key] = hit

    return sorted(
        best_by_position.values(),
        key=lambda h: (
            h.bin_id,
            h.contig_id,
            h.gene_index,
            h.protein_id,
            h.accession,
        ),
    )


def _index_hits_by_accession(
    hits: list[HitRecord],
) -> dict[str, list[HitRecord]]:
    """Group hits by accession (profile name)."""
    by_acc: dict[str, list[HitRecord]] = defaultdict(list)
    for h in hits:
        by_acc[h.accession].append(h)
    return dict(by_acc)


# ------------------------------------------------------------------ #
# HMM profile name mapping (Bug 1 fix)
# ------------------------------------------------------------------ #


def _build_hmm_name_mapping(profiles_dir: Path) -> dict[str, str]:
    """Build a mapping from HMM internal NAME to filename stem.

    The annotations table stores HMM internal NAME fields (e.g. 'flgB', 'abc'),
    but XML model <gene name="..."> uses filename stems (e.g. 'Flg_flgB',
    'T1SS_abc'). This function reads all .hmm files in profiles_dir and builds
    a mapping {internal_name: filename_stem} for cases where they differ.

    Args:
        profiles_dir: Directory containing .hmm profile files.

    Ambiguous internal names are mapped to a reserved non-profile value. Astra's
    flattened annotation rows no longer retain which source HMM file produced
    such a hit, so allowing either filename would turn an ambiguous observation
    into a named component call.

    Returns:
        Dict mapping internal_name -> filename_stem for unambiguous mismatches,
        or to a reserved value for ambiguous names.
    """
    stems_by_name: dict[str, list[str]] = defaultdict(list)
    if not profiles_dir.is_dir():
        logger.warning(f"Profiles directory not found: {profiles_dir}")
        return {}

    name_re = re.compile(r"^NAME\s+(\S+)")
    for hmm_path in sorted(profiles_dir.glob("*.hmm")):
        stem = hmm_path.stem
        try:
            with open(hmm_path) as f:
                for line in f:
                    m = name_re.match(line)
                    if m:
                        internal_name = m.group(1)
                        stems_by_name[internal_name].append(stem)
                        break
        except (OSError, UnicodeDecodeError):
            continue

    mapping: dict[str, str] = {}
    ambiguous: dict[str, list[str]] = {}
    for internal_name, stems in stems_by_name.items():
        unique_stems = sorted(set(stems))
        if len(unique_stems) > 1:
            ambiguous[internal_name] = unique_stems
            mapping[internal_name] = f"{_AMBIGUOUS_HMM_NAME_PREFIX}{internal_name}"
        elif internal_name != unique_stems[0]:
            mapping[internal_name] = unique_stems[0]

    if ambiguous:
        examples = "; ".join(
            f"{name} -> {','.join(stems)}" for name, stems in sorted(ambiguous.items())[:5]
        )
        logger.warning(
            "Quarantining %d ambiguous HMM internal names that map to multiple "
            "profile files (%s). The raw annotation remains available, but it "
            "cannot support a named system component.",
            len(ambiguous),
            examples,
        )

    return mapping


def _translate_hit_accessions(
    hits: list[HitRecord], name_to_stem: dict[str, str]
) -> list[HitRecord]:
    """Translate hit accession values from internal HMM names to filename stems.

    Modifies hits in-place for efficiency. Only translates accessions that exist
    in the mapping (i.e. where internal_name != filename_stem).

    Args:
        hits: HitRecords loaded from DuckDB.
        name_to_stem: Mapping from internal_name -> filename_stem.

    Returns:
        The same list (mutated in-place).
    """
    if not name_to_stem:
        return hits

    translated = 0
    for h in hits:
        if h.accession in name_to_stem:
            h.accession = name_to_stem[h.accession]
            translated += 1

    if translated > 0:
        logger.info(f"  Translated {translated} hit accessions (internal NAME -> filename stem)")

    return hits


def _find_profiles_dir(family_dir: Path) -> Path | None:
    """Locate the profiles/ directory for a model family.

    Handles both DefenseFinder (profiles/ at top level) and TXSScan
    (profiles/ at top level) directory structures.
    """
    profiles_dir = family_dir / "profiles"
    if profiles_dir.is_dir():
        return profiles_dir

    # Fallback: search for any directory containing .hmm files
    for d in family_dir.rglob("*.hmm"):
        return d.parent

    return None


# ------------------------------------------------------------------ #
# Clustering
# ------------------------------------------------------------------ #


def _cluster_contig(hits: list[HitRecord], max_gap: int) -> list[list[HitRecord]]:
    """Group hits on the same contig within max_gap gene positions.

    Hits must all be from the same contig. They are sorted by gene_index
    and grouped greedily. ``max_gap`` is the number of intervening genes,
    matching MacSyFinder's ``inter_gene_max_space`` semantics.

    Args:
        hits: HitRecords from a single contig.
        max_gap: Maximum gene-index distance between consecutive hits
                 in the same cluster.

    Returns:
        List of clusters, each a list of HitRecords.
    """
    if not hits:
        return []

    replicons = {hit.replicon_key for hit in hits}
    if len(replicons) != 1:
        raise ValueError(
            f"_cluster_contig requires one replicon, received {sorted(replicons)}"
        )

    sorted_hits = sorted(
        _select_best_hits_per_position(hits),
        key=lambda h: (h.gene_index, h.protein_id, h.accession),
    )
    clusters: list[list[HitRecord]] = [[sorted_hits[0]]]

    for hit in sorted_hits[1:]:
        intervening_genes = hit.gene_index - clusters[-1][-1].gene_index - 1
        if 0 <= intervening_genes <= max_gap:
            clusters[-1].append(hit)
        else:
            clusters.append([hit])

    return clusters


def _effective_max_gap(
    left: HitRecord,
    right: HitRecord,
    model: SystemModel,
) -> int:
    """Return the model/gene-specific gap allowed between two hits."""
    left_slot = model.profile_to_gene_slot[left.accession]
    right_slot = model.profile_to_gene_slot[right.accession]
    left_gap = left_slot.inter_gene_max_space
    right_gap = right_slot.inter_gene_max_space

    if left_gap is None and right_gap is None:
        return model.inter_gene_max_space
    if left_gap is None:
        return int(right_gap)
    if right_gap is None:
        return int(left_gap)
    return min(left_gap, right_gap)


def _hits_colocate(
    left: HitRecord,
    right: HitRecord,
    model: SystemModel,
) -> bool:
    """Return whether two ordered hits co-localize under model semantics."""
    if left.replicon_key != right.replicon_key:
        return False
    intervening_genes = right.gene_index - left.gene_index - 1
    return 0 <= intervening_genes <= _effective_max_gap(left, right, model)


def _loner_sort_key(hit: HitRecord, model: SystemModel) -> tuple:
    """Prefer a primary loner profile, then the deterministic best HMM hit."""
    slot = model.profile_to_gene_slot[hit.accession]
    return (hit.accession != slot.name, *_best_hit_sort_key(hit))


def _build_model_clusters(
    hits: list[HitRecord],
    model: SystemModel,
) -> tuple[list[list[HitRecord]], list[HitRecord]]:
    """Build regular clusters and true loners on one contig.

    This mirrors the ordered-replicon behavior in ``macsypy.cluster``:

    - cluster by the number of intervening genes;
    - retain regular clusters only when they encode multiple functions
      (except valid one-gene models);
    - treat a ``loner=True`` hit as a true loner only when it is outside a
      regular multi-function cluster;
    - keep one deterministic best true loner per canonical function.
    """
    if not hits:
        return [], []

    relevant = [
        hit for hit in _select_best_hits_per_position(hits) if hit.accession in model.all_profiles
    ]
    if not relevant:
        return [], []

    replicons = {hit.replicon_key for hit in relevant}
    if len(replicons) != 1:
        raise ValueError(
            f"_build_model_clusters requires one replicon, received {sorted(replicons)}"
        )

    relevant.sort(
        key=lambda hit: (
            hit.gene_index,
            *_best_hit_sort_key(hit),
        )
    )
    scaffolds: list[list[HitRecord]] = [[relevant[0]]]
    for hit in relevant[1:]:
        if _hits_colocate(scaffolds[-1][-1], hit, model):
            scaffolds[-1].append(hit)
        else:
            scaffolds.append([hit])

    regular_clusters: list[list[HitRecord]] = []
    loners_by_function: dict[str, HitRecord] = {}

    for scaffold in scaffolds:
        slots = [model.profile_to_gene_slot[hit.accession] for hit in scaffold]
        gene_refs = {hit.accession for hit in scaffold}

        if len(gene_refs) > 1:
            if not all(slot.presence == "neutral" for slot in slots):
                regular_clusters.append(scaffold)
            continue

        slot = slots[0]
        if slot.loner:
            best = min(scaffold, key=lambda hit: _loner_sort_key(hit, model))
            current = loners_by_function.get(slot.name)
            if current is None or _loner_sort_key(best, model) < _loner_sort_key(current, model):
                loners_by_function[slot.name] = best
        elif len(scaffold) == 1 and model.min_genes_required == 1:
            regular_clusters.append(scaffold)

    regular_clusters.sort(
        key=lambda cluster: (
            cluster[0].gene_index,
            cluster[-1].gene_index,
            tuple(hit.protein_id for hit in cluster),
        )
    )
    true_loners = sorted(
        loners_by_function.values(),
        key=lambda hit: (hit.gene_index, hit.protein_id, hit.accession),
    )
    return regular_clusters, true_loners


# ------------------------------------------------------------------ #
# Validation
# ------------------------------------------------------------------ #


def _validate_cluster(
    cluster: list[HitRecord],
    model: SystemModel,
    *,
    true_loner_proteins: set[str] | None = None,
    rescued_multi_system_proteins: set[str] | None = None,
    regular_clusters: list[list[HitRecord]] | None = None,
) -> SystemHit | None:
    """Check if a cluster of hits satisfies a model's quorum constraints.

    Performs:
      1. Map each hit's accession to a gene slot via the model's lookup tables
      2. Enforce one protein position and one contig
      3. Reject if any forbidden profile is present
      4. Count satisfied mandatory and accessory gene slots (by canonical name)
      5. Check quorum: min_mandatory_genes_required, min_genes_required

    Args:
        cluster: HitRecords from one contig.
        model: The SystemModel to validate against.
        true_loner_proteins: Protein IDs represented by out-of-cluster loners.
        rescued_multi_system_proteins: Protein IDs copied from another valid
            candidate under MacSyFinder's ``multi_system`` rule.
        regular_clusters: Original regular cluster partition, used for
            MacSyFinder-style model scoring.

    Returns:
        A SystemHit if the cluster is valid, None otherwise.
    """
    true_loner_proteins = set(true_loner_proteins or ())
    rescued_multi_system_proteins = set(rescued_multi_system_proteins or ())
    out_of_cluster_proteins = (
        true_loner_proteins | rescued_multi_system_proteins
    )

    # Repeat-domain and competing-profile rows cannot create multiple genes.
    unique_hits = _select_best_hits_per_position(cluster)
    relevant_hits = [hit for hit in unique_hits if hit.accession in model.profile_to_gene_slot]
    if not relevant_hits:
        return None

    replicons = {hit.replicon_key for hit in relevant_hits}
    if len(replicons) != 1:
        return None

    # 1. Classify hits by gene slot
    satisfied_mandatory: set[str] = set()  # canonical gene slot names
    satisfied_accessory: set[str] = set()
    satisfied_neutral: set[str] = set()
    primary_mandatory_slots: set[str] = set()  # Bug 2: slots hit by their PRIMARY profile
    multi_system_pids: set[str] = set()
    multi_model_pids: set[str] = set()

    for hit in relevant_hits:
        gene_slot = model.profile_to_gene_slot[hit.accession]
        if gene_slot.presence == "forbidden":
            return None

        slot_name = gene_slot.name  # canonical parent name

        if gene_slot.presence == "mandatory":
            satisfied_mandatory.add(slot_name)
            # Bug 2: track if this was a PRIMARY profile hit (not exchangeable)
            if hit.accession == gene_slot.name:
                primary_mandatory_slots.add(slot_name)
        elif gene_slot.presence == "accessory":
            satisfied_accessory.add(slot_name)
        elif gene_slot.presence == "neutral":
            satisfied_neutral.add(slot_name)

        if gene_slot.multi_system:
            multi_system_pids.add(hit.protein_id)
        if gene_slot.multi_model:
            multi_model_pids.add(hit.protein_id)

    # 2. Check quorum
    n_mandatory = len(satisfied_mandatory)
    n_total = n_mandatory + len(satisfied_accessory)

    if n_mandatory < model.min_mandatory_genes_required:
        return None
    if n_total < model.min_genes_required:
        return None

    # 3. Build the SystemHit.
    relevant_sorted = sorted(
        relevant_hits,
        key=lambda h: (h.gene_index, h.protein_id, h.accession),
    )
    genome_id = relevant_hits[0].bin_id

    presences: list[str] = []
    for hit in relevant_sorted:
        gene_slot = model.profile_to_gene_slot[hit.accession]
        if hit.protein_id in true_loner_proteins:
            presences.append("loner")
        elif hit.protein_id in rescued_multi_system_proteins:
            presences.append("multi_system")
        else:
            presences.append(gene_slot.presence)

    model_score = _score_candidate(
        relevant_sorted,
        model,
        out_of_cluster_proteins=out_of_cluster_proteins,
        regular_clusters=regular_clusters,
    )

    system_hit = SystemHit(
        model_name=model.name,
        model_family=model.family,
        genome_id=genome_id,
        contig_ids=[relevant_sorted[0].contig_id],
        protein_ids=[h.protein_id for h in relevant_sorted],
        gene_indices=[h.gene_index for h in relevant_sorted],
        profile_names=[h.accession for h in relevant_sorted],
        scores=[h.score for h in relevant_sorted],
        gene_presences=presences,
        mandatory_satisfied=n_mandatory,
        total_satisfied=n_total,
        score=sum(h.score for h in relevant_sorted),
        primary_mandatory_hits=len(primary_mandatory_slots),
        n_mandatory_total=model.n_mandatory_genes,
        max_genes_total=model.max_nb_genes,
        model_score=model_score,
    )
    system_hit.multi_system_proteins = multi_system_pids
    system_hit.multi_model_proteins = multi_model_pids
    system_hit.true_loner_proteins = true_loner_proteins
    system_hit.rescued_multi_system_proteins = rescued_multi_system_proteins
    return system_hit


def _profile_role_weight(hit: HitRecord, model: SystemModel) -> float:
    """Return MacSyFinder's default role/exchangeable weight for one hit."""
    slot = model.profile_to_gene_slot[hit.accession]
    role_weight = {
        "mandatory": 1.0,
        "accessory": 0.5,
        "neutral": 0.0,
        "forbidden": 0.0,
    }[slot.presence]
    if hit.accession != slot.name:
        role_weight *= 0.8
    return role_weight


def _score_candidate(
    hits: list[HitRecord],
    model: SystemModel,
    *,
    out_of_cluster_proteins: set[str],
    regular_clusters: list[list[HitRecord]] | None,
) -> float:
    """Compute the default MacSyFinder model score for conflict selection."""
    if regular_clusters is None:
        regular_clusters = [
            [hit for hit in hits if hit.protein_id not in out_of_cluster_proteins]
        ]

    score = 0.0
    cluster_functions: list[set[str]] = []
    for regular in regular_clusters:
        best_by_function: dict[str, float] = {}
        for hit in regular:
            if hit.protein_id in out_of_cluster_proteins:
                continue
            slot = model.profile_to_gene_slot[hit.accession]
            if slot.presence == "forbidden":
                continue
            function = slot.name
            hit_weight = _profile_role_weight(hit, model)
            if len(regular) == 1 and slot.multi_system:
                hit_weight *= 0.7
            best_by_function[function] = max(
                best_by_function.get(function, 0.0),
                hit_weight,
            )
        score += sum(best_by_function.values())
        cluster_functions.append(set(best_by_function))

    for gene in model.genes:
        if gene.presence not in {"mandatory", "accessory"}:
            continue
        occurrences = sum(gene.name in functions for functions in cluster_functions)
        if occurrences > 1:
            score -= (occurrences - 1) * 1.5

    regular_functions = set().union(*cluster_functions) if cluster_functions else set()
    for hit in hits:
        if hit.protein_id not in out_of_cluster_proteins:
            continue
        slot = model.profile_to_gene_slot[hit.accession]
        if slot.name not in regular_functions:
            score += _profile_role_weight(hit, model) * 0.7

    return score


def _cluster_functions(
    cluster: list[HitRecord],
    model: SystemModel,
) -> set[str]:
    return {
        model.profile_to_gene_slot[hit.accession].name
        for hit in cluster
        if hit.accession in model.profile_to_gene_slot
    }


def _candidate_cluster_combinations(
    regular_clusters: list[list[HitRecord]],
    true_loners: list[HitRecord],
    model: SystemModel,
) -> list[SystemHit]:
    """Build all model candidates within one contig.

    This mirrors MacSyFinder's two ordered-replicon passes:

    1. test every permitted regular-cluster/true-loner combination;
    2. copy the best ``multi_system`` hit per function from a valid candidate
       into otherwise rejected combinations and test them again.

    All combinations stay replicon-local. Conflict resolution later chooses
    the maximum-score compatible set, so optional loci are not discarded
    prematurely.
    """
    if model.multi_loci:
        regular_index_combinations = [
            combo
            for size in range(1, len(regular_clusters) + 1)
            for combo in combinations(range(len(regular_clusters)), size)
        ]
    else:
        regular_index_combinations = [(index,) for index in range(len(regular_clusters))]

    # An empty regular selection permits true-loner-only models.
    regular_index_combinations.append(())

    candidates: list[SystemHit] = []
    rejected_units: list[tuple[list[list[HitRecord]], list[HitRecord]]] = []
    seen_memberships: set[tuple[str, ...]] = set()

    for regular_indices in regular_index_combinations:
        selected_regular = [regular_clusters[index] for index in regular_indices]
        regular_functions = (
            set().union(*(_cluster_functions(cluster, model) for cluster in selected_regular))
            if selected_regular
            else set()
        )

        eligible_loner_indices = [
            index
            for index, loner in enumerate(true_loners)
            if model.profile_to_gene_slot[loner.accession].name not in regular_functions
        ]

        for loner_count in range(len(eligible_loner_indices) + 1):
            for loner_indices in combinations(eligible_loner_indices, loner_count):
                if not regular_indices and not loner_indices:
                    continue

                selected_loners = [true_loners[index] for index in loner_indices]
                hits = [hit for regular in selected_regular for hit in regular] + selected_loners
                loner_pids = {hit.protein_id for hit in selected_loners}
                candidate = _validate_cluster(
                    hits,
                    model,
                    true_loner_proteins=loner_pids,
                    regular_clusters=selected_regular,
                )
                if candidate is None:
                    rejected_units.append((selected_regular, selected_loners))
                    continue

                membership = tuple(candidate.protein_ids)
                if membership in seen_memberships:
                    continue
                seen_memberships.add(membership)
                candidates.append(candidate)

    # MacSyFinder's second pass can reuse a gene explicitly marked
    # multi_system to rescue a different occurrence of the same model.
    hit_by_protein = {
        hit.protein_id: hit
        for cluster in regular_clusters
        for hit in cluster
    }
    hit_by_protein.update({hit.protein_id: hit for hit in true_loners})

    best_multi_system_by_function: dict[str, HitRecord] = {}
    for candidate in candidates:
        for protein_id in candidate.multi_system_proteins:
            hit = hit_by_protein[protein_id]
            function = model.profile_to_gene_slot[hit.accession].name
            current = best_multi_system_by_function.get(function)
            if current is None or _loner_sort_key(hit, model) < _loner_sort_key(
                current, model
            ):
                best_multi_system_by_function[function] = hit

    best_multi_system_hits = [
        best_multi_system_by_function[function]
        for function in sorted(best_multi_system_by_function)
    ]
    for selected_regular, selected_loners in rejected_units:
        base_hits = [
            hit for regular in selected_regular for hit in regular
        ] + selected_loners
        base_functions = _cluster_functions(base_hits, model)
        eligible_multi_system_hits = [
            hit
            for hit in best_multi_system_hits
            if model.profile_to_gene_slot[hit.accession].name not in base_functions
        ]

        for rescue_count in range(1, len(eligible_multi_system_hits) + 1):
            for rescue_hits in combinations(eligible_multi_system_hits, rescue_count):
                candidate = _validate_cluster(
                    [*base_hits, *rescue_hits],
                    model,
                    true_loner_proteins={
                        hit.protein_id for hit in selected_loners
                    },
                    rescued_multi_system_proteins={
                        hit.protein_id for hit in rescue_hits
                    },
                    regular_clusters=selected_regular,
                )
                if candidate is None:
                    continue
                membership = tuple(candidate.protein_ids)
                if membership in seen_memberships:
                    continue
                seen_memberships.add(membership)
                candidates.append(candidate)

    return sorted(
        candidates,
        key=lambda hit: (
            hit.contig_id,
            min(hit.gene_indices),
            hit.model_name,
            tuple(hit.gene_indices),
            tuple(hit.profile_names),
        ),
    )


def _detect_model_systems(
    hits: list[HitRecord],
    model: SystemModel,
) -> list[SystemHit]:
    """Detect model instances independently on each contig."""
    relevant_by_replicon: dict[tuple[str, str], list[HitRecord]] = defaultdict(list)
    for hit in hits:
        if hit.accession in model.all_profiles:
            relevant_by_replicon[hit.replicon_key].append(hit)

    systems: list[SystemHit] = []
    for replicon_key in sorted(relevant_by_replicon):
        regular_clusters, true_loners = _build_model_clusters(
            relevant_by_replicon[replicon_key],
            model,
        )
        systems.extend(
            _candidate_cluster_combinations(
                regular_clusters,
                true_loners,
                model,
            )
        )
    return systems


# ------------------------------------------------------------------ #
# Conflict resolution
# ------------------------------------------------------------------ #


def resolve_conflicts(hits: list[SystemHit]) -> list[SystemHit]:
    """Select a deterministic maximum-score set of compatible systems.

    Compatibility follows MacSyFinder's model:

    - systems on different contigs are independent;
    - same-model candidates cannot share regular hits (true loners and copied
      ``multi_system`` hits may be shared);
    - different models may share a hit only when that hit is ``multi_model``
      in both models.

    Each connected conflict component is solved exactly as a weighted
    independent-set problem. Components are normally tiny because conflicts
    require an identical protein.

    Args:
        hits: All candidate SystemHits across all models and genomes.

    Returns:
        The selected non-conflicting subset.
    """
    if not hits:
        return []

    selected: list[SystemHit] = []
    by_replicon: dict[tuple[str, str], list[SystemHit]] = defaultdict(list)
    for hit in hits:
        by_replicon[(hit.genome_id, hit.contig_id)].append(hit)

    for replicon_key in sorted(by_replicon):
        replicon_hits = sorted(
            by_replicon[replicon_key],
            key=_system_sort_key,
        )
        conflicts = _conflict_masks(replicon_hits)
        remaining = set(range(len(replicon_hits)))

        while remaining:
            seed = min(remaining)
            component = {seed}
            frontier = [seed]
            remaining.remove(seed)
            while frontier:
                current = frontier.pop()
                neighbors = {index for index in remaining if conflicts[current] & (1 << index)}
                if neighbors:
                    component.update(neighbors)
                    remaining.difference_update(neighbors)
                    frontier.extend(sorted(neighbors))

            component_indices = sorted(component)
            selected_indices = _solve_conflict_component(
                replicon_hits,
                conflicts,
                component_indices,
            )
            selected.extend(replicon_hits[index] for index in selected_indices)

    return sorted(selected, key=_system_sort_key)


def _system_sort_key(hit: SystemHit) -> tuple:
    return (
        hit.genome_id,
        hit.contig_id,
        min(hit.gene_indices),
        hit.model_name,
        -hit.model_score,
        tuple(hit.gene_indices),
        tuple(hit.protein_ids),
        tuple(hit.profile_names),
    )


def _systems_compatible(left: SystemHit, right: SystemHit) -> bool:
    if left.genome_id != right.genome_id or left.contig_id != right.contig_id:
        return True

    shared = left.protein_set & right.protein_set
    if not shared:
        return True

    if left.model_name == right.model_name:
        shareable = (
            left.true_loner_proteins
            | right.true_loner_proteins
            | left.rescued_multi_system_proteins
            | right.rescued_multi_system_proteins
        )
        regular_shared = shared - shareable
        return not regular_shared

    allowed = left.multi_model_proteins & right.multi_model_proteins
    return shared <= allowed


def _conflict_masks(hits: list[SystemHit]) -> list[int]:
    masks = [0] * len(hits)
    for left_index, right_index in combinations(range(len(hits)), 2):
        if not _systems_compatible(hits[left_index], hits[right_index]):
            masks[left_index] |= 1 << right_index
            masks[right_index] |= 1 << left_index
    return masks


def _solution_metrics(
    hits: list[SystemHit],
    selected: tuple[int, ...],
) -> tuple[float, int, int, float, tuple[int, ...]]:
    if not selected:
        return (0.0, 0, 0, 0.0, ())

    ordered_systems = sorted(
        (hits[index] for index in selected),
        key=lambda hit: (
            tuple(hit.gene_indices),
            hit.model_name,
            -hit.model_score,
        ),
    )
    wholeness = sum(
        hits[index].total_satisfied / max(hits[index].max_genes_total, 1) for index in selected
    ) / len(selected)
    return (
        sum(hits[index].model_score for index in selected),
        sum(len(hits[index].protein_ids) for index in selected),
        len(selected),
        wholeness,
        tuple(
            position
            for system in ordered_systems
            for position in system.gene_indices
        ),
    )


def _better_solution(
    hits: list[SystemHit],
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    left_metrics = _solution_metrics(hits, left)
    right_metrics = _solution_metrics(hits, right)
    if left_metrics != right_metrics:
        return left if left_metrics > right_metrics else right

    left_signature = tuple(_system_sort_key(hits[index]) for index in left)
    right_signature = tuple(_system_sort_key(hits[index]) for index in right)
    return left if left_signature < right_signature else right


def _solve_conflict_component(
    hits: list[SystemHit],
    conflict_masks: list[int],
    component_indices: list[int],
) -> tuple[int, ...]:
    """Solve one conflict component exactly with memoized branching."""

    component_mask = sum(1 << index for index in component_indices)

    @functools.cache
    def solve(mask: int) -> tuple[int, ...]:
        if not mask:
            return ()

        active = [index for index in component_indices if mask & (1 << index)]
        pivot = max(
            active,
            key=lambda index: (
                (conflict_masks[index] & mask).bit_count(),
                -index,
            ),
        )
        pivot_bit = 1 << pivot

        without_pivot = solve(mask & ~pivot_bit)
        with_mask = mask & ~pivot_bit & ~conflict_masks[pivot]
        with_pivot = tuple(sorted((pivot, *solve(with_mask))))
        return _better_solution(hits, with_pivot, without_pivot)

    return solve(component_mask)


# ------------------------------------------------------------------ #
# Main orchestrator
# ------------------------------------------------------------------ #


def _find_models_dir() -> Path | None:
    """Locate the MacSyFinder models directory."""
    candidates = [
        Path.home() / ".macsyfinder" / "models",
        Path.home() / ".mdmlab" / "macsyfinder" / "models",
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return None


def _model_family_for_source(source: str) -> str:
    """Map annotation source name to model family directory name."""
    mapping = {
        "defensefinder": "defense-finder-models",
        "txsscan": "TXSScan",
    }
    return mapping.get(source, source)


def validate_systems(
    db_path: str | Path,
    models_dir: str | Path | None = None,
    source: str = "defensefinder",
    model_family: str | None = None,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run co-location validation on HMM hits from DuckDB.

    This is the top-level entry point that replaces the MacSyFinder subprocess
    calls in validate_defense_systems.py and validate_secretion_systems.py.

    Args:
        db_path: Path to sharur.duckdb.
        models_dir: Path to MacSyFinder models root directory
                    (e.g. ~/.macsyfinder/models). Auto-detected if None.
        source: Annotation source in the DuckDB annotations table
                ('defensefinder' or 'txsscan').
        model_family: Model family directory name (e.g. 'defense-finder-models',
                      'TXSScan'). Auto-detected from source if None.
        verbose: Log progress messages.

    Returns:
        (systems_df, genes_df) — DataFrames matching the output schema of
        validate_defense_systems.py / validate_secretion_systems.py.

        systems_df columns:
            system_id, genome_id, contig_id, type, subtype, protein_ids,
            genes_count, profile_names, score, hit_id

        genes_df columns:
            protein_id, system_id, accession, score, presence
    """
    t0 = time.time()

    db_path = Path(db_path)

    # Resolve models directory
    if models_dir is None:
        models_dir = _find_models_dir()
        if models_dir is None:
            raise FileNotFoundError(
                "MacSyFinder models directory not found. Expected at ~/.macsyfinder/models/"
            )
    else:
        models_dir = Path(models_dir)

    # Resolve model family
    if model_family is None:
        model_family = _model_family_for_source(source)

    # Phase 1: Parse models
    if verbose:
        logger.info(f"Parsing {model_family} models from {models_dir / model_family}...")
    models = parse_models_for_family(models_dir, model_family)
    if not models:
        raise RuntimeError(
            f"No parseable models found for family {model_family!r} under {models_dir}"
        )
    if verbose:
        logger.info(f"  Loaded {len(models)} system models")

    # Phase 2: Load hits from DuckDB
    if verbose:
        logger.info(f"Loading HMM hits (source='{source}') from {db_path}...")
    conn = duckdb.connect(str(db_path), read_only=True)
    all_hits = load_hits(conn, source)
    conn.close()

    if not all_hits:
        if verbose:
            logger.info("  No HMM hits found — nothing to validate")
        return pd.DataFrame(), pd.DataFrame()
    if verbose:
        logger.info(f"  Loaded {len(all_hits)} HMM hits")

    # Bug 1 fix: translate hit accessions from internal HMM NAMEs to
    # filename stems so they match the XML <gene name="..."> values
    family_dir = models_dir / model_family
    profiles_dir = _find_profiles_dir(family_dir)
    if profiles_dir is not None:
        name_to_stem = _build_hmm_name_mapping(profiles_dir)
        if name_to_stem:
            if verbose:
                logger.info(
                    f"  Found {len(name_to_stem)} profile name mismatches "
                    f"(internal NAME != filename stem) in {profiles_dir}"
                )
            _translate_hit_accessions(all_hits, name_to_stem)
    else:
        raise FileNotFoundError(
            f"Could not find the HMM profiles directory for {model_family!r} "
            f"under {family_dir}"
        )

    # MacSyFinder resolves competing profiles globally before any model sees
    # the hits. This also collapses repeated-domain rows on the same protein.
    raw_hit_count = len(all_hits)
    all_hits = _select_best_hits_per_position(all_hits)
    if verbose:
        logger.info(
            "  Selected %d best protein-position hits from %d raw/domain rows",
            len(all_hits),
            raw_hit_count,
        )

    # Build deterministic indices after best-hit selection.
    hits_by_accession = _index_hits_by_accession(all_hits)

    # Build a quick set of all accessions present in the data
    all_accessions_in_data = set(hits_by_accession.keys())

    # Phase 3+4: Cluster and validate per model
    if verbose:
        logger.info("Clustering and validating...")

    candidate_hits: list[SystemHit] = []
    models_with_hits = 0

    for model_name in sorted(models):
        model = models[model_name]
        # Quick check: does this model have ANY relevant profiles in the data?
        relevant_profiles = model.all_profiles & all_accessions_in_data
        if not relevant_profiles:
            continue
        models_with_hits += 1

        model_hits = [
            hit
            for profile in sorted(relevant_profiles)
            for hit in hits_by_accession.get(profile, ())
        ]
        candidate_hits.extend(_detect_model_systems(model_hits, model))

    if verbose:
        logger.info(
            f"  {models_with_hits} models matched hits; "
            f"{len(candidate_hits)} candidate systems before conflict resolution"
        )

    # Phase 5: Conflict resolution
    selected_hits = resolve_conflicts(candidate_hits)

    if verbose:
        logger.info(f"  {len(selected_hits)} systems after conflict resolution")

    # Build output DataFrames
    systems_df, genes_df = _build_output_dataframes(selected_hits, source)

    elapsed = time.time() - t0
    if verbose:
        n_genomes = systems_df["genome_id"].nunique() if not systems_df.empty else 0
        logger.info(
            f"\nCo-location validation complete in {elapsed:.1f}s: "
            f"{len(systems_df)} systems, {len(genes_df)} gene assignments "
            f"across {n_genomes} genomes"
        )

        if not systems_df.empty:
            type_counts = systems_df.groupby("type").size().sort_values(ascending=False)
            logger.info("\nSystem type summary:")
            for stype, count in type_counts.head(20).items():
                logger.info(f"  {stype}: {count}")
            if len(type_counts) > 20:
                logger.info(f"  ... and {len(type_counts) - 20} more types")

    return systems_df, genes_df


def _build_output_dataframes(
    hits: list[SystemHit], _source: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert SystemHits into the output DataFrames expected by downstream code.

    systems_df matches the schema from validate_defense_systems.py:
        system_id, genome_id, contig_id, type, subtype, protein_ids, genes_count,
        profile_names, score, hit_id

    genes_df provides one row per gene in a validated system:
        protein_id, system_id, accession, score, presence

    The system_id format is: {genome_id}_{model_name}_{N} where N is the
    instance counter for that model+genome combination.
    """
    if not hits:
        return pd.DataFrame(), pd.DataFrame()

    # Generate unique system IDs
    # Count instances per (genome_id, model_name)
    instance_counter: dict[tuple[str, str], int] = defaultdict(int)

    systems_rows: list[dict] = []
    genes_rows: list[dict] = []

    for hit in hits:
        if len(hit.protein_ids) != len(set(hit.protein_ids)):
            raise ValueError(
                f"System {hit.model_name!r} on {hit.genome_id!r} contains duplicate protein members"
            )
        key = (hit.genome_id, hit.model_name)
        instance_counter[key] += 1
        instance_num = instance_counter[key]

        system_id = f"{hit.genome_id}_{hit.model_name}_{instance_num}"

        systems_rows.append(
            {
                "system_id": system_id,
                "genome_id": hit.genome_id,
                "contig_id": hit.contig_id,
                "type": hit.model_name,
                "subtype": hit.model_name,
                "protein_ids": ",".join(hit.protein_ids),
                "genes_count": len(hit.protein_ids),
                "profile_names": ",".join(hit.profile_names),
                "score": hit.score,
                "hit_id": system_id,
            }
        )

        for i, pid in enumerate(hit.protein_ids):
            genes_rows.append(
                {
                    "protein_id": pid,
                    "system_id": system_id,
                    "accession": hit.profile_names[i],
                    "score": hit.scores[i],
                    "presence": hit.gene_presences[i],
                    "type": hit.model_name,
                    "subtype": hit.model_name,
                    "genome_id": hit.genome_id,
                    "contig_id": hit.contig_id,
                }
            )

    systems_df = pd.DataFrame(systems_rows)
    genes_df = pd.DataFrame(genes_rows)

    return systems_df, genes_df


def _assert_system_dataframe_invariants(
    conn: duckdb.DuckDBPyConnection,
    systems_df: pd.DataFrame,
    genes_df: pd.DataFrame,
) -> None:
    """Fail closed before persisting malformed structured-system calls."""
    if systems_df.empty:
        if not genes_df.empty:
            raise ValueError("Gene assignments exist without structured-system rows")
        return
    if genes_df.empty:
        raise ValueError("Structured-system rows exist without gene assignments")

    required = {
        "system_id",
        "genome_id",
        "contig_id",
        "protein_ids",
        "profile_names",
        "genes_count",
    }
    missing_columns = required - set(systems_df.columns)
    if missing_columns:
        raise ValueError(
            f"Structured-system output is missing required columns: {sorted(missing_columns)}"
        )
    if systems_df["system_id"].duplicated().any():
        raise ValueError("Structured-system output contains duplicate system_id values")

    expected_rows: list[dict[str, str]] = []
    for _, row in systems_df.iterrows():
        system_id = str(row["system_id"])
        genome_id = str(row["genome_id"])
        contig_id = str(row["contig_id"])
        protein_ids = [value for value in str(row["protein_ids"]).split(",") if value]
        profile_names = [value for value in str(row["profile_names"]).split(",") if value]
        genes_count = int(row["genes_count"])

        if not contig_id:
            raise ValueError(f"Structured system {system_id!r} has no contig_id")
        if not protein_ids:
            raise ValueError(f"Structured system {system_id!r} has no protein members")
        if len(protein_ids) != len(set(protein_ids)):
            raise ValueError(f"Structured system {system_id!r} repeats a protein member")
        if len(protein_ids) != len(profile_names) or genes_count != len(protein_ids):
            raise ValueError(f"Structured system {system_id!r} has inconsistent member arrays")

        expected_rows.extend(
            {
                "system_id": system_id,
                "protein_id": protein_id,
                "accession": accession,
                "expected_genome_id": genome_id,
                "expected_contig_id": contig_id,
            }
            for protein_id, accession in zip(
                protein_ids, profile_names, strict=True
            )
        )

    expected_df = pd.DataFrame(expected_rows)
    required_gene_columns = {"system_id", "protein_id", "accession"}
    missing_gene_columns = required_gene_columns - set(genes_df.columns)
    if missing_gene_columns:
        raise ValueError(
            "Structured-system gene output is missing required columns: "
            f"{sorted(missing_gene_columns)}"
        )
    if genes_df.duplicated(["system_id", "protein_id"]).any():
        raise ValueError("Structured-system gene output repeats a system/protein pair")

    expected_members = {
        (str(row.system_id), str(row.protein_id), str(row.accession))
        for row in expected_df.itertuples()
    }
    observed_members = {
        (str(row.system_id), str(row.protein_id), str(row.accession))
        for row in genes_df.itertuples()
    }
    if expected_members != observed_members:
        missing = sorted(expected_members - observed_members)[:10]
        extra = sorted(observed_members - expected_members)[:10]
        raise ValueError(
            "Structured-system rows and gene assignments disagree; "
            f"missing={missing}, extra={extra}"
        )

    conn.register("tmp_expected_system_members", expected_df)
    try:
        invalid = conn.execute(
            """
            SELECT
                e.system_id,
                e.protein_id,
                e.expected_genome_id,
                p.bin_id AS observed_genome_id,
                e.expected_contig_id,
                p.contig_id AS observed_contig_id
            FROM tmp_expected_system_members e
            LEFT JOIN proteins p ON p.protein_id = e.protein_id
            WHERE p.protein_id IS NULL
               OR p.bin_id != e.expected_genome_id
               OR p.contig_id != e.expected_contig_id
            LIMIT 10
            """
        ).fetchall()
    finally:
        conn.unregister("tmp_expected_system_members")

    if invalid:
        raise ValueError(
            "Structured-system membership failed protein/genome/contig "
            f"validation; examples: {invalid}"
        )


# ------------------------------------------------------------------ #
# DB integration (mirrors the integrate_results in validate_defense_systems.py)
# ------------------------------------------------------------------ #


def integrate_defense_results(
    db_path: str | Path,
    systems_df: pd.DataFrame,
    genes_df: pd.DataFrame,
) -> set[str]:
    """Load defense system validation results into the Sharur database.

    Creates/replaces the defense_systems table and inserts
    source='defensefinder_system' annotations.
    """
    return _integrate_results(
        db_path=db_path,
        systems_df=systems_df,
        genes_df=genes_df,
        table_name="defense_systems",
        source_label="defensefinder_system",
        raw_source="defensefinder",
        activity="Defense",
    )


def integrate_secretion_results(
    db_path: str | Path,
    systems_df: pd.DataFrame,
    genes_df: pd.DataFrame,
) -> set[str]:
    """Load secretion system validation results into the Sharur database.

    Creates/replaces the secretion_systems table and inserts
    source='txsscan_system' annotations.
    """
    return _integrate_results(
        db_path=db_path,
        systems_df=systems_df,
        genes_df=genes_df,
        table_name="secretion_systems",
        source_label="txsscan_system",
        raw_source="txsscan",
        activity=None,
    )


def _integrate_results(
    db_path: str | Path,
    systems_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    table_name: str,
    source_label: str,
    raw_source: str,
    activity: str | None,
) -> set[str]:
    """Generic integration: create/replace system table + annotations.

    This mirrors the integrate_results() functions in the existing validation
    scripts, producing the same table schemas for backward compatibility.
    """
    db_path = Path(db_path)
    conn = duckdb.connect(str(db_path))

    _assert_system_dataframe_invariants(conn, systems_df, genes_df)
    existing_tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    capture_sources = {source_label}
    capture_tables = {table_name}
    if (
        "schema_version" in existing_tables
        and get_current_version(conn) < 5
    ):
        # Migration 5 quarantines both legacy callers. Include every protein it
        # invalidates in the returned refresh set, even when this integration
        # invocation is rebuilding only one source.
        capture_sources.update({"defensefinder_system", "txsscan_system"})
        capture_tables.update({"defense_systems", "secretion_systems"})

    capture_source_values = sorted(capture_sources)
    capture_placeholders = ", ".join("?" for _ in capture_source_values)
    affected_protein_ids: set[str] = set()
    if "annotations" in existing_tables:
        affected_protein_ids.update(
            str(row[0])
            for row in conn.execute(
                "SELECT DISTINCT protein_id FROM annotations "
                f"WHERE source IN ({capture_placeholders})",
                capture_source_values,
            ).fetchall()
        )
    if "system_proteins" in existing_tables:
        affected_protein_ids.update(
            str(row[0])
            for row in conn.execute(
                f"""
                SELECT DISTINCT protein_id
                FROM system_proteins
                WHERE system_source IN ({capture_placeholders})
                """,
                capture_source_values,
            ).fetchall()
        )
    for capture_table in sorted(capture_tables & existing_tables):
        for (protein_text,) in conn.execute(
            f"SELECT protein_ids FROM {capture_table} "
            "WHERE protein_ids IS NOT NULL AND protein_ids != ''"
        ).fetchall():
            affected_protein_ids.update(
                cleaned
                for value in re.split(r"[,;\s]+", str(protein_text).strip("[](){}"))
                if (cleaned := value.strip().strip("'\""))
            )
    if "semantic_atoms" in existing_tables:
        affected_protein_ids.update(
            str(row[0])
            for row in conn.execute(
                "SELECT DISTINCT protein_id FROM semantic_atoms "
                f"WHERE source_db IN ({capture_placeholders})",
                capture_source_values,
            ).fetchall()
        )
    if "semantic_terms" in existing_tables:
        affected_protein_ids.update(
            str(row[0])
            for row in conn.execute(
                "SELECT DISTINCT protein_id FROM semantic_terms "
                f"WHERE source_db IN ({capture_placeholders})",
                capture_source_values,
            ).fetchall()
        )
    if not genes_df.empty and "protein_id" in genes_df:
        affected_protein_ids.update(str(value) for value in genes_df["protein_id"])

    if "schema_version" in existing_tables:
        run_migrations(conn)

    conn.execute("BEGIN TRANSACTION")

    tables = [r[0] for r in conn.execute("SHOW TABLES").fetchall()]

    # Create or migrate the systems table.
    before = 0
    if table_name == "defense_systems":
        if table_name not in tables:
            conn.execute("""
                CREATE TABLE defense_systems (
                    system_id VARCHAR PRIMARY KEY,
                    genome_id VARCHAR,
                    contig_id VARCHAR,
                    system_type VARCHAR,
                    system_subtype VARCHAR,
                    activity VARCHAR,
                    genes_count INTEGER,
                    protein_ids VARCHAR,
                    profile_names VARCHAR,
                    sys_beg VARCHAR,
                    sys_end VARCHAR,
                    created_at TIMESTAMP
                )
            """)
            logger.info(f"Created {table_name} table")
        else:
            columns = {row[0] for row in conn.execute("DESCRIBE defense_systems").fetchall()}
            if "contig_id" not in columns:
                conn.execute("ALTER TABLE defense_systems ADD COLUMN contig_id VARCHAR")
            before = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
    elif table_name == "secretion_systems":
        if table_name not in tables:
            conn.execute("""
                CREATE TABLE secretion_systems (
                    system_id VARCHAR PRIMARY KEY,
                    genome_id VARCHAR,
                    contig_id VARCHAR,
                    system_type VARCHAR,
                    system_subtype VARCHAR,
                    genes_count INTEGER,
                    protein_ids VARCHAR,
                    profile_names VARCHAR,
                    sys_beg VARCHAR,
                    sys_end VARCHAR,
                    created_at TIMESTAMP
                )
            """)
            logger.info(f"Created {table_name} table")
        else:
            columns = {row[0] for row in conn.execute("DESCRIBE secretion_systems").fetchall()}
            if "contig_id" not in columns:
                conn.execute("ALTER TABLE secretion_systems ADD COLUMN contig_id VARCHAR")
            before = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]

    # Insert systems
    if not systems_df.empty:
        now = datetime.now(timezone.utc)
        rows: list[dict] = []
        for _, row in systems_df.iterrows():
            entry: dict = {
                "system_id": row["system_id"],
                "genome_id": row.get("genome_id", ""),
                "contig_id": row.get("contig_id", ""),
                "system_type": row.get("type", ""),
                "system_subtype": row.get("subtype", ""),
                "genes_count": int(row.get("genes_count", 0)),
                "protein_ids": row.get("protein_ids", ""),
                "profile_names": row.get("profile_names", ""),
                "sys_beg": "",
                "sys_end": "",
                "created_at": now,
            }
            if activity is not None:
                entry["activity"] = activity
            rows.append(entry)

        insert_df = pd.DataFrame(rows)

        if table_name == "defense_systems":
            cols = [
                "system_id",
                "genome_id",
                "contig_id",
                "system_type",
                "system_subtype",
                "activity",
                "genes_count",
                "protein_ids",
                "profile_names",
                "sys_beg",
                "sys_end",
                "created_at",
            ]
        else:
            cols = [
                "system_id",
                "genome_id",
                "contig_id",
                "system_type",
                "system_subtype",
                "genes_count",
                "protein_ids",
                "profile_names",
                "sys_beg",
                "sys_end",
                "created_at",
            ]
        insert_df = insert_df.reindex(columns=cols, fill_value="")
        conn.register("tmp_systems", insert_df)
        column_sql = ", ".join(cols)
        conn.execute(
            f"INSERT OR REPLACE INTO {table_name} ({column_sql}) "
            f"SELECT {column_sql} FROM tmp_systems"
        )
        conn.execute(
            f"DELETE FROM {table_name} "
            "WHERE system_id NOT IN (SELECT system_id FROM tmp_systems)"
        )
        conn.unregister("tmp_systems")
        logger.info(
            "Replaced %d existing %s rows with %d validated rows",
            before,
            table_name,
            len(insert_df),
        )
    else:
        conn.execute(f"DELETE FROM {table_name}")
        logger.info(f"Cleared {before} existing {table_name} rows")

    # Keep normalized membership synchronized with the authoritative summary
    # table. V2 and downstream operators prefer this join table when present.
    conn.execute("""
        CREATE TABLE IF NOT EXISTS system_proteins (
            system_id VARCHAR NOT NULL,
            protein_id VARCHAR NOT NULL,
            system_source VARCHAR NOT NULL,
            position INTEGER,
            profile_name VARCHAR,
            score DOUBLE,
            PRIMARY KEY (system_id, protein_id, system_source)
        )
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_protein "
        "ON system_proteins(protein_id)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_system "
        "ON system_proteins(system_id)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_source "
        "ON system_proteins(system_source)"
    )
    if not genes_df.empty:
        positions: dict[str, int] = defaultdict(int)
        system_scores = {
            str(row.system_id): float(row.genes_count)
            for row in systems_df.itertuples()
        }
        member_rows: list[dict] = []
        for row in genes_df.itertuples():
            system_id = str(row.system_id)
            positions[system_id] += 1
            member_rows.append(
                {
                    "system_id": system_id,
                    "protein_id": str(row.protein_id),
                    "system_source": source_label,
                    "position": positions[system_id],
                    "profile_name": str(row.accession),
                    "score": system_scores[system_id],
                }
            )
        member_df = pd.DataFrame(member_rows)
        conn.register("tmp_system_proteins", member_df)
        conn.execute("""
            INSERT OR REPLACE INTO system_proteins (
                system_id, protein_id, system_source, position, profile_name, score
            )
            SELECT system_id, protein_id, system_source, position, profile_name, score
            FROM tmp_system_proteins
        """)
        conn.execute(
            """
            DELETE FROM system_proteins
            WHERE system_source = ?
              AND (system_id, protein_id) NOT IN (
                  SELECT system_id, protein_id
                  FROM tmp_system_proteins
              )
            """,
            [source_label],
        )
        conn.unregister("tmp_system_proteins")
    else:
        conn.execute(
            "DELETE FROM system_proteins WHERE system_source = ?",
            [source_label],
        )

    # Replace system-validated annotations even when the new call set is empty.
    next_id = conn.execute(
        "SELECT COALESCE(MAX(annotation_id), 0) FROM annotations"
    ).fetchone()[0]
    existing = conn.execute(
        f"SELECT COUNT(*) FROM annotations WHERE source = '{source_label}'"
    ).fetchone()[0]
    if existing > 0:
        conn.execute(f"DELETE FROM annotations WHERE source = '{source_label}'")
        logger.info(f"Cleared {existing} existing {source_label} annotations")

    if not genes_df.empty:
        ann_rows: list[dict] = []
        for _, row in genes_df.iterrows():
            next_id += 1
            ann_rows.append(
                {
                    "annotation_id": next_id,
                    "protein_id": str(row.get("protein_id", "")),
                    "source": source_label,
                    "accession": row.get("accession", ""),
                    "name": f"{row.get('type', '')}/{row.get('subtype', '')}",
                    "description": f"System: {row.get('system_id', '')}",
                    "evalue": None,
                    "score": row.get("score"),
                    "start_aa": None,
                    "end_aa": None,
                }
            )

        if ann_rows:
            ann_df = pd.DataFrame(ann_rows)
            keep_cols = [
                "annotation_id",
                "protein_id",
                "source",
                "accession",
                "name",
                "description",
                "evalue",
                "score",
                "start_aa",
                "end_aa",
            ]
            ann_df = ann_df.reindex(columns=keep_cols, fill_value=None)
            conn.register("tmp_ann", ann_df)
            conn.execute("""
                INSERT INTO annotations (annotation_id, protein_id, source,
                                         accession, name, description,
                                         evalue, score, start_aa, end_aa)
                SELECT * FROM tmp_ann
            """)
            conn.unregister("tmp_ann")
            logger.info(f"Inserted {len(ann_df)} {source_label} annotations")

    # Structured caller replacement invalidates every derived semantic surface
    # for both former and current members. Raw atoms for other sources remain,
    # but state, search terms, and the legacy compatibility row must be rebuilt
    # before they can be trusted.
    current_tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    if "semantic_atoms" in current_tables:
        conn.execute(
            "DELETE FROM semantic_atoms WHERE source_db = ?",
            [source_label],
        )
    if affected_protein_ids:
        affected_df = pd.DataFrame(
            {"protein_id": sorted(affected_protein_ids)}
        )
        conn.register("tmp_affected_system_proteins", affected_df)
        for semantic_table in (
            "semantic_state",
            "semantic_terms",
            "protein_predicates",
        ):
            if semantic_table in current_tables:
                conn.execute(
                    f"""
                    DELETE FROM {semantic_table}
                    WHERE protein_id IN (
                        SELECT protein_id
                        FROM tmp_affected_system_proteins
                    )
                    """
                )
        conn.unregister("tmp_affected_system_proteins")
        logger.info(
            "Invalidated derived semantic rows for %d former/current %s members",
            len(affected_protein_ids),
            source_label,
        )

    # Summary
    astra_count = conn.execute(
        f"SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = '{raw_source}'"
    ).fetchone()[0]
    system_count = conn.execute(
        f"SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = '{source_label}'"
    ).fetchone()[0]
    logger.info("\nComparison:")
    logger.info(f"  Astra HMM-only (source='{raw_source}'): {astra_count} unique proteins")
    logger.info(f"  System-validated (source='{source_label}'): {system_count} unique proteins")
    if astra_count > 0:
        fp_rate = (astra_count - system_count) / astra_count * 100
        logger.info(f"  FP reduction: {astra_count - system_count} proteins ({fp_rate:.1f}%)")

    conn.commit()
    # DuckDB 1.1 cannot create an index with outstanding updates and does not
    # update a newly added indexed column during INSERT OR REPLACE. Build this
    # nonessential lookup index only after the replacement transaction commits.
    conn.execute(
        f"CREATE INDEX IF NOT EXISTS idx_{table_name}_replicon "
        f"ON {table_name}(genome_id, contig_id)"
    )
    conn.close()
    return affected_protein_ids
