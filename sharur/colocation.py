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
  5. Resolve conflicts (greedy score-ranked non-overlapping selection)

Output: (systems_df, genes_df) DataFrames matching the schema produced by
validate_defense_systems.py and validate_secretion_systems.py.
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Data classes
# ------------------------------------------------------------------ #

@dataclass
class GeneSpec:
    """One gene slot in a system model definition."""

    name: str                     # HMM profile name (= annotations.accession)
    presence: str                 # mandatory | accessory | neutral | forbidden
    loner: bool = False
    multi_system: bool = False
    multi_model: bool = False
    exchangeables: list[str] = field(default_factory=list)


@dataclass
class SystemModel:
    """A system definition parsed from a MacSyFinder XML model file."""

    name: str                     # e.g. "Gabija", "T3SS"
    family: str                   # e.g. "defense-finder-models", "TXSScan"
    inter_gene_max_space: int
    min_mandatory_genes_required: int
    min_genes_required: int
    multi_loci: bool
    genes: list[GeneSpec]

    # Derived lookup tables — built by _build_lookups() after parsing
    n_mandatory_genes: int = 0        # total mandatory gene slots in model
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
        self.n_mandatory_genes = sum(
            1 for g in self.genes if g.presence == "mandatory"
        )
        for gene in self.genes:
            # The parent profile itself
            all_names = [gene.name] + gene.exchangeables
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
    contig_ids: list[str]          # may span contigs if multi_loci
    protein_ids: list[str]
    profile_names: list[str]
    scores: list[float]
    gene_presences: list[str]      # per-gene: mandatory/accessory/neutral/loner
    mandatory_satisfied: int
    total_satisfied: int           # mandatory + accessory (not neutral)
    score: float                   # sum of HMM scores
    primary_mandatory_hits: int = 0  # Bug 2 fix: mandatory slots hit by their PRIMARY profile
    n_mandatory_total: int = 0       # Bug 2 fix: total mandatory slots in the model

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


# ------------------------------------------------------------------ #
# XML model parser
# ------------------------------------------------------------------ #

def _parse_bool(val: Optional[str]) -> bool:
    """Parse XML boolean attribute (1, True, true → True; everything else → False)."""
    if val is None:
        return False
    return val.strip().lower() in ("1", "true")


def _parse_model_xml(xml_path: Path, family: str) -> Optional[SystemModel]:
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
    inter_gene_max_space = int(root.get("inter_gene_max_space", "5"))
    min_mandatory = int(root.get("min_mandatory_genes_required", "1"))
    min_genes = int(root.get("min_genes_required", "1"))
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

        # Parse exchangeables
        exchangeables: list[str] = []
        exch_elem = gene_elem.find("exchangeables")
        if exch_elem is not None:
            for exch_gene in exch_elem.findall("gene"):
                exch_name = exch_gene.get("name", "")
                if exch_name:
                    exchangeables.append(exch_name)

        genes.append(GeneSpec(
            name=name,
            presence=presence,
            loner=loner,
            multi_system=multi_system,
            multi_model=multi_model,
            exchangeables=exchangeables,
        ))

    model = SystemModel(
        name=model_name,
        family=family,
        inter_gene_max_space=inter_gene_max_space,
        min_mandatory_genes_required=min_mandatory,
        min_genes_required=min_genes,
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


def parse_models_for_family(
    models_dir: Path, model_family: str
) -> dict[str, SystemModel]:
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


def load_hits(
    conn: duckdb.DuckDBPyConnection, source: str
) -> list[HitRecord]:
    """Load HMM hits from the DuckDB annotations + proteins tables.

    Args:
        conn: Open DuckDB connection (read-only is fine).
        source: Annotation source filter (e.g. 'defensefinder', 'txsscan').

    Returns:
        List of HitRecord, sorted by contig_id then gene_index.
    """
    rows = conn.execute("""
        SELECT a.protein_id, a.accession, a.score,
               p.contig_id, p.bin_id, p.gene_index
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE a.source = ?
        ORDER BY p.contig_id, p.gene_index
    """, [source]).fetchall()

    return [
        HitRecord(
            protein_id=r[0],
            accession=r[1],
            score=float(r[2]) if r[2] is not None else 0.0,
            contig_id=r[3],
            bin_id=r[4],
            gene_index=int(r[5]) if r[5] is not None else 0,
        )
        for r in rows
    ]


# ------------------------------------------------------------------ #
# Indexing helpers
# ------------------------------------------------------------------ #

def _index_hits_by_contig(
    hits: list[HitRecord],
) -> dict[str, list[HitRecord]]:
    """Group hits by contig_id."""
    by_contig: dict[str, list[HitRecord]] = defaultdict(list)
    for h in hits:
        by_contig[h.contig_id].append(h)
    return dict(by_contig)


def _index_hits_by_genome(
    hits: list[HitRecord],
) -> dict[str, list[HitRecord]]:
    """Group hits by bin_id (genome)."""
    by_genome: dict[str, list[HitRecord]] = defaultdict(list)
    for h in hits:
        by_genome[h.bin_id].append(h)
    return dict(by_genome)


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

    Returns:
        Dict mapping internal_name -> filename_stem (only for mismatches).
    """
    mapping: dict[str, str] = {}
    if not profiles_dir.is_dir():
        logger.warning(f"Profiles directory not found: {profiles_dir}")
        return mapping

    name_re = re.compile(r"^NAME\s+(\S+)")
    for hmm_path in profiles_dir.glob("*.hmm"):
        stem = hmm_path.stem
        try:
            with open(hmm_path, "r") as f:
                for line in f:
                    m = name_re.match(line)
                    if m:
                        internal_name = m.group(1)
                        if internal_name != stem:
                            mapping[internal_name] = stem
                        break
        except (OSError, UnicodeDecodeError):
            continue

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


def _find_profiles_dir(family_dir: Path) -> Optional[Path]:
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
# Forbidden proximity check (Bug 4 fix)
# ------------------------------------------------------------------ #

def _cluster_has_nearby_forbidden(
    cluster: list[HitRecord],
    forbidden_by_contig: dict[str, list[HitRecord]],
    max_gap: int,
) -> bool:
    """Check if any forbidden-profile hit is within max_gap of a cluster gene.

    Args:
        cluster: HitRecords forming a candidate cluster (regular hits only).
        forbidden_by_contig: Forbidden-profile hits indexed by contig_id.
        max_gap: Maximum gene-index distance to consider "nearby".

    Returns:
        True if a forbidden hit is within range (cluster should be rejected).
    """
    for hit in cluster:
        for fh in forbidden_by_contig.get(hit.contig_id, []):
            if abs(hit.gene_index - fh.gene_index) <= max_gap:
                return True
    return False


# ------------------------------------------------------------------ #
# Clustering
# ------------------------------------------------------------------ #

def _cluster_contig(
    hits: list[HitRecord], max_gap: int
) -> list[list[HitRecord]]:
    """Group hits on the same contig within max_gap gene positions.

    Hits must all be from the same contig. They are sorted by gene_index
    and grouped greedily: a new cluster starts whenever the gap between
    consecutive hits exceeds max_gap.

    Args:
        hits: HitRecords from a single contig.
        max_gap: Maximum gene-index distance between consecutive hits
                 in the same cluster.

    Returns:
        List of clusters, each a list of HitRecords.
    """
    if not hits:
        return []

    sorted_hits = sorted(hits, key=lambda h: h.gene_index)
    clusters: list[list[HitRecord]] = [[sorted_hits[0]]]

    for hit in sorted_hits[1:]:
        if hit.gene_index - clusters[-1][-1].gene_index <= max_gap:
            clusters[-1].append(hit)
        else:
            clusters.append([hit])

    return clusters


# ------------------------------------------------------------------ #
# Validation
# ------------------------------------------------------------------ #

def _validate_cluster(
    cluster: list[HitRecord], model: SystemModel
) -> Optional[SystemHit]:
    """Check if a cluster of hits satisfies a model's quorum constraints.

    Performs:
      1. Map each hit's accession to a gene slot via the model's lookup tables
      2. Reject if any forbidden profile is present
      3. Count satisfied mandatory and accessory gene slots (by canonical name)
      4. Check quorum: min_mandatory_genes_required, min_genes_required
      5. Verify no internal gap exceeds inter_gene_max_space

    Args:
        cluster: HitRecords from one contig (or merged across contigs for multi_loci).
        model: The SystemModel to validate against.

    Returns:
        A SystemHit if the cluster is valid, None otherwise.
    """
    # 1. Classify hits by gene slot
    satisfied_mandatory: set[str] = set()   # canonical gene slot names
    satisfied_accessory: set[str] = set()
    satisfied_neutral: set[str] = set()
    primary_mandatory_slots: set[str] = set()  # Bug 2: slots hit by their PRIMARY profile
    multi_system_pids: set[str] = set()

    # Track the actual hits that map to this model (exclude irrelevant ones)
    relevant_hits: list[HitRecord] = []

    for hit in cluster:
        gene_slot = model.profile_to_gene_slot.get(hit.accession)
        if gene_slot is None:
            continue  # Hit doesn't belong to this model

        # 2. Forbidden check
        if gene_slot.presence == "forbidden":
            return None

        relevant_hits.append(hit)
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

    if not relevant_hits:
        return None

    # 3. Check quorum
    n_mandatory = len(satisfied_mandatory)
    n_total = n_mandatory + len(satisfied_accessory)

    if n_mandatory < model.min_mandatory_genes_required:
        return None
    if n_total < model.min_genes_required:
        return None

    # 4. Verify gap constraint on relevant hits only
    #    Sort by gene_index and check consecutive gaps
    relevant_sorted = sorted(relevant_hits, key=lambda h: (h.contig_id, h.gene_index))

    # For single-contig clusters, check internal gaps
    if not model.multi_loci:
        gene_indices = [h.gene_index for h in relevant_sorted]
        for i in range(1, len(gene_indices)):
            if gene_indices[i] - gene_indices[i - 1] > model.inter_gene_max_space:
                # Gap too large — try to split and validate sub-clusters
                # For simplicity, reject this cluster; the per-model clustering
                # in the main loop should have already produced correct clusters
                return None

    # 5. Build the SystemHit
    # Determine genome_id from the hits
    genome_id = relevant_hits[0].bin_id

    # Determine presence label for each gene
    presences: list[str] = []
    for hit in relevant_sorted:
        gene_slot = model.profile_to_gene_slot.get(hit.accession)
        if gene_slot:
            presences.append(gene_slot.presence)
        else:
            presences.append("unknown")

    system_hit = SystemHit(
        model_name=model.name,
        model_family=model.family,
        genome_id=genome_id,
        contig_ids=list(dict.fromkeys(h.contig_id for h in relevant_sorted)),
        protein_ids=[h.protein_id for h in relevant_sorted],
        profile_names=[h.accession for h in relevant_sorted],
        scores=[h.score for h in relevant_sorted],
        gene_presences=presences,
        mandatory_satisfied=n_mandatory,
        total_satisfied=n_total,
        score=sum(h.score for h in relevant_sorted),
        primary_mandatory_hits=len(primary_mandatory_slots),
        n_mandatory_total=model.n_mandatory_genes,
    )
    system_hit.multi_system_proteins = multi_system_pids
    return system_hit


# ------------------------------------------------------------------ #
# Loner support (TXSScan only)
# ------------------------------------------------------------------ #

def _add_loner_genes(
    hit: SystemHit,
    model: SystemModel,
    genome_hits: list[HitRecord],
) -> SystemHit:
    """Add loner-eligible genes from anywhere in the genome to a system hit.

    Loner genes are those marked loner=True in the model. They can exist
    outside the inter_gene_max_space of the main cluster and still count
    toward quorum.

    Only adds loner genes that are NOT already in the cluster and whose
    profile matches the model's loner-eligible profiles.

    Args:
        hit: An existing validated SystemHit (the core cluster).
        model: The SystemModel being validated.
        genome_hits: All HMM hits in this genome.

    Returns:
        A new SystemHit with loner genes appended (or the same hit if
        no loners found).
    """
    if not model.loner_profiles:
        return hit

    existing_pids = hit.protein_set
    # Track which mandatory/accessory slots are already satisfied
    satisfied_mandatory: set[str] = set()
    satisfied_accessory: set[str] = set()
    for profile in hit.profile_names:
        gene_slot = model.profile_to_gene_slot.get(profile)
        if gene_slot:
            if gene_slot.presence == "mandatory":
                satisfied_mandatory.add(gene_slot.name)
            elif gene_slot.presence == "accessory":
                satisfied_accessory.add(gene_slot.name)

    loner_additions: list[HitRecord] = []
    loner_presences: list[str] = []
    new_multi_system: set[str] = set(hit.multi_system_proteins)
    # Bug 2: track primary mandatory hits (propagate from parent + add loner primaries)
    primary_mandatory_slots: set[str] = set()
    for profile in hit.profile_names:
        gene_slot = model.profile_to_gene_slot.get(profile)
        if gene_slot and gene_slot.presence == "mandatory" and profile == gene_slot.name:
            primary_mandatory_slots.add(gene_slot.name)

    for h in genome_hits:
        if h.protein_id in existing_pids:
            continue
        if h.accession not in model.loner_profiles:
            continue

        gene_slot = model.profile_to_gene_slot.get(h.accession)
        if gene_slot is None:
            continue

        # Check this hit is actually outside the cluster's reach
        # (i.e., it's not within inter_gene_max_space of any cluster gene
        #  on the same contig)
        is_near_cluster = False
        if h.contig_id in hit.contig_ids:
            for pid, existing_profile in zip(hit.protein_ids, hit.profile_names):
                # Find the gene_index of the existing hit — we need to look it up
                # from the genome_hits list
                for gh in genome_hits:
                    if gh.protein_id == pid and gh.contig_id == h.contig_id:
                        if abs(h.gene_index - gh.gene_index) <= model.inter_gene_max_space:
                            is_near_cluster = True
                            break
                if is_near_cluster:
                    break

        if is_near_cluster:
            # This gene is close enough to be in the cluster already;
            # it should have been picked up by clustering. Skip to avoid dupes.
            continue

        slot_name = gene_slot.name
        # Only add if this slot isn't already satisfied
        if gene_slot.presence == "mandatory" and slot_name not in satisfied_mandatory:
            loner_additions.append(h)
            loner_presences.append("loner")
            satisfied_mandatory.add(slot_name)
            if h.accession == gene_slot.name:
                primary_mandatory_slots.add(slot_name)
            if gene_slot.multi_system:
                new_multi_system.add(h.protein_id)
        elif gene_slot.presence == "accessory" and slot_name not in satisfied_accessory:
            loner_additions.append(h)
            loner_presences.append("loner")
            satisfied_accessory.add(slot_name)
            if gene_slot.multi_system:
                new_multi_system.add(h.protein_id)

    if not loner_additions:
        return hit

    # Build a new SystemHit with the loner genes included
    new_hit = SystemHit(
        model_name=hit.model_name,
        model_family=hit.model_family,
        genome_id=hit.genome_id,
        contig_ids=list(dict.fromkeys(
            hit.contig_ids + [h.contig_id for h in loner_additions]
        )),
        protein_ids=hit.protein_ids + [h.protein_id for h in loner_additions],
        profile_names=hit.profile_names + [h.accession for h in loner_additions],
        scores=hit.scores + [h.score for h in loner_additions],
        gene_presences=hit.gene_presences + loner_presences,
        mandatory_satisfied=len(satisfied_mandatory),
        total_satisfied=len(satisfied_mandatory) + len(satisfied_accessory),
        score=hit.score + sum(h.score for h in loner_additions),
        primary_mandatory_hits=len(primary_mandatory_slots),
        n_mandatory_total=hit.n_mandatory_total,
    )
    new_hit.multi_system_proteins = new_multi_system
    return new_hit


# ------------------------------------------------------------------ #
# Multi-loci support (TXSScan only)
# ------------------------------------------------------------------ #

def _merge_multi_loci(
    clusters: list[SystemHit], model: SystemModel
) -> Optional[SystemHit]:
    """Merge validated clusters across contigs for multi_loci models.

    For models with multi_loci=True, genes from the same system can be
    spread across multiple contigs. This function merges all valid clusters
    for a model within a genome into a single system hit, then re-checks
    quorum on the merged set.

    Args:
        clusters: All validated SystemHits for this model in one genome.
        model: The SystemModel being validated.

    Returns:
        A merged SystemHit if quorum is met, None otherwise.
    """
    if not clusters:
        return None
    if len(clusters) == 1:
        return clusters[0]

    # Merge all genes from all clusters
    all_pids: list[str] = []
    all_profiles: list[str] = []
    all_scores: list[float] = []
    all_presences: list[str] = []
    all_contigs: list[str] = []
    multi_sys_pids: set[str] = set()
    seen_pids: set[str] = set()

    for c in clusters:
        for i, pid in enumerate(c.protein_ids):
            if pid not in seen_pids:
                seen_pids.add(pid)
                all_pids.append(pid)
                all_profiles.append(c.profile_names[i])
                all_scores.append(c.scores[i])
                all_presences.append(c.gene_presences[i])
        all_contigs.extend(c.contig_ids)
        multi_sys_pids |= c.multi_system_proteins

    # Re-check quorum on the merged set
    satisfied_mandatory: set[str] = set()
    satisfied_accessory: set[str] = set()
    primary_mandatory_slots: set[str] = set()
    for profile in all_profiles:
        gene_slot = model.profile_to_gene_slot.get(profile)
        if gene_slot:
            if gene_slot.presence == "mandatory":
                satisfied_mandatory.add(gene_slot.name)
                if profile == gene_slot.name:
                    primary_mandatory_slots.add(gene_slot.name)
            elif gene_slot.presence == "accessory":
                satisfied_accessory.add(gene_slot.name)

    n_mandatory = len(satisfied_mandatory)
    n_total = n_mandatory + len(satisfied_accessory)

    if n_mandatory < model.min_mandatory_genes_required:
        return None
    if n_total < model.min_genes_required:
        return None

    merged = SystemHit(
        model_name=model.name,
        model_family=model.family,
        genome_id=clusters[0].genome_id,
        contig_ids=list(dict.fromkeys(all_contigs)),
        protein_ids=all_pids,
        profile_names=all_profiles,
        scores=all_scores,
        gene_presences=all_presences,
        mandatory_satisfied=n_mandatory,
        total_satisfied=n_total,
        score=sum(all_scores),
        primary_mandatory_hits=len(primary_mandatory_slots),
        n_mandatory_total=model.n_mandatory_genes,
    )
    merged.multi_system_proteins = multi_sys_pids
    return merged


# ------------------------------------------------------------------ #
# Conflict resolution
# ------------------------------------------------------------------ #

def resolve_conflicts(hits: list[SystemHit]) -> list[SystemHit]:
    """Select best non-overlapping system assignments via greedy selection.

    Sorts hits by score descending, then greedily selects hits whose
    non-multi_system proteins don't overlap with already-selected hits.

    Proteins flagged as multi_system in their model are exempt from the
    overlap check — they can be assigned to multiple system instances.

    Args:
        hits: All candidate SystemHits across all models and genomes.

    Returns:
        The selected non-conflicting subset.
    """
    if not hits:
        return []

    # Group by genome for efficiency — conflicts only happen within a genome
    by_genome: dict[str, list[SystemHit]] = defaultdict(list)
    for h in hits:
        by_genome[h.genome_id].append(h)

    selected: list[SystemHit] = []

    for genome_id, genome_hits in by_genome.items():
        # Bug 2 fix: sort by (mandatory_fraction DESC, primary_mandatory_hits DESC, score DESC)
        # A model satisfying 3/3 mandatory (1.0) beats 3/4 (0.75) — this prevents
        # models with lots of exchangeables (e.g. PrrC) from stealing loci that
        # belong to simpler models (e.g. RM_Type_I) when the defining gene is absent.
        genome_hits_sorted = sorted(
            genome_hits,
            key=lambda h: (
                h.mandatory_satisfied / max(h.n_mandatory_total, 1),
                h.primary_mandatory_hits,
                h.score,
            ),
            reverse=True,
        )
        used_proteins: set[str] = set()

        for hit in genome_hits_sorted:
            # Proteins that participate in conflict checking
            # (exclude multi_system proteins — they can be shared)
            conflict_proteins = hit.protein_set - hit.multi_system_proteins

            if conflict_proteins & used_proteins:
                continue  # Conflicts with an already-selected system

            selected.append(hit)
            used_proteins |= conflict_proteins

    return selected


# ------------------------------------------------------------------ #
# Main orchestrator
# ------------------------------------------------------------------ #

def _find_models_dir() -> Optional[Path]:
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
            system_id, genome_id, type, subtype, protein_ids,
            genes_count, profile_names, score, hit_id

        genes_df columns:
            protein_id, system_id, accession, score, presence
    """
    import time as _time
    t0 = _time.time()

    db_path = Path(db_path)

    # Resolve models directory
    if models_dir is None:
        models_dir = _find_models_dir()
        if models_dir is None:
            logger.error(
                "MacSyFinder models directory not found. "
                "Expected at ~/.macsyfinder/models/"
            )
            return pd.DataFrame(), pd.DataFrame()
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
        logger.error(f"No models found for family '{model_family}'")
        return pd.DataFrame(), pd.DataFrame()
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
        logger.warning(f"  Could not find profiles directory for {model_family}")

    # Build indices
    hits_by_contig = _index_hits_by_contig(all_hits)
    hits_by_genome = _index_hits_by_genome(all_hits)
    hits_by_accession = _index_hits_by_accession(all_hits)

    # Build a quick set of all accessions present in the data
    all_accessions_in_data = set(hits_by_accession.keys())

    # Phase 3+4: Cluster and validate per model
    if verbose:
        logger.info("Clustering and validating...")

    candidate_hits: list[SystemHit] = []
    models_with_hits = 0

    for model_name, model in models.items():
        # Quick check: does this model have ANY relevant profiles in the data?
        relevant_profiles = model.all_profiles & all_accessions_in_data
        if not relevant_profiles:
            continue
        models_with_hits += 1

        # Bug 4 fix: split hits into regular (mandatory/accessory/neutral) and
        # forbidden. Cluster only regular hits; check forbidden proximity after.
        regular_hits: list[HitRecord] = []
        forbidden_hits: list[HitRecord] = []
        for profile in relevant_profiles:
            for h in hits_by_accession.get(profile, []):
                if profile in model.forbidden_profiles:
                    forbidden_hits.append(h)
                else:
                    regular_hits.append(h)

        if not regular_hits:
            continue

        # Index forbidden hits by (contig_id) for proximity checks
        forbidden_by_contig: dict[str, list[HitRecord]] = defaultdict(list)
        for h in forbidden_hits:
            forbidden_by_contig[h.contig_id].append(h)

        # Group regular hits by genome, then by contig
        model_by_genome: dict[str, dict[str, list[HitRecord]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for h in regular_hits:
            model_by_genome[h.bin_id][h.contig_id].append(h)

        # Process each genome
        for genome_id, contigs in model_by_genome.items():
            per_genome_clusters: list[SystemHit] = []

            if model.multi_loci:
                # Bug 3 fix: for multi_loci models, collect ALL clusters
                # genome-wide (without quorum validation), merge into one
                # candidate, THEN validate the merged result.
                all_raw_clusters: list[list[HitRecord]] = []
                for contig_id, contig_hits in contigs.items():
                    clusters = _cluster_contig(contig_hits, model.inter_gene_max_space)
                    # Bug 4: reject clusters near forbidden hits
                    for cluster in clusters:
                        if not _cluster_has_nearby_forbidden(
                            cluster, forbidden_by_contig, model.inter_gene_max_space
                        ):
                            all_raw_clusters.append(cluster)

                if all_raw_clusters:
                    # Merge all raw hits into a single candidate cluster
                    merged_hits: list[HitRecord] = []
                    seen_pids: set[str] = set()
                    for cluster in all_raw_clusters:
                        for h in cluster:
                            if h.protein_id not in seen_pids:
                                seen_pids.add(h.protein_id)
                                merged_hits.append(h)

                    # Validate the merged genome-wide cluster
                    system_hit = _validate_cluster(merged_hits, model)
                    if system_hit is not None:
                        per_genome_clusters.append(system_hit)

                    # Also try individual clusters that pass on their own
                    # (for genomes with multiple independent instances)
                    if not per_genome_clusters:
                        for cluster in all_raw_clusters:
                            system_hit = _validate_cluster(cluster, model)
                            if system_hit is not None:
                                per_genome_clusters.append(system_hit)
            else:
                for contig_id, contig_hits in contigs.items():
                    # Cluster hits on this contig for this model's max_gap
                    clusters = _cluster_contig(contig_hits, model.inter_gene_max_space)

                    for cluster in clusters:
                        # Bug 4: skip clusters near forbidden-profile hits
                        if _cluster_has_nearby_forbidden(
                            cluster, forbidden_by_contig, model.inter_gene_max_space
                        ):
                            continue
                        system_hit = _validate_cluster(cluster, model)
                        if system_hit is not None:
                            per_genome_clusters.append(system_hit)

            # Phase 4c: Add loner genes (TXSScan only)
            if model.loner_profiles:
                genome_all_hits = hits_by_genome.get(genome_id, [])
                per_genome_clusters = [
                    _add_loner_genes(hit, model, genome_all_hits)
                    for hit in per_genome_clusters
                ]

                # Re-validate quorum after adding loners (for clusters that
                # previously didn't meet quorum, they were already rejected).
                # Loners can only help, never hurt, so no re-rejection needed.

            candidate_hits.extend(per_genome_clusters)

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

    elapsed = _time.time() - t0
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
    hits: list[SystemHit], source: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert SystemHits into the output DataFrames expected by downstream code.

    systems_df matches the schema from validate_defense_systems.py:
        system_id, genome_id, type, subtype, protein_ids, genes_count,
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

    # Determine the system source annotation label
    system_source = f"{source}_system"

    for hit in hits:
        key = (hit.genome_id, hit.model_name)
        instance_counter[key] += 1
        instance_num = instance_counter[key]

        system_id = f"{hit.genome_id}_{hit.model_name}_{instance_num}"

        systems_rows.append({
            "system_id": system_id,
            "genome_id": hit.genome_id,
            "type": hit.model_name,
            "subtype": hit.model_name,
            "protein_ids": ",".join(hit.protein_ids),
            "genes_count": len(hit.protein_ids),
            "profile_names": ",".join(hit.profile_names),
            "score": hit.score,
            "hit_id": system_id,
        })

        for i, pid in enumerate(hit.protein_ids):
            genes_rows.append({
                "protein_id": pid,
                "system_id": system_id,
                "accession": hit.profile_names[i],
                "score": hit.scores[i],
                "presence": hit.gene_presences[i],
                "type": hit.model_name,
                "subtype": hit.model_name,
                "genome_id": hit.genome_id,
            })

    systems_df = pd.DataFrame(systems_rows)
    genes_df = pd.DataFrame(genes_rows)

    return systems_df, genes_df


# ------------------------------------------------------------------ #
# DB integration (mirrors the integrate_results in validate_defense_systems.py)
# ------------------------------------------------------------------ #

def integrate_defense_results(
    db_path: str | Path,
    systems_df: pd.DataFrame,
    genes_df: pd.DataFrame,
) -> None:
    """Load defense system validation results into the Sharur database.

    Creates/replaces the defense_systems table and inserts
    source='defensefinder_system' annotations.
    """
    _integrate_results(
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
) -> None:
    """Load secretion system validation results into the Sharur database.

    Creates/replaces the secretion_systems table and inserts
    source='txsscan_system' annotations.
    """
    _integrate_results(
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
    activity: Optional[str],
) -> None:
    """Generic integration: create/replace system table + annotations.

    This mirrors the integrate_results() functions in the existing validation
    scripts, producing the same table schemas for backward compatibility.
    """
    from datetime import datetime, timezone

    db_path = Path(db_path)
    conn = duckdb.connect(str(db_path))

    tables = [r[0] for r in conn.execute("SHOW TABLES").fetchall()]

    # Create or clear the systems table
    if table_name == "defense_systems":
        if table_name not in tables:
            conn.execute("""
                CREATE TABLE defense_systems (
                    system_id VARCHAR PRIMARY KEY,
                    genome_id VARCHAR,
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
            before = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            conn.execute(f"DELETE FROM {table_name}")
            logger.info(f"Cleared {before} existing {table_name} rows")
    elif table_name == "secretion_systems":
        if table_name not in tables:
            conn.execute("""
                CREATE TABLE secretion_systems (
                    system_id VARCHAR PRIMARY KEY,
                    genome_id VARCHAR,
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
            before = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            conn.execute(f"DELETE FROM {table_name}")
            logger.info(f"Cleared {before} existing {table_name} rows")

    # Insert systems
    if not systems_df.empty:
        now = datetime.now(timezone.utc)
        rows: list[dict] = []
        for _, row in systems_df.iterrows():
            entry: dict = {
                "system_id": row["system_id"],
                "genome_id": row.get("genome_id", ""),
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
                "system_id", "genome_id", "system_type", "system_subtype",
                "activity", "genes_count", "protein_ids", "profile_names",
                "sys_beg", "sys_end", "created_at",
            ]
        else:
            cols = [
                "system_id", "genome_id", "system_type", "system_subtype",
                "genes_count", "protein_ids", "profile_names",
                "sys_beg", "sys_end", "created_at",
            ]
        insert_df = insert_df.reindex(columns=cols, fill_value="")
        conn.register("tmp_systems", insert_df)
        conn.execute(f"INSERT INTO {table_name} SELECT * FROM tmp_systems")
        conn.unregister("tmp_systems")
        logger.info(f"Inserted {len(insert_df)} {table_name}")

    # Add system-validated annotations
    if not genes_df.empty:
        existing = conn.execute(
            f"SELECT COUNT(*) FROM annotations WHERE source = '{source_label}'"
        ).fetchone()[0]
        if existing > 0:
            conn.execute(
                f"DELETE FROM annotations WHERE source = '{source_label}'"
            )
            logger.info(f"Cleared {existing} existing {source_label} annotations")

        next_id = conn.execute(
            "SELECT COALESCE(MAX(annotation_id), 0) FROM annotations"
        ).fetchone()[0]

        existing_proteins = set(
            r[0]
            for r in conn.execute(
                "SELECT DISTINCT protein_id FROM proteins"
            ).fetchall()
        )

        ann_rows: list[dict] = []
        for _, row in genes_df.iterrows():
            pid = row.get("protein_id", "")
            if pid not in existing_proteins:
                continue
            next_id += 1
            ann_rows.append({
                "annotation_id": next_id,
                "protein_id": pid,
                "source": source_label,
                "accession": row.get("accession", ""),
                "name": f"{row.get('type', '')}/{row.get('subtype', '')}",
                "description": f"System: {row.get('system_id', '')}",
                "evalue": None,
                "score": row.get("score"),
                "start_aa": None,
                "end_aa": None,
            })

        if ann_rows:
            ann_df = pd.DataFrame(ann_rows)
            keep_cols = [
                "annotation_id", "protein_id", "source", "accession", "name",
                "description", "evalue", "score", "start_aa", "end_aa",
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

    # Summary
    astra_count = conn.execute(
        f"SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = '{raw_source}'"
    ).fetchone()[0]
    system_count = conn.execute(
        f"SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = '{source_label}'"
    ).fetchone()[0]
    logger.info(f"\nComparison:")
    logger.info(f"  Astra HMM-only (source='{raw_source}'): {astra_count} unique proteins")
    logger.info(
        f"  System-validated (source='{source_label}'): {system_count} unique proteins"
    )
    if astra_count > 0:
        fp_rate = (astra_count - system_count) / astra_count * 100
        logger.info(
            f"  FP reduction: {astra_count - system_count} proteins ({fp_rate:.1f}%)"
        )

    conn.commit()
    conn.close()
