"""
Predicate V2 data model: typed semantic atoms and aggregated state.

This module defines the core types for the V2 predicate system:
- SemanticFacet: what kind of claim an atom makes (activity, role, etc.)
- ClaimRelation: how strong the evidence-to-claim link is
- SemanticAtom: a single typed claim about a protein
- SemanticState: aggregated per-protein view across all atoms
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class SemanticFacet(str, Enum):
    """What kind of claim a semantic atom makes."""

    activity = "activity"  # What the protein catalyzes (hydrogenase, kinase)
    role = "role"  # Functional role (defense_system, transporter, regulator)
    architecture = "architecture"  # Structural features (multi_domain, beta_barrel)
    localization = "localization"  # Where the protein is (membrane, periplasmic)
    topology = "topology"  # TM helix count, signal peptide, orientation
    size_class = "size_class"  # Size category (tiny, small, medium, large, giant)
    quality_flag = "quality_flag"  # Annotation quality (unannotated, well_annotated)


class ClaimRelation(str, Enum):
    """How confidently source evidence supports an atom.

    These are deterministic based on rules, never LLM-generated.
    """

    implies = "implies"  # Strong: KEGG ortholog, HydDB classification
    supports = "supports"  # Moderate: PFAM domain consistent with function
    flags = "flags"  # Weak/ambiguous: pattern match, DUF, superfamily
    excludes = "excludes"  # Counter-evidence (Complex I excludes hydrogenase)
    unresolved = "unresolved"  # Annotation exists but no curated mapping


@dataclass
class SemanticAtom:
    """A single typed claim about a protein, derived from one annotation.

    One protein can have many atoms. Each atom is a typed claim
    backed by a specific annotation source.
    """

    protein_id: str
    atom_id: str  # The semantic tag, e.g. "hydrogenase", "membrane"
    facet: SemanticFacet  # What kind of claim
    relation: ClaimRelation  # How strong the evidence is
    source_accession: str  # Annotation accession that produced this atom
    source_db: str  # pfam, kegg, cazy, hyddb, defensefinder, vogdb, ...
    evidence_evalue: float | None = None  # Pass-through from HMM/DIAMOND
    evidence_score: float | None = None  # Pass-through bitscore

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary for storage."""
        return {
            "protein_id": self.protein_id,
            "atom_id": self.atom_id,
            "facet": self.facet.value,
            "relation": self.relation.value,
            "source_accession": self.source_accession,
            "source_db": self.source_db,
            "evidence_evalue": self.evidence_evalue,
            "evidence_score": self.evidence_score,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> SemanticAtom:
        """Deserialize from dictionary."""
        return cls(
            protein_id=d["protein_id"],
            atom_id=d["atom_id"],
            facet=SemanticFacet(d["facet"]),
            relation=ClaimRelation(d["relation"]),
            source_accession=d["source_accession"],
            source_db=d["source_db"],
            evidence_evalue=d.get("evidence_evalue"),
            evidence_score=d.get("evidence_score"),
        )


@dataclass
class SemanticState:
    """Per-protein aggregated view after resolving all atoms.

    This is the resolved state per facet — conflicts are handled,
    duplicates removed, and unresolved accessions collected.
    """

    protein_id: str
    activities: list[str] = field(default_factory=list)
    roles: list[str] = field(default_factory=list)
    architecture: list[str] = field(default_factory=list)
    localization: list[str] = field(default_factory=list)
    topology: dict[str, Any] = field(default_factory=dict)
    size_class: str = ""
    quality_flags: list[str] = field(default_factory=list)
    composite_predicates: list[str] = field(default_factory=list)
    unresolved_accessions: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary for storage."""
        return {
            "protein_id": self.protein_id,
            "activities": self.activities,
            "roles": self.roles,
            "architecture": self.architecture,
            "localization": self.localization,
            "topology": self.topology,
            "size_class": self.size_class,
            "quality_flags": self.quality_flags,
            "composite_predicates": self.composite_predicates,
            "unresolved_accessions": self.unresolved_accessions,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> SemanticState:
        """Deserialize from dictionary."""
        return cls(
            protein_id=d["protein_id"],
            activities=d.get("activities", []),
            roles=d.get("roles", []),
            architecture=d.get("architecture", []),
            localization=d.get("localization", []),
            topology=d.get("topology", {}),
            size_class=d.get("size_class", ""),
            quality_flags=d.get("quality_flags", []),
            composite_predicates=d.get("composite_predicates", []),
            unresolved_accessions=d.get("unresolved_accessions", []),
        )


# ---------------------------------------------------------------------------
# DuckDB table creation helper (standalone, does NOT modify schema.py)
# ---------------------------------------------------------------------------

V2_SCHEMA = """
-- Per-atom storage (the raw evidence layer)
CREATE TABLE IF NOT EXISTS semantic_atoms (
    protein_id VARCHAR NOT NULL,
    atom_id VARCHAR NOT NULL,
    facet VARCHAR NOT NULL,
    relation VARCHAR NOT NULL,
    source_accession VARCHAR,
    source_db VARCHAR,
    evidence_evalue DOUBLE,
    evidence_score DOUBLE,
    PRIMARY KEY (protein_id, atom_id, source_accession)
);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_atom ON semantic_atoms(atom_id);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_facet_atom ON semantic_atoms(facet, atom_id);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_source ON semantic_atoms(source_db, source_accession);

-- Per-protein aggregated state (resolved view)
CREATE TABLE IF NOT EXISTS semantic_state (
    protein_id VARCHAR PRIMARY KEY,
    activities VARCHAR[],
    roles VARCHAR[],
    architecture VARCHAR[],
    localization VARCHAR[],
    topology JSON,
    size_class VARCHAR,
    quality_flags VARCHAR[],
    composite_predicates VARCHAR[],
    unresolved_count INTEGER,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
-- Unified materialized search terms
CREATE TABLE IF NOT EXISTS semantic_terms (
    protein_id VARCHAR NOT NULL,
    term_id VARCHAR NOT NULL,
    term_kind VARCHAR NOT NULL,
    facet VARCHAR,
    relation VARCHAR,
    source_db VARCHAR NOT NULL DEFAULT '',
    source_accession VARCHAR NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_semantic_terms_term ON semantic_terms(term_id, protein_id);
CREATE INDEX IF NOT EXISTS idx_semantic_terms_protein ON semantic_terms(protein_id);

-- Constraint-free append targets for resumable full refreshes. Completed
-- generations are promoted into the canonical constrained/indexed tables in
-- one bulk operation, avoiding random index maintenance on every chunk.
CREATE TABLE IF NOT EXISTS v2_generation_atoms (
    protein_id VARCHAR NOT NULL,
    atom_id VARCHAR NOT NULL,
    facet VARCHAR NOT NULL,
    relation VARCHAR NOT NULL,
    source_accession VARCHAR,
    source_db VARCHAR,
    evidence_evalue DOUBLE,
    evidence_score DOUBLE
);
CREATE TABLE IF NOT EXISTS v2_generation_state (
    protein_id VARCHAR NOT NULL,
    activities VARCHAR[],
    roles VARCHAR[],
    architecture VARCHAR[],
    localization VARCHAR[],
    topology JSON,
    size_class VARCHAR,
    quality_flags VARCHAR[],
    composite_predicates VARCHAR[],
    unresolved_count INTEGER,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE IF NOT EXISTS v2_generation_terms (
    protein_id VARCHAR NOT NULL,
    term_id VARCHAR NOT NULL,
    term_kind VARCHAR NOT NULL,
    facet VARCHAR,
    relation VARCHAR,
    source_db VARCHAR NOT NULL DEFAULT '',
    source_accession VARCHAR NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS v2_generation_legacy (
    protein_id VARCHAR NOT NULL,
    predicates VARCHAR[] NOT NULL,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Atomic full-refresh resume boundary. One row records the latest committed
-- protein prefix for the active semantic generation contract.
CREATE TABLE IF NOT EXISTS v2_generation_checkpoint (
    generation_key VARCHAR PRIMARY KEY,
    semantic_fingerprint VARCHAR NOT NULL,
    input_signature VARCHAR NOT NULL,
    status VARCHAR NOT NULL,
    last_protein_id VARCHAR,
    processed_count BIGINT NOT NULL DEFAULT 0,
    total_count BIGINT NOT NULL,
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP
);
"""


def _v2_schema_statements() -> list[str]:
    """Return normalized V2 DDL statements."""
    return [
        statement.strip() + ";"
        for statement in V2_SCHEMA.strip().split(";")
        if statement.strip()
    ]


V2_INDEX_STATEMENTS = tuple(
    statement
    for statement in _v2_schema_statements()
    if statement.lstrip().upper().startswith("CREATE INDEX")
)
V2_INDEX_NAMES = tuple(
    statement.split(" ON ", 1)[0].split()[-1]
    for statement in V2_INDEX_STATEMENTS
)
V2_RETIRED_INDEX_NAMES = (
    "idx_semantic_atoms_protein",
    "idx_semantic_atoms_relation",
    "idx_semantic_state_size",
    "idx_semantic_state_updated",
    "idx_semantic_terms_facet_term",
    "idx_semantic_terms_source",
    "idx_semantic_terms_kind",
)
V2_ALL_INDEX_NAMES = tuple(dict.fromkeys((*V2_INDEX_NAMES, *V2_RETIRED_INDEX_NAMES)))


def create_v2_tables(store: Any, *, create_indexes: bool = True) -> None:
    """Create V2 semantic tables in an existing DuckDB database.

    Args:
        store: A DuckDBStore or any object with an execute() method.
        create_indexes: Build secondary query indexes. Full refreshes disable
            these during append and restore them after bulk promotion.
    """
    for statement in _v2_schema_statements():
        if not create_indexes and statement in V2_INDEX_STATEMENTS:
            continue
        store.execute(statement)


def create_v2_indexes(store: Any) -> None:
    """Create the canonical V2 secondary indexes."""
    for statement in V2_INDEX_STATEMENTS:
        store.execute(statement)


def drop_v2_indexes(store: Any) -> None:
    """Drop secondary V2 indexes before a bulk full refresh."""
    for index_name in V2_ALL_INDEX_NAMES:
        store.execute(f"DROP INDEX IF EXISTS {index_name};")


__all__ = [
    "V2_ALL_INDEX_NAMES",
    "V2_INDEX_NAMES",
    "V2_INDEX_STATEMENTS",
    "V2_SCHEMA",
    "ClaimRelation",
    "SemanticAtom",
    "SemanticFacet",
    "SemanticState",
    "create_v2_indexes",
    "create_v2_tables",
    "drop_v2_indexes",
]
