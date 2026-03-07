"""
Aggregator: atoms -> SemanticState per protein.

Takes a list of SemanticAtom objects for a single protein and produces
a resolved SemanticState. Handles:
- Grouping atoms by facet
- Conflict resolution (excludes beats supports/flags)
- Deduplication (multiple sources supporting the same atom)
- Collecting unresolved accessions
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
)


def aggregate_atoms(
    protein_id: str,
    atoms: list[SemanticAtom],
) -> SemanticState:
    """Aggregate a list of atoms into a resolved SemanticState.

    Resolution rules:
    - excludes beats supports and flags (but NOT implies — that's a conflict)
    - implies + excludes on the same atom_id = conflict -> surfaced as unresolved
    - Multiple supports for the same atom = keep (convergent evidence)
    - unresolved atoms don't appear in resolved state — they go to review queue

    Args:
        protein_id: The protein this state belongs to.
        atoms: All atoms for this protein.

    Returns:
        Resolved SemanticState.
    """
    # Group atoms by (facet, atom_id) for resolution
    groups: dict[tuple[str, str], list[SemanticAtom]] = defaultdict(list)
    unresolved_accessions: list[str] = []

    for atom in atoms:
        if atom.relation == ClaimRelation.unresolved:
            # Collect unresolved accessions — don't put them in facet groups
            acc = atom.source_accession
            if acc and acc not in unresolved_accessions:
                unresolved_accessions.append(acc)
            continue

        key = (atom.facet.value, atom.atom_id)
        groups[key].append(atom)

    # Resolve each group
    resolved_by_facet: dict[str, list[str]] = defaultdict(list)
    conflicts: list[str] = []

    for (facet_str, atom_id), group_atoms in groups.items():
        relations = {a.relation for a in group_atoms}

        # Check for implies + excludes conflict
        if ClaimRelation.implies in relations and ClaimRelation.excludes in relations:
            conflicts.append(atom_id)
            continue

        # excludes beats supports and flags
        if ClaimRelation.excludes in relations:
            # Atom is excluded — don't add to resolved state
            continue

        # Atom survives — add to resolved facet
        if atom_id not in resolved_by_facet[facet_str]:
            resolved_by_facet[facet_str].append(atom_id)

    # Build topology dict from topology atoms
    topology: dict[str, Any] = {}
    topology_atoms = resolved_by_facet.get("topology", [])
    if "transmembrane_predicted" in topology_atoms:
        # Count TM helices from specific atoms
        tm_count = 0
        if "polytopic_membrane" in topology_atoms:
            tm_count = 4  # Minimum for polytopic
        elif "multi_pass_membrane" in topology_atoms:
            tm_count = 2  # Minimum for multi-pass
        elif "single_pass_membrane" in topology_atoms:
            tm_count = 1
        topology["tm_count"] = tm_count
        topology["signal_peptide"] = False  # TODO: integrate signal peptide prediction

        # N-terminus orientation
        if "n_in_topology" in topology_atoms:
            topology["n_terminus"] = "inside"
        elif "n_out_topology" in topology_atoms:
            topology["n_terminus"] = "outside"
    elif "soluble_predicted" in topology_atoms:
        topology["tm_count"] = 0
        topology["signal_peptide"] = False

    # Determine size_class (exactly one)
    size_atoms = resolved_by_facet.get("size_class", [])
    size_class = ""
    # Priority order: massive > giant > large > medium > small > tiny
    size_priority = ["massive", "giant", "large", "medium", "small", "tiny"]
    for sz in size_priority:
        if sz in size_atoms:
            size_class = sz
            break

    # Add conflicts to unresolved
    for conflict_atom in conflicts:
        conflict_label = f"conflict:{conflict_atom}"
        if conflict_label not in unresolved_accessions:
            unresolved_accessions.append(conflict_label)

    return SemanticState(
        protein_id=protein_id,
        activities=sorted(resolved_by_facet.get("activity", [])),
        roles=sorted(resolved_by_facet.get("role", [])),
        architecture=sorted(resolved_by_facet.get("architecture", [])),
        localization=sorted(resolved_by_facet.get("localization", [])),
        topology=topology,
        size_class=size_class,
        quality_flags=sorted(resolved_by_facet.get("quality_flag", [])),
        unresolved_accessions=sorted(unresolved_accessions),
    )


__all__ = ["aggregate_atoms"]
