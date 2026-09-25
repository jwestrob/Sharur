"""Protein cards and predicate explanations.

:func:`why` traces one predicate on one protein back to its evidence: the
annotation hit that produced each V2 atom, the mapping evidence that ties that
annotation to the predicate (the per-pair evidence strings of the Pfam and KEGG
maps, the HydDB classification, or the validated system call), and, for
predicates reached by hierarchy expansion, the mapped child predicate and the
is-a chain. :func:`card` summarizes one protein in a bounded view for reading
and for agents: location and genome context, annotations by source,
predicates with their best evidence, classifications, system membership, a
neighborhood strip and the predicate-map status. Neither returns sequences.
"""

from __future__ import annotations

from typing import Any

from sharur.operators.predicates_v2 import get_atoms, get_semantic_state
from sharur.predicates.mappings import kegg_map
from sharur.predicates.mappings.cazy_map import cazy_evidence
from sharur.predicates.mappings.pfam_map import PFAM_EVIDENCE
from sharur.predicates.pfam_identity import normalize_pfam_accession
from sharur.predicates.provenance import map_status
from sharur.predicates.vocabulary import PREDICATE_BY_ID, get_hierarchy
from sharur.predicates_v2.composites import explain_composites


COMPUTED_SOURCES = {
    "_property": "computed from protein properties (length, GC)",
    "_computed": "computed from annotation status",
    "_topology": "transmembrane topology prediction",
    "_validation": "cross-annotation validation rule",
}
RULE_SOURCES = {
    "vog": "VOG rule (sharur.predicates.mappings.vog_map)",
    "vogdb": "VOG rule (sharur.predicates.mappings.vog_map)",
    "hyddb": "HydDB HMM class rule",
    "defensefinder": "DefenseFinder profile rule (component level)",
}
SYSTEM_CALLER_SOURCES = {"defensefinder_system", "txsscan_system"}


def _tables(store) -> set[str]:
    return {r[0] for r in store.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'")}


def _map_evidence(source_db: str, accession: str) -> dict[str, str] | None:
    """Per-pair evidence of the Pfam or KEGG map entry an annotation maps through."""
    if source_db == "pfam":
        return PFAM_EVIDENCE.get(normalize_pfam_accession(accession))
    if source_db in ("kegg", "kofam"):
        return kegg_map.KEGG_EVIDENCE.get(accession)
    if source_db == "cazy":
        return cazy_evidence(accession) or None
    return None


def mapping_evidence(source_db: str, accession: str, predicate: str) -> dict[str, Any]:
    """Why an annotation (source, accession) yields ``predicate``."""
    source_db = (source_db or "").lower()
    if source_db in COMPUTED_SOURCES:
        return {"kind": "computed", "evidence": COMPUTED_SOURCES[source_db]}
    if source_db in SYSTEM_CALLER_SOURCES:
        return {"kind": "system_call", "evidence": f"validated system call ({source_db})"}
    if source_db == "hyddb_subgroup":
        return {"kind": "classification", "evidence": "HydDB nearest-reference classification (see classification)"}
    if source_db in ("kegg", "kofam") and not kegg_map.kegg_map_available():
        return {"kind": "rule", "evidence": "no local KEGG map; KEGG-stated EC numbers only (run `sharur setup-kegg`)"}
    if source_db in ("pfam", "kegg", "kofam", "cazy"):
        evidence = _map_evidence(source_db, accession)
        if evidence is None:
            return {"kind": "unsupported", "evidence": "accession absent from the current map (regenerate predicates)"}
        if predicate in evidence:
            return {"kind": "direct", "evidence": evidence[predicate]}
        via = sorted(child for child in evidence if predicate in get_hierarchy(child)[1:])
        if via:
            chain = get_hierarchy(via[0])
            return {"kind": "expanded", "via": via[0], "evidence": evidence[via[0]],
                    "chain": chain[: chain.index(predicate) + 1]}
        return {"kind": "unsupported",
                "evidence": "not supported by the current map (atom from an earlier map; regenerate predicates)"}
    if source_db in RULE_SOURCES:
        return {"kind": "rule", "evidence": RULE_SOURCES[source_db]}
    return {"kind": "rule", "evidence": f"{source_db} mapping rule"}


def _annotation_names(store, protein_id: str) -> dict[tuple[str, str], tuple[str, str]]:
    rows = store.execute(
        "SELECT LOWER(source), accession, name, description FROM annotations WHERE protein_id = ?",
        [protein_id])
    return {(src, acc): (name or "", desc or "") for src, acc, name, desc in rows}


def _hydrogenase_row(store, protein_id: str, tables: set[str]) -> dict[str, Any] | None:
    if "hydrogenase_classifications" not in tables:
        return None
    columns = [r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'hydrogenase_classifications' "
        "ORDER BY ordinal_position")]
    rows = store.execute("SELECT * FROM hydrogenase_classifications WHERE protein_id = ?", [protein_id])
    if not rows:
        return None
    keep = ("outcome", "reference_label", "interpretation_status", "pident", "bitscore", "curation_status",
            "curation_reason", "ko_support", "ko_support_detail", "classifier_version")
    row = dict(zip(columns, rows[0], strict=True))
    return {k: row[k] for k in keep if k in row}


def _map_status(store) -> dict[str, Any]:
    status = map_status(store)
    labels = {"pfam_map_sha256": "Pfam map", "kegg_map_sha256": "KEGG map", "kegg_rules_sha256": "KEGG rules",
              "vocabulary_sha256": "vocabulary", "v2_config_sha256": "V2 config"}
    return {"state": status.state, "changed": [labels[f] for f in status.changed],
            "generated_at": str((status.stamp or {}).get("generated_at")) if status.stamp else None}


# --------------------------------------------------------------------------- #
# why
# --------------------------------------------------------------------------- #


def why(store, protein_id: str, predicate: str) -> dict[str, Any]:
    """Every evidence path from this protein's annotations to ``predicate``."""
    tables = _tables(store)
    result: dict[str, Any] = {
        "protein_id": protein_id,
        "predicate": predicate,
        "definition": None,
        "present": False,
        "paths": [],
        "composite": None,
        "classification": None,
        "map_status": _map_status(store) if "predicate_provenance" in tables else {"state": "unstamped"},
    }
    vocab = PREDICATE_BY_ID.get(predicate)
    if vocab is not None:
        result["definition"] = {"description": vocab.description, "level": vocab.level,
                                "is_a": get_hierarchy(predicate)[1:], "part_of": vocab.part_of}
    atoms = [a for a in get_atoms(store, protein_id) if a.atom_id == predicate] \
        if "semantic_atoms" in tables else []
    names = _annotation_names(store, protein_id)
    for atom in sorted(atoms, key=lambda a: (a.source_db, a.source_accession)):
        name, desc = names.get((atom.source_db, atom.source_accession), ("", ""))
        result["paths"].append({
            "source_db": atom.source_db,
            "accession": atom.source_accession,
            "annotation": name or desc,
            "relation": atom.relation.value,
            "facet": atom.facet.value,
            "evalue": atom.evidence_evalue,
            "score": atom.evidence_score,
            "mapping": mapping_evidence(atom.source_db, atom.source_accession, predicate),
        })
        if atom.source_db == "hyddb_subgroup":
            result["classification"] = _hydrogenase_row(store, protein_id, tables)
    state = get_semantic_state(store, protein_id) if "semantic_state" in tables else None
    if state is not None and predicate in state.composite_predicates:
        result["composite"] = explain_composites(get_atoms(store, protein_id), topology=state.topology,
                                                 only=[predicate]).get(predicate)
    result["present"] = bool(result["paths"] or result["composite"])
    return result


def why_markdown(result: dict[str, Any]) -> str:
    lines = [f"# Why `{result['predicate']}` on {result['protein_id']}"]
    definition = result["definition"]
    if definition:
        lines.append(f"{definition['description']} (level: {definition['level']}"
                     + (f"; is-a: {' > '.join(definition['is_a'])}" if definition["is_a"] else "")
                     + (f"; part of: {definition['part_of']}" if definition["part_of"] else "") + ")")
    if not result["present"]:
        lines.append("\nNot present on this protein.")
    for i, path in enumerate(result["paths"], 1):
        mapping = path["mapping"]
        hit = f"{path['source_db']}:{path['accession']}"
        if path["annotation"]:
            hit += f" ({path['annotation']})"
        score = f"e={path['evalue']:.1e}" if path["evalue"] is not None else ""
        lines.append(f"\n{i}. {hit} {score} -> {path['relation']} [{path['facet']}]")
        if mapping["kind"] == "expanded":
            lines.append(f"   via `{mapping['via']}` ({mapping['evidence']}); is-a {' > '.join(mapping['chain'])}")
        else:
            lines.append(f"   {mapping['kind']}: {mapping['evidence']}")
    if result["classification"]:
        c = result["classification"]
        lines.append(f"\nClassification: {c.get('reference_label')} ({c.get('interpretation_status')}), "
                     f"{c.get('curation_status')}; KOfam support: {c.get('ko_support')}"
                     + (f" [{c['ko_support_detail']}]" if c.get("ko_support_detail") else ""))
    if result["composite"]:
        lines.append(f"\nComposite: {result['composite']}")
    status = result["map_status"]
    lines.append(f"\nPredicate maps: {status['state']}"
                 + (f" (changed: {', '.join(status['changed'])})" if status.get("changed") else ""))
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# card
# --------------------------------------------------------------------------- #

FACETS = ("activities", "roles", "architecture", "localization", "topology", "quality_flags")


def _best_evidence(atoms, predicate: str) -> str:
    """Shortest useful evidence string for a predicate: strongest atom's source and mapping."""
    rank = {"implies": 0, "supports": 1, "flags": 2, "excludes": 3, "unresolved": 4}
    candidates = sorted((a for a in atoms if a.atom_id == predicate),
                        key=lambda a: (rank.get(a.relation.value, 9), a.evidence_evalue or 1.0))
    if not candidates:
        return ""
    atom = candidates[0]
    mapping = mapping_evidence(atom.source_db, atom.source_accession, predicate)
    where = atom.source_db if atom.source_db.startswith("_") else f"{atom.source_db}:{atom.source_accession}"
    how = mapping["evidence"] if mapping["kind"] != "expanded" else f"via {mapping['via']}"
    more = f" +{len(candidates) - 1}" if len(candidates) > 1 else ""
    return f"{where} ({atom.relation.value}; {how}){more}"


def card(store, protein_id: str, window: int = 5, max_per_source: int = 8) -> dict[str, Any]:
    """A bounded summary of one protein (no sequences)."""
    tables = _tables(store)
    rows = store.execute(
        "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length, gc_content "
        "FROM proteins WHERE protein_id = ?", [protein_id])
    if not rows:
        return {"protein_id": protein_id, "found": False}
    pid, contig, bin_id, start, end, strand, gene_index, length, gc = rows[0]
    genome: dict[str, Any] = {"bin_id": bin_id}
    bin_rows = store.execute("SELECT * FROM bins WHERE bin_id = ?", [bin_id]) if bin_id else []
    if bin_rows:
        columns = [r[0] for r in store.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name = 'bins' ORDER BY ordinal_position")]
        row = dict(zip(columns, bin_rows[0], strict=True))
        genome.update({k: row.get(k) for k in ("taxonomy", "completeness", "contamination") if k in row})
    contig_rows = store.execute("SELECT length FROM contigs WHERE contig_id = ?", [contig])

    annotations: dict[str, list[dict[str, Any]]] = {}
    for source, accession, name, desc, evalue, score in store.execute(
            "SELECT LOWER(source), accession, name, description, evalue, score FROM annotations "
            "WHERE protein_id = ? ORDER BY source, evalue NULLS LAST, score DESC NULLS LAST", [protein_id]):
        annotations.setdefault(source, []).append(
            {"accession": accession, "name": name, "description": desc, "evalue": evalue, "score": score})
    truncated = {src: len(hits) - max_per_source for src, hits in annotations.items() if len(hits) > max_per_source}
    annotations = {src: hits[:max_per_source] for src, hits in annotations.items()}

    atoms = get_atoms(store, protein_id) if "semantic_atoms" in tables else []
    state = get_semantic_state(store, protein_id) if "semantic_state" in tables else None
    predicates: dict[str, dict[str, str]] = {}
    if state is not None:
        state_dict = state.to_dict()
        for facet in FACETS:
            values = state_dict.get(facet) or []
            values = values if isinstance(values, list) else list(values)
            predicates[facet] = {p: _best_evidence(atoms, p) for p in sorted(values)}
        predicates["composites"] = dict.fromkeys(sorted(state.composite_predicates), "composite")

    systems = []
    if "system_proteins" in tables:
        systems = [{"system_id": s, "source": src, "profile": prof}
                   for s, src, prof in store.execute(
                       "SELECT system_id, system_source, profile_name FROM system_proteins WHERE protein_id = ?",
                       [protein_id])]

    neighborhood = []
    if gene_index is not None:
        neighbors = store.execute(
            "SELECT protein_id, gene_index, strand, sequence_length FROM proteins "
            "WHERE contig_id = ? AND gene_index BETWEEN ? AND ? ORDER BY gene_index",
            [contig, gene_index - window, gene_index + window])
        ids = [n[0] for n in neighbors]
        best: dict[str, str] = {}
        if ids:
            placeholders = ",".join("?" * len(ids))
            for nid, source, accession, name in store.execute(
                    f"SELECT protein_id, LOWER(source), accession, COALESCE(NULLIF(name, ''), description) "
                    f"FROM annotations WHERE protein_id IN ({placeholders}) "
                    f"AND LOWER(source) IN ('pfam', 'kofam', 'kegg', 'hyddb', 'cazy', 'defensefinder') "
                    f"ORDER BY evalue NULLS LAST", ids):
                best.setdefault(nid, f"{source}:{name or accession}")
        neighborhood = [{"offset": gi - gene_index, "protein_id": nid, "strand": st, "length": ln,
                         "top_annotation": best.get(nid, "")} for nid, gi, st, ln in neighbors]

    return {
        "protein_id": pid,
        "found": True,
        "location": {"contig_id": contig, "contig_length": contig_rows[0][0] if contig_rows else None,
                     "start": start, "end": end, "strand": strand, "gene_index": gene_index,
                     "length_aa": length, "gc": gc},
        "genome": genome,
        "annotations": annotations,
        "annotations_truncated": truncated,
        "predicates": predicates,
        "hydrogenase_classification": _hydrogenase_row(store, protein_id, tables),
        "validated_systems": systems,
        "neighborhood": neighborhood,
        "map_status": _map_status(store) if "predicate_provenance" in tables else {"state": "unstamped"},
    }


def card_markdown(c: dict[str, Any]) -> str:
    if not c.get("found"):
        return f"Protein not found: {c['protein_id']}"
    loc, genome = c["location"], c["genome"]
    lines = [f"# {c['protein_id']}",
             f"{loc['length_aa']} aa, {loc['contig_id']}:{loc['start']}-{loc['end']} ({loc['strand']}), "
             f"gene {loc['gene_index']}" + (f" of a {loc['contig_length']:,} bp contig" if loc["contig_length"] else ""),
             f"Genome {genome.get('bin_id')}"
             + (f": {genome['taxonomy']}" if genome.get("taxonomy") else "")
             + (f" ({genome['completeness']}% complete, {genome['contamination']}% contamination)"
                if genome.get("completeness") is not None else "")]
    lines.append("\n## Annotations")
    for source, hits in c["annotations"].items():
        shown = "; ".join(
            f"{h['name'] or h['accession']}" + (f" [{h['accession']}]" if h["name"] and h["name"] != h["accession"] else "")
            + (f" e={h['evalue']:.0e}" if h["evalue"] is not None else "") for h in hits)
        extra = c["annotations_truncated"].get(source)
        lines.append(f"- {source}: {shown}" + (f" (+{extra} more)" if extra else ""))
    if c["predicates"]:
        lines.append("\n## Predicates")
        for facet, preds in c["predicates"].items():
            if preds:
                lines.append(f"- {facet}: " + "; ".join(f"{p}" + (f" <- {e}" if e and e != 'composite' else "")
                                                         for p, e in preds.items()))
    if c["hydrogenase_classification"]:
        h = c["hydrogenase_classification"]
        lines.append(f"\n## Hydrogenase\n{h.get('outcome')}: {h.get('reference_label')} "
                     f"({h.get('interpretation_status')}), pident {h.get('pident')}, {h.get('curation_status')}; "
                     f"KOfam support {h.get('ko_support')}")
    if c["validated_systems"]:
        lines.append("\n## Validated systems\n" + "; ".join(
            f"{s['system_id']} ({s['source']}, {s['profile']})" for s in c["validated_systems"]))
    if c["neighborhood"]:
        lines.append("\n## Neighborhood")
        for n in c["neighborhood"]:
            mark = ">>" if n["offset"] == 0 else "  "
            lines.append(f"{mark} {n['offset']:+d} {n['strand']} {n['length']} aa {n['top_annotation'] or '-'}")
    status = c["map_status"]
    lines.append(f"\nPredicate maps: {status['state']}"
                 + (f" (changed: {', '.join(status['changed'])})" if status.get("changed") else ""))
    return "\n".join(lines)
