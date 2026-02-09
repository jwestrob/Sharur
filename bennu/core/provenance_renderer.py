"""
Render provenance DAGs as Mermaid diagrams for publication figures.

Produces Mermaid-format graph definitions showing:
- Provenance entries as nodes (labeled with truncated query text)
- Parent linkages as directed edges
- Hypothesis nodes with status indicators
- Evidence links from provenance entries to hypotheses (support/refute styling)

Usage:
    from bennu.core.provenance_renderer import render_provenance_mermaid
    from bennu.core.session import ExplorationSession

    session = ExplorationSession.load(Path("session.json"))
    mermaid_str = render_provenance_mermaid(session)

    # Write to file for rendering
    Path("provenance.mermaid").write_text(mermaid_str)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from bennu.core.hypothesis_registry import HypothesisRegistry
    from bennu.core.session import ExplorationSession


STATUS_ICONS = {
    "proposed": "?",
    "supported": "+",
    "refuted": "x",
    "uncertain": "~",
}


def _sanitize_label(text: str, max_len: int = 40) -> str:
    """Sanitize text for use as a Mermaid node label."""
    clean = text.replace('"', "'").replace("\n", " ").strip()
    if len(clean) > max_len:
        clean = clean[: max_len - 3] + "..."
    return clean


def render_provenance_mermaid(
    session: ExplorationSession,
    registry: Optional[HypothesisRegistry] = None,
    title: Optional[str] = None,
) -> str:
    """Render the provenance DAG from a session as a Mermaid diagram.

    Args:
        session: ExplorationSession with provenance entries and hypotheses.
        registry: Optional HypothesisRegistry for cross-session hypotheses.
            If provided, its hypotheses are included in addition to
            session-scoped hypotheses.
        title: Optional title displayed at top of diagram.

    Returns:
        Mermaid-format string suitable for rendering to SVG/PNG.
    """
    lines = ["graph TD"]

    if title:
        lines.append(f'    title_node["{_sanitize_label(title, 60)}"]:::titleStyle')
        lines.append("    style title_node fill:none,stroke:none,font-size:16px")

    # --- Provenance entry nodes ---
    provenance = session.get_provenance()
    for entry in provenance:
        node_id = entry.entry_id.hex[:8]
        label = _sanitize_label(entry.query)
        if entry.error:
            lines.append(f'    {node_id}["{label}"]:::errorNode')
        else:
            lines.append(f'    {node_id}["{label}"]')

        # Parent edges
        for pid in entry.parent_ids:
            parent_node = pid.hex[:8]
            lines.append(f"    {parent_node} --> {node_id}")

    # --- Hypothesis nodes (session-scoped) ---
    all_hypotheses = list(session.list_hypotheses())

    # Add registry hypotheses if provided
    if registry is not None:
        registry_ids = {h.hypothesis_id for h in all_hypotheses}
        for h in registry.list_all():
            if h.hypothesis_id not in registry_ids:
                all_hypotheses.append(h)

    for hypo in all_hypotheses:
        hid = hypo.hypothesis_id.hex[:8]
        icon = STATUS_ICONS.get(hypo.status.value, "")
        label = _sanitize_label(hypo.statement, 35)
        lines.append(f'    {hid}["{icon} {label}"]:::hypothesis')

        # Evidence edges from provenance entries to hypotheses
        for ev in hypo.evidence:
            if ev.provenance_id:
                eid = ev.provenance_id.hex[:8]
                if ev.supports:
                    lines.append(f"    {eid} -->|supports| {hid}")
                else:
                    lines.append(f"    {eid} -.->|against| {hid}")

    # --- Styles ---
    lines.append("")
    lines.append("    classDef hypothesis fill:#e8d5f5,stroke:#7b2d8e,stroke-width:2px")
    lines.append("    classDef errorNode fill:#fde8e8,stroke:#c0392b,stroke-width:2px")

    return "\n".join(lines)


def render_provenance_summary(session: ExplorationSession) -> str:
    """Render a text summary of the provenance DAG structure.

    Useful for quick inspection without rendering the Mermaid diagram.
    """
    provenance = session.get_provenance()
    if not provenance:
        return "No provenance entries."

    lines = [f"Provenance DAG: {len(provenance)} entries"]

    # Find roots (entries with no parents)
    all_ids = {e.entry_id for e in provenance}
    roots = [e for e in provenance if not e.parent_ids or not any(p in all_ids for p in e.parent_ids)]
    lines.append(f"Root entries: {len(roots)}")

    # Find leaves (entries not referenced as parent by any other entry)
    referenced_as_parent: set = set()
    for e in provenance:
        referenced_as_parent.update(e.parent_ids)
    leaves = [e for e in provenance if e.entry_id not in referenced_as_parent]
    lines.append(f"Leaf entries: {len(leaves)}")

    # Max depth
    parent_map = {e.entry_id: e.parent_ids for e in provenance}

    def depth(eid: object, seen: set | None = None) -> int:
        if seen is None:
            seen = set()
        if eid in seen:
            return 0  # cycle protection
        seen.add(eid)
        parents = parent_map.get(eid, [])
        if not parents:
            return 0
        return 1 + max(depth(p, seen) for p in parents if p in parent_map)

    max_depth = max(depth(e.entry_id) for e in provenance) if provenance else 0
    lines.append(f"Max chain depth: {max_depth}")

    # Hypotheses linked
    hypotheses = session.list_hypotheses()
    linked = sum(
        1 for h in hypotheses
        if any(ev.provenance_id for ev in h.evidence)
    )
    lines.append(f"Hypotheses: {len(hypotheses)} ({linked} with provenance links)")

    return "\n".join(lines)
