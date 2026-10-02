"""Annotation hits for display: names for bare accessions, external links, domain lanes.

KOfam rows store only the KO number (``name``) and how the hit passed
(``description``: ``GA`` for KOfam's adaptive score threshold, ``evalue_1e-15``
for KOs without one). KO symbols and definitions come from the local KEGG build
(``sharur setup-kegg``), falling back to the dataset's own ``kegg`` rows.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from markupsafe import Markup, escape

KO = re.compile(r"^K\d{5}$")
EC = re.compile(r"\[EC:([0-9.\- ]+)\]")
PFAM = re.compile(r"^PF\d{5}")
CAZY = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+")

SOURCE_LABELS = {"pfam": "Pfam", "pfam_relaxed": "Pfam (relaxed)", "kofam": "KOfam", "kegg": "KEGG", "defensefinder": "DefenseFinder",
                 "txsscan": "TXSScan", "hyddb": "HydDB", "hyddb_subgroup": "HydDB subgroup", "cazy": "CAZy",
                 "vogdb": "VOGdb", "vog": "VOGdb", "defensefinder_system": "DefenseFinder system",
                 "txsscan_system": "TXSScan system", "cctyper": "CCTyper"}


def _rows(store, sql: str, params: list | None = None) -> list:
    """Rows from a Sharur store or a raw DuckDB connection."""
    result = store.execute(sql, params or [])
    return result.fetchall() if hasattr(result, "fetchall") else list(result)


class KoNames:
    """KO -> (symbols, definition), loaded once."""

    def __init__(self, store, kegg_dir: Path | None = None):
        self._store, self._dir, self._names = store, kegg_dir, None

    def _load(self) -> dict[str, tuple[str, str]]:
        names: dict[str, tuple[str, str]] = {}
        try:
            for ko, desc in _rows(self._store,
                                  "SELECT accession, MIN(description) FROM annotations "
                                  "WHERE LOWER(source) = 'kegg' AND description IS NOT NULL GROUP BY 1"):
                names[ko] = ("", desc)
        except Exception:  # noqa: BLE001 - datasets without annotations still render
            pass
        kegg_dir = self._dir
        if kegg_dir is None:
            try:
                from sharur.predicates.mappings.kegg_map import default_kegg_dir

                kegg_dir = default_kegg_dir()
            except Exception:  # noqa: BLE001
                kegg_dir = None
        snapshot = Path(kegg_dir) / "kegg_evidence_snapshot.tsv" if kegg_dir else None
        if snapshot and snapshot.is_file():
            with open(snapshot) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    cells = line.rstrip("\n").split("\t")
                    if len(cells) >= 3 and KO.match(cells[0]):
                        names[cells[0]] = (cells[1], cells[2])
        return names

    def get(self, ko: str) -> tuple[str, str] | None:
        if self._names is None:
            self._names = self._load()
        return self._names.get(ko)


def external_url(source: str, accession: str) -> str | None:
    acc = (accession or "").split(".")[0]
    if KO.match(acc):
        return f"https://www.kegg.jp/entry/{acc}"
    if source == "pfam" and PFAM.match(acc):
        return f"https://www.ebi.ac.uk/interpro/entry/pfam/{acc}/"
    if source == "cazy" and CAZY.match(acc):
        return f"https://www.cazy.org/{acc}.html"
    return None


def link_ec(text: str) -> Markup:
    """Escape ``text`` and link its ``[EC:...]`` numbers to KEGG ENZYME."""
    out, last = [], 0
    for m in EC.finditer(text):
        out.append(escape(text[last:m.start()]))
        ecs = [e for e in m.group(1).split() if e]
        links = " ".join(f'<a href="https://www.kegg.jp/entry/ec:{escape(e)}" target="_blank" rel="noopener">{escape(e)}</a>'
                         if "-" not in e else str(escape(e)) for e in ecs)
        out.append(Markup(f"[EC:{links}]"))
        last = m.end()
    out.append(escape(text[last:]))
    return Markup("").join(out)


def threshold_text(description: str | None) -> str:
    if description == "GA":
        return "passes KOfam score threshold"
    if description and description.startswith("evalue_"):
        return f"E ≤ {description[len('evalue_'):]} (KO has no KOfam threshold)"
    return description or ""


def hit_view(source: str, hit: dict[str, Any], ko_names: KoNames) -> dict[str, Any]:
    """Display fields for one annotation row: name, detail and external link."""
    acc, name, desc = hit["accession"], hit.get("name") or "", hit.get("description") or ""
    view = {"accession": acc, "external": external_url(source, acc), "name": name, "detail": Markup(escape(desc)),
            "symbol": ""}
    if KO.match(acc.split(".")[0]):
        known = ko_names.get(acc.split(".")[0])
        if known:
            view["symbol"], definition = known
            view["name"] = link_ec(definition)
        if source == "kofam":
            view["detail"] = Markup(escape(threshold_text(desc)))
        elif known and desc == known[1]:
            view["detail"] = Markup("")
    return view


def describe_ko(label: str, ko_names: KoNames) -> str:
    """'kofam:K03737' -> 'kofam:K03737 por (pyruvate-ferredoxin/flavodoxin oxidoreductase)'."""
    m = re.search(r"\bK\d{5}\b", label)
    if not m:
        return label
    known = ko_names.get(m.group(0))
    if not known:
        return label
    symbols, definition = known
    definition = EC.sub("", definition).strip()
    first = symbols.split(",")[0].strip()
    return f"{label} {first} ({definition})" if first else f"{label} ({definition})"


def domain_lanes(store, protein_id: str) -> list[tuple[str, list[dict[str, Any]]]]:
    """Placed hits per source, overlaps resolved within each source; Pfam first."""
    from sharur.architecture import Domain, resolve

    by_source: dict[str, list] = {}
    for name, accession, start, end, evalue, source in _rows(
            store, """SELECT COALESCE(NULLIF(name, ''), accession), accession, start_aa, end_aa, evalue, LOWER(source)
               FROM annotations WHERE protein_id = ? AND start_aa IS NOT NULL AND end_aa IS NOT NULL
                 AND end_aa >= start_aa""", [protein_id]):
        by_source.setdefault(source, []).append(Domain(name, accession, start, end, evalue, source))
    order = sorted(by_source, key=lambda s: (s != "pfam", s))
    return [(s, [d.to_dict() for d in resolve(by_source[s])]) for s in order]
