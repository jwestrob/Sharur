"""Toy defense-finder-models directory and hit tables for DefenseFinder tests."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd


if TYPE_CHECKING:
    from pathlib import Path


VERSION = "9.9.9"

# stem -> (internal NAME, LENG, GA or None)
PROFILES = {
    "Toy__ToyA": ("Toy__ToyA", 100, (25.0, 25.0)),
    "Toy__ToyB": ("ToyB_internal", 200, None),
    "Solo__SoloA": ("Solo__SoloA", 100, (30.0, 30.0)),
    "Anti__AntiX": ("Anti__AntiX", 100, (20.0, 20.0)),
    "RM_Type_X__MTases_FAM_1": ("RM_Type_X__MTases_FAM_1", 150, (20.0, 20.0)),
    "RM_Type_X__MTases_FAM_2": ("RM_Type_X__MTases_FAM_2", 150, (20.0, 20.0)),
}

MODELS = {
    ("DefenseFinder", "Toy", "Toy"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="2" min_genes_required="2" vers="2.0">\n'
        '<gene name="Toy__ToyA" presence="mandatory"/>\n'
        '<gene name="Toy__ToyB" presence="mandatory"/>\n'
        "</model>\n"
    ),
    ("DefenseFinder", "Solo", "Solo"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="1" min_genes_required="1" vers="2.0">\n'
        '<gene name="Solo__SoloA" presence="mandatory"/>\n'
        "</model>\n"
    ),
    ("RM", "RM_Type_X", "RM_Type_X"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="1" min_genes_required="1" vers="2.0">\n'
        '<gene name="RM_Type_X__MTases_FAM_1" presence="mandatory">\n'
        '<exchangeables><gene name="RM_Type_X__MTases_FAM_2"/></exchangeables>\n'
        "</gene>\n"
        "</model>\n"
    ),
    ("ADF", "Anti", "AntiX"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="1" min_genes_required="1" vers="2.0">\n'
        '<gene name="Anti__AntiX" presence="mandatory"/>\n'
        "</model>\n"
    ),
}

HIT_COLUMNS = [
    "sequence_id", "hmm_name", "bitscore", "evalue", "c_evalue", "i_evalue", "env_from",
    "env_to", "dom_bitscore", "ali_from", "ali_to", "hmm_from", "hmm_to",
]


def write_model_set(root: Path, version: str = VERSION) -> Path:
    """Write a defense-finder-models directory under ``root``; return its path."""
    family = root / "defense-finder-models"
    profiles = family / "profiles"
    profiles.mkdir(parents=True)
    (family / "metadata.yml").write_text(f"short_desc: toy\nvers: {version}\n")
    for stem, (name, length, ga) in PROFILES.items():
        lines = ["HMMER3/f [3.3 | toy]", f"NAME  {name}", f"LENG  {length}", "ALPH  amino"]
        if ga is not None:
            lines.append(f"GA    {ga[0]:g} {ga[1]:g};")
        lines += ["STATS LOCAL MSV      -9.0  0.70", "HMM          A        C", "//"]
        (profiles / f"{stem}.hmm").write_text("\n".join(lines) + "\n")
    for (subfamily, group, model), xml in MODELS.items():
        path = family / "definitions" / subfamily / group / f"{model}.xml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(xml)
    return family


def hit(
    protein: str,
    profile: str,
    *,
    score: float = 80.0,
    dom_score: float | None = None,
    evalue: float = 1e-20,
    i_evalue: float | None = None,
    hmm_from: int = 1,
    hmm_to: int = 90,
    ali: tuple[int, int] = (5, 95),
) -> dict:
    """One per-domain hit row; ``profile`` is the HMM NAME the search reports."""
    return {
        "sequence_id": protein, "hmm_name": profile, "bitscore": score,
        "evalue": evalue, "c_evalue": evalue, "i_evalue": evalue if i_evalue is None else i_evalue,
        "env_from": ali[0], "env_to": ali[1],
        "dom_bitscore": score if dom_score is None else dom_score,
        "ali_from": ali[0], "ali_to": ali[1], "hmm_from": hmm_from, "hmm_to": hmm_to,
    }


def hits_frame(*rows: dict) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=HIT_COLUMNS)


# Gene order on one contig: p1-p2 a two-component Toy locus, p3 a Solo locus, p4 an
# anti-defense (ADF) gene, p5 a former Solo hit, p6 a Pfam-only protein, p7 a lone ToyA.
PROTEINS = ["p1", "p2", "p3", "p4", "p5", "p6", "p7"]
GENE_INDEX = {"p1": 1, "p2": 2, "p3": 10, "p4": 20, "p5": 30, "p6": 40, "p7": 60}


def build_toy_db(path: Path, family: Path, *, call_systems: bool = True,
                 semantics: bool = True) -> Path:
    """A one-genome database with Pfam rows and one former DefenseFinder hit (p5, Solo)."""
    from sharur.colocation import integrate_defense_results, validate_systems
    from sharur.predicates_v2.persistence import generate_and_persist_v2
    from sharur.storage.duckdb_store import DuckDBStore

    path.parent.mkdir(parents=True, exist_ok=True)
    store = DuckDBStore(path)
    store.execute("INSERT INTO bins (bin_id) VALUES ('bin_1')")
    store.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'bin_1', 100000)")
    for pid in PROTEINS:
        i = GENE_INDEX[pid]
        store.execute(
            "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, "
            "gene_index, sequence_length) VALUES (?, 'c1', 'bin_1', ?, ?, '+', ?, 300)",
            [pid, i * 1000, i * 1000 + 900, i],
        )
    store.conn.executemany(
        "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, "
        "evalue, score) VALUES (?, ?, ?, ?, ?, '', ?, ?)",
        [
            (1, "p6", "pfam", "PF00005", "ABC_tran", 1e-40, 110.0),
            (2, "p5", "defensefinder", "Solo__SoloA", "Solo__SoloA", 1e-30, 90.0),
            (3, "p1", "pfam", "PF00004", "AAA", 1e-30, 95.0),
        ],
    )
    store.close()
    if call_systems:
        systems, genes = validate_systems(path, models_dir=family.parent, source="defensefinder",
                                          verbose=False)
        integrate_defense_results(path, systems, genes)
    if semantics:
        store = DuckDBStore(path)
        generate_and_persist_v2(store, update_legacy_predicates=True, return_states=False)
        store.close()
    return path
