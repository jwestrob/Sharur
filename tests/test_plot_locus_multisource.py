from __future__ import annotations

import importlib.util
from pathlib import Path

import duckdb


SCRIPT_PATH = Path(__file__).parents[1] / "scripts" / "plot_locus_multisource.py"
SPEC = importlib.util.spec_from_file_location("plot_locus_multisource", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
PLOT_LOCUS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT_LOCUS)


def _annotation_db() -> duckdb.DuckDBPyConnection:
    connection = duckdb.connect(":memory:")
    connection.execute(
        """
        CREATE TABLE annotations (
            protein_id VARCHAR,
            source VARCHAR,
            accession VARCHAR,
            name VARCHAR,
            evalue DOUBLE,
            score DOUBLE
        )
        """
    )
    return connection


def test_structured_system_call_precedes_raw_profile() -> None:
    connection = _annotation_db()
    connection.executemany(
        "INSERT INTO annotations VALUES (?, ?, ?, ?, ?, ?)",
        [
            (
                "p1",
                "defensefinder",
                "pAgo__pAgo_S2B",
                "pAgo__pAgo_S2B",
                0.0,
                350.0,
            ),
            (
                "p1",
                "defensefinder_system",
                "pAgo__pAgo_SPARTA",
                "pAgo_SPARTA/pAgo_SPARTA",
                None,
                558.0,
            ),
        ],
    )

    label, source = PLOT_LOCUS.get_best_annotation(
        "p1", connection, {}, {}, {}
    )

    assert label == "pAgo_SPARTA: pAgo"
    assert source == "DefenseFinder call"


def test_raw_system_profile_does_not_displace_pfam_observation() -> None:
    connection = _annotation_db()
    connection.executemany(
        "INSERT INTO annotations VALUES (?, ?, ?, ?, ?, ?)",
        [
            (
                "p2",
                "defensefinder",
                "Dnd_ABCDEFGH__DptH",
                "Dnd_ABCDEFGH__DptH",
                0.0,
                80.0,
            ),
            ("p2", "pfam", "PF01850", "PIN", 1e-12, 42.0),
        ],
    )

    label, source = PLOT_LOCUS.get_best_annotation(
        "p2", connection, {}, {}, {}
    )

    assert label == "PIN"
    assert source == "PFAM"


def test_unopposed_raw_system_profile_is_visibly_observational() -> None:
    connection = _annotation_db()
    connection.execute(
        "INSERT INTO annotations VALUES (?, ?, ?, ?, ?, ?)",
        (
            "p3",
            "defensefinder",
            "Dnd_ABCDEFGH__DptH",
            "Dnd_ABCDEFGH__DptH",
            0.0,
            80.0,
        ),
    )

    label, source = PLOT_LOCUS.get_best_annotation(
        "p3", connection, {}, {}, {}
    )

    assert label == "profile: Dnd_ABCDEFGH__DptH"
    assert source == "DefenseFinder profile"
