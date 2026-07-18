"""Shared focused fixtures for cross-surface case/evidence tests."""

from __future__ import annotations

import pytest

from sharur.storage.duckdb_store import DuckDBStore


def build_case_database(path):
    """Build a tiny multi-replicon caller dataset with known context counts."""
    store = DuckDBStore(path)
    bins = [
        ("fg_plus", "d__Archaea;p__Testarchaeota"),
        ("fg_minus", "d__Archaea;p__Testarchaeota"),
        ("bg_positive", "d__Archaea;p__Testarchaeota"),
        ("bg_negative", "d__Archaea;p__Testarchaeota"),
        ("bg_edge", "d__Archaea;p__Testarchaeota"),
    ]
    store.conn.executemany(
        """
        INSERT INTO bins(bin_id, taxonomy, n_contigs, total_length)
        VALUES (?, ?, 1, 10000)
        """,
        bins,
    )
    store.conn.executemany(
        """
        INSERT INTO contigs(contig_id, bin_id, length, gc_content)
        VALUES (?, ?, 10000, 0.5)
        """,
        [(f"{bin_id}_contig", bin_id) for bin_id, _ in bins],
    )

    protein_rows = []
    for bin_id, _taxonomy in bins:
        gene_count = 4 if bin_id == "bg_edge" else 9
        for gene_index in range(gene_count):
            strand = "+"
            if bin_id == "fg_minus" and gene_index in {3, 4}:
                strand = "-1"
            protein_rows.append(
                (
                    f"{bin_id}_p{gene_index}",
                    f"{bin_id}_contig",
                    bin_id,
                    gene_index * 300 + 1,
                    gene_index * 300 + 270,
                    strand,
                    gene_index,
                    "M" + "A" * 89,
                    90,
                )
            )
    store.conn.executemany(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence, sequence_length
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        protein_rows,
    )
    store.conn.executemany(
        "INSERT INTO protein_predicates(protein_id, predicates) VALUES (?, [])",
        [(row[0],) for row in protein_rows],
    )

    systems = [
        ("target_plus", "fg_plus", "Target", 3, 4),
        ("target_minus", "fg_minus", "Target", 3, 4),
        ("control_positive", "bg_positive", "Control", 3, 4),
        ("control_positive_same_replicon", "bg_positive", "Control", 5, 6),
        ("control_negative", "bg_negative", "Control", 3, 4),
        ("control_edge", "bg_edge", "Control", 0, 1),
    ]
    store.conn.executemany(
        """
        INSERT INTO defense_systems(
            system_id, genome_id, contig_id, system_type, system_subtype,
            genes_count, protein_ids
        ) VALUES (?, ?, ?, ?, ?, 2, ?)
        """,
        [
            (
                system_id,
                bin_id,
                f"{bin_id}_contig",
                system_type,
                system_type,
                f"{bin_id}_p{left},{bin_id}_p{right}",
            )
            for system_id, bin_id, system_type, left, right in systems
        ],
    )
    store.conn.executemany(
        """
        INSERT INTO system_proteins(
            system_id, protein_id, system_source, position, profile_name, score
        ) VALUES (?, ?, 'defensefinder_system', ?, ?, 100)
        """,
        [
            (
                system_id,
                f"{bin_id}_p{gene_index}",
                position,
                f"{system_type}__component_{position}",
            )
            for system_id, bin_id, system_type, left, right in systems
            for position, gene_index in enumerate((left, right), start=1)
        ],
    )

    annotations = [
        ("fg_plus_p1", "pfam", "PFTEST", "test_domain"),
        ("fg_minus_p6", "pfam", "PFTEST", "test_domain"),
        ("bg_positive_p1", "pfam", "PFTEST", "test_domain"),
        ("bg_positive_p7", "pfam", "PFTEST", "test_domain"),
        ("fg_plus_p3", "defensefinder", "raw_profile", "raw_component_profile"),
        (
            "fg_plus_p3",
            "defensefinder_system",
            "Target",
            "structured_projection_row",
        ),
    ]
    store.conn.executemany(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name, description,
            evalue, score
        ) VALUES (?, ?, ?, ?, ?, ?, 1e-20, 80)
        """,
        [
            (index, protein_id, source, accession, name, name)
            for index, (protein_id, source, accession, name) in enumerate(
                annotations,
                start=1,
            )
        ],
    )
    return store


@pytest.fixture
def case_database(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = build_case_database(path)
    store.close()
    return path
