"""Tests for V2 predicate persistence."""

import pytest

import sharur.predicates_v2.persistence as persistence
from sharur.predicates_v2.model import (
    V2_ALL_INDEX_NAMES,
    V2_INDEX_NAMES,
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
)
from sharur.predicates_v2.persistence import (
    generate_and_persist_v2,
    materialize_semantic_terms_from_v2,
    refresh_composite_predicates_from_v2,
)
from sharur.predicates_v2.validated_systems import materialize_system_proteins
from sharur.storage.duckdb_store import DuckDBStore


def _seed_store() -> DuckDBStore:
    """Create a tiny in-memory store with two proteins and one annotation."""
    store = DuckDBStore()
    conn = store.conn
    conn.execute("""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy)
        VALUES ('bin1', 90.0, 2.0, 'd__Bacteria');

        INSERT INTO contigs (contig_id, bin_id, length, gc_content)
        VALUES ('contig1', 'bin1', 10000, 0.5);

        INSERT INTO proteins (
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence_length, gc_content
        )
        VALUES
            ('p1', 'contig1', 'bin1', 1, 300, '+', 1, 100, 0.5),
            ('p2', 'contig1', 'bin1', 301, 900, '+', 2, 200, 0.5);

        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name, description,
            evalue, score
        )
        VALUES
            (1, 'p1', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter',
             1e-30, 100.0);
    """)
    return store


def test_subset_refresh_preserves_other_v2_rows():
    """Regenerating a subset must not clear semantic rows for other proteins."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO semantic_atoms (
            protein_id, atom_id, facet, relation, source_accession, source_db
        )
        VALUES ('p2', 'existing_atom', 'role', 'implies', 'manual', 'manual');

        INSERT INTO semantic_state (
            protein_id, activities, roles, architecture, localization, topology,
            size_class, quality_flags, composite_predicates, unresolved_count
        )
        VALUES (
            'p2', [], ['existing_atom'], [], [], '{}', 'medium', [], [], 0
        );
    """)

    states = generate_and_persist_v2(
        store,
        protein_ids=["p1"],
        chunk_size=1,
        return_states=True,
    )

    assert set(states) == {"p1"}
    assert store.conn.execute(
        "SELECT COUNT(*) FROM semantic_state WHERE protein_id = 'p2'"
    ).fetchone()[0] == 1
    assert store.conn.execute(
        "SELECT COUNT(*) FROM semantic_atoms WHERE protein_id = 'p2'"
    ).fetchone()[0] == 1


def test_full_refresh_clears_stale_v2_rows():
    """Full regeneration should replace stale semantic rows."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO semantic_state (
            protein_id, activities, roles, architecture, localization, topology,
            size_class, quality_flags, composite_predicates, unresolved_count
        )
        VALUES (
            'stale', ['obsolete'], [], [], [], '{}', 'small', [], [], 0
        );
    """)

    states = generate_and_persist_v2(
        store,
        chunk_size=1,
        return_states=False,
    )

    assert states == {}
    assert store.conn.execute(
        "SELECT COUNT(*) FROM semantic_state WHERE protein_id = 'stale'"
    ).fetchone()[0] == 0
    assert store.conn.execute("SELECT COUNT(*) FROM semantic_state").fetchone()[0] == 2


def test_update_legacy_predicates_materializes_direct_access():
    """Opt-in legacy cache writes V2-derived flat predicates."""
    store = _seed_store()

    generate_and_persist_v2(
        store,
        protein_ids=["p1"],
        chunk_size=1,
        update_legacy_predicates=True,
    )

    row = store.conn.execute(
        "SELECT predicates FROM protein_predicates WHERE protein_id = 'p1'"
    ).fetchone()

    assert row is not None
    assert "pfam:PF00005" in row[0]
    assert "pfam_annotated" in row[0]


def test_generation_materializes_semantic_terms():
    """Generation should write atom, direct-access, and composite search terms."""
    store = _seed_store()

    generate_and_persist_v2(
        store,
        protein_ids=["p1"],
        chunk_size=1,
        return_states=True,
    )

    rows = store.conn.execute("""
        SELECT term_id, term_kind, facet, relation, source_db, source_accession
        FROM semantic_terms
        WHERE protein_id = 'p1'
        ORDER BY term_kind, term_id
    """).fetchall()

    assert ("abc_transporter", "atom", "role", "supports", "pfam", "PF00005") in rows
    assert ("pfam:PF00005", "direct_access", "role", "supports", "pfam", "PF00005") in rows
    assert ("abc_transporter_complete", "composite", None, "implies", "_composite", "abc_transporter_complete") in rows


def test_materialize_semantic_terms_from_existing_v2_state():
    """Backfill should rebuild semantic_terms without recomputing atoms."""
    store = _seed_store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    store.execute("DELETE FROM semantic_terms;")

    written = materialize_semantic_terms_from_v2(store, chunk_size=1)

    assert written > 0
    assert store.conn.execute("""
        SELECT COUNT(*)
        FROM semantic_terms
        WHERE protein_id = 'p1' AND term_id = 'pfam:PF00005'
    """).fetchone()[0] == 1


def test_refresh_composites_updates_state_terms_and_legacy_cache():
    """Composite config fixes should apply without regenerating atoms."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO semantic_atoms (
            protein_id, atom_id, facet, relation, source_accession, source_db
        )
        VALUES
            (
                'p2', 'defense_system', 'role', 'implies',
                'sys_def_1', 'defensefinder_system'
            ),
            (
                'p2', 'restriction_modification', 'role', 'implies',
                'sys_def_1', 'defensefinder_system'
            );

        INSERT INTO semantic_state (
            protein_id, activities, roles, architecture, localization, topology,
            size_class, quality_flags, composite_predicates, unresolved_count
        )
        VALUES (
            'p2', [], ['defense_system', 'restriction_modification'], [], [],
            '{}', 'medium', [], ['crispr_validated'], 0
        );

        INSERT INTO semantic_terms (
            protein_id, term_id, term_kind, facet, relation,
            source_db, source_accession
        )
        VALUES (
            'p2', 'crispr_validated', 'composite', NULL, 'implies',
            '_composite', 'crispr_validated'
        );
    """)

    refreshed = refresh_composite_predicates_from_v2(
        store,
        protein_ids=["p2"],
        chunk_size=1,
        update_semantic_terms=True,
        update_legacy_predicates=True,
    )

    assert refreshed == 1
    composites = store.conn.execute("""
        SELECT composite_predicates
        FROM semantic_state
        WHERE protein_id = 'p2'
    """).fetchone()[0]
    assert "defense_system_validated" in composites
    assert "restriction_modification_validated" in composites
    assert "crispr_validated" not in composites

    terms = {
        row[0]
        for row in store.conn.execute("""
            SELECT term_id
            FROM semantic_terms
            WHERE protein_id = 'p2' AND term_kind = 'composite'
        """).fetchall()
    }
    assert "defense_system_validated" in terms
    assert "restriction_modification_validated" in terms
    assert "crispr_validated" not in terms

    legacy = store.conn.execute("""
        SELECT predicates FROM protein_predicates WHERE protein_id = 'p2'
    """).fetchone()[0]
    assert "defense_system_validated" in legacy
    assert "restriction_modification_validated" in legacy
    assert "crispr_validated" not in legacy


def test_crispr_locus_overlap_emits_v2_quality_atom():
    """ORFs overlapping CRISPR arrays should be represented in V2 state."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO loci (
            locus_id, locus_type, contig_id, start, end_coord, confidence,
            metadata
        )
        VALUES ('crispr_1', 'crispr', 'contig1', 10, 40, 1.0, '{}');
    """)

    states = generate_and_persist_v2(
        store,
        protein_ids=["p1"],
        chunk_size=1,
        update_legacy_predicates=True,
        return_states=True,
    )

    assert "in_crispr_array" in states["p1"].quality_flags

    legacy = store.conn.execute(
        "SELECT predicates FROM protein_predicates WHERE protein_id = 'p1'"
    ).fetchone()[0]
    assert "in_crispr_array" in legacy

    rows = store.conn.execute("""
        SELECT atom_id, relation, source_db, source_accession
        FROM semantic_atoms
        WHERE protein_id = 'p1' AND atom_id = 'in_crispr_array'
    """).fetchall()
    assert rows == [("in_crispr_array", "implies", "_locus", "crispr_1")]


def test_validated_defense_system_table_emits_atoms_without_annotations():
    """Validated defense_systems rows should feed V2 directly."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO defense_systems (
            system_id, genome_id, system_type, system_subtype, activity,
            genes_count, protein_ids, profile_names
        )
        VALUES (
            'sys_def_1', 'bin1', 'RM_Type_II', 'RM_Type_II', 'Defense',
            2, 'p2', 'HsdR'
        );
    """)

    states = generate_and_persist_v2(
        store,
        protein_ids=["p2"],
        chunk_size=1,
        return_states=True,
    )

    state = states["p2"]
    assert "defense_system" in state.roles
    assert "restriction_modification" in state.roles
    assert "rm_type_ii" in state.roles
    assert "unannotated" not in state.quality_flags

    rows = store.conn.execute("""
        SELECT atom_id, relation, source_db, source_accession
        FROM semantic_atoms
        WHERE protein_id = 'p2' AND source_db = 'defensefinder_system'
        ORDER BY atom_id
    """).fetchall()
    assert ("defense_system", "implies", "defensefinder_system", "sys_def_1") in rows

    members = store.conn.execute("""
        SELECT system_id, protein_id, system_source, position, profile_name, score
        FROM system_proteins
        WHERE system_id = 'sys_def_1'
    """).fetchall()
    assert members == [("sys_def_1", "p2", "defensefinder_system", 1, "HsdR", 2.0)]


def test_validated_secretion_system_table_emits_type_atoms_without_annotations():
    """Validated secretion_systems rows should produce type-specific atoms."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO secretion_systems (
            system_id, genome_id, system_type, system_subtype,
            genes_count, protein_ids, profile_names
        )
        VALUES (
            'sys_sec_1', 'bin1', 'T5aSS', 'T5aSS',
            1, 'p2', 'T5aSS_profile'
        );
    """)

    states = generate_and_persist_v2(
        store,
        protein_ids=["p2"],
        chunk_size=1,
        return_states=True,
    )

    state = states["p2"]
    assert "secretion_system" in state.roles
    assert "type_v_secretion" in state.roles
    assert "unannotated" not in state.quality_flags

    rows = store.conn.execute("""
        SELECT atom_id, relation, source_db, source_accession
        FROM semantic_atoms
        WHERE protein_id = 'p2' AND source_db = 'txsscan_system'
        ORDER BY atom_id
    """).fetchall()
    assert ("secretion_system", "implies", "txsscan_system", "sys_sec_1") in rows

    members = store.conn.execute("""
        SELECT system_id, protein_id, system_source, position, profile_name, score
        FROM system_proteins
        WHERE system_id = 'sys_sec_1'
    """).fetchall()
    assert members == [("sys_sec_1", "p2", "txsscan_system", 1, "T5aSS_profile", 1.0)]


def test_materialize_system_proteins_parses_multiple_members():
    """Validated system string membership should normalize to one row per protein."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO defense_systems (
            system_id, genome_id, system_type, system_subtype, activity,
            genes_count, protein_ids, profile_names
        )
        VALUES (
            'sys_def_2', 'bin1', 'RM_Type_I', 'RM_Type_I', 'Defense',
            2, 'p1,p2', 'HsdR,HsdM'
        );
    """)

    written = materialize_system_proteins(store)

    assert written == 2
    assert store.conn.execute("""
        SELECT protein_id, position, profile_name
        FROM system_proteins
        WHERE system_id = 'sys_def_2'
        ORDER BY position
    """).fetchall() == [("p1", 1, "HsdR"), ("p2", 2, "HsdM")]


def _semantic_snapshot(store: DuckDBStore) -> dict[str, list[tuple]]:
    """Return stable semantic outputs while excluding timestamps."""
    return {
        "atoms": store.conn.execute("""
            SELECT protein_id, atom_id, facet, relation, source_accession,
                   source_db, evidence_evalue, evidence_score
            FROM semantic_atoms
            ORDER BY ALL
        """).fetchall(),
        "states": store.conn.execute("""
            SELECT protein_id, activities, roles, architecture, localization,
                   topology, size_class, quality_flags, composite_predicates,
                   unresolved_count
            FROM semantic_state
            ORDER BY protein_id
        """).fetchall(),
        "terms": store.conn.execute("""
            SELECT protein_id, term_id, term_kind, facet, relation,
                   source_db, source_accession
            FROM semantic_terms
            ORDER BY ALL
        """).fetchall(),
        "legacy": store.conn.execute("""
            SELECT protein_id, predicates
            FROM protein_predicates
            ORDER BY protein_id
        """).fetchall(),
    }


def _explicit_index_names(store: DuckDBStore) -> set[str]:
    """Return user-created DuckDB indexes, excluding implicit constraints."""
    return {
        row[0]
        for row in store.conn.execute(
            "SELECT index_name FROM duckdb_indexes()"
        ).fetchall()
    }


def _add_parallel_test_proteins(store: DuckDBStore, count: int = 32) -> None:
    """Add deterministic, annotation-free records for multi-chunk coverage."""
    rows = []
    for index in range(3, count + 3):
        start = index * 300
        rows.append((
            f"z{index:03d}",
            "contig1",
            "bin1",
            start,
            start + 299,
            "+",
            index,
            100,
            0.5,
        ))
    store.conn.executemany(
        """
        INSERT INTO proteins (
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence_length, gc_content
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def test_process_parallel_generation_matches_serial_output():
    """Spawned workers should produce byte-equivalent semantic content."""
    serial_store = _seed_store()
    parallel_store = _seed_store()
    _add_parallel_test_proteins(serial_store)
    _add_parallel_test_proteins(parallel_store)

    generate_and_persist_v2(
        serial_store,
        chunk_size=7,
        workers=1,
        update_legacy_predicates=True,
        return_states=False,
    )
    generate_and_persist_v2(
        parallel_store,
        chunk_size=7,
        workers=4,
        worker_batch_size=2,
        update_legacy_predicates=True,
        return_states=False,
    )

    assert _semantic_snapshot(parallel_store) == _semantic_snapshot(serial_store)
    assert parallel_store.conn.execute("""
        SELECT status, last_protein_id, processed_count, total_count
        FROM v2_generation_checkpoint
        WHERE generation_key = 'full_v2'
    """).fetchone() == ("complete", "z034", 34, 34)


def test_full_generation_defers_indexes_and_releases_build_tables(monkeypatch):
    """Full refreshes should append index-free and publish indexed tables once."""
    store = _seed_store()
    original = persistence._persist_generation_chunk
    observed_during_append: set[str] | None = None

    def inspect_then_persist(*args, **kwargs):
        nonlocal observed_during_append
        observed_during_append = _explicit_index_names(store)
        return original(*args, **kwargs)

    monkeypatch.setattr(persistence, "_persist_generation_chunk", inspect_then_persist)
    generate_and_persist_v2(
        store,
        chunk_size=1,
        workers=1,
        update_legacy_predicates=True,
        return_states=False,
    )

    assert observed_during_append is not None
    assert observed_during_append.isdisjoint(V2_ALL_INDEX_NAMES)
    completed_indexes = _explicit_index_names(store)
    assert set(V2_INDEX_NAMES).issubset(completed_indexes)
    assert completed_indexes.isdisjoint(set(V2_ALL_INDEX_NAMES) - set(V2_INDEX_NAMES))
    tables = {row[0] for row in store.conn.execute("SHOW TABLES").fetchall()}
    assert not tables.intersection(persistence._GENERATION_TABLES)


def test_completed_generation_can_be_resumed_after_build_cleanup():
    """A completed checkpoint should remain resumable after scratch tables drop."""
    store = _seed_store()
    generate_and_persist_v2(
        store,
        chunk_size=1,
        workers=1,
        update_legacy_predicates=True,
        return_states=False,
    )
    before = _semantic_snapshot(store)

    generate_and_persist_v2(
        store,
        chunk_size=1,
        workers=1,
        resume=True,
        update_legacy_predicates=True,
        return_states=False,
    )

    assert _semantic_snapshot(store) == before


def test_generation_preserves_raw_duplicate_evidence_in_terms_and_legacy():
    """Set-oriented persistence must retain pre-dedup confidence evidence."""
    store = _seed_store()
    persistence.create_v2_tables(store, create_indexes=False)
    persistence._clear_generation_tables(store)
    atoms = [
        SemanticAtom(
            protein_id="p1",
            atom_id="duplicate_atom",
            facet=SemanticFacet.role,
            relation=ClaimRelation.supports,
            source_accession="PF00005",
            source_db="pfam",
        ),
        SemanticAtom(
            protein_id="p1",
            atom_id="duplicate_atom",
            facet=SemanticFacet.role,
            relation=ClaimRelation.flags,
            source_accession="PF00005",
            source_db="pfam",
        ),
    ]
    states = {
        "p1": SemanticState(
            protein_id="p1",
            roles=["duplicate_atom"],
            size_class="small",
        )
    }
    exact_legacy = [
        ("p1", ["duplicate_atom", "confident_hit", "weak_hit"]),
    ]

    persistence._persist_generation_chunk(
        store,
        atoms,
        states,
        legacy_rows=exact_legacy,
        update_legacy_predicates=True,
        checkpoint_processed=None,
        checkpoint_last_protein_id=None,
    )
    persistence._promote_full_generation(
        store,
        update_legacy_predicates=True,
        expected_count=1,
    )

    assert store.conn.execute("""
        SELECT relation FROM semantic_atoms
        WHERE protein_id = 'p1' AND atom_id = 'duplicate_atom'
    """).fetchall() == [("flags",)]
    assert store.conn.execute("""
        SELECT term_id, term_kind, relation
        FROM semantic_terms
        WHERE protein_id = 'p1'
          AND term_id IN ('duplicate_atom', 'pfam:PF00005')
        ORDER BY term_kind, relation
    """).fetchall() == [
        ("duplicate_atom", "atom", "flags"),
        ("duplicate_atom", "atom", "supports"),
        ("pfam:PF00005", "direct_access", "flags"),
        ("pfam:PF00005", "direct_access", "supports"),
    ]
    assert store.conn.execute("""
        SELECT predicates FROM protein_predicates WHERE protein_id = 'p1'
    """).fetchone()[0] == exact_legacy[0][1]


def test_full_generation_resumes_at_last_atomic_chunk(monkeypatch):
    """A failed run should retain and resume from its committed protein prefix."""
    store = _seed_store()
    original_persist_chunk = persistence._persist_chunk
    calls = 0

    def fail_before_second_commit(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected generation failure")
        return original_persist_chunk(*args, **kwargs)

    monkeypatch.setattr(persistence, "_persist_chunk", fail_before_second_commit)
    with pytest.raises(RuntimeError, match="injected generation failure"):
        generate_and_persist_v2(
            store,
            chunk_size=1,
            workers=1,
            return_states=False,
        )

    assert store.conn.execute("""
        SELECT status, last_protein_id, processed_count, total_count
        FROM v2_generation_checkpoint
        WHERE generation_key = 'full_v2'
    """).fetchone() == ("failed", "p1", 1, 2)
    assert store.conn.execute(
        "SELECT protein_id FROM semantic_state ORDER BY protein_id"
    ).fetchall() == []
    assert store.conn.execute(
        "SELECT protein_id FROM v2_generation_state ORDER BY protein_id"
    ).fetchall() == [("p1",)]

    monkeypatch.setattr(persistence, "_persist_chunk", original_persist_chunk)
    generate_and_persist_v2(
        store,
        chunk_size=1,
        workers=1,
        resume=True,
        return_states=False,
    )

    assert store.conn.execute(
        "SELECT protein_id FROM semantic_state ORDER BY protein_id"
    ).fetchall() == [("p1",), ("p2",)]
    assert store.conn.execute("""
        SELECT status, last_protein_id, processed_count, total_count
        FROM v2_generation_checkpoint
        WHERE generation_key = 'full_v2'
    """).fetchone() == ("complete", "p2", 2, 2)


def test_full_generation_resumes_after_pre_promotion_failure(monkeypatch):
    """Completed scratch output should survive a failed canonical promotion."""
    store = _seed_store()
    original_promote = persistence._promote_full_generation

    def fail_promotion(*args, **kwargs):
        raise RuntimeError("injected promotion failure")

    monkeypatch.setattr(persistence, "_promote_full_generation", fail_promotion)
    with pytest.raises(RuntimeError, match="injected promotion failure"):
        generate_and_persist_v2(
            store,
            chunk_size=1,
            workers=1,
            update_legacy_predicates=True,
            return_states=False,
        )

    assert store.conn.execute("""
        SELECT status, processed_count, total_count
        FROM v2_generation_checkpoint
        WHERE generation_key = 'full_v2'
    """).fetchone() == ("failed", 2, 2)
    assert store.conn.execute(
        "SELECT COUNT(*) FROM v2_generation_state"
    ).fetchone()[0] == 2
    assert store.conn.execute(
        "SELECT COUNT(*) FROM semantic_state"
    ).fetchone()[0] == 0

    monkeypatch.setattr(persistence, "_promote_full_generation", original_promote)
    generate_and_persist_v2(
        store,
        chunk_size=1,
        workers=1,
        resume=True,
        update_legacy_predicates=True,
        return_states=False,
    )

    assert store.conn.execute(
        "SELECT COUNT(*) FROM semantic_state"
    ).fetchone()[0] == 2


def test_resume_rejects_source_table_drift(monkeypatch):
    """Resume should fail closed when upstream inputs changed after a checkpoint."""
    store = _seed_store()
    original_persist_chunk = persistence._persist_chunk
    calls = 0

    def fail_before_second_commit(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected generation failure")
        return original_persist_chunk(*args, **kwargs)

    monkeypatch.setattr(persistence, "_persist_chunk", fail_before_second_commit)
    with pytest.raises(RuntimeError):
        generate_and_persist_v2(
            store,
            chunk_size=1,
            workers=1,
            return_states=False,
        )

    store.conn.execute("""
        UPDATE proteins SET sequence_length = sequence_length + 1
        WHERE protein_id = 'p2'
    """)
    monkeypatch.setattr(persistence, "_persist_chunk", original_persist_chunk)

    with pytest.raises(ValueError, match="source tables changed"):
        generate_and_persist_v2(
            store,
            chunk_size=1,
            workers=1,
            resume=True,
            return_states=False,
        )


def test_review_queue_is_aggregated_from_persisted_atoms(tmp_path):
    """Review output should retain exact protein and genome aggregation."""
    store = _seed_store()
    store.conn.execute("""
        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name, description,
            evalue, score
        ) VALUES
            (2, 'p1', 'pfam', 'UNMAPPED_X', 'unmapped', 'unmapped', 1e-5, 20),
            (3, 'p2', 'pfam', 'UNMAPPED_X', 'unmapped', 'unmapped', 1e-5, 20)
    """)
    queue_path = tmp_path / "review.tsv"

    generate_and_persist_v2(
        store,
        output_review_queue=str(queue_path),
        chunk_size=1,
        workers=1,
        return_states=False,
    )

    lines = queue_path.read_text().splitlines()
    row = next(line for line in lines if "UNMAPPED_X" in line).split("\t")
    assert row[:4] == ["UNMAPPED_X", "pfam", "2", "1"]
    assert row[5] == "p1;p2"
