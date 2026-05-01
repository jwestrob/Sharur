"""Tests for V2 predicate persistence."""

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
