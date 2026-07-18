"""Regression tests for contig-local MacSyFinder-style system calling."""

from __future__ import annotations

from typing import TYPE_CHECKING

import duckdb
import pandas as pd
import pytest

from sharur.colocation import (
    GeneSpec,
    HitRecord,
    SystemModel,
    _build_hmm_name_mapping,
    _build_model_clusters,
    _cluster_contig,
    _detect_model_systems,
    _parse_model_xml,
    _select_best_hits_per_position,
    _translate_hit_accessions,
    integrate_secretion_results,
    resolve_conflicts,
    validate_systems,
)


if TYPE_CHECKING:
    from pathlib import Path


def _model(
    *genes: GeneSpec,
    name: str = "test-system",
    gap: int = 5,
    min_mandatory: int | None = None,
    min_genes: int | None = None,
    multi_loci: bool = False,
) -> SystemModel:
    quorum = sum(gene.presence == "mandatory" for gene in genes)
    model = SystemModel(
        name=name,
        family="test-models",
        inter_gene_max_space=gap,
        min_mandatory_genes_required=(quorum if min_mandatory is None else min_mandatory),
        min_genes_required=quorum if min_genes is None else min_genes,
        max_nb_genes=sum(gene.presence in {"mandatory", "accessory"} for gene in genes),
        multi_loci=multi_loci,
        genes=list(genes),
    )
    model._build_lookups()
    return model


def _hit(
    protein_id: str,
    accession: str,
    gene_index: int,
    *,
    contig: str = "contig-1",
    genome: str = "genome-1",
    score: float = 100.0,
    evalue: float = 1e-20,
    annotation_id: int = 1,
) -> HitRecord:
    return HitRecord(
        protein_id=protein_id,
        accession=accession,
        score=score,
        contig_id=contig,
        bin_id=genome,
        gene_index=gene_index,
        evalue=evalue,
        annotation_id=annotation_id,
    )


def _seed_integration_db(path: Path) -> None:
    conn = duckdb.connect(str(path))
    conn.execute("""
        CREATE TABLE proteins (
            protein_id VARCHAR,
            contig_id VARCHAR,
            bin_id VARCHAR,
            gene_index INTEGER
        );
        INSERT INTO proteins VALUES
            ('p1', 'contig-1', 'genome-1', 1),
            ('p2', 'contig-1', 'genome-1', 2),
            ('p3', 'contig-2', 'genome-1', 1);

        CREATE TABLE annotations (
            annotation_id BIGINT,
            protein_id VARCHAR,
            source VARCHAR,
            accession VARCHAR,
            name VARCHAR,
            description VARCHAR,
            evalue DOUBLE,
            score DOUBLE,
            start_aa INTEGER,
            end_aa INTEGER
        );
        INSERT INTO annotations VALUES
            (1, 'p1', 'txsscan', 'a', 'a', '', 1e-20, 100, NULL, NULL),
            (2, 'p1', 'txsscan_system', 'stale', 'stale', '', NULL, 1, NULL, NULL);

        CREATE TABLE secretion_systems (
            system_id VARCHAR PRIMARY KEY,
            genome_id VARCHAR,
            system_type VARCHAR,
            system_subtype VARCHAR,
            genes_count INTEGER,
            protein_ids VARCHAR,
            profile_names VARCHAR,
            sys_beg VARCHAR,
            sys_end VARCHAR,
            created_at TIMESTAMP
        );
        INSERT INTO secretion_systems VALUES
            ('stale', 'genome-1', 'stale', 'stale', 1, 'p1', 'stale', '', '', CURRENT_TIMESTAMP);
    """)
    conn.close()


def _system_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    systems = pd.DataFrame(
        [
            {
                "system_id": "genome-1_test_1",
                "genome_id": "genome-1",
                "contig_id": "contig-1",
                "type": "test",
                "subtype": "test",
                "protein_ids": "p1,p2",
                "genes_count": 2,
                "profile_names": "a,b",
                "score": 200.0,
                "hit_id": "genome-1_test_1",
            }
        ]
    )
    genes = pd.DataFrame(
        [
            {
                "protein_id": protein_id,
                "system_id": "genome-1_test_1",
                "accession": accession,
                "score": 100.0,
                "presence": "mandatory",
                "type": "test",
                "subtype": "test",
                "genome_id": "genome-1",
                "contig_id": "contig-1",
            }
            for protein_id, accession in (("p1", "a"), ("p2", "b"))
        ]
    )
    return systems, genes


def test_xml_parser_preserves_model_gap_when_genes_have_overrides(
    tmp_path: Path,
) -> None:
    definition = tmp_path / "spacing.xml"
    definition.write_text(
        """
        <model inter_gene_max_space="9"
               min_mandatory_genes_required="2"
               min_genes_required="2">
          <gene name="a" presence="mandatory" inter_gene_max_space="2" />
          <gene name="b" presence="mandatory" />
        </model>
        """
    )

    model = _parse_model_xml(definition, "test-models")

    assert model is not None
    assert model.inter_gene_max_space == 9
    assert model.genes[0].inter_gene_max_space == 2
    assert model.genes[1].inter_gene_max_space is None


def test_best_profile_is_selected_once_per_protein_position() -> None:
    hits = [
        _hit("p1", "lower", 4, score=80.0, annotation_id=1),
        _hit("p1", "best", 4, score=100.0, annotation_id=2),
        _hit("p1", "duplicate-domain", 4, score=90.0, annotation_id=3),
    ]

    selected = _select_best_hits_per_position(hits)

    assert [(hit.protein_id, hit.accession) for hit in selected] == [("p1", "best")]


def test_best_profile_ties_are_deterministic() -> None:
    left = _hit("p1", "z-profile", 4, score=100.0, evalue=1e-20, annotation_id=2)
    right = _hit("p1", "a-profile", 4, score=100.0, evalue=1e-20, annotation_id=1)

    forward = _select_best_hits_per_position([left, right])
    reverse = _select_best_hits_per_position([right, left])

    assert [hit.accession for hit in forward] == ["a-profile"]
    assert [hit.accession for hit in reverse] == ["a-profile"]


def test_best_hit_position_key_includes_genome_for_reused_contig_labels() -> None:
    hits = [
        _hit("p1", "a", 4, contig="shared", genome="genome-1"),
        _hit("p2", "b", 4, contig="shared", genome="genome-2"),
    ]

    selected = _select_best_hits_per_position(hits)

    assert [hit.protein_id for hit in selected] == ["p1", "p2"]


def test_inter_gene_max_space_counts_intervening_genes() -> None:
    within = [_hit("p1", "a", 0), _hit("p2", "b", 6)]
    beyond = [_hit("p1", "a", 0), _hit("p2", "b", 7)]

    assert len(_cluster_contig(within, max_gap=5)) == 1
    assert len(_cluster_contig(beyond, max_gap=5)) == 2


def test_multi_locus_quorum_never_combines_different_contigs() -> None:
    model = _model(
        *(GeneSpec(name, "mandatory") for name in ("a", "b", "c", "d")),
        multi_loci=True,
    )
    hits = [
        _hit("p1", "a", 1, contig="contig-1"),
        _hit("p2", "b", 2, contig="contig-1"),
        _hit("p3", "c", 1, contig="contig-2"),
        _hit("p4", "d", 2, contig="contig-2"),
    ]

    assert _detect_model_systems(hits, model) == []


def test_reused_contig_labels_are_called_independently_per_genome() -> None:
    model = _model(
        GeneSpec("a", "mandatory"),
        GeneSpec("b", "mandatory"),
    )
    hits = [
        _hit("p1", "a", 1, contig="shared", genome="genome-1"),
        _hit("p2", "b", 2, contig="shared", genome="genome-1"),
        _hit("p3", "a", 1, contig="shared", genome="genome-2"),
        _hit("p4", "b", 2, contig="shared", genome="genome-2"),
    ]

    systems = _detect_model_systems(hits, model)

    assert [(system.genome_id, system.protein_ids) for system in systems] == [
        ("genome-1", ["p1", "p2"]),
        ("genome-2", ["p3", "p4"]),
    ]


def test_multi_locus_quorum_can_combine_loci_on_one_contig() -> None:
    model = _model(
        *(GeneSpec(name, "mandatory") for name in ("a", "b", "c", "d")),
        gap=1,
        multi_loci=True,
    )
    hits = [
        _hit("p1", "a", 1),
        _hit("p2", "b", 2),
        _hit("p3", "c", 20),
        _hit("p4", "d", 21),
    ]

    systems = _detect_model_systems(hits, model)

    assert len(systems) == 1
    assert systems[0].contig_ids == ["contig-1"]
    assert systems[0].protein_ids == ["p1", "p2", "p3", "p4"]


def test_two_complete_loci_remain_two_system_instances() -> None:
    model = _model(
        GeneSpec("a", "mandatory"),
        GeneSpec("b", "mandatory"),
        gap=1,
        multi_loci=True,
    )
    hits = [
        _hit("p1", "a", 1),
        _hit("p2", "b", 2),
        _hit("p3", "a", 20),
        _hit("p4", "b", 21),
    ]

    selected = resolve_conflicts(_detect_model_systems(hits, model))

    assert [system.protein_ids for system in selected] == [
        ["p1", "p2"],
        ["p3", "p4"],
    ]


def test_loner_can_rescue_only_on_the_same_contig() -> None:
    model = _model(
        GeneSpec("a", "mandatory"),
        GeneSpec("b", "mandatory"),
        GeneSpec("c", "mandatory", loner=True),
        gap=1,
    )
    split_hits = [
        _hit("p1", "a", 1, contig="contig-1"),
        _hit("p2", "b", 2, contig="contig-1"),
        _hit("p3", "c", 50, contig="contig-2"),
    ]
    same_contig_hits = [
        _hit("p1", "a", 1),
        _hit("p2", "b", 2),
        _hit("p3", "c", 50),
    ]

    assert _detect_model_systems(split_hits, model) == []
    systems = _detect_model_systems(same_contig_hits, model)
    assert len(systems) == 1
    assert systems[0].gene_presences == ["mandatory", "mandatory", "loner"]


def test_forbidden_hit_in_cluster_rejects_candidate() -> None:
    model = _model(
        GeneSpec("a", "mandatory"),
        GeneSpec("b", "mandatory"),
        GeneSpec("x", "forbidden"),
        min_mandatory=2,
        min_genes=2,
    )
    hits = [
        _hit("p1", "a", 1),
        _hit("p2", "x", 2),
        _hit("p3", "b", 3),
    ]

    assert _detect_model_systems(hits, model) == []


def test_multi_system_second_pass_can_rescue_an_independent_instance() -> None:
    model = _model(
        GeneSpec("shared", "mandatory", multi_system=True),
        GeneSpec("left", "mandatory"),
        GeneSpec("right", "mandatory"),
        GeneSpec("neutral", "neutral"),
        gap=1,
        min_mandatory=2,
        min_genes=2,
    )
    hits = [
        _hit("p1", "shared", 1),
        _hit("p2", "left", 2),
        _hit("p3", "right", 20),
        _hit("p4", "neutral", 21),
    ]

    selected = resolve_conflicts(_detect_model_systems(hits, model))

    assert [system.protein_ids for system in selected] == [
        ["p1", "p2"],
        ["p1", "p3", "p4"],
    ]
    assert selected[1].rescued_multi_system_proteins == {"p1"}


def test_ambiguous_hmm_internal_names_are_quarantined(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles"
    profiles.mkdir()
    (profiles / "model-a.hmm").write_text("HMMER3/f\nNAME shared-name\n")
    (profiles / "model-b.hmm").write_text("HMMER3/f\nNAME shared-name\n")
    hits = [_hit("p1", "shared-name", 1)]

    mapping = _build_hmm_name_mapping(profiles)
    _translate_hit_accessions(hits, mapping)

    assert mapping["shared-name"].startswith("__sharur_ambiguous_hmm_name__:")
    assert hits[0].accession not in {"model-a", "model-b"}


def test_validation_fails_closed_when_model_family_is_unavailable(
    tmp_path: Path,
) -> None:
    with pytest.raises(RuntimeError, match="No parseable models"):
        validate_systems(
            db_path=tmp_path / "unused.duckdb",
            models_dir=tmp_path,
            model_family="missing-model-family",
            verbose=False,
        )


def test_same_function_loner_stretch_keeps_one_best_hit() -> None:
    model = _model(
        GeneSpec("a", "mandatory", loner=True),
        min_mandatory=1,
        min_genes=1,
    )
    hits = [
        _hit("p1", "a", 1, score=50.0),
        _hit("p2", "a", 2, score=100.0),
    ]

    regular, loners = _build_model_clusters(hits, model)

    assert regular == []
    assert [hit.protein_id for hit in loners] == ["p2"]


def test_integration_migrates_contig_id_and_replaces_stale_rows(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)
    systems, genes = _system_frames()

    integrate_secretion_results(db_path, systems, genes)

    conn = duckdb.connect(str(db_path), read_only=True)
    columns = {
        row[0] for row in conn.execute("DESCRIBE secretion_systems").fetchall()
    }
    stored = conn.execute(
        "SELECT system_id, genome_id, contig_id, protein_ids FROM secretion_systems"
    ).fetchall()
    validated = conn.execute(
        "SELECT protein_id, accession FROM annotations "
        "WHERE source = 'txsscan_system' ORDER BY protein_id"
    ).fetchall()
    normalized = conn.execute(
        "SELECT system_id, protein_id, system_source, position, profile_name "
        "FROM system_proteins ORDER BY position"
    ).fetchall()
    conn.close()

    assert "contig_id" in columns
    assert stored == [
        ("genome-1_test_1", "genome-1", "contig-1", "p1,p2")
    ]
    assert validated == [("p1", "a"), ("p2", "b")]
    assert normalized == [
        ("genome-1_test_1", "p1", "txsscan_system", 1, "a"),
        ("genome-1_test_1", "p2", "txsscan_system", 2, "b"),
    ]


def test_integration_invalidates_derived_semantics_for_old_and_new_members(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)
    conn = duckdb.connect(str(db_path))
    conn.execute("""
        CREATE TABLE semantic_atoms (
            protein_id VARCHAR,
            source_db VARCHAR
        );
        INSERT INTO semantic_atoms VALUES
            ('p1', 'txsscan_system'),
            ('p3', 'pfam');
        CREATE TABLE semantic_state (protein_id VARCHAR);
        INSERT INTO semantic_state VALUES ('p1'), ('p3');
        CREATE TABLE semantic_terms (
            protein_id VARCHAR,
            source_db VARCHAR
        );
        INSERT INTO semantic_terms VALUES
            ('p1', 'txsscan_system'),
            ('p3', 'pfam');
        CREATE TABLE protein_predicates (protein_id VARCHAR);
        INSERT INTO protein_predicates VALUES ('p1'), ('p3');
    """)
    conn.close()
    systems, genes = _system_frames()

    affected = integrate_secretion_results(db_path, systems, genes)

    conn = duckdb.connect(str(db_path), read_only=True)
    atoms = conn.execute(
        "SELECT protein_id, source_db FROM semantic_atoms ORDER BY protein_id"
    ).fetchall()
    states = conn.execute(
        "SELECT protein_id FROM semantic_state ORDER BY protein_id"
    ).fetchall()
    terms = conn.execute(
        "SELECT protein_id, source_db FROM semantic_terms ORDER BY protein_id"
    ).fetchall()
    predicates = conn.execute(
        "SELECT protein_id FROM protein_predicates ORDER BY protein_id"
    ).fetchall()
    conn.close()

    assert affected == {"p1", "p2"}
    assert atoms == [("p3", "pfam")]
    assert states == [("p3",)]
    assert terms == [("p3", "pfam")]
    assert predicates == [("p3",)]


def test_integration_returns_every_member_quarantined_by_v5_migration(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)
    conn = duckdb.connect(str(db_path))
    conn.execute("""
        CREATE TABLE schema_version (
            version INTEGER PRIMARY KEY,
            applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            description TEXT
        );
        INSERT INTO schema_version (version, description)
            VALUES (4, 'legacy caller state');
        CREATE TABLE defense_systems (
            system_id VARCHAR PRIMARY KEY,
            genome_id VARCHAR,
            system_type VARCHAR,
            system_subtype VARCHAR,
            protein_ids VARCHAR,
            profile_names VARCHAR
        );
        INSERT INTO defense_systems VALUES
            ('legacy-defense', 'genome-1', 'legacy', 'legacy', 'p3', 'd');
        INSERT INTO annotations VALUES
            (3, 'p3', 'defensefinder_system', 'd', 'legacy', '',
             NULL, 1, NULL, NULL);
    """)
    conn.close()
    systems, genes = _system_frames()

    affected = integrate_secretion_results(db_path, systems, genes)

    conn = duckdb.connect(str(db_path), read_only=True)
    defense_count = conn.execute(
        "SELECT COUNT(*) FROM defense_systems"
    ).fetchone()[0]
    current_sources = conn.execute(
        "SELECT source, COUNT(*) FROM annotations GROUP BY source ORDER BY source"
    ).fetchall()
    conn.close()

    assert affected == {"p1", "p2", "p3"}
    assert defense_count == 0
    assert current_sources == [
        ("txsscan", 1),
        ("txsscan_system", 2),
    ]


def test_empty_integration_clears_stale_structured_calls(tmp_path: Path) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)

    integrate_secretion_results(db_path, pd.DataFrame(), pd.DataFrame())

    conn = duckdb.connect(str(db_path), read_only=True)
    system_count = conn.execute("SELECT COUNT(*) FROM secretion_systems").fetchone()[0]
    annotation_count = conn.execute(
        "SELECT COUNT(*) FROM annotations WHERE source = 'txsscan_system'"
    ).fetchone()[0]
    normalized_count = conn.execute(
        "SELECT COUNT(*) FROM system_proteins "
        "WHERE system_source = 'txsscan_system'"
    ).fetchone()[0]
    conn.close()

    assert system_count == 0
    assert annotation_count == 0
    assert normalized_count == 0


def test_integration_replaces_an_existing_primary_key_atomically(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)
    conn = duckdb.connect(str(db_path))
    conn.execute(
        "UPDATE secretion_systems "
        "SET system_id = 'genome-1_test_1' WHERE system_id = 'stale'"
    )
    conn.close()
    systems, genes = _system_frames()

    integrate_secretion_results(db_path, systems, genes)

    conn = duckdb.connect(str(db_path), read_only=True)
    stored = conn.execute(
        "SELECT system_id, contig_id, protein_ids FROM secretion_systems"
    ).fetchall()
    conn.close()
    assert stored == [("genome-1_test_1", "contig-1", "p1,p2")]


def test_integration_rejects_cross_contig_membership_before_mutation(
    tmp_path: Path,
) -> None:
    db_path = tmp_path / "systems.duckdb"
    _seed_integration_db(db_path)
    systems, genes = _system_frames()
    systems.loc[0, "protein_ids"] = "p1,p3"
    genes.loc[1, "protein_id"] = "p3"

    with pytest.raises(ValueError, match="protein/genome/contig"):
        integrate_secretion_results(db_path, systems, genes)

    conn = duckdb.connect(str(db_path), read_only=True)
    stored_systems = conn.execute(
        "SELECT system_id FROM secretion_systems"
    ).fetchall()
    stored_annotations = conn.execute(
        "SELECT accession FROM annotations WHERE source = 'txsscan_system'"
    ).fetchall()
    conn.close()

    assert stored_systems == [("stale",)]
    assert stored_annotations == [("stale",)]
