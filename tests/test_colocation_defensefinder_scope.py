"""The co-location caller's DefenseFinder scope: sub-families run and component quorum."""

from __future__ import annotations

import duckdb
import pytest

from sharur.colocation import MACSYFINDER_EXCHANGEABLE_WEIGHT, _model_passes, validate_systems
from tests.defensefinder_fixtures import build_toy_db, write_model_set


@pytest.fixture
def family(tmp_path):
    return write_model_set(tmp_path / "models")


@pytest.fixture
def db_path(tmp_path, family):
    return build_toy_db(tmp_path / "sharur.duckdb", family, call_systems=False, semantics=False)


class TestCallerScope:
    def test_antidefense_models_are_not_called_as_defense_systems(self, db_path, family):
        conn = duckdb.connect(str(db_path))
        conn.execute("DELETE FROM annotations WHERE source = 'defensefinder'")
        conn.execute(
            "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, score) "
            "VALUES (10, 'p4', 'defensefinder', 'Anti__AntiX', 'Anti__AntiX', 80.0), "
            "(11, 'p3', 'defensefinder', 'Solo__SoloA', 'Solo__SoloA', 60.0)"
        )
        conn.close()
        default, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert set(default["type"]) == {"Solo"}
        every, _ = validate_systems(db_path, models_dir=family.parent, verbose=False,
                                    subfamilies=("DefenseFinder", "RM", "ADF"))
        assert set(every["type"]) == {"Solo", "AntiX"}

    def test_two_component_model_requires_both_mandatory_genes(self, db_path, family):
        conn = duckdb.connect(str(db_path))
        conn.execute("DELETE FROM annotations WHERE source = 'defensefinder'")
        conn.execute(
            "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, score) "
            "VALUES (10, 'p7', 'defensefinder', 'Toy__ToyA', 'Toy__ToyA', 50.0)"
        )
        conn.close()
        systems, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert systems.empty
        conn = duckdb.connect(str(db_path))
        conn.execute(
            "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, score) "
            "VALUES (11, 'p1', 'defensefinder', 'Toy__ToyA', 'Toy__ToyA', 40.0), "
            "(12, 'p2', 'defensefinder', 'ToyB_internal', 'ToyB_internal', 35.0)"
        )
        conn.close()
        systems, genes = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert list(systems["type"]) == ["Toy"]
        assert sorted(genes["protein_id"]) == ["p1", "p2"]


def _replace_hits(db_path, *rows):
    """Replace the defensefinder rows with (annotation_id, protein, profile NAME, score)."""
    conn = duckdb.connect(str(db_path))
    conn.execute("DELETE FROM annotations WHERE source = 'defensefinder'")
    conn.executemany(
        "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, score, evalue) "
        "VALUES (?, ?, 'defensefinder', ?, ?, ?, 1e-20)",
        [(aid, protein, profile, profile, score) for aid, protein, profile, score in rows],
    )
    conn.close()


class TestSubfamilyPasses:
    """DefenseFinder runs one MacSyFinder invocation per sub-family (DefenseFinder, RM):
    each searches only its own models' profiles, keeps the best hit per protein among
    those and picks its own solution; the post-treatment concatenates the runs."""

    def test_rm_hit_does_not_displace_a_defensefinder_hit_on_the_same_protein(self, db_path, family):
        _replace_hits(db_path, (10, "p3", "Solo__SoloA", 60.0))
        alone, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
        _replace_hits(db_path, (10, "p3", "Solo__SoloA", 60.0),
                      (11, "p3", "RM_Type_X__MTases_FAM_1", 90.0))
        both, genes = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert sorted(both["type"]) == ["RM_Type_X", "Solo"]
        assert sorted(genes["protein_id"]) == ["p3", "p3"]
        # The Solo call keeps the system ID it has without the RM hit.
        assert both.loc[both["type"] == "Solo", "system_id"].tolist() == alone["system_id"].tolist()

    def test_profiles_of_uncalled_models_do_not_compete(self, db_path, family):
        # Anti__AntiX belongs only to an ADF model, which DefenseFinder does not run by default.
        _replace_hits(db_path, (10, "p3", "Solo__SoloA", 60.0), (11, "p3", "Anti__AntiX", 80.0))
        systems, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert list(systems["type"]) == ["Solo"]

    def test_best_hit_competition_still_applies_within_a_pass(self, db_path, family):
        # Toy__ToyA outscores Solo__SoloA on p3 in the DefenseFinder pass; a lone ToyA
        # cannot form a Toy system, so nothing is called.
        _replace_hits(db_path, (10, "p3", "Solo__SoloA", 60.0), (11, "p3", "Toy__ToyA", 80.0))
        systems, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
        assert systems.empty


# Two DefenseFinder models claim the protein carrying Shared__S. ExchA reaches its other
# three mandatory genes only through exchangeable profiles: score 1 + 3w. PlainB scores
# 3 mandatory + 1 accessory = 3.5. MacSyFinder's default w = 0.8 gives 3.4 (PlainB wins);
# DefenseFinder passes --exchangeable-weight 1, giving 4.0 (ExchA wins).
_EXCH_PROFILES = {stem: (stem, 100, (20.0, 20.0)) for stem in (
    "Shared__S", "ExchA__A1", "ExchA__A1x", "ExchA__A2", "ExchA__A2x", "ExchA__A3", "ExchA__A3x",
    "PlainB__B1", "PlainB__B2", "PlainB__B3",
)}
_EXCH_MODELS = {
    ("DefenseFinder", "ExchA", "ExchA"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="4" min_genes_required="4" vers="2.0">\n'
        '<gene name="Shared__S" presence="mandatory"/>\n'
        + "".join(
            f'<gene name="ExchA__A{i}" presence="mandatory">'
            f'<exchangeables><gene name="ExchA__A{i}x"/></exchangeables></gene>\n'
            for i in (1, 2, 3)
        )
        + "</model>\n"
    ),
    ("DefenseFinder", "PlainB", "PlainB"): (
        '<model inter_gene_max_space="5" min_mandatory_genes_required="3" min_genes_required="3" vers="2.0">\n'
        '<gene name="Shared__S" presence="mandatory"/>\n'
        '<gene name="PlainB__B1" presence="mandatory"/>\n'
        '<gene name="PlainB__B2" presence="mandatory"/>\n'
        '<gene name="PlainB__B3" presence="accessory"/>\n'
        "</model>\n"
    ),
}


def test_defensefinder_scores_exchangeable_hits_at_full_weight(tmp_path):
    family = write_model_set(tmp_path / "models", extra_profiles=_EXCH_PROFILES,
                             extra_models=_EXCH_MODELS)
    db_path = build_toy_db(tmp_path / "sharur.duckdb", family, call_systems=False, semantics=False)
    conn = duckdb.connect(str(db_path))
    conn.execute("DELETE FROM annotations WHERE source = 'defensefinder'")
    conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c2', 'bin_1', 50000)")
    layout = [("x1", "ExchA__A1x"), ("x2", "ExchA__A2x"), ("x3", "ExchA__A3x"), ("x4", "Shared__S"),
              ("x5", "PlainB__B1"), ("x6", "PlainB__B2"), ("x7", "PlainB__B3")]
    for index, (protein, profile) in enumerate(layout, 1):
        conn.execute(
            "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
            "sequence_length) VALUES (?, 'c2', 'bin_1', ?, ?, '+', ?, 300)",
            [protein, index * 1000, index * 1000 + 900, index],
        )
        conn.execute(
            "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, score, evalue) "
            "VALUES (?, ?, 'defensefinder', ?, ?, 50.0, 1e-20)",
            [100 + index, protein, profile, profile],
        )
    conn.close()

    systems, _ = validate_systems(db_path, models_dir=family.parent, verbose=False)
    assert list(systems["type"]) == ["ExchA"]
    # Other callers keep MacSyFinder's default.
    (_name, models), *_ = _model_passes(family.parent, family.name, None, separate=False)
    assert {model.exchangeable_weight for model in models.values()} == {MACSYFINDER_EXCHANGEABLE_WEIGHT}
