"""The co-location caller's DefenseFinder scope: sub-families run and component quorum."""

from __future__ import annotations

import duckdb
import pytest

from sharur.colocation import validate_systems
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
