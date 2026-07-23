"""Tests for Sharur facade class."""

from sharur.operators import Sharur, SharurResult


class TestSharurFacade:
    """Tests for Sharur facade class."""

    def test_init_with_path(self, tmp_path):
        """Sharur should accept db_path."""
        db_file = tmp_path / "test.duckdb"
        b = Sharur(db_file)
        assert b._db_path == db_file

    def test_init_read_only(self, tmp_path):
        """Sharur should pass read-only mode through to the store."""
        from sharur.storage.duckdb_store import DuckDBStore

        db_file = tmp_path / "test.duckdb"
        store = DuckDBStore(db_file)
        store.execute("INSERT INTO bins (bin_id) VALUES ('bin1')")
        store.conn.close()

        b = Sharur(db_file, read_only=True)

        assert b._read_only is True
        assert b.store.read_only is True
        assert b.store.execute("SELECT COUNT(*) FROM bins")[0][0] == 1

    def test_init_without_path(self):
        """Sharur should work without db_path (in-memory)."""
        b = Sharur()
        assert b._db_path is None

    def test_overview(self, sharur):
        """Sharur.overview() should return SharurResult."""
        result = sharur.overview()
        assert isinstance(result, SharurResult)
        assert result.meta.rows > 0

    def test_list_genomes(self, sharur):
        """Sharur.list_genomes() should return SharurResult."""
        result = sharur.list_genomes()
        assert isinstance(result, SharurResult)
        assert result.meta.rows == 3

    def test_list_genomes_with_filter(self, sharur):
        """Sharur.list_genomes() should accept filters."""
        result = sharur.list_genomes(taxonomy_filter="Archaea")
        assert result.meta.rows == 1

    def test_list_proteins(self, sharur):
        """Sharur.list_proteins() should return SharurResult."""
        result = sharur.list_proteins()
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows == 10

    def test_list_proteins_with_filter(self, sharur):
        """Sharur.list_proteins() should accept filters."""
        result = sharur.list_proteins(genome_id="bin_001")
        assert result.meta.total_rows == 5

    def test_get_genome(self, sharur):
        """Sharur.get_genome() should return SharurResult."""
        result = sharur.get_genome("bin_001")
        assert isinstance(result, SharurResult)
        assert "bin_001" in result.data

    def test_get_protein(self, sharur):
        """Sharur.get_protein() should return SharurResult."""
        result = sharur.get_protein("prot_001")
        assert isinstance(result, SharurResult)
        assert "prot_001" in result.data

    def test_get_neighborhood(self, sharur):
        """Sharur.get_neighborhood() should return SharurResult."""
        result = sharur.get_neighborhood("prot_003", window=2)
        assert isinstance(result, SharurResult)
        assert result.meta.rows > 0

    def test_search_by_predicates(self, sharur):
        """Sharur.search_by_predicates() should return SharurResult."""
        result = sharur.search_by_predicates(has=["giant"])
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows >= 1

    def test_search_proteins(self, sharur):
        """Sharur.search_proteins() should return SharurResult."""
        result = sharur.search_proteins(annotation_pattern="hydrogenase")
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows >= 1

    def test_v2_search_by_atoms_direct_access(self, sharur):
        """V2 atom search should include V1-style direct accession keys."""
        sharur.generate_v2(
            protein_ids=["prot_008"],
            chunk_size=1,
            pipeline_depth=1,
        )

        result = sharur.search_by_atoms(has=["pfam:PF00005"])

        assert result == ["prot_008"]

    def test_v2_search_atoms_filters_source_and_relation(self, sharur):
        """Atom evidence search should expose relation/source filters."""
        sharur.generate_v2(protein_ids=["prot_008"], chunk_size=1)

        rows = sharur.search_atoms(
            atom_id="abc_transporter",
            relation="supports",
            source_db="pfam",
            source_accession="PF00005",
        )

        assert rows == [
            {
                "protein_id": "prot_008",
                "atom_id": "abc_transporter",
                "facet": "role",
                "relation": "supports",
                "source_db": "pfam",
                "source_accession": "PF00005",
            }
        ]

    def test_explain_v2_protein(self, sharur):
        """explain() should expose resolved state, terms, and raw atoms."""
        sharur.generate_v2(protein_ids=["prot_008"], chunk_size=1)

        explanation = sharur.explain("prot_008")

        assert explanation["protein_id"] == "prot_008"
        assert explanation["protein"]["bin_id"] == "bin_002"
        assert explanation["semantic_state"]["protein_id"] == "prot_008"
        assert "abc_transporter" in explanation["semantic_state"]["roles"]
        assert "pfam:PF00005" in explanation["direct_access_terms"]
        assert "abc_transporter_complete" in explanation["composite_terms"]
        assert [
            witness["atom_id"]
            for witness in explanation["composite_explanations"][
                "abc_transporter_complete"
            ]
        ] == ["abc_transporter", "atp_binding"]
        assert any(atom["atom_id"] == "abc_transporter" for atom in explanation["atoms"])
        assert explanation["conflicts"] == []

    def test_explain_includes_validated_system_membership(self, sharur):
        """explain() should include normalized validated-system membership."""
        sharur.store.execute("""
            INSERT INTO defense_systems (
                system_id, genome_id, system_type, system_subtype, activity,
                genes_count, protein_ids, profile_names
            )
            VALUES (
                'sys_def_explain', 'bin_001', 'RM_Type_II', 'RM_Type_II',
                'Defense', 1, 'prot_002', 'HsdR'
            );
        """)
        sharur.generate_v2(protein_ids=["prot_002"], chunk_size=1)

        explanation = sharur.explain("prot_002")

        assert explanation["validated_systems"] == [
            {
                "system_id": "sys_def_explain",
                "system_source": "defensefinder_system",
                "position": 1,
                "profile_name": "HsdR",
                "score": 1.0,
            }
        ]
        assert "defense_system" in explanation["semantic_state"]["roles"]

    def test_explain_degrades_without_v2_tables(self, sharur):
        """explain() should still inspect legacy-only databases."""
        sharur.store.execute("""
            INSERT OR REPLACE INTO protein_predicates (protein_id, predicates)
            VALUES ('prot_008', ['abc_transporter']);
        """)
        for table in [
            "semantic_terms",
            "semantic_state",
            "semantic_atoms",
            "system_proteins",
        ]:
            sharur.store.execute(f"DROP TABLE IF EXISTS {table};")

        explanation = sharur.explain("prot_008")

        assert explanation["protein"]["protein_id"] == "prot_008"
        assert explanation["semantic_state"] is None
        assert explanation["atoms"] == []
        assert explanation["terms"] == []
        assert explanation["composite_explanations"] == {}
        assert explanation["validated_systems"] == []
        assert explanation["legacy_predicates"] == ["abc_transporter"]

    def test_v2_review_queue_filters_and_writes_tsv(self, sharur, tmp_path):
        """V2 review queue should support filters without loading all atoms."""
        sharur.store.execute("""
            INSERT INTO semantic_atoms (
                protein_id, atom_id, facet, relation, source_accession,
                source_db
            )
            VALUES
                (
                    'prot_001', 'unresolved:PF_HIGH', 'quality_flag',
                    'unresolved', 'PF_HIGH', 'pfam'
                ),
                (
                    'prot_006', 'unresolved:PF_HIGH', 'quality_flag',
                    'unresolved', 'PF_HIGH', 'pfam'
                ),
                (
                    'prot_008', 'unresolved:PF_LOW', 'quality_flag',
                    'unresolved', 'PF_LOW', 'pfam'
                ),
                (
                    'prot_002', 'unresolved:TXS_RAW', 'quality_flag',
                    'unresolved', 'TXS_RAW', 'txsscan'
                ),
                (
                    'prot_003', 'unresolved:SYS_VALIDATED', 'quality_flag',
                    'unresolved', 'SYS_VALIDATED', 'txsscan_system'
                );
        """)

        output = tmp_path / "queue.tsv"
        queue = sharur.v2_review_queue(
            limit=10,
            source=["pfam", "txsscan", "txsscan_system"],
            min_proteins=1,
            exclude_raw_system_profiles=True,
            output_tsv=output,
        )

        by_accession = {row["accession"]: row for row in queue}
        assert "TXS_RAW" not in by_accession
        assert by_accession["PF_HIGH"]["n_proteins"] == 2
        assert by_accession["PF_HIGH"]["n_genomes"] == 2
        assert by_accession["SYS_VALIDATED"]["source_db"] == "txsscan_system"
        assert output.read_text().startswith("accession\tsource_db")


class TestSharurResultStr:
    """Tests for SharurResult string representation."""

    def test_str_returns_data(self, sharur):
        """str(SharurResult) should return data."""
        result = sharur.overview()
        assert str(result) == result.data

    def test_repr(self, sharur):
        """repr(SharurResult) should be informative."""
        result = sharur.list_genomes()
        repr_str = repr(result)
        assert "SharurResult" in repr_str
        assert "rows=" in repr_str
