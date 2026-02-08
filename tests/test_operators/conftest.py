"""
Pytest fixtures for operator tests.

Provides seeded in-memory DuckDB databases for testing.
"""

import pytest

from bennu.storage.duckdb_store import DuckDBStore
from bennu.operators import Bennu


@pytest.fixture
def store():
    """Create an in-memory DuckDB store with test data."""
    store = DuckDBStore(db_path=None)

    # Insert test bins (genomes)
    store.execute(
        """
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length)
        VALUES
            ('bin_001', 95.5, 2.1, 'd__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales', 15, 2500000),
            ('bin_002', 88.2, 4.5, 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales', 28, 4800000),
            ('bin_003', 72.0, 8.2, 'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales', 45, 3200000)
        """
    )

    # Insert test contigs
    store.execute(
        """
        INSERT INTO contigs (contig_id, bin_id, length, gc_content, is_circular, taxonomy)
        VALUES
            ('contig_001', 'bin_001', 250000, 0.45, FALSE, 'd__Archaea;p__Euryarchaeota'),
            ('contig_002', 'bin_001', 180000, 0.42, FALSE, 'd__Archaea;p__Euryarchaeota'),
            ('contig_003', 'bin_002', 500000, 0.52, FALSE, 'd__Bacteria;p__Proteobacteria'),
            ('contig_004', 'bin_003', 200000, 0.38, FALSE, 'd__Bacteria;p__Firmicutes')
        """
    )

    # Insert test proteins
    store.execute(
        """
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length, gc_content)
        VALUES
            -- bin_001, contig_001: 5 proteins
            ('prot_001', 'contig_001', 'bin_001', 1000, 2500, '+', 1, 500, 0.45),
            ('prot_002', 'contig_001', 'bin_001', 3000, 4200, '+', 2, 400, 0.44),
            ('prot_003', 'contig_001', 'bin_001', 4500, 5100, '-', 3, 200, 0.46),
            ('prot_004', 'contig_001', 'bin_001', 5500, 12500, '+', 4, 2333, 0.43),  -- giant
            ('prot_005', 'contig_001', 'bin_001', 13000, 13120, '+', 5, 40, 0.42),   -- tiny

            -- bin_002, contig_003: 3 proteins
            ('prot_006', 'contig_003', 'bin_002', 1000, 2800, '+', 1, 600, 0.52),
            ('prot_007', 'contig_003', 'bin_002', 3200, 19200, '-', 2, 5333, 0.51),  -- massive
            ('prot_008', 'contig_003', 'bin_002', 20000, 21500, '+', 3, 500, 0.53),

            -- bin_003, contig_004: 2 proteins
            ('prot_009', 'contig_004', 'bin_003', 500, 1700, '+', 1, 400, 0.38),
            ('prot_010', 'contig_004', 'bin_003', 2000, 3500, '-', 2, 500, 0.37)
        """
    )

    # Insert test annotations
    store.execute(
        """
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score)
        VALUES
            (1, 'prot_001', 'pfam', 'PF00142', 'NiFe-hydrogenase', 'Nickel-iron hydrogenase large subunit', 1e-50, 180.5),
            (2, 'prot_002', 'pfam', 'PF01512', 'ATP_synthase', 'ATP synthase subunit alpha', 1e-45, 165.2),
            (3, 'prot_002', 'kegg', 'K02111', 'atpA', 'ATP synthase F1 subunit alpha', 1e-40, 155.0),
            (4, 'prot_003', 'pfam', 'PF00001', 'hypothetical', 'Hypothetical protein with unknown function', 1e-5, 45.0),
            (5, 'prot_006', 'pfam', 'PF00072', 'Response_reg', 'Response regulator receiver domain', 1e-35, 120.0),
            (6, 'prot_008', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-30, 110.0),
            (7, 'prot_009', 'pfam', 'PF00009', 'GTP_EFTU', 'Elongation factor Tu GTP binding domain', 1e-55, 190.0),
            (8, 'prot_009', 'pfam', 'PF03144', 'GTP_EFTU_D2', 'Elongation factor Tu domain 2', 1e-40, 150.0),
            (9, 'prot_009', 'pfam', 'PF00009', 'GTP_EFTU_D3', 'Elongation factor Tu domain 3', 1e-35, 130.0)
        """
    )

    return store


@pytest.fixture
def bennu(store):
    """Create Bennu instance with test store."""
    b = Bennu()
    # Force session creation and inject our test store
    _ = b.session  # Trigger lazy creation
    b._session._db = store
    return b
