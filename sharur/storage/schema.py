"""DuckDB schema for Sharur."""

SCHEMA_VERSION = 1

SCHEMA = """
CREATE TABLE IF NOT EXISTS bins (
    bin_id VARCHAR PRIMARY KEY,
    completeness FLOAT,
    contamination FLOAT,
    taxonomy VARCHAR,
    n_contigs INTEGER,
    total_length INTEGER
);

CREATE TABLE IF NOT EXISTS contigs (
    contig_id VARCHAR PRIMARY KEY,
    bin_id VARCHAR,
    length INTEGER NOT NULL,
    gc_content FLOAT,
    is_circular BOOLEAN DEFAULT FALSE,
    taxonomy VARCHAR,
    
    FOREIGN KEY (bin_id) REFERENCES bins(bin_id)
);

CREATE TABLE IF NOT EXISTS proteins (
    protein_id VARCHAR PRIMARY KEY,
    contig_id VARCHAR NOT NULL,
    bin_id VARCHAR,
    start INTEGER NOT NULL,
    end_coord INTEGER NOT NULL,  -- 'end' is reserved
    strand VARCHAR(1) NOT NULL,
    gene_index INTEGER,
    sequence TEXT,
    sequence_length INTEGER,
    gc_content FLOAT,
    
    -- Indexes
    FOREIGN KEY (contig_id) REFERENCES contigs(contig_id)
);
CREATE INDEX IF NOT EXISTS idx_proteins_contig ON proteins(contig_id);
CREATE INDEX IF NOT EXISTS idx_proteins_bin ON proteins(bin_id);
CREATE INDEX IF NOT EXISTS idx_proteins_coords ON proteins(contig_id, start, end_coord);

CREATE TABLE IF NOT EXISTS annotations (
    annotation_id INTEGER PRIMARY KEY,
    protein_id VARCHAR NOT NULL,
    source VARCHAR NOT NULL,  -- 'pfam', 'kegg', 'cazy', etc.
    accession VARCHAR NOT NULL,
    name VARCHAR,
    description TEXT,
    evalue FLOAT,
    score FLOAT,
    start_aa INTEGER,
    end_aa INTEGER,
    
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);
CREATE INDEX IF NOT EXISTS idx_annotations_protein ON annotations(protein_id);
CREATE INDEX IF NOT EXISTS idx_annotations_source_acc ON annotations(source, accession);

CREATE TABLE IF NOT EXISTS loci (
    locus_id VARCHAR PRIMARY KEY,
    locus_type VARCHAR NOT NULL,
    contig_id VARCHAR NOT NULL,
    start INTEGER NOT NULL,
    end_coord INTEGER NOT NULL,
    confidence FLOAT,
    metadata JSON,
    
    FOREIGN KEY (contig_id) REFERENCES contigs(contig_id)
);

CREATE TABLE IF NOT EXISTS locus_proteins (
    locus_id VARCHAR NOT NULL,
    protein_id VARCHAR NOT NULL,
    position INTEGER NOT NULL,  -- Order within locus
    
    PRIMARY KEY (locus_id, protein_id),
    FOREIGN KEY (locus_id) REFERENCES loci(locus_id),
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);

-- Feature store
CREATE TABLE IF NOT EXISTS feature_store (
    protein_id VARCHAR NOT NULL,
    metric_name VARCHAR NOT NULL,
    metric_value FLOAT NOT NULL,
    computed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    version VARCHAR DEFAULT '1.0',
    
    PRIMARY KEY (protein_id, metric_name)
);
CREATE INDEX IF NOT EXISTS idx_features_metric ON feature_store(metric_name);

-- Helpful secondary indexes
CREATE INDEX IF NOT EXISTS idx_proteins_contig ON proteins(contig_id);
CREATE INDEX IF NOT EXISTS idx_proteins_bin ON proteins(bin_id);
CREATE INDEX IF NOT EXISTS idx_proteins_coords ON proteins(contig_id, start, end_coord);
CREATE INDEX IF NOT EXISTS idx_annotations_protein ON annotations(protein_id);
CREATE INDEX IF NOT EXISTS idx_annotations_source_acc ON annotations(source, accession);
CREATE INDEX IF NOT EXISTS idx_annotations_accession ON annotations(accession);
CREATE INDEX IF NOT EXISTS idx_loci_contig ON loci(contig_id, start, end_coord);
CREATE INDEX IF NOT EXISTS idx_loci_type ON loci(locus_type);

-- Predicate definitions (operator system)
CREATE TABLE IF NOT EXISTS predicate_definitions (
    predicate_id VARCHAR PRIMARY KEY,
    name VARCHAR NOT NULL,
    description TEXT NOT NULL,
    category VARCHAR NOT NULL,
    entity_type VARCHAR NOT NULL,
    eval_query TEXT NOT NULL,
    window_bp INTEGER,
    ref_urls VARCHAR[],
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Pre-computed predicates per protein (DuckDB array syntax)
CREATE TABLE IF NOT EXISTS protein_predicates (
    protein_id VARCHAR PRIMARY KEY,
    predicates VARCHAR[] NOT NULL DEFAULT [],
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);

-- Interest scores for ranking
CREATE TABLE IF NOT EXISTS protein_scores (
    protein_id VARCHAR PRIMARY KEY,
    interest_score FLOAT NOT NULL DEFAULT 0,
    size_zscore FLOAT,
    annotation_gap FLOAT,
    gc_deviation FLOAT,
    computed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);
CREATE INDEX IF NOT EXISTS idx_protein_scores_interest ON protein_scores(interest_score DESC);

-- Refs for expand() pagination
CREATE TABLE IF NOT EXISTS refs (
    ref_id VARCHAR PRIMARY KEY,
    query_hash VARCHAR NOT NULL,
    query_params JSON NOT NULL,
    result_count INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP
);
"""

__all__ = ["SCHEMA", "SCHEMA_VERSION"]
