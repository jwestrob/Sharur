"""DuckDB schema for Sharur."""

SCHEMA_VERSION = 6

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
CREATE INDEX IF NOT EXISTS idx_contigs_bin ON contigs(bin_id, contig_id);

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

-- Predicate V2 semantic atom storage
CREATE TABLE IF NOT EXISTS semantic_atoms (
    protein_id VARCHAR NOT NULL,
    atom_id VARCHAR NOT NULL,
    facet VARCHAR NOT NULL,
    relation VARCHAR NOT NULL,
    source_accession VARCHAR,
    source_db VARCHAR,
    evidence_evalue DOUBLE,
    evidence_score DOUBLE,

    PRIMARY KEY (protein_id, atom_id, source_accession)
);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_protein ON semantic_atoms(protein_id);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_atom ON semantic_atoms(atom_id);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_facet ON semantic_atoms(facet);
CREATE INDEX IF NOT EXISTS idx_semantic_atoms_source ON semantic_atoms(source_db, source_accession);

-- Predicate V2 per-protein resolved state
CREATE TABLE IF NOT EXISTS semantic_state (
    protein_id VARCHAR PRIMARY KEY,
    activities VARCHAR[],
    roles VARCHAR[],
    architecture VARCHAR[],
    localization VARCHAR[],
    topology JSON,
    size_class VARCHAR,
    quality_flags VARCHAR[],
    composite_predicates VARCHAR[],
    unresolved_count INTEGER,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
-- Predicate V2 unified search terms. This materializes atom IDs, V1-style
-- direct accession keys, and composite predicates so search does not rebuild
-- that union for every query.
CREATE TABLE IF NOT EXISTS semantic_terms (
    protein_id VARCHAR NOT NULL,
    term_id VARCHAR NOT NULL,
    term_kind VARCHAR NOT NULL,
    facet VARCHAR,
    relation VARCHAR,
    source_db VARCHAR NOT NULL DEFAULT '',
    source_accession VARCHAR NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_semantic_terms_term ON semantic_terms(term_id);
CREATE INDEX IF NOT EXISTS idx_semantic_terms_protein ON semantic_terms(protein_id);

-- Constraint-free append targets used while a resumable full V2 refresh is
-- running. They are bulk-promoted into the canonical tables at completion.
CREATE TABLE IF NOT EXISTS v2_generation_atoms (
    protein_id VARCHAR NOT NULL,
    atom_id VARCHAR NOT NULL,
    facet VARCHAR NOT NULL,
    relation VARCHAR NOT NULL,
    source_accession VARCHAR,
    source_db VARCHAR,
    evidence_evalue DOUBLE,
    evidence_score DOUBLE
);
CREATE TABLE IF NOT EXISTS v2_generation_state (
    protein_id VARCHAR NOT NULL,
    activities VARCHAR[],
    roles VARCHAR[],
    architecture VARCHAR[],
    localization VARCHAR[],
    topology JSON,
    size_class VARCHAR,
    quality_flags VARCHAR[],
    composite_predicates VARCHAR[],
    unresolved_count INTEGER,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE IF NOT EXISTS v2_generation_terms (
    protein_id VARCHAR NOT NULL,
    term_id VARCHAR NOT NULL,
    term_kind VARCHAR NOT NULL,
    facet VARCHAR,
    relation VARCHAR,
    source_db VARCHAR NOT NULL DEFAULT '',
    source_accession VARCHAR NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS v2_generation_legacy (
    protein_id VARCHAR NOT NULL,
    predicates VARCHAR[] NOT NULL,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Predicate V2 full-refresh checkpoint. Chunk data and this boundary commit
-- in one transaction so interrupted generations can resume safely.
CREATE TABLE IF NOT EXISTS v2_generation_checkpoint (
    generation_key VARCHAR PRIMARY KEY,
    semantic_fingerprint VARCHAR NOT NULL,
    input_signature VARCHAR NOT NULL,
    status VARCHAR NOT NULL,
    last_protein_id VARCHAR,
    processed_count BIGINT NOT NULL DEFAULT 0,
    total_count BIGINT NOT NULL,
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP
);

-- Defense system validation (MacSyFinder co-localization)
CREATE TABLE IF NOT EXISTS defense_systems (
    system_id VARCHAR PRIMARY KEY,
    genome_id VARCHAR,
    contig_id VARCHAR,
    system_type VARCHAR,
    system_subtype VARCHAR,
    activity VARCHAR,
    genes_count INTEGER,
    protein_ids VARCHAR,
    profile_names VARCHAR,
    sys_beg VARCHAR,
    sys_end VARCHAR,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (genome_id) REFERENCES bins(bin_id)
);
CREATE INDEX IF NOT EXISTS idx_defense_systems_genome ON defense_systems(genome_id);
CREATE INDEX IF NOT EXISTS idx_defense_systems_replicon
    ON defense_systems(genome_id, contig_id);
CREATE INDEX IF NOT EXISTS idx_defense_systems_type ON defense_systems(system_type);

-- Secretion system validation (MacSyFinder/TXSScan co-localization)
CREATE TABLE IF NOT EXISTS secretion_systems (
    system_id VARCHAR PRIMARY KEY,
    genome_id VARCHAR,
    contig_id VARCHAR,
    system_type VARCHAR,
    system_subtype VARCHAR,
    genes_count INTEGER,
    protein_ids VARCHAR,
    profile_names VARCHAR,
    sys_beg VARCHAR,
    sys_end VARCHAR,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (genome_id) REFERENCES bins(bin_id)
);
CREATE INDEX IF NOT EXISTS idx_secretion_systems_genome ON secretion_systems(genome_id);
CREATE INDEX IF NOT EXISTS idx_secretion_systems_replicon
    ON secretion_systems(genome_id, contig_id);
CREATE INDEX IF NOT EXISTS idx_secretion_systems_type ON secretion_systems(system_type);

-- Normalized validated-system membership. The summary tables above retain
-- their string protein_ids fields for compatibility; this table is the query
-- substrate used by V2 and downstream operators.
CREATE TABLE IF NOT EXISTS system_proteins (
    system_id VARCHAR NOT NULL,
    protein_id VARCHAR NOT NULL,
    system_source VARCHAR NOT NULL,
    position INTEGER,
    profile_name VARCHAR,
    score DOUBLE,

    PRIMARY KEY (system_id, protein_id, system_source),
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);
CREATE INDEX IF NOT EXISTS idx_system_proteins_protein ON system_proteins(protein_id);
CREATE INDEX IF NOT EXISTS idx_system_proteins_system ON system_proteins(system_id);
CREATE INDEX IF NOT EXISTS idx_system_proteins_source ON system_proteins(system_source);

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
