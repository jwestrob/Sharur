-- @sharur-sql-version: 1
-- @requires: proteins, annotations, bins
--
-- Sharur SQL surface: per-connection macros over a Sharur dataset database.
--
-- Every object here is a TEMP macro, so the prelude needs no write access to the
-- database file, adds nothing to information_schema, and leaves dataset seals unchanged.
-- Sharur's own connection helper (DuckDBStore, so `Sharur(path, read_only=True)` too)
-- installs it on every connection. In a generic DuckDB client, attach the database,
-- make it the default catalog, then run this file:
--
--     ATTACH 'data/DATASET/sharur.duckdb' AS ds (READ_ONLY);
--     USE ds;
--     .read sharur/sql/prelude.sql          -- or paste the file, or `sharur sql prelude`
--
-- Table macros are called in FROM: `FROM genes()`, `FROM neighbors('PF00069', k := 3)`.
-- Feature specs (a comma-separated string or a list):
--     PF00069, K00001, Pkinase  accession or name in any annotation source
--     pfam:Pkinase, kofam:K00001  that source only (accession or name); `pfam:*` = any Pfam hit
--     ko:K00001                 KofamScan or KEGG
--     predicate:hydrogenase     protein_predicates
--     defense:CBASS             a curated defense_systems call (type or subtype)
--     secretion:TYPE            a curated secretion_systems call (type or subtype)
--     crispr:I-B                a crispr_cas_systems call with status crispr_cas or cas
--     a protein ID              (anchor arguments only) that protein itself
-- Results carry an `evidence` column: observed (per-domain hits), caller (curated system
-- calls and their annotation projections), predicate, or resembles (remote-homology rows:
-- CDD per protein from protein_resemblance, family-tier HHsearch hits; an observed
-- similarity, never a name). Name a system only from `caller` rows. Resembles-tier sources
-- (`cdd`) match only when named (`cdd:cd00009`, `cdd:*`); an unqualified spec matches the
-- annotations table. CDD rows that a release before protein_resemblance wrote to annotations
-- are left out everywhere (`sharur remote-homology migrate` moves them).
--
-- Sequences stay out of every result here: `genes()` is the protein table without its
-- sequence column, and FASTA leaves the database only as a file (`write_fasta`, installed by
-- the Python helper, or the COPY idiom in docs/guides/sql-surface.md), which reports a path
-- and a count.
--
-- Each `-- @object NAME KIND: DOC` line starts one statement. The Python installer runs the
-- statements one by one and records what installed (`FROM sharur_sql()` lists them). A
-- statement marked `-- @if-absent: TABLE` is the stand-in for the regular statement of the
-- same name that follows it; the installer and `sharur sql prelude --db` pick the one that
-- matches the database.


-- @object strand_sign scalar_macro: +1 or -1 for any strand encoding ('+', '-', '1', '-1', '+1'); NULL otherwise
CREATE OR REPLACE TEMP MACRO strand_sign(s) AS
    CASE WHEN CAST(s AS VARCHAR) IN ('+', '1', '+1') THEN 1
         WHEN CAST(s AS VARCHAR) IN ('-', '-1') THEN -1 END;


-- @object id_list scalar_macro: a comma-separated string or a list of strings as VARCHAR[] (items trimmed, empty items dropped)
CREATE OR REPLACE TEMP MACRO id_list(x) AS
    list_filter(
        list_transform(string_split(regexp_replace(CAST(x AS VARCHAR), '^\[|\]$', '', 'g'), ','), s -> trim(s)),
        s -> s <> '');


-- @object rank_of scalar_macro: one GTDB rank (domain, phylum, class, order, family, genus, species) from a taxonomy string
CREATE OR REPLACE TEMP MACRO rank_of(taxonomy, rank) AS
    nullif(regexp_extract(taxonomy,
        CASE lower(rank) WHEN 'domain' THEN 'd' WHEN 'phylum' THEN 'p' WHEN 'class' THEN 'c'
                         WHEN 'order' THEN 'o' WHEN 'family' THEN 'f' WHEN 'genus' THEN 'g'
                         WHEN 'species' THEN 's' END || '__([^;]*)', 1), '');


-- @object genes table_macro: one row per called gene; the proteins table without its sequence column, plus strand_sign
CREATE OR REPLACE TEMP MACRO genes() AS TABLE
    SELECT * EXCLUDE (sequence), strand_sign(strand) AS strand_sign
    FROM proteins;


-- @object genomes table_macro: bins with GTDB ranks split out (domain, phylum, class, "order", family, genus, species)
CREATE OR REPLACE TEMP MACRO genomes() AS TABLE
    SELECT b.*,
           rank_of(b.taxonomy, 'domain') AS domain,
           rank_of(b.taxonomy, 'phylum') AS phylum,
           rank_of(b.taxonomy, 'class') AS class,
           rank_of(b.taxonomy, 'order') AS "order",
           rank_of(b.taxonomy, 'family') AS family,
           rank_of(b.taxonomy, 'genus') AS genus,
           rank_of(b.taxonomy, 'species') AS species
    FROM bins b;


-- @object protein_resemblance_rows table_macro: per-protein resemblance rows (remote homology: CDD per protein): protein_id, source, search_run, db_name, db_version, accession, name, description, evalue, score, start_aa, end_aa; empty on a database without protein_resemblance
-- @if-absent: protein_resemblance
CREATE OR REPLACE TEMP MACRO protein_resemblance_rows() AS TABLE
    SELECT NULL::VARCHAR AS protein_id, NULL::VARCHAR AS source, NULL::VARCHAR AS search_run,
           NULL::VARCHAR AS db_name, NULL::VARCHAR AS db_version, NULL::VARCHAR AS accession,
           NULL::VARCHAR AS name, NULL::VARCHAR AS description, NULL::FLOAT AS evalue, NULL::FLOAT AS score,
           NULL::INTEGER AS start_aa, NULL::INTEGER AS end_aa
    WHERE FALSE;


-- @object protein_resemblance_rows table_macro: per-protein resemblance rows (remote homology: CDD per protein): protein_id, source, search_run, db_name, db_version, accession, name, description, evalue, score, start_aa, end_aa; empty on a database without protein_resemblance
CREATE OR REPLACE TEMP MACRO protein_resemblance_rows() AS TABLE
    SELECT protein_id, source, search_run, db_name, db_version, accession, name, description, evalue, score,
           start_aa, end_aa
    FROM protein_resemblance;


-- @object feature_proteins table_macro: proteins carrying any of the feature specs: feature, protein_id, bin_id, contig_id, gene_index, strand, source, accession, name, evidence
CREATE OR REPLACE TEMP MACRO feature_proteins(features) AS TABLE
    WITH spec AS MATERIALIZED (
        SELECT DISTINCT f AS feature,
               CASE WHEN strpos(f, ':') > 0 THEN lower(split_part(f, ':', 1)) END AS ns,
               CASE WHEN strpos(f, ':') > 0 THEN substr(f, strpos(f, ':') + 1) ELSE f END AS key
        FROM (SELECT unnest(id_list(features)) AS f)
    ),
    spec_src AS MATERIALIZED (
        SELECT feature, key, unnest(CASE WHEN ns = 'ko' THEN ['kofam', 'kegg'] ELSE [ns] END) AS spec_source
        FROM spec
        WHERE ns IS NOT NULL AND ns NOT IN ('predicate', 'defense', 'secretion', 'crispr')
    ),
    projection AS MATERIALIZED (
        SELECT DISTINCT system_source FROM system_proteins
    ),
    ann AS MATERIALIZED (
        SELECT a.protein_id, a.source, a.accession, a.name, FALSE AS resembles
        FROM annotations a
        WHERE (a.accession IN (SELECT key FROM spec)
               OR a.name IN (SELECT key FROM spec)
               OR a.source IN (SELECT spec_source FROM spec_src WHERE key = '*'))
          AND a.source NOT IN ('cdd')  -- resembles-tier rows left in annotations by an earlier release
        UNION ALL
        -- resembles-tier rows match only when their source is named
        SELECT r.protein_id, r.source, r.accession, r.name, TRUE AS resembles
        FROM protein_resemblance_rows() r
        WHERE r.source IN (SELECT spec_source FROM spec_src)
          AND (r.accession IN (SELECT key FROM spec)
               OR r.name IN (SELECT key FROM spec)
               OR r.source IN (SELECT spec_source FROM spec_src WHERE key = '*'))
    ),
    ann_hits AS (
        SELECT s.feature, a.protein_id, a.source, a.accession, a.name, a.resembles
        FROM ann a JOIN spec s ON s.ns IS NULL AND (a.accession = s.key OR a.name = s.key)
        WHERE NOT a.resembles
        UNION ALL
        SELECT s.feature, a.protein_id, a.source, a.accession, a.name, a.resembles
        FROM ann a JOIN spec_src s
          ON a.source = s.spec_source AND (a.accession = s.key OR a.name = s.key OR s.key = '*')
    ),
    pred_hits AS (
        SELECT s.feature, pp.protein_id, 'predicate' AS source, s.key AS accession, s.key AS name
        FROM protein_predicates pp
        JOIN (SELECT * FROM spec WHERE ns = 'predicate') s ON list_contains(pp.predicates, s.key)
    ),
    systems AS (
        SELECT system_id, 'defense' AS ns, system_type, system_subtype FROM defense_systems
        UNION ALL
        SELECT system_id, 'secretion', system_type, system_subtype FROM secretion_systems
        UNION ALL
        SELECT system_id, 'crispr', best_type, best_type FROM crispr_cas_systems
        WHERE status IN ('crispr_cas', 'cas')
    ),
    sys_hits AS (
        SELECT s.feature, sp.protein_id, sp.system_source AS source, y.system_type AS accession,
               coalesce(y.system_subtype, y.system_type) AS name
        FROM spec s
        JOIN systems y ON y.ns = s.ns AND (s.key = y.system_type OR s.key = y.system_subtype OR s.key = '*')
        JOIN system_proteins sp ON sp.system_id = y.system_id
    ),
    hits AS (
        SELECT feature, protein_id, source, accession, name,
               CASE WHEN resembles THEN 'resembles'
                    WHEN source IN (SELECT system_source FROM projection) THEN 'caller' ELSE 'observed' END AS evidence
        FROM ann_hits
        UNION ALL
        SELECT *, 'predicate' FROM pred_hits
        UNION ALL
        SELECT *, 'caller' FROM sys_hits
    )
    SELECT DISTINCT h.feature, h.protein_id, p.bin_id, p.contig_id, p.gene_index, p.strand,
           h.source, h.accession, h.name, h.evidence
    FROM hits h JOIN proteins p ON p.protein_id = h.protein_id;


-- @object feature_presence table_macro: every genome x feature, with present, n_proteins, sources, evidence, completeness, contamination and GTDB ranks (one uniform shape across annotation sources and callers)
CREATE OR REPLACE TEMP MACRO feature_presence(features, min_completeness := 0, max_contamination := 100) AS TABLE
    WITH hit AS MATERIALIZED (
        SELECT feature, bin_id, count(DISTINCT protein_id) AS n_proteins,
               list(DISTINCT source ORDER BY source) AS sources,
               list(DISTINCT evidence ORDER BY evidence) AS evidence
        FROM feature_proteins(features)
        GROUP BY feature, bin_id
    ),
    g AS MATERIALIZED (
        SELECT * FROM genomes()
        WHERE coalesce(completeness, 0) >= min_completeness
          AND coalesce(contamination, 0) <= max_contamination
    ),
    f AS (SELECT DISTINCT unnest(id_list(features)) AS feature)
    SELECT g.bin_id, f.feature, h.bin_id IS NOT NULL AS present, coalesce(h.n_proteins, 0) AS n_proteins,
           h.sources, h.evidence, g.completeness, g.contamination, g.n_contigs,
           g.domain, g.phylum, g.class, g."order", g.family, g.genus, g.species
    FROM g CROSS JOIN f
    LEFT JOIN hit h ON h.bin_id = g.bin_id AND h.feature = f.feature;


-- @object anchor_proteins table_macro: proteins named directly by ID or matched by a feature spec, with the features they matched
CREATE OR REPLACE TEMP MACRO anchor_proteins(anchor) AS TABLE
    WITH hits AS (
        SELECT protein_id, NULL::VARCHAR AS feature
        FROM proteins WHERE protein_id IN (SELECT unnest(id_list(anchor)))
        UNION ALL
        SELECT protein_id, feature FROM feature_proteins(anchor)
    )
    SELECT h.protein_id, list(DISTINCT h.feature) FILTER (WHERE h.feature IS NOT NULL) AS anchor_features
    FROM hits h
    GROUP BY h.protein_id;


-- @object protein_summary table_macro: per protein (a list or comma-separated IDs): pfam names in domain order, KO accessions, KO names, curated system calls, other observed hits, n_annotations (rows of annotations), resembles (protein_resemblance rows, source:accession)
CREATE OR REPLACE TEMP MACRO protein_summary(ids) AS TABLE
    WITH want AS MATERIALIZED (SELECT DISTINCT unnest(id_list(ids)) AS protein_id),
    projection AS MATERIALIZED (SELECT DISTINCT system_source FROM system_proteins),
    ann AS MATERIALIZED (
        SELECT a.* FROM annotations a WHERE a.protein_id IN (SELECT protein_id FROM want)
          AND a.source NOT IN ('cdd')  -- resembles-tier rows left in annotations by an earlier release
    ),
    a_sum AS (
        SELECT protein_id,
               list(name ORDER BY start_aa NULLS LAST, name) FILTER (WHERE source = 'pfam') AS pfam,
               list(DISTINCT accession ORDER BY accession) FILTER (WHERE source IN ('kofam', 'kegg')) AS ko,
               min(CASE WHEN source = 'kegg' THEN name END) AS ko_name,
               list(DISTINCT source || ':' || accession ORDER BY source || ':' || accession)
                   FILTER (WHERE source NOT IN ('pfam', 'kofam', 'kegg')
                             AND source NOT IN (SELECT system_source FROM projection)) AS other,
               count(*) AS n_annotations
        FROM ann GROUP BY protein_id
    ),
    r_sum AS (
        SELECT protein_id,
               list(DISTINCT source || ':' || accession ORDER BY source || ':' || accession) AS resembles
        FROM protein_resemblance_rows()
        WHERE protein_id IN (SELECT protein_id FROM want)
        GROUP BY protein_id
    ),
    systems AS (
        SELECT system_id, 'defense' AS ns, coalesce(system_subtype, system_type) AS label FROM defense_systems
        UNION ALL
        SELECT system_id, 'secretion', coalesce(system_subtype, system_type) FROM secretion_systems
        UNION ALL
        SELECT system_id, 'crispr', best_type FROM crispr_cas_systems WHERE status IN ('crispr_cas', 'cas')
    ),
    s_sum AS (
        SELECT sp.protein_id, list(DISTINCT y.ns || ':' || y.label ORDER BY y.ns || ':' || y.label) AS systems
        FROM system_proteins sp JOIN systems y USING (system_id)
        WHERE sp.protein_id IN (SELECT protein_id FROM want)
        GROUP BY sp.protein_id
    )
    SELECT w.protein_id,
           coalesce(a.pfam, []::VARCHAR[]) AS pfam,
           coalesce(a.ko, []::VARCHAR[]) AS ko,
           a.ko_name,
           coalesce(s.systems, []::VARCHAR[]) AS systems,
           coalesce(a.other, []::VARCHAR[]) AS other,
           coalesce(a.n_annotations, 0) AS n_annotations,
           coalesce(r.resembles, []::VARCHAR[]) AS resembles
    FROM want w
    LEFT JOIN a_sum a USING (protein_id)
    LEFT JOIN r_sum r USING (protein_id)
    LEFT JOIN s_sum s USING (protein_id);


-- @object neighbors table_macro: genes within k positions of each anchor (protein IDs or feature specs), offsets in the anchor's orientation (negative = upstream), each neighbor's annotations summarized
CREATE OR REPLACE TEMP MACRO neighbors(anchor, k := 5) AS TABLE
    WITH a AS MATERIALIZED (
        SELECT ap.protein_id, ap.anchor_features, p.contig_id
        FROM anchor_proteins(anchor) ap JOIN proteins p USING (protein_id)
    ),
    cg AS MATERIALIZED (
        SELECT p.protein_id, p.bin_id, p.contig_id, p.gene_index, p.start, p.end_coord, p.strand,
               p.sequence_length,
               row_number() OVER (PARTITION BY p.contig_id
                                  ORDER BY p.gene_index NULLS LAST, p.start, p.protein_id) AS pos,
               count(*) OVER (PARTITION BY p.contig_id) AS n_genes
        FROM proteins p
        WHERE p.contig_id IN (SELECT contig_id FROM a)
    ),
    ap AS MATERIALIZED (
        SELECT a.protein_id AS anchor_id, a.anchor_features, c.bin_id, c.contig_id, c.pos AS a_pos,
               c.strand AS anchor_strand, coalesce(strand_sign(c.strand), 1) AS a_sign, c.n_genes
        FROM a JOIN cg c USING (protein_id)
    ),
    nb AS MATERIALIZED (
        SELECT ap.*, n.protein_id, n.gene_index, n.start, n.end_coord, n.strand, n.sequence_length,
               (n.pos - ap.a_pos) * ap.a_sign AS "offset"
        FROM ap JOIN cg n ON n.contig_id = ap.contig_id AND abs(n.pos - ap.a_pos) <= k
    ),
    summ AS (
        FROM protein_summary((SELECT list(DISTINCT protein_id) FROM nb))
    )
    SELECT nb.anchor_id, nb.anchor_features, nb.bin_id, nb.contig_id, nb.anchor_strand,
           CASE WHEN nb.a_sign = 1 THEN nb.a_pos - 1 ELSE nb.n_genes - nb.a_pos END AS genes_upstream_to_edge,
           CASE WHEN nb.a_sign = 1 THEN nb.n_genes - nb.a_pos ELSE nb.a_pos - 1 END AS genes_downstream_to_edge,
           nb."offset", nb.protein_id, nb.gene_index, nb.start, nb.end_coord, nb.strand,
           strand_sign(nb.strand) = nb.a_sign AS same_strand,
           nb.sequence_length AS length_aa,
           s.pfam, s.ko, s.ko_name, s.systems, s.other, s.n_annotations, s.resembles
    FROM nb LEFT JOIN summ s USING (protein_id)
    ORDER BY nb.bin_id, nb.contig_id, nb.anchor_id, nb."offset";


-- @object neighbor_features table_macro: how often each feature of one annotation source (src := 'pfam'; '*' for all) sits at each oriented offset around the anchors: offset, same_strand, feature, n_genomes, n_anchors
CREATE OR REPLACE TEMP MACRO neighbor_features(anchor, k := 5, src := 'pfam') AS TABLE
    WITH n AS MATERIALIZED (
        SELECT * FROM neighbors(anchor, k := k) WHERE "offset" <> 0
    ),
    rows AS (
        SELECT protein_id, source, accession, name FROM annotations
        WHERE source NOT IN ('cdd')  -- resembles-tier rows left in annotations by an earlier release
        UNION ALL
        SELECT protein_id, source, accession, name FROM protein_resemblance_rows()
    ),
    f AS (
        SELECT n.bin_id, n.anchor_id, n."offset", n.same_strand, a.source || ':' || a.accession AS feature,
               any_value(a.name) AS name
        FROM n JOIN rows a ON a.protein_id = n.protein_id AND (a.source = src OR src = '*')
        GROUP BY n.bin_id, n.anchor_id, n."offset", n.same_strand, a.source, a.accession
    )
    SELECT "offset", same_strand, feature, any_value(name) AS name,
           count(DISTINCT bin_id) AS n_genomes, count(DISTINCT anchor_id) AS n_anchors
    FROM f GROUP BY "offset", same_strand, feature
    ORDER BY n_genomes DESC, feature, "offset";


-- @object slot table_macro: what sits beside one anchor (side := 'upstream' or 'downstream', in the anchor's orientation) or between an anchor and its nearest partner (partner := spec, side := 'upstream', 'downstream' or 'either'), per anchor protein, among called genes, with genome completeness and GTDB ranks
CREATE OR REPLACE TEMP MACRO slot(anchor, side := 'upstream', partner := NULL, max_span := 10,
                                  min_completeness := 0, max_contamination := 100) AS TABLE
    WITH g AS MATERIALIZED (
        SELECT * FROM genomes()
        WHERE coalesce(completeness, 0) >= min_completeness
          AND coalesce(contamination, 0) <= max_contamination
    ),
    a AS MATERIALIZED (
        SELECT ap.protein_id, ap.anchor_features, p.contig_id
        FROM anchor_proteins(anchor) ap
        JOIN proteins p USING (protein_id)
        WHERE p.bin_id IN (SELECT bin_id FROM g)
    ),
    cg AS MATERIALIZED (
        SELECT p.protein_id, p.bin_id, p.contig_id, p.strand, p.sequence_length,
               row_number() OVER (PARTITION BY p.contig_id
                                  ORDER BY p.gene_index NULLS LAST, p.start, p.protein_id) AS pos
        FROM proteins p
        WHERE p.contig_id IN (SELECT contig_id FROM a)
    ),
    ap AS MATERIALIZED (
        SELECT a.protein_id AS anchor_id, a.anchor_features, c.bin_id, c.contig_id, c.pos AS a_pos,
               c.strand AS anchor_strand, coalesce(strand_sign(c.strand), 1) AS a_sign
        FROM a JOIN cg c USING (protein_id)
    ),
    -- single-anchor mode: the gene at offset -1 (upstream) or +1 (downstream)
    single AS (
        SELECT ap.*, NULL::VARCHAR AS partner_id, NULL::BIGINT AS partner_offset,
               CASE WHEN o.protein_id IS NULL THEN 0 ELSE 1 END AS n_between,
               o.protein_id AS occupant_id,
               CASE WHEN o.protein_id IS NULL THEN NULL
                    ELSE (o.pos - ap.a_pos) * ap.a_sign END AS occupant_offset,
               o.protein_id IS NULL AS at_contig_edge
        FROM ap
        LEFT JOIN cg o ON o.contig_id = ap.contig_id
                      AND o.pos = ap.a_pos + (CASE WHEN lower(side) = 'downstream' THEN 1 ELSE -1 END) * ap.a_sign
        WHERE partner IS NULL
    ),
    -- two-anchor mode: the nearest partner within max_span genes on the requested side
    pp AS (
        SELECT protein_id FROM anchor_proteins(partner)
    ),
    pair_all AS (
        SELECT ap.anchor_id, b.protein_id AS partner_id, (b.pos - ap.a_pos) * ap.a_sign AS partner_offset,
               b.pos AS b_pos
        FROM ap JOIN cg b ON b.contig_id = ap.contig_id AND b.protein_id <> ap.anchor_id
                         AND abs(b.pos - ap.a_pos) BETWEEN 1 AND max_span + 1
        WHERE b.protein_id IN (SELECT protein_id FROM pp)
    ),
    pair_side AS (
        SELECT * FROM pair_all
        WHERE (lower(side) = 'upstream' AND partner_offset < 0)
           OR (lower(side) = 'downstream' AND partner_offset > 0)
           OR lower(side) = 'either'
        QUALIFY row_number() OVER (PARTITION BY anchor_id ORDER BY abs(partner_offset), partner_offset) = 1
    ),
    -- Keep `pairs` and `between_genes` inlined: in DuckDB 1.2.2, materializing a CTE whose
    -- filter is the constant `partner IS NOT NULL` empties the single-anchor branch as well
    -- (tests/test_sql_surface.py covers both modes).
    pairs AS (
        SELECT ap.*, ps.partner_id, ps.partner_offset, ps.b_pos
        FROM ap LEFT JOIN pair_side ps USING (anchor_id)
        WHERE partner IS NOT NULL
    ),
    between_genes AS (
        SELECT pr.anchor_id, o.protein_id AS occupant_id, (o.pos - pr.a_pos) * pr.a_sign AS occupant_offset
        FROM pairs pr JOIN cg o ON o.contig_id = pr.contig_id
                               AND o.pos > least(pr.a_pos, pr.b_pos) AND o.pos < greatest(pr.a_pos, pr.b_pos)
        WHERE pr.partner_id IS NOT NULL
    ),
    double AS (
        SELECT pr.* EXCLUDE (b_pos),
               coalesce(abs(pr.partner_offset) - 1, 0) AS n_between,
               bg.occupant_id, bg.occupant_offset, FALSE AS at_contig_edge
        FROM pairs pr LEFT JOIN between_genes bg USING (anchor_id)
    ),
    s AS MATERIALIZED (
        SELECT * FROM single
        UNION ALL BY NAME
        SELECT * FROM double
    ),
    summ AS (
        FROM protein_summary((SELECT list(DISTINCT occupant_id) FROM s WHERE occupant_id IS NOT NULL))
    )
    SELECT s.anchor_id, s.anchor_features, s.bin_id, s.contig_id, s.anchor_strand, lower(side) AS side,
           s.partner_id, s.partner_offset, s.n_between,
           CASE WHEN s.occupant_offset IS NULL THEN NULL ELSE abs(s.occupant_offset) END AS slot_position,
           s.occupant_id, s.occupant_offset, o.strand AS occupant_strand,
           strand_sign(o.strand) = s.a_sign AS occupant_same_strand,
           o.sequence_length AS occupant_length_aa,
           s.at_contig_edge,
           sm.pfam, sm.ko, sm.ko_name, sm.systems, sm.other, sm.n_annotations, sm.resembles,
           g.completeness, g.contamination, g.taxonomy,
           g.domain, g.phylum, g.class, g."order", g.family, g.genus, g.species
    FROM s
    JOIN g ON g.bin_id = s.bin_id
    LEFT JOIN cg o ON o.protein_id = s.occupant_id
    LEFT JOIN summ sm ON sm.protein_id = s.occupant_id
    ORDER BY s.bin_id, s.contig_id, s.anchor_id, slot_position NULLS FIRST;


-- @object slot_tally table_macro: per taxon (rank := 'order' by default) and, with split := spec, by whether the genome carries that feature: genomes with an anchor, with an occupant matching `occupant` (same strand / any strand), with an occupant lacking any Pfam hit, at a contig edge, and without a partner in range
CREATE OR REPLACE TEMP MACRO slot_tally(anchor, occupant, side := 'upstream', partner := NULL, rank := 'order',
                                        split := NULL, max_span := 10,
                                        min_completeness := 0, max_contamination := 100) AS TABLE
    WITH s AS MATERIALIZED (
        FROM slot(anchor, side := side, partner := partner, max_span := max_span,
                  min_completeness := min_completeness, max_contamination := max_contamination)
    ),
    occ AS MATERIALIZED (SELECT DISTINCT protein_id FROM feature_proteins(occupant)),
    spl AS MATERIALIZED (SELECT DISTINCT bin_id FROM feature_proteins(split)),
    x AS (
        SELECT s.*,
               coalesce(s.occupant_id IN (SELECT protein_id FROM occ), FALSE) AS occ_match,
               CASE WHEN split IS NULL THEN NULL ELSE s.bin_id IN (SELECT bin_id FROM spl) END AS genome_has_split,
               rank_of(s.taxonomy, rank) AS taxon
        FROM s
    )
    SELECT lower(rank) AS rank,
           CASE WHEN lower(rank) IN ('domain', 'phylum') THEN NULL ELSE phylum END AS phylum,
           taxon,
           CAST(split AS VARCHAR) AS split_feature,
           genome_has_split,
           count(DISTINCT bin_id) AS genomes_with_anchor,
           count(DISTINCT bin_id) FILTER (WHERE occ_match AND occupant_same_strand) AS genomes_occupant_same_strand,
           count(DISTINCT bin_id) FILTER (WHERE occ_match) AS genomes_occupant_any_strand,
           count(DISTINCT bin_id) FILTER (WHERE occupant_id IS NOT NULL AND len(pfam) = 0) AS genomes_occupant_no_pfam,
           count(DISTINCT bin_id) FILTER (WHERE at_contig_edge) AS genomes_at_contig_edge,
           count(DISTINCT bin_id) FILTER (WHERE partner IS NOT NULL AND partner_id IS NULL) AS genomes_no_partner_in_span,
           count(DISTINCT anchor_id) AS anchor_proteins,
           count(DISTINCT occupant_id) FILTER (WHERE occ_match AND occupant_same_strand) AS occupant_proteins_same_strand
    FROM x
    GROUP BY ALL
    ORDER BY phylum NULLS FIRST, genomes_with_anchor DESC, taxon, genome_has_split;


-- @object family_proteins table_macro: members of protein families or subfamilies (family_id values of any clustering run) in the feature_proteins shape: feature, protein_id, bin_id, contig_id, gene_index, strand, source (the clustering run), accession (the family id), name (the run's level), evidence; needs the protein-family layer (docs/guides/protein-families.md)
CREATE OR REPLACE TEMP MACRO family_proteins(families) AS TABLE
    SELECT f.family_id AS feature, f.protein_id, p.bin_id, p.contig_id, p.gene_index, p.strand,
           f.clustering_run AS source, f.family_id AS accession, r.level AS name, 'observed' AS evidence
    FROM protein_families f
    JOIN protein_family_runs r USING (clustering_run)
    JOIN proteins p USING (protein_id)
    WHERE f.family_id IN (SELECT unnest(id_list(families)));


-- @object family_presence table_macro: every genome x family (or subfamily) in the feature_presence shape: present, n_proteins, sources, evidence, completeness, contamination and GTDB ranks; needs the protein-family layer
CREATE OR REPLACE TEMP MACRO family_presence(families, min_completeness := 0, max_contamination := 100) AS TABLE
    WITH hit AS MATERIALIZED (
        SELECT feature, bin_id, count(DISTINCT protein_id) AS n_proteins,
               list(DISTINCT source ORDER BY source) AS sources,
               list(DISTINCT evidence ORDER BY evidence) AS evidence
        FROM family_proteins(families)
        GROUP BY feature, bin_id
    ),
    g AS MATERIALIZED (
        SELECT * FROM genomes()
        WHERE coalesce(completeness, 0) >= min_completeness
          AND coalesce(contamination, 0) <= max_contamination
    ),
    f AS (SELECT DISTINCT unnest(id_list(families)) AS feature)
    SELECT g.bin_id, f.feature, h.bin_id IS NOT NULL AS present, coalesce(h.n_proteins, 0) AS n_proteins,
           h.sources, h.evidence, g.completeness, g.contamination, g.n_contigs,
           g.domain, g.phylum, g.class, g."order", g.family, g.genus, g.species
    FROM g CROSS JOIN f
    LEFT JOIN hit h ON h.bin_id = g.bin_id AND h.feature = f.feature;


-- @object families_of table_macro: per protein (a list or comma-separated IDs), its family in every clustering run: clustering_run, level, is_default, family_id, is_representative, parent_family_id, n_members, n_genomes, n_orders, pfam_dark, dominant_pfam_architecture, dominant_ko; needs the protein-family layer
CREATE OR REPLACE TEMP MACRO families_of(ids) AS TABLE
    SELECT f.protein_id, f.clustering_run, r.level, r.is_default, f.family_id, f.is_representative,
           s.parent_family_id, s.n_members, s.n_genomes, s.n_orders, s.pfam_dark,
           s.dominant_pfam_architecture, s.dominant_ko
    FROM protein_families f
    JOIN protein_family_runs r USING (clustering_run)
    LEFT JOIN family_stats s ON s.clustering_run = f.clustering_run AND s.family_id = f.family_id
    WHERE f.protein_id IN (SELECT unnest(id_list(ids)))
    ORDER BY f.protein_id, r.level, f.clustering_run;


-- @object resembles_of table_macro: per protein (a list or comma-separated IDs), what it resembles: hits of its families in every clustering run that meet the rule (Pfam or PDB by HHsearch or HHblits, probability >= min_probability over >= min_cols aligned columns; tier 'family') and its own protein_resemblance rows (CDD; tier 'protein'), with db_name, db_version, msa_kind, probability, evalue, aligned columns, query range, unknown_function, hit_class (characterized or unknown_function) and evidence 'resembles'; needs the remote-homology layer (docs/guides/remote-homology.md)
CREATE OR REPLACE TEMP MACRO resembles_of(ids, min_probability := 90, min_cols := 30) AS TABLE
    WITH want AS MATERIALIZED (SELECT DISTINCT unnest(id_list(ids)) AS protein_id),
    fam AS (
        SELECT f.protein_id, 'family' AS tier, h.clustering_run, h.family_id, h.source, h.db_name, h.db_version,
               h.msa_kind, h.hit_accession, h.hit_name, h.hit_description, h.probability, h.evalue,
               h.aligned_cols, h.q_start, h.q_end, h.unknown_function
        FROM protein_families f
        JOIN family_annotations h ON h.family_id = f.family_id AND h.clustering_run = f.clustering_run
        WHERE f.protein_id IN (SELECT protein_id FROM want)
          AND h.source IN ('hhsearch_pfam', 'hhsearch_pdb', 'hhblits_pdb')
          AND h.probability >= min_probability AND h.aligned_cols >= min_cols
    ),
    prot AS (
        SELECT r.protein_id, 'protein' AS tier, NULL::VARCHAR AS clustering_run, NULL::VARCHAR AS family_id,
               r.source, r.db_name, r.db_version, NULL::VARCHAR AS msa_kind,
               r.accession AS hit_accession, r.name AS hit_name, CAST(r.description AS VARCHAR) AS hit_description,
               NULL::DOUBLE AS probability, CAST(r.evalue AS DOUBLE) AS evalue, NULL::INTEGER AS aligned_cols,
               r.start_aa AS q_start, r.end_aa AS q_end,
               regexp_matches(coalesce(r.name, '') || ' ' || coalesce(CAST(r.description AS VARCHAR), ''),
                              '(?i)DUF\d|UPF\d|unknown function|uncharacteri|hypothetical|function unknown|unknown protein')
                   AS unknown_function
        FROM protein_resemblance_rows() r
        WHERE r.protein_id IN (SELECT protein_id FROM want)
    )
    SELECT x.*, CASE WHEN x.unknown_function THEN 'unknown_function' ELSE 'characterized' END AS hit_class,
           'resembles' AS evidence
    FROM (SELECT * FROM fam UNION ALL SELECT * FROM prot) x
    ORDER BY x.protein_id, x.tier, x.probability DESC NULLS LAST, x.evalue;


-- @object resembles_proteins table_macro: proteins that resemble any of the targets (Pfam or PDB accessions or names reached through their family's hits at the rule, or accessions or names of their own protein_resemblance rows, e.g. CDD) in the feature_proteins shape, evidence 'resembles'; needs the remote-homology layer
CREATE OR REPLACE TEMP MACRO resembles_proteins(targets, min_probability := 90, min_cols := 30) AS TABLE
    WITH spec AS MATERIALIZED (SELECT DISTINCT unnest(id_list(targets)) AS key),
    fam AS (
        SELECT s.key AS feature, f.protein_id, h.source, h.hit_accession AS accession, h.hit_name AS name
        FROM spec s
        JOIN family_annotations h ON h.hit_accession = s.key OR h.hit_name = s.key
        JOIN protein_families f ON f.family_id = h.family_id AND f.clustering_run = h.clustering_run
        WHERE h.source IN ('hhsearch_pfam', 'hhsearch_pdb', 'hhblits_pdb')
          AND h.probability >= min_probability AND h.aligned_cols >= min_cols
    ),
    prot AS (
        SELECT s.key AS feature, r.protein_id, r.source, r.accession, r.name
        FROM spec s JOIN protein_resemblance_rows() r ON r.accession = s.key OR r.name = s.key
    )
    SELECT DISTINCT h.feature, h.protein_id, p.bin_id, p.contig_id, p.gene_index, p.strand,
           h.source, h.accession, h.name, 'resembles' AS evidence
    FROM (SELECT * FROM fam UNION ALL SELECT * FROM prot) h JOIN proteins p USING (protein_id);


-- @object resembles_presence table_macro: every genome x target in the feature_presence shape (present, n_proteins, sources, evidence, completeness, contamination, GTDB ranks) from resembles_proteins; needs the remote-homology layer
CREATE OR REPLACE TEMP MACRO resembles_presence(targets, min_probability := 90, min_cols := 30,
                                                min_completeness := 0, max_contamination := 100) AS TABLE
    WITH hit AS MATERIALIZED (
        SELECT feature, bin_id, count(DISTINCT protein_id) AS n_proteins,
               list(DISTINCT source ORDER BY source) AS sources,
               list(DISTINCT evidence ORDER BY evidence) AS evidence
        FROM resembles_proteins(targets, min_probability := min_probability, min_cols := min_cols)
        GROUP BY feature, bin_id
    ),
    g AS MATERIALIZED (
        SELECT * FROM genomes()
        WHERE coalesce(completeness, 0) >= min_completeness
          AND coalesce(contamination, 0) <= max_contamination
    ),
    f AS (SELECT DISTINCT unnest(id_list(targets)) AS feature)
    SELECT g.bin_id, f.feature, h.bin_id IS NOT NULL AS present, coalesce(h.n_proteins, 0) AS n_proteins,
           h.sources, h.evidence, g.completeness, g.contamination, g.n_contigs,
           g.domain, g.phylum, g.class, g."order", g.family, g.genus, g.species
    FROM g CROSS JOIN f
    LEFT JOIN hit h ON h.bin_id = g.bin_id AND h.feature = f.feature;


-- @object family_annotations_via table_macro: family-tier hits of other clustering runs carried to the families of target_run that share the searched family's representative (family_stats.rep_protein_id; the stored view family_annotations_carried holds every target run): family_id and clustering_run of the target family, source_run, source_family_id, source_members, target_members, then source, db_name, db_version, msa_kind, hit columns, unknown_function and evidence 'resembles'; meaningful when target_run only adds members to the searched families (verify with docs/guides/remote-homology.md); needs the remote-homology layer
CREATE OR REPLACE TEMP MACRO family_annotations_via(target_run) AS TABLE
    SELECT t.family_id, t.clustering_run, h.clustering_run AS source_run, h.family_id AS source_family_id,
           s.n_members AS source_members, t.n_members AS target_members,
           h.source, h.db_name, h.db_version, h.msa_kind, h.hit_rank, h.hit_accession, h.hit_name,
           h.hit_description, h.probability, h.evalue, h.aligned_cols, h.q_start, h.q_end, h.t_start, h.t_end,
           h.unknown_function, 'resembles' AS evidence
    FROM family_annotations h
    JOIN family_stats s ON s.family_id = h.family_id AND s.clustering_run = h.clustering_run
    JOIN family_stats t ON t.rep_protein_id = s.rep_protein_id AND t.clustering_run = target_run
    WHERE h.clustering_run <> target_run;
