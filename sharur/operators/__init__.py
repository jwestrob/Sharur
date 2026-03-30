"""
Sharur Operator System.

The Sharur class provides a high-level API for exploring metagenomic datasets.
All methods return SharurResult objects with formatted output and metadata.

Example usage:
    from sharur.operators import Sharur

    b = Sharur("data/sharur.duckdb")
    print(b.overview())
    print(b.search_by_predicates(has=["giant", "unannotated"]))
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from sharur.core.session import ExplorationSession
from sharur.operators.base import SharurResult, OperatorContext, OperatorTrace, ResultMeta
from sharur.operators.introspection import describe_schema, overview
from sharur.operators.navigation import (
    get_genome,
    get_neighborhood,
    get_protein,
    list_genomes,
    list_proteins,
)
from sharur.operators.search import search_by_predicates, search_proteins
from sharur.operators.export import export_fasta, export_neighborhood_fasta, get_sequence
from sharur.operators.similarity import find_similar, find_similar_to_set
from sharur.operators.visualization import (
    visualize_neighborhood,
    visualize_domain_architecture,
    get_kegg_pathway_context,
)
from sharur.operators.structure import (
    predict_structure,
    predict_structure_from_sequence,
    batch_predict_structures,
)
from sharur.operators.foldseek import (
    search_foldseek,
    search_foldseek_for_protein,
    list_foldseek_databases,
)
from sharur.operators.validation import (
    PROBLEM_DOMAINS,
    validate_annotation,
    validate_context,
    analyze_crispr_systems,
    detect_annotation_errors,
)
from sharur.operators.manifest import AnalysisManifest


class Sharur:
    """
    Main facade for Sharur operator system.

    Provides a clean API for dataset exploration with automatic
    session and database management.
    """

    def __init__(self, db_path: Optional[Path | str] = None):
        """
        Initialize Sharur instance.

        Args:
            db_path: Path to DuckDB database file. If None, uses in-memory DB.
        """
        self._db_path = Path(db_path) if db_path else None
        self._session: Optional[ExplorationSession] = None
        self._manifest: Optional[AnalysisManifest] = None
        self._hypothesis_registry = None

    @property
    def session(self) -> ExplorationSession:
        """Lazy-load exploration session."""
        if self._session is None:
            self._session = ExplorationSession(db_path=self._db_path)
        return self._session

    @property
    def store(self):
        """Direct access to DuckDB store."""
        return self.session.db

    @property
    def manifest(self) -> AnalysisManifest:
        """Lazy-load analysis manifest."""
        if self._manifest is None:
            self._manifest = AnalysisManifest(self._db_path, store=self.store)
        return self._manifest

    def resume(self) -> str:
        """
        Generate a summary of current analysis state for resuming a session.

        Call this at the start of a new session to understand what has been
        done and what remains to explore.

        Returns:
            Formatted markdown summary of analysis state
        """
        return self.manifest.get_status_summary()

    def save_manifest(self) -> None:
        """Save the manifest to disk."""
        self.manifest.save()

    # ------------------------------------------------------------------ #
    # Introspection operators
    # ------------------------------------------------------------------ #

    def overview(self) -> SharurResult:
        """
        Get dataset overview (~400-600 tokens).

        Returns summary statistics including genome/protein counts,
        annotation coverage, taxonomy distribution, and predicate summary.
        """
        return overview(self.store)

    def describe_schema(self) -> SharurResult:
        """
        Describe database schema and available predicates.

        Useful for understanding what data and operations are available.
        """
        return describe_schema(self.store)

    # ------------------------------------------------------------------ #
    # Navigation operators
    # ------------------------------------------------------------------ #

    def list_genomes(
        self,
        taxonomy_filter: Optional[str] = None,
        min_completeness: Optional[float] = None,
        max_contamination: Optional[float] = None,
        limit: int = 20,
        offset: int = 0,
    ) -> SharurResult:
        """
        List genomes (MAGs) with optional filtering.

        Args:
            taxonomy_filter: Filter by taxonomy substring (e.g., "Archaea")
            min_completeness: Minimum completeness percentage
            max_contamination: Maximum contamination percentage
            limit: Maximum results to return
            offset: Pagination offset
        """
        return list_genomes(
            self.store,
            taxonomy_filter=taxonomy_filter,
            min_completeness=min_completeness,
            max_contamination=max_contamination,
            limit=limit,
            offset=offset,
        )

    def list_proteins(
        self,
        genome_id: Optional[str] = None,
        contig_id: Optional[str] = None,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        has_annotation: Optional[bool] = None,
        limit: int = 50,
        offset: int = 0,
    ) -> SharurResult:
        """
        List proteins with optional filtering.

        Args:
            genome_id: Filter to specific genome (bin_id)
            contig_id: Filter to specific contig
            min_length: Minimum protein length (aa)
            max_length: Maximum protein length (aa)
            has_annotation: Filter by annotation status
            limit: Maximum results to return
            offset: Pagination offset
        """
        return list_proteins(
            self.store,
            genome_id=genome_id,
            contig_id=contig_id,
            min_length=min_length,
            max_length=max_length,
            has_annotation=has_annotation,
            limit=limit,
            offset=offset,
        )

    def get_genome(self, genome_id: str, verbosity: int = 1) -> SharurResult:
        """
        Get detailed information about a specific genome.

        Args:
            genome_id: Bin ID to retrieve
            verbosity: 0=minimal, 1=standard, 2=detailed
        """
        return get_genome(self.store, genome_id=genome_id, verbosity=verbosity)

    def get_protein(self, protein_id: str, verbosity: int = 1) -> SharurResult:
        """
        Get detailed information about a specific protein.

        Args:
            protein_id: Protein ID to retrieve
            verbosity: 0=minimal, 1=standard, 2=detailed
        """
        return get_protein(self.store, protein_id=protein_id, verbosity=verbosity)

    def get_neighborhood(
        self,
        entity_id: str,
        window: int = 10,
        verbosity: int = 1,
        all_annotations: bool = False,
    ) -> SharurResult:
        """
        Get genomic neighborhood around a protein.

        Args:
            entity_id: Protein ID as anchor
            window: Number of genes on each side
            verbosity: 0=minimal, 1=standard, 2=detailed
            all_annotations: If True, return all annotation sources per
                gene (PFAM, KEGG, VOGdb, DefenseFinder, HydDB, CAZy).
                Use this for context-based functional interpretation.
        """
        return get_neighborhood(
            self.store,
            entity_id=entity_id,
            window=window,
            verbosity=verbosity,
            all_annotations=all_annotations,
        )

    # ------------------------------------------------------------------ #
    # Search operators
    # ------------------------------------------------------------------ #

    def search_by_predicates(
        self,
        has: Optional[list[str]] = None,
        lacks: Optional[list[str]] = None,
        limit: int = 50,
        offset: int = 0,
    ) -> SharurResult:
        """
        Search proteins by predicate membership.

        Uses set intersection/difference logic:
        - has: Protein must have ALL of these predicates (AND)
        - lacks: Protein must have NONE of these predicates (AND NOT)

        Example:
            b.search_by_predicates(has=["giant", "unannotated"])
            b.search_by_predicates(has=["multi_domain"], lacks=["hypothetical"])

        Args:
            has: List of predicates that must be true
            lacks: List of predicates that must be false
            limit: Maximum results
            offset: Pagination offset
        """
        return search_by_predicates(
            self.store,
            has=has,
            lacks=lacks,
            limit=limit,
            offset=offset,
        )

    def search_proteins(
        self,
        annotation_pattern: Optional[str] = None,
        accession: Optional[str] = None,
        taxonomy_filter: Optional[str] = None,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        limit: int = 50,
        offset: int = 0,
    ) -> SharurResult:
        """
        Search proteins by annotation, accession, or taxonomy.

        Args:
            annotation_pattern: Pattern to match in annotation name/description
            accession: Exact accession to match (e.g., "PF00142")
            taxonomy_filter: Filter by genome taxonomy
            min_length: Minimum protein length (aa)
            max_length: Maximum protein length (aa)
            limit: Maximum results
            offset: Pagination offset
        """
        return search_proteins(
            self.store,
            annotation_pattern=annotation_pattern,
            accession=accession,
            taxonomy_filter=taxonomy_filter,
            min_length=min_length,
            max_length=max_length,
            limit=limit,
            offset=offset,
        )

    # ------------------------------------------------------------------ #
    # Export operators
    # ------------------------------------------------------------------ #

    def export_fasta(
        self,
        protein_ids: list[str],
        output_path: Optional[str] = None,
        include_annotations: bool = False,
    ) -> SharurResult:
        """
        Export proteins as FASTA format.

        Args:
            protein_ids: List of protein IDs to export
            output_path: Optional file path to write FASTA
            include_annotations: Include top annotation in header
        """
        return export_fasta(
            self.store,
            protein_ids=protein_ids,
            output_path=output_path,
            include_annotations=include_annotations,
        )

    def export_neighborhood_fasta(
        self,
        protein_id: str,
        window: int = 10,
        output_path: Optional[str] = None,
    ) -> SharurResult:
        """
        Export genomic neighborhood as multi-FASTA.

        Args:
            protein_id: Center protein ID
            window: Number of genes on each side
            output_path: Optional file path
        """
        return export_neighborhood_fasta(
            self.store,
            protein_id=protein_id,
            window=window,
            output_path=output_path,
        )

    def get_sequence(self, protein_id: str) -> SharurResult:
        """
        Get raw sequence for a single protein.

        Args:
            protein_id: Protein ID
        """
        return get_sequence(self.store, protein_id=protein_id)

    # ------------------------------------------------------------------ #
    # Similarity operators
    # ------------------------------------------------------------------ #

    def find_similar(
        self,
        protein_id: str,
        k: int = 10,
        threshold: Optional[float] = None,
        include_self_genome: bool = True,
    ) -> SharurResult:
        """
        Find proteins with similar ESM2 embeddings.

        Uses LanceDB kNN search to find structurally/functionally similar proteins.

        Args:
            protein_id: Query protein ID
            k: Number of similar proteins to return
            threshold: Minimum similarity score (0-1)
            include_self_genome: Include hits from the same genome
        """
        return find_similar(
            self.store,
            self.session.vector_store,
            protein_id=protein_id,
            k=k,
            threshold=threshold,
            include_self_genome=include_self_genome,
        )

    def find_similar_to_set(
        self,
        protein_ids: list[str],
        k: int = 10,
        threshold: Optional[float] = 0.7,
    ) -> SharurResult:
        """
        Find proteins similar to ANY in a set.

        Args:
            protein_ids: Set of query protein IDs
            k: Number of similar proteins per query
            threshold: Minimum similarity score
        """
        return find_similar_to_set(
            self.store,
            self.session.vector_store,
            protein_ids=protein_ids,
            k=k,
            threshold=threshold,
        )

    # ------------------------------------------------------------------ #
    # Visualization operators
    # ------------------------------------------------------------------ #

    def visualize_neighborhood(
        self,
        protein_id: str,
        window: int = 10,
        output_path: Optional[str] = None,
        figure_width: int = 14,
        title: Optional[str] = None,
        legend: Optional[str] = None,
        finding_id: Optional[int] = None,
        require_legend: bool = True,
    ) -> SharurResult:
        """
        Generate gene arrow diagram for a genomic neighborhood.

        Automatically records the figure in the manifest.

        Args:
            protein_id: Center protein ID
            window: Number of genes on each side
            output_path: Path to save PNG (None for temp file)
            figure_width: Width of figure in inches
            title: Figure title for manifest (recommended)
            legend: Figure legend/caption for manifest (REQUIRED by default)
            finding_id: Optional finding ID to associate with figure
            require_legend: If True (default), warn when legend is missing

        Note:
            Figure legends are essential for interpretable reports. Each legend
            should include: what the locus shows, color key, annotation sources,
            and the key biological observation.
        """
        # Warn if legend is missing (figures without legends are hard to interpret)
        if require_legend and not legend:
            import warnings
            warnings.warn(
                f"visualize_neighborhood({protein_id}): No legend provided. "
                "Figures without legends are difficult to interpret in reports. "
                "Please provide a legend explaining what the reader should observe.",
                UserWarning
            )

        result = visualize_neighborhood(
            self.store,
            protein_id=protein_id,
            window=window,
            output_path=output_path,
            figure_width=figure_width,
        )

        # Auto-update manifest with figure
        if result._raw and result._raw.get("output_path"):
            self.manifest.add_figure(
                path=result._raw["output_path"],
                figure_type="neighborhood",
                title=title,
                legend=legend,
                center_protein=protein_id,
                finding_id=finding_id,
                window=window,
                gene_count=result._raw.get("gene_count"),
            )
            self.manifest.save()

        return result

    def visualize_domains(
        self,
        protein_id: str,
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        legend: Optional[str] = None,
    ) -> SharurResult:
        """
        Generate domain architecture diagram for a protein.

        Automatically records the figure in the manifest.

        Args:
            protein_id: Protein ID
            output_path: Path to save image
            title: Optional figure title for manifest
            legend: Optional figure legend/caption for manifest
        """
        result = visualize_domain_architecture(
            self.store,
            protein_id=protein_id,
            output_path=output_path,
        )

        # Auto-update manifest with figure
        if result._raw and result._raw.get("output_path"):
            self.manifest.add_figure(
                path=result._raw["output_path"],
                figure_type="domain",
                title=title,
                legend=legend,
                center_protein=protein_id,
                domain_count=result._raw.get("domain_count"),
            )
            self.manifest.save()

        return result

    def get_pathway_context(self, protein_id: str) -> SharurResult:
        """
        Get KEGG pathway context for a protein.

        Args:
            protein_id: Protein ID
        """
        return get_kegg_pathway_context(self.store, protein_id=protein_id)

    # ------------------------------------------------------------------ #
    # Structure prediction operators
    # ------------------------------------------------------------------ #

    def predict_structure(
        self,
        protein_id: str,
        output_path: Optional[str] = None,
    ) -> SharurResult:
        """
        Predict protein structure using ESM3.

        Requires ESM_API_KEY environment variable.
        Limited to proteins <= 1024 aa for the open model.
        Automatically records the prediction in the manifest.

        Args:
            protein_id: Protein ID to predict structure for
            output_path: Optional path to save PDB file
        """
        result = predict_structure(
            self.store,
            protein_id=protein_id,
            output_path=output_path,
        )

        # Auto-update manifest with structure prediction
        if result._raw and result._raw.get("pdb_path"):
            raw = result._raw
            self.manifest.add_structure(
                protein_id=protein_id,
                pdb_path=raw["pdb_path"],
                metrics={
                    "length": raw.get("length"),
                    "plddt_mean": raw.get("plddt_mean"),
                    "ptm": raw.get("ptm"),
                },
            )
            self.manifest.save()

        return result

    def batch_predict_structures(
        self,
        protein_ids: list[str],
        output_dir: Optional[str] = None,
        max_length: int = 1024,
    ) -> SharurResult:
        """
        Predict structures for multiple proteins.

        Automatically records all predictions in the manifest.

        Args:
            protein_ids: List of protein IDs
            output_dir: Directory to save PDB files
            max_length: Skip proteins longer than this
        """
        result = batch_predict_structures(
            self.store,
            protein_ids=protein_ids,
            output_dir=output_dir,
            max_length=max_length,
        )

        # Auto-update manifest with all predictions
        if result._raw and result._raw.get("predictions"):
            for pred in result._raw["predictions"]:
                self.manifest.add_structure(
                    protein_id=pred["protein_id"],
                    pdb_path=pred["pdb_path"],
                    metrics={
                        "length": pred.get("length"),
                        "plddt_mean": pred.get("plddt_mean"),
                    },
                )
            self.manifest.save()

        return result

    # ------------------------------------------------------------------ #
    # Foldseek structure search operators
    # ------------------------------------------------------------------ #

    def search_foldseek(
        self,
        pdb_path: str,
        databases: Optional[list[str]] = None,
        top_k: int = 10,
    ) -> SharurResult:
        """
        Search a PDB structure against Foldseek databases.

        Args:
            pdb_path: Path to PDB file
            databases: Databases to search (default: afdb50, afdb-swissprot, pdb100)
            top_k: Number of top hits to return
        """
        return search_foldseek(pdb_path, databases, top_k)

    def search_foldseek_for_protein(
        self,
        protein_id: str,
        pdb_path: Optional[str] = None,
        databases: Optional[list[str]] = None,
        top_k: int = 10,
    ) -> SharurResult:
        """
        Search Foldseek for structural homologs of a protein.

        Requires an existing PDB file (run predict_structure first).

        Args:
            protein_id: Protein ID
            pdb_path: Optional path to PDB (auto-detects from /tmp/sharur_structures/)
            databases: Databases to search
            top_k: Number of hits to return
        """
        return search_foldseek_for_protein(
            self.store,
            protein_id=protein_id,
            pdb_path=pdb_path,
            databases=databases,
            top_k=top_k,
        )

    def list_foldseek_databases(self) -> SharurResult:
        """List available Foldseek databases."""
        return list_foldseek_databases()

    # ------------------------------------------------------------------ #
    # Validation operators
    # ------------------------------------------------------------------ #

    def validate_annotation(
        self,
        protein_ids: list[str],
        expected_function: str,
        domain_name: Optional[str] = None,
    ) -> SharurResult:
        """
        Validate that proteins match expected function.

        Cross-checks domain hits against protein annotations to detect
        potential annotation errors or domain misinterpretations.

        Args:
            protein_ids: List of protein IDs to validate
            expected_function: Expected function keyword (e.g., "CRISPR", "transposase")
            domain_name: Optional domain name that led to the function inference

        Example:
            # Found proteins with Cas12f1-like domain - validate they're actually CRISPR
            b.validate_annotation(cas12f_proteins, "CRISPR", domain_name="Cas12f1-like_TNB")
        """
        return validate_annotation(
            self.store,
            protein_ids=protein_ids,
            expected_function=expected_function,
            domain_name=domain_name,
        )

    def validate_context(
        self,
        protein_id: str,
        expected_neighbors: list[str],
        window: int = 10,
    ) -> SharurResult:
        """
        Check if a protein is in expected genomic context.

        Validates whether nearby genes match expected functional context
        (e.g., CRISPR proteins should be near cas genes and CRISPR arrays).

        Args:
            protein_id: Target protein ID
            expected_neighbors: List of keywords/predicates expected nearby
            window: Number of genes to check on each side

        Example:
            b.validate_context("protein_123", ["cas1", "cas2", "crispr_array"])
        """
        return validate_context(
            self.store,
            protein_id=protein_id,
            expected_neighbors=expected_neighbors,
            window=window,
        )

    def analyze_crispr_systems(self) -> SharurResult:
        """
        Analyze all CRISPR-Cas systems in the database.

        Identifies complete vs incomplete systems, orphan arrays,
        and systems affected by assembly fragmentation.

        Returns analysis of:
        - Total arrays and system types
        - Complete vs incomplete systems
        - Orphan arrays (no associated cas genes)
        - Arrays at contig edges (potential fragmentation)
        """
        return analyze_crispr_systems(self.store)

    def detect_annotation_errors(self, limit: int = 50) -> SharurResult:
        """
        Scan for likely annotation errors in the database.

        Checks for:
        1. Domain-annotation mismatches (especially for known problem domains)
        2. Proteins in wrong genomic context
        3. Spurious ORFs (overlapping CRISPR arrays)

        Args:
            limit: Maximum errors to report per category
        """
        return detect_annotation_errors(self.store, limit=limit)

    # ------------------------------------------------------------------ #
    # Hypothesis & Provenance operators
    # ------------------------------------------------------------------ #

    @property
    def hypothesis_registry(self):
        """Persistent cross-session hypothesis store.

        Hypotheses registered here survive across sessions and subagent runs.
        Storage: {dataset_dir}/exploration/hypotheses.json
        """
        if self._hypothesis_registry is None:
            from sharur.core.hypothesis_registry import HypothesisRegistry

            path = self._db_path.parent / "exploration" / "hypotheses.json"
            self._hypothesis_registry = HypothesisRegistry(path)
        return self._hypothesis_registry

    def propose_hypothesis(self, statement: str) -> "Hypothesis":
        """Propose a new hypothesis and register it in the persistent registry.

        Returns the Hypothesis object (use .hypothesis_id for later reference).

        Example:
            h = b.propose_hypothesis("Group 4 NiFe are energy-conserving")
            # later: b.add_evidence(h.hypothesis_id, ...)
        """
        from sharur.core.types import Hypothesis

        h = Hypothesis(statement=statement)
        self.hypothesis_registry.register(h)
        return h

    def add_evidence(
        self,
        hypothesis_id,
        query: str,
        result_summary: str,
        supports: bool,
        confidence: float,
    ) -> None:
        """Add evidence to an existing hypothesis.

        Args:
            hypothesis_id: UUID or string UUID of the hypothesis
            query: Description of the analytical step
            result_summary: What was found
            supports: True if evidence supports the hypothesis
            confidence: Confidence in this evidence (0-1)

        Example:
            b.add_evidence(
                h.hypothesis_id,
                query="Count NiFe Group 4 across genomes",
                result_summary="Found in 12/41 genomes",
                supports=True,
                confidence=0.8,
            )
        """
        import uuid as _uuid

        from sharur.core.types import Evidence

        if isinstance(hypothesis_id, str):
            hypothesis_id = _uuid.UUID(hypothesis_id)
        ev = Evidence(
            query=query,
            result_summary=result_summary,
            supports=supports,
            confidence=confidence,
        )
        self.hypothesis_registry.add_evidence(hypothesis_id, ev)

    def list_hypotheses(self) -> list:
        """List all hypotheses from the persistent registry."""
        return self.hypothesis_registry.list_all()

    def hypothesis_summary(self) -> str:
        """Human-readable summary of all hypotheses with evidence counts."""
        return self.hypothesis_registry.summary()

    def render_provenance(
        self,
        title: Optional[str] = None,
        output_path: Optional[str] = None,
    ) -> str:
        """Render the provenance DAG as a Mermaid diagram.

        Combines session provenance entries with persistent hypothesis
        registry to produce a publication-ready figure.

        Args:
            title: Optional title for the diagram
            output_path: Optional file path to write the Mermaid output

        Returns:
            Mermaid-format string
        """
        from sharur.core.provenance_renderer import render_provenance_mermaid

        mermaid = render_provenance_mermaid(
            session=self.session,
            registry=self._hypothesis_registry,
            title=title,
        )

        if output_path is not None:
            p = Path(output_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(mermaid)

        return mermaid

    def log_provenance(
        self,
        query: str,
        result_summary: str,
        parent_ids: Optional[list] = None,
    ) -> "ProvenanceEntry":
        """Log a provenance entry manually.

        Most operators log provenance automatically. Use this for custom
        analytical steps or when chaining provenance explicitly.

        Args:
            query: Description of what was done
            result_summary: What was found
            parent_ids: Optional list of parent entry UUIDs (UUID or str)
                for DAG chaining

        Returns:
            The ProvenanceEntry (use .entry_id for parent chaining)

        Example:
            e1 = b.log_provenance("Count hydrogenases", "42 found")
            e2 = b.log_provenance("Check neighborhoods", "12 with maturation",
                                  parent_ids=[e1.entry_id])
        """
        import uuid as _uuid

        if parent_ids:
            parent_ids = [
                _uuid.UUID(p) if isinstance(p, str) else p for p in parent_ids
            ]
        return self.session.log_query(
            query=query,
            tool_calls=[],
            results_summary=result_summary,
            duration_ms=0,
            parent_ids=parent_ids,
        )

    # ------------------------------------------------------------------ #
    # V2 Semantic Atom operators
    # ------------------------------------------------------------------ #

    def generate_v2(
        self,
        protein_ids: Optional[list[str]] = None,
        output_review_queue: Optional[str] = None,
    ) -> dict:
        """Generate V2 semantic atoms + states and persist to DuckDB.

        Full pipeline: reads proteins + annotations, generates typed
        SemanticAtoms, aggregates into SemanticState per protein,
        evaluates composite predicates, writes to DB tables.

        Args:
            protein_ids: Optional subset (None = all proteins).
            output_review_queue: Path to write unmapped-accession TSV.

        Returns:
            Dict mapping protein_id -> SemanticState.

        Example:
            states = b.generate_v2()
            states = b.generate_v2(output_review_queue="review_queue.tsv")
        """
        from sharur.predicates_v2.persistence import generate_and_persist_v2

        return generate_and_persist_v2(
            self.store,
            protein_ids=protein_ids,
            output_review_queue=output_review_queue,
        )

    def get_semantic_state(self, protein_id: str) -> Optional["SemanticState"]:
        """Get the V2 semantic state for a protein.

        Reads from the semantic_state table. Requires generate_v2() first.

        Args:
            protein_id: Protein ID.

        Returns:
            SemanticState dataclass, or None if not found.

        Example:
            state = b.get_semantic_state("protein_123")
            print(state.activities)   # ["hydrogenase", "nife_hydrogenase"]
            print(state.roles)        # ["defense_system"]
            print(state.size_class)   # "giant"
            print(state.composite_predicates)  # ["energy_conserving_hydrogenase"]
        """
        import json

        from sharur.predicates_v2.model import SemanticState

        rows = self.store.execute(
            "SELECT activities, roles, architecture, localization, "
            "topology, size_class, quality_flags, composite_predicates, "
            "unresolved_count FROM semantic_state WHERE protein_id = ?",
            [protein_id],
        )
        if not rows:
            return None
        row = rows[0]
        topo = json.loads(row[4]) if isinstance(row[4], str) else (row[4] or {})
        return SemanticState(
            protein_id=protein_id,
            activities=list(row[0] or []),
            roles=list(row[1] or []),
            architecture=list(row[2] or []),
            localization=list(row[3] or []),
            topology=topo,
            size_class=row[5] or "",
            quality_flags=list(row[6] or []),
            composite_predicates=list(row[7] or []),
        )

    def get_atoms(self, protein_id: str) -> list:
        """Get all V2 semantic atoms for a protein.

        Reads from the semantic_atoms table. Requires generate_v2() first.

        Args:
            protein_id: Protein ID.

        Returns:
            List of SemanticAtom dataclasses.

        Example:
            atoms = b.get_atoms("protein_123")
            for a in atoms:
                print(f"{a.atom_id} ({a.facet.value}): {a.relation.value} via {a.source_db}")
        """
        from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet

        rows = self.store.execute(
            "SELECT atom_id, facet, relation, source_accession, "
            "source_db, evidence_evalue, evidence_score "
            "FROM semantic_atoms WHERE protein_id = ?",
            [protein_id],
        )
        return [
            SemanticAtom(
                protein_id=protein_id,
                atom_id=row[0],
                facet=SemanticFacet(row[1]),
                relation=ClaimRelation(row[2]),
                source_accession=row[3],
                source_db=row[4],
                evidence_evalue=row[5],
                evidence_score=row[6],
            )
            for row in rows
        ]

    def search_by_facet(
        self,
        facet: str,
        atom_ids: Optional[list[str]] = None,
        relation: Optional[str] = None,
        limit: int = 50,
    ) -> list:
        """Search proteins by V2 semantic facet.

        Queries the semantic_atoms table for proteins with atoms in
        the given facet, optionally filtered by specific atom IDs
        and/or relation strength.

        Args:
            facet: Facet name (activity, role, architecture, localization,
                   topology, size_class, quality_flag).
            atom_ids: Optional list of specific atoms to filter on.
            relation: Optional relation filter (implies, supports, flags,
                      excludes).
            limit: Maximum results.

        Returns:
            List of (protein_id, atom_id, relation) tuples.

        Example:
            # All proteins with any activity atom
            b.search_by_facet("activity")

            # Proteins with hydrogenase activity, strong evidence only
            b.search_by_facet("activity", atom_ids=["hydrogenase"], relation="implies")

            # All defense-role proteins
            b.search_by_facet("role", atom_ids=["defense_system"])
        """
        params: list = [facet]
        clauses = ["facet = ?"]

        if atom_ids:
            placeholders = ",".join(["?"] * len(atom_ids))
            clauses.append(f"atom_id IN ({placeholders})")
            params.extend(atom_ids)

        if relation:
            clauses.append("relation = ?")
            params.append(relation)

        where = " AND ".join(clauses)
        rows = self.store.execute(
            f"SELECT DISTINCT protein_id, atom_id, relation "
            f"FROM semantic_atoms WHERE {where} "
            f"ORDER BY protein_id LIMIT ?",
            params + [limit],
        )
        return [tuple(row) for row in rows]

    def search_by_atoms(
        self,
        has: Optional[list[str]] = None,
        lacks: Optional[list[str]] = None,
        limit: int = 50,
    ) -> list[str]:
        """Search proteins by V2 atom or composite predicate presence/absence.

        V2-native analog of search_by_predicates. Searches both
        semantic_atoms and composite_predicates in semantic_state.

        Args:
            has: Atoms/composites that must ALL be present (AND logic).
            lacks: Atoms/composites that must ALL be absent (AND NOT logic).
            limit: Maximum results.

        Returns:
            List of protein_id strings.

        Example:
            # Giant unannotated proteins
            b.search_by_atoms(has=["giant", "unannotated"])

            # Hydrogenases that are NOT membrane proteins
            b.search_by_atoms(has=["hydrogenase"], lacks=["membrane"])

            # Composite predicates work too
            b.search_by_atoms(has=["nife_hydrogenase_validated"])
        """
        has = has or []
        lacks = lacks or []
        if not has and not lacks:
            return []

        params: list = []

        # Search both semantic_atoms and composite_predicates
        # UNION gives us a unified view of atom_ids + composite names
        unified = (
            "SELECT protein_id, atom_id FROM semantic_atoms "
            "UNION ALL "
            "SELECT protein_id, UNNEST(composite_predicates) as atom_id "
            "FROM semantic_state WHERE LEN(composite_predicates) > 0"
        )

        if has:
            placeholders = ",".join(["?"] * len(has))
            base = (
                f"SELECT protein_id FROM ({unified}) u "
                f"WHERE atom_id IN ({placeholders}) "
                f"GROUP BY protein_id "
                f"HAVING COUNT(DISTINCT atom_id) = ?"
            )
            params.extend(has)
            params.append(len(has))
        else:
            base = f"SELECT DISTINCT protein_id FROM ({unified}) u"

        if lacks:
            placeholders = ",".join(["?"] * len(lacks))
            query = (
                f"WITH base AS ({base}) "
                f"SELECT protein_id FROM base "
                f"WHERE protein_id NOT IN ("
                f"  SELECT DISTINCT protein_id FROM ({unified}) u2 "
                f"  WHERE atom_id IN ({placeholders})"
                f") LIMIT ?"
            )
            params.extend(lacks)
            params.append(limit)
        else:
            query = f"{base} LIMIT ?"
            params.append(limit)

        rows = self.store.execute(query, params)
        return [row[0] for row in rows]

    def list_composites(self) -> list:
        """List available composite predicate definitions.

        Returns composite definitions from config/predicates_v2/composites.yaml.
        Each has name, description, facet, and requires (DSL condition).

        Returns:
            List of CompositeDefinition objects.

        Example:
            for comp in b.list_composites():
                print(f"{comp.name}: {comp.description}")
        """
        from sharur.predicates_v2.composites import load_composites

        return load_composites()

    def v2_review_queue(self, limit: int = 50) -> list[dict]:
        """Get unmapped accession review queue.

        Queries semantic_atoms for unresolved atoms and ranks by
        frequency. Use this to identify annotation accessions that
        need curation rules added to the mapping files.

        Requires generate_v2() to have been run first.

        Args:
            limit: Maximum entries.

        Returns:
            List of dicts with accession, source_db, n_proteins.

        Example:
            queue = b.v2_review_queue(limit=20)
            for entry in queue:
                print(f"{entry['accession']} ({entry['source_db']}): "
                      f"{entry['n_proteins']} proteins")
        """
        rows = self.store.execute(
            "SELECT source_accession, source_db, "
            "COUNT(DISTINCT protein_id) as n_proteins "
            "FROM semantic_atoms WHERE relation = 'unresolved' "
            "GROUP BY source_accession, source_db "
            "ORDER BY n_proteins DESC LIMIT ?",
            [limit],
        )
        return [
            {"accession": row[0], "source_db": row[1], "n_proteins": row[2]}
            for row in rows
        ]

    def run_shadow_diff(
        self,
        protein_ids: Optional[list[str]] = None,
        output_report: Optional[str] = None,
        output_jsonl: Optional[str] = None,
    ) -> dict:
        """Run V1 vs V2 shadow comparison.

        Generates both V1 flat predicates and V2 atoms for the same
        proteins and compares the outputs. Useful for validating V2
        coverage before switching.

        Args:
            protein_ids: Optional subset (None = all).
            output_report: Path to write Markdown report.
            output_jsonl: Path to write per-protein JSONL diffs.

        Returns:
            Dict with 'per_protein' diffs and 'summary' statistics.

        Example:
            diff = b.run_shadow_diff(output_report="shadow_report.md")
            print(f"Match rate: {diff['summary']['match_rate_pct']}%")
        """
        from sharur.predicates_v2.shadow import run_shadow_diff_on_store

        return run_shadow_diff_on_store(
            self.store,
            protein_ids=protein_ids,
            output_report=output_report,
            output_jsonl=output_jsonl,
        )


__all__ = [
    # Main facade
    "Sharur",
    # Base classes
    "SharurResult",
    "ResultMeta",
    "OperatorTrace",
    "OperatorContext",
    # Introspection
    "overview",
    "describe_schema",
    # Navigation
    "list_genomes",
    "list_proteins",
    "get_genome",
    "get_protein",
    "get_neighborhood",
    # Search
    "search_by_predicates",
    "search_proteins",
    # Export
    "export_fasta",
    "export_neighborhood_fasta",
    "get_sequence",
    # Similarity
    "find_similar",
    "find_similar_to_set",
    # Visualization
    "visualize_neighborhood",
    "visualize_domain_architecture",
    "get_kegg_pathway_context",
    # Structure prediction
    "predict_structure",
    "predict_structure_from_sequence",
    "batch_predict_structures",
    # Foldseek
    "search_foldseek",
    "search_foldseek_for_protein",
    "list_foldseek_databases",
    # Validation
    "PROBLEM_DOMAINS",
    "validate_annotation",
    "validate_context",
    "analyze_crispr_systems",
    "detect_annotation_errors",
    # Manifest
    "AnalysisManifest",
]
