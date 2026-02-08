"""
DSPy signatures for Bennu agent.

These signatures define the LLM interfaces for:
- Query routing (understanding what the user wants)
- Parameter extraction (converting natural language to tool parameters)
- Result synthesis (explaining results in natural language)

Each signature is a contract: inputs -> outputs with typed fields.
DSPy optimizers can then tune prompts to maximize performance.
"""

from typing import Any, Optional
import json

import dspy


# =============================================================================
# QUERY ANALYSIS AND ROUTING
# =============================================================================


class QueryAnalysis(dspy.Signature):
    """Analyze a user query to decide intent and required tools."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Summary of session state")
    available_tools: str = dspy.InputField(desc="Tools and their descriptions")

    intent: str = dspy.OutputField(desc="What the user wants")
    tools_needed: list[str] = dspy.OutputField(
        desc=(
            'Respond with tool names from available_tools in a plaintext comma-separated list, no brackets. No prose.'
        ),
        format=lambda x: x if isinstance(x, str) else json.dumps(x),
        parser=lambda x: (
            x
            if isinstance(x, list)
            else (
                (json.loads(x).get("value", [])) if isinstance(json.loads(x), dict) else json.loads(x)
            )
        ),
    )
    reasoning: str = dspy.OutputField(desc="Why those tools fit (one short sentence)")


# =============================================================================
# TOOL PARAMETER EXTRACTION
# =============================================================================


class FindProteinsParams(dspy.Signature):
    """Extract parameters for the find_proteins tool. Must include a non-empty rationale."""

    query: str = dspy.InputField(desc="User wording of the search")
    session_context: str = dspy.InputField(desc="Session state summary")
    focus_hint: str = dspy.InputField(desc="Type of current focus entity, if any")

    domains: list[str] = dspy.OutputField(desc="Pfam accessions only (PFxxxxx); leave empty if none implied")
    functions: list[str] = dspy.OutputField(
        desc=(
            "Query tokens for text matching; include ALL alphanumeric tokens from the query that could identify a protein—"
            "enzyme names, subunit codes, complex names, descriptive terms; keep each as a SEPARATE list item; "
            "exclude ONLY request verbs (find, show, get, list); empty if none implied"
        )
    )
    kegg: list[str] = dspy.OutputField(desc="KEGG ortholog IDs only; empty if none implied")
    taxonomy: str = dspy.OutputField(desc="GTDB-style taxonomic pattern; empty if none implied")
    similar_to: str = dspy.OutputField(desc="Anchor protein ID; empty if none")
    similarity_level: str = dspy.OutputField(desc="protein or locus")
    near_protein: str = dspy.OutputField(desc="Protein to center context on; empty if none")
    window_genes: int = dspy.OutputField(desc="Genes in context window")
    in_set: str = dspy.OutputField(desc="Working set to include; empty if none")
    exclude_set: str = dspy.OutputField(desc="Working set to exclude; empty if none")
    limit: int = dspy.OutputField(desc="Max results (0 or empty means no cap)")
    rationale: str = dspy.OutputField(desc="Required: one-sentence rationale for chosen parameters; do not leave empty")


class GetContextParams(dspy.Signature):
    """Parameters for get_context (neighborhood lookup)."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Session state summary")
    focus_hint: str = dspy.InputField(desc="Type of current focus entity, if any")

    protein_id: str = dspy.OutputField(desc="Protein ID to center on")
    window_genes: int = dspy.OutputField(desc="Genes on each side")
    window_bp: int = dspy.OutputField(desc="Window size in bp (0 if unused)")
    include_annotations: bool = dspy.OutputField(desc="Include annotation details")
    expand_operon: bool = dspy.OutputField(desc="Expand to operon if True")


class DetectLociParams(dspy.Signature):
    """Parameters for detect_loci."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Session state summary")

    locus_type: str = dspy.OutputField(desc="Requested locus type")
    scope: str = dspy.OutputField(desc="genome | bin | contig")
    scope_id: str = dspy.OutputField(desc="ID for scoped search")
    required_domains: list[str] = dspy.OutputField(desc="Required domains (custom loci)")
    optional_domains: list[str] = dspy.OutputField(desc="Optional domains (custom loci)")


class FindSimilarParams(dspy.Signature):
    """Parameters for find_similar."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Session state summary")
    focus_hint: str = dspy.InputField(desc="Type of current focus entity, if any")

    query_id: str = dspy.OutputField(desc="Anchor protein ID")
    level: str = dspy.OutputField(desc="protein or locus")
    n: int = dspy.OutputField(desc="How many neighbors")
    threshold: float = dspy.OutputField(desc="Similarity threshold")
    exclude_same_bin: bool = dspy.OutputField(desc="Exclude same bin if True")
    taxonomy_filter: str = dspy.OutputField(desc="GTDB filter string")


class FindAnomaliesParams(dspy.Signature):
    """Parameters for find_anomalies."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Session state summary")

    scope: str = dspy.OutputField(desc="genome | bin | set")
    scope_id: str = dspy.OutputField(desc="ID for scoped search")
    anomaly_types: list[str] = dspy.OutputField(desc="Types to include")
    signals: list[str] = dspy.OutputField(desc="Signals of interest")
    min_score: float = dspy.OutputField(desc="Score threshold")
    limit: int = dspy.OutputField(desc="Max results")


class CompareAcrossParams(dspy.Signature):
    """Parameters for compare_across."""

    query: str = dspy.InputField(desc="User query text")
    session_context: str = dspy.InputField(desc="Session state summary")

    feature: str = dspy.OutputField(desc="Feature type to compare")
    feature_id: str = dspy.OutputField(desc="Specific feature id (optional)")
    group_by: str = dspy.OutputField(desc="Grouping mode")
    taxonomy_level: int = dspy.OutputField(desc="Taxonomy depth level")
    metric: str = dspy.OutputField(desc="Comparison metric")


class ManageSetsParams(dspy.Signature):
    """
    Extract parameters for the manage_sets tool.
    
    This tool manages working sets for iterative analysis.
    """
    
    query: str = dspy.InputField(
        desc="The user's query about managing sets"
    )
    session_context: str = dspy.InputField(
        desc="Session state including current working sets"
    )
    focus_hint: str = dspy.InputField(
        desc="Type of current focus entity, if any"
    )
    
    action: str = dspy.OutputField(
        desc="Action: 'create', 'add', 'remove', 'delete', 'list', 'info', 'union', 'intersect', 'difference'"
    )
    set_name: str = dspy.OutputField(
        desc="Name of the set to operate on. Empty string for 'list' action."
    )
    protein_ids: list[str] = dspy.OutputField(
        desc="Protein IDs for create/add/remove. Usually resolved from focus stack."
    )
    description: str = dspy.OutputField(
        desc="Description for new set. Empty string if not specified."
    )
    other_set: str = dspy.OutputField(
        desc="Second set for union/intersect/difference. Empty string if not applicable."
    )


class ExportParams(dspy.Signature):
    """
    Extract parameters for the export tool.
    
    This tool exports data in various formats.
    """
    
    query: str = dspy.InputField(
        desc="The user's query about exporting"
    )
    session_context: str = dspy.InputField(
        desc="Session state"
    )
    focus_hint: str = dspy.InputField(
        desc="Type of current focus entity, if any"
    )
    
    format: str = dspy.OutputField(
        desc="Export format: 'fasta', 'gff', 'tsv', 'json', 'notebook'"
    )
    source: str = dspy.OutputField(
        desc="What to export: 'set', 'query', 'loci'"
    )
    source_id: str = dspy.OutputField(
        desc="Set name or query ID. May use focus stack if not specified."
    )
    include_sequences: bool = dspy.OutputField(
        desc="Whether to include protein sequences. Default True."
    )
    include_annotations: bool = dspy.OutputField(
        desc="Whether to include annotations. Default True."
    )


# =============================================================================
# RESULT SYNTHESIS
# =============================================================================


class ResultSynthesis(dspy.Signature):
    """
    Synthesize tool results into a natural language response.
    
    This is the final step: converting structured results into
    a helpful, conversational response for the user.
    """
    
    query: str = dspy.InputField(
        desc="The original user query"
    )
    tool_results: str = dspy.InputField(
        desc="JSON string of tool execution results"
    )
    session_context: str = dspy.InputField(
        desc="Current session state"
    )
    
    response: str = dspy.OutputField(
        desc="Natural language response to the user. Should directly answer their question."
    )
    key_findings: list[str] = dspy.OutputField(
        desc="Bullet points of most important findings (2-5 items)"
    )
    follow_ups: list[str] = dspy.OutputField(
        desc="Suggested follow-up questions (2-3 items)"
    )
    warnings: list[str] = dspy.OutputField(
        desc="Any caveats or warnings about the results"
    )


class AnomalyExplanation(dspy.Signature):
    """
    Generate human-readable explanations for anomalies.
    
    Takes raw anomaly scores and genomic context, produces
    a clear explanation of why the protein is unusual.
    """
    
    protein_id: str = dspy.InputField(
        desc="The anomalous protein ID"
    )
    anomaly_signals: str = dspy.InputField(
        desc="JSON list of anomaly signals with scores"
    )
    genomic_context: str = dspy.InputField(
        desc="Description of surrounding genes and annotations"
    )
    
    explanation: str = dspy.OutputField(
        desc="Plain-English explanation of why this protein is anomalous"
    )
    biological_significance: str = dspy.OutputField(
        desc="Why this might be interesting biologically"
    )
    suggested_investigations: list[str] = dspy.OutputField(
        desc="What to look at next (2-3 suggestions)"
    )


# =============================================================================
# HYPOTHESIS MANAGEMENT
# =============================================================================


class HypothesisExtraction(dspy.Signature):
    """
    Extract implicit hypotheses from user queries and results.
    
    Users often make implicit claims ("This looks like a predator")
    that should be tracked as hypotheses.
    """
    
    query: str = dspy.InputField(
        desc="The user's query"
    )
    response: str = dspy.InputField(
        desc="The system's response"
    )
    session_context: str = dspy.InputField(
        desc="Current session state including existing hypotheses"
    )
    
    new_hypotheses: list[str] = dspy.OutputField(
        desc="New hypothesis statements extracted from the conversation"
    )
    evidence_updates: list[dict] = dspy.OutputField(
        desc="Updates to existing hypotheses: [{'hypothesis_id': str, 'supports': bool, 'reason': str}]"
    )


# =============================================================================
# ERROR HANDLING
# =============================================================================


class ErrorRecovery(dspy.Signature):
    """
    Generate helpful responses when tool execution fails.
    
    Should explain what went wrong and suggest alternatives.
    """
    
    query: str = dspy.InputField(
        desc="The user's original query"
    )
    error_message: str = dspy.InputField(
        desc="The error that occurred"
    )
    tool_name: str = dspy.InputField(
        desc="Which tool failed"
    )
    session_context: str = dspy.InputField(
        desc="Session state"
    )
    
    explanation: str = dspy.OutputField(
        desc="User-friendly explanation of what went wrong"
    )
    suggestions: list[str] = dspy.OutputField(
        desc="Alternative approaches the user could try"
    )
    can_retry: bool = dspy.OutputField(
        desc="Whether a retry with different parameters might help"
    )


# =============================================================================
# MODULE DEFINITIONS (for DSPy optimizers)
# =============================================================================


class BennuRouter(dspy.Module):
    """
    Router module that analyzes queries and dispatches to tools.
    
    Uses ChainOfThought for explicit reasoning about tool selection.
    """
    
    def __init__(self):
        super().__init__()
        self.analyze = dspy.ChainOfThought(QueryAnalysis)
    
    def forward(self, query: str, session_context: str, available_tools: str):
        return self.analyze(
            query=query,
            session_context=session_context,
            available_tools=available_tools
        )


class BennuSynthesizer(dspy.Module):
    """
    Synthesizer module that converts results to natural language.
    
    Uses ChainOfThought for reasoning about how to present results.
    """
    
    def __init__(self):
        super().__init__()
        self.synthesize = dspy.ChainOfThought(ResultSynthesis)
    
    def forward(self, query: str, tool_results: str, session_context: str):
        return self.synthesize(
            query=query,
            tool_results=tool_results,
            session_context=session_context
        )


# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    # Routing
    "QueryAnalysis",
    
    # Parameter extraction
    "FindProteinsParams",
    "GetContextParams",
    "DetectLociParams",
    "FindSimilarParams",
    "FindAnomaliesParams",
    "CompareAcrossParams",
    "ManageSetsParams",
    "ExportParams",
    
    # Synthesis
    "ResultSynthesis",
    "AnomalyExplanation",
    
    # Hypothesis
    "HypothesisExtraction",
    
    # Error handling
    "ErrorRecovery",
    
    # Modules
    "BennuRouter",
    "BennuSynthesizer",
]
