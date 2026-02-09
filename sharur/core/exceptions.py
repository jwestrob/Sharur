"""
Custom exceptions for Sharur.

All Sharur-specific exceptions inherit from SharurError.
"""


class SharurError(Exception):
    """Base exception for all Sharur errors."""

    pass


# =============================================================================
# DATA ERRORS
# =============================================================================


class DataError(SharurError):
    """Base class for data-related errors."""

    pass


class SchemaError(DataError):
    """Database schema is invalid or missing."""

    pass


class DataNotFoundError(DataError):
    """Requested data does not exist."""

    def __init__(self, entity_type: str, entity_id: str):
        self.entity_type = entity_type
        self.entity_id = entity_id
        super().__init__(f"{entity_type} not found: {entity_id}")


class DataValidationError(DataError):
    """Data fails validation constraints."""

    pass


# =============================================================================
# QUERY ERRORS
# =============================================================================


class QueryError(SharurError):
    """Base class for query-related errors."""

    pass


class InvalidQueryError(QueryError):
    """Query is malformed or invalid."""

    pass


class QueryTimeoutError(QueryError):
    """Query exceeded time limit."""

    def __init__(self, query: str, timeout_seconds: float):
        self.query = query
        self.timeout_seconds = timeout_seconds
        super().__init__(
            f"Query timed out after {timeout_seconds}s: {query[:100]}..."
        )


class QueryTooExpensiveError(QueryError):
    """Query estimated cost exceeds threshold."""

    def __init__(self, estimated_rows: int, threshold: int):
        self.estimated_rows = estimated_rows
        self.threshold = threshold
        super().__init__(
            f"Query would scan {estimated_rows:,} rows (threshold: {threshold:,})"
        )


# =============================================================================
# SESSION ERRORS
# =============================================================================


class SessionError(SharurError):
    """Base class for session-related errors."""

    pass


class WorkingSetError(SessionError):
    """Error with working set operations."""

    pass


class WorkingSetNotFoundError(WorkingSetError):
    """Referenced working set does not exist."""

    def __init__(self, name: str):
        self.name = name
        super().__init__(f"Working set not found: {name}")


class WorkingSetExistsError(WorkingSetError):
    """Working set with this name already exists."""

    def __init__(self, name: str):
        self.name = name
        super().__init__(f"Working set already exists: {name}")


class ReferenceResolutionError(SessionError):
    """Could not resolve a natural language reference."""

    def __init__(self, reference: str):
        self.reference = reference
        super().__init__(f"Could not resolve reference: '{reference}'")


# =============================================================================
# TOOL ERRORS
# =============================================================================


class ToolError(SharurError):
    """Base class for tool execution errors."""

    pass


class ToolNotFoundError(ToolError):
    """Referenced tool does not exist."""

    def __init__(self, tool_name: str, available_tools: list[str]):
        self.tool_name = tool_name
        self.available_tools = available_tools
        super().__init__(
            f"Unknown tool: {tool_name}. Available: {', '.join(available_tools)}"
        )


class ToolExecutionError(ToolError):
    """Tool failed during execution."""

    def __init__(self, tool_name: str, message: str, cause: Exception | None = None):
        self.tool_name = tool_name
        self.cause = cause
        super().__init__(f"Tool '{tool_name}' failed: {message}")


class ToolValidationError(ToolError):
    """Tool parameters failed validation."""

    def __init__(self, tool_name: str, param_name: str, message: str):
        self.tool_name = tool_name
        self.param_name = param_name
        super().__init__(
            f"Invalid parameter '{param_name}' for tool '{tool_name}': {message}"
        )


# =============================================================================
# VECTOR STORE ERRORS
# =============================================================================


class VectorStoreError(SharurError):
    """Base class for vector store errors."""

    pass


class EmbeddingNotFoundError(VectorStoreError):
    """Requested embedding does not exist."""

    def __init__(self, protein_id: str, index_name: str):
        self.protein_id = protein_id
        self.index_name = index_name
        super().__init__(
            f"Embedding not found for {protein_id} in index {index_name}"
        )


class IndexNotFoundError(VectorStoreError):
    """Vector index does not exist."""

    def __init__(self, index_name: str):
        self.index_name = index_name
        super().__init__(f"Vector index not found: {index_name}")


# =============================================================================
# AGENT ERRORS
# =============================================================================


class AgentError(SharurError):
    """Base class for agent errors."""

    pass


class RoutingError(AgentError):
    """Could not determine appropriate tool for query."""

    def __init__(self, query: str, reason: str):
        self.query = query
        self.reason = reason
        super().__init__(f"Could not route query: {reason}")


class SynthesisError(AgentError):
    """Failed to synthesize results into response."""

    pass


class ClarificationNeededError(AgentError):
    """Query is ambiguous and requires clarification."""

    def __init__(self, query: str, question: str):
        self.query = query
        self.question = question
        super().__init__(f"Clarification needed: {question}")


# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    # Base
    "SharurError",
    # Data
    "DataError",
    "SchemaError",
    "DataNotFoundError",
    "DataValidationError",
    # Query
    "QueryError",
    "InvalidQueryError",
    "QueryTimeoutError",
    "QueryTooExpensiveError",
    # Session
    "SessionError",
    "WorkingSetError",
    "WorkingSetNotFoundError",
    "WorkingSetExistsError",
    "ReferenceResolutionError",
    # Tool
    "ToolError",
    "ToolNotFoundError",
    "ToolExecutionError",
    "ToolValidationError",
    # Vector
    "VectorStoreError",
    "EmbeddingNotFoundError",
    "IndexNotFoundError",
    # Agent
    "AgentError",
    "RoutingError",
    "SynthesisError",
    "ClarificationNeededError",
]
