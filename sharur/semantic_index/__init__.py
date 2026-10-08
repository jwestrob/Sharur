"""Compact semantic-term index: membership postings, original rich rows and genome scopes.

An opt-in alternative to scanning ``semantic_terms`` for V2 term search and per-protein
term explanations. See :mod:`sharur.semantic_index.provider` for the shared contract,
:mod:`sharur.semantic_index.generation` for generation binding and selection, and
:mod:`sharur.semantic_index.build` for builders. The compact readers need the optional
``pyroaring`` dependency (``pip install 'sharur[compact]'``); the SQL provider does not.
"""

from sharur.semantic_index.generation import (
    Generation,
    GenerationError,
    StaleGenerationError,
    assemble,
    resolve,
    select,
    verify,
)
from sharur.semantic_index.provider import (
    CompactSemanticProvider,
    SearchPage,
    SearchRequest,
    SqlSemanticProvider,
    attach,
    attached,
    provider_for,
)


__all__ = [
    "CompactSemanticProvider",
    "Generation",
    "GenerationError",
    "SearchPage",
    "SearchRequest",
    "SqlSemanticProvider",
    "StaleGenerationError",
    "assemble",
    "attach",
    "attached",
    "provider_for",
    "resolve",
    "select",
    "verify",
]
