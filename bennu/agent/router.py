"""DSPy-powered BennuRouter with heuristic fallback (Part 6.2)."""

from __future__ import annotations

from typing import List
import typer

import dspy

from bennu.agent.signatures import QueryAnalysis
from bennu.core.session import ExplorationSession
from bennu.tools.registry import ToolRegistry


class BennuRouter:
    """Routes queries to appropriate tools."""

    def __init__(self, use_llm: bool = True, debug: bool = False):
        self.use_llm = use_llm
        self.debug = debug
        self.analyzer = dspy.ChainOfThought(QueryAnalysis) if use_llm else None

    def forward(self, query: str, session: ExplorationSession, registry: ToolRegistry) -> QueryAnalysis:
        session_context = self._summarize_session(session)
        available_tools = self._format_tools(registry)

        if self.use_llm and self.analyzer:
            try:
                qa = self.analyzer(query=query, session_context=session_context, available_tools=available_tools)
                if self.debug:
                    typer.echo(f"[debug] Router LLM payload: {qa.__dict__}", err=True)
                return qa
            except Exception:
                # fall back to heuristic
                pass

        return self._heuristic(query, available_tools)

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _summarize_session(session: ExplorationSession) -> str:
        parts = []
        sets = getattr(session, "list_sets", lambda: [])()
        if sets:
            parts.append(f"Working sets: {', '.join(s.name for s in sets)}")
        focus = session.get_focus() if hasattr(session, "get_focus") else None
        if focus:
            parts.append(f"Current focus: {focus.entity_type} {focus.entity_id}")
        prov = session.get_provenance() if hasattr(session, "get_provenance") else []
        if prov:
            parts.append("Recent queries: " + "; ".join(p.query for p in prov[-3:]))
        return "\n".join(parts) if parts else "New session."

    @staticmethod
    def _format_tools(registry: ToolRegistry) -> str:
        return "\n".join(f"- {t['name']}: {t['description']}" for t in registry.list_tools())

    @staticmethod
    def _heuristic(query: str, available_tools: str) -> QueryAnalysis:
        q_lower = query.lower()
        tools_needed: List[str] = []

        if any(k in q_lower for k in ["near", "context", "around", "neighbors"]):
            tools_needed.append("get_context")
        elif "similar" in q_lower:
            tools_needed.append("find_similar")
        elif "anomal" in q_lower or "weird" in q_lower:
            tools_needed.append("find_anomalies")
        elif "compare" in q_lower:
            tools_needed.append("compare_across")
        elif "export" in q_lower or "fasta" in q_lower or "tsv" in q_lower or "notebook" in q_lower:
            tools_needed.append("export")
        elif "locus" in q_lower or "bgc" in q_lower or "prophage" in q_lower:
            tools_needed.append("detect_loci")
        else:
            tools_needed.append("find_proteins")

        return QueryAnalysis(
            intent=f"Handle user query: {query}",
            tools_needed=tools_needed,
            reasoning="Heuristic routing based on keywords (LLM unavailable).",
            query=query,
            session_context="",
            available_tools=available_tools,
        )


__all__ = ["BennuRouter"]
