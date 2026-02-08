"""BennuAgent orchestrator with DSPy param extraction and safe fallback."""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import os

import dspy
import typer

from bennu.agent import signatures as sig
from bennu.agent.router import BennuRouter
from bennu.core.session import ExplorationSession
from bennu.tools import (
    CompareAcrossTool,
    DetectLociTool,
    ExportTool,
    FindAnomaliesTool,
    FindProteinsTool,
    FindSimilarTool,
    GetContextTool,
    ManageSetsTool,
    ToolRegistry,
)
from bennu.tools.compare_across import CompareAcrossParams as PydCompareAcrossParams
from bennu.tools.detect_loci import DetectLociParams as PydDetectLociParams
from bennu.tools.export import ExportParams as PydExportParams
from bennu.tools.find_anomalies import FindAnomaliesParams as PydFindAnomaliesParams
from bennu.tools.find_proteins import FindProteinsParams as PydFindProteinsParams
from bennu.tools.find_similar import FindSimilarParams as PydFindSimilarParams
from bennu.tools.get_context import GetContextParams as PydGetContextParams
from bennu.tools.manage_sets import ManageSetsParams as PydManageSetsParams
from bennu.tools.registry import ToolResult


class BennuAgent:
    """
    Main agent orchestrating query processing.

    Flow:
    1. Router analyzes query (DSPy if available)
    2. For each tool needed:
       a. Extract parameters via TypedPredictor (fallback to defaults)
       b. Execute tool
    3. Return synthesized textual summary (LLM synthesis can be added later)
    """

    def __init__(self, session: ExplorationSession, lm: Optional[dspy.LM] = None, debug_lm: bool = False):
        self.session = session
        self.debug_lm = debug_lm
        if lm:
            dspy.configure(lm=lm)
        self.use_llm = lm is not None

        self.registry = ToolRegistry()
        self.router = BennuRouter(use_llm=self.use_llm, debug=self.debug_lm)
        self._predictor_factory = None  # injectable for tests
        self._register_default_tools()

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def process(self, query: str) -> str:
        routing: sig.QueryAnalysis = self.router.forward(query, self.session, self.registry)
        tools_needed = self._coerce_tools_needed(routing.tools_needed)
        self._debug(f"Route: tools={tools_needed}")

        tool_results: List[ToolResult] = []
        for tool_name in tools_needed:
            tool = self.registry.get(tool_name)
            params = self._extract_params(tool_name, query)
            self._debug(f"Tool={tool_name} params={params}")
            if hasattr(params, "rationale") and getattr(params, "rationale", None):
                self._debug(f"Tool={tool_name} rationale={getattr(params, 'rationale')}")
            result = tool.execute(params, self.session)
            tool_results.append(result)

        return self._synthesize_response(query, tool_results)

    # ------------------------------------------------------------------ #
    # Internals
    # ------------------------------------------------------------------ #
    def _register_default_tools(self) -> None:
        self.registry.register(FindProteinsTool())
        self.registry.register(GetContextTool())
        self.registry.register(ManageSetsTool())
        self.registry.register(FindSimilarTool())
        self.registry.register(FindAnomaliesTool())
        self.registry.register(DetectLociTool())
        self.registry.register(CompareAcrossTool())
        self.registry.register(ExportTool())

    def _extract_params(self, tool_name: str, query: str):
        """
        Use DSPy TypedPredictor (signatures) to extract, then coerce into
        the tool's Pydantic params. Falls back to defaults on failure.
        """
        sig_cls, pyd_cls = self._signature_mapping(tool_name)

        # Build keyword args expected by signatures
        kwargs: Dict = {"query": query, "session_context": self._session_context()}
        annotations = getattr(sig_cls, "__annotations__", {})
        if "available_tools" in annotations:
            kwargs["available_tools"] = "\n".join(t["name"] for t in self.registry.list_tools())
        if "tool_name" in annotations:
            kwargs["tool_name"] = tool_name
        if "tool_schema" in annotations:
            kwargs["tool_schema"] = pyd_cls.model_json_schema()
        # Provide focus/set hints without concrete examples to reduce bias
        if "focus_hint" in annotations:
            kwargs["focus_hint"] = self._focus_hint()

        if self.use_llm:
            predictor = None
            if self._predictor_factory:
                predictor = self._predictor_factory(sig_cls)
            elif hasattr(dspy, "TypedPredictor"):
                predictor = dspy.TypedPredictor(sig_cls)  # type: ignore[attr-defined]

            if predictor:
                self._debug(f"LLM call → {tool_name}")
                sig_obj = predictor(**kwargs)
                payload = sig_obj.__dict__
                if self.debug_lm:
                    self._debug(f"LLM raw payload: {payload}")
                filtered = {k: v for k, v in payload.items() if k in pyd_cls.model_fields}
                if self.debug_lm:
                    self._debug(f"LLM filtered payload: {filtered}")
                params = pyd_cls(**filtered)
                return self._postprocess_params(tool_name, params, query)
            else:
                raise RuntimeError("LLM predictor not available")

        # No LLM configured: raise instead of silently falling back
        raise RuntimeError("LLM is required but not configured")

    @staticmethod
    def _signature_mapping(tool_name: str) -> Tuple[Optional[type], type]:
        """Map tool name to (DSPy signature, Pydantic params)."""
        mapping = {
            "find_proteins": (sig.FindProteinsParams, PydFindProteinsParams),
            "get_context": (sig.GetContextParams, PydGetContextParams),
            "manage_sets": (sig.ManageSetsParams, PydManageSetsParams),
            "find_similar": (sig.FindSimilarParams, PydFindSimilarParams),
            "find_anomalies": (sig.FindAnomaliesParams, PydFindAnomaliesParams),
            "detect_loci": (sig.DetectLociParams, PydDetectLociParams),
            "compare_across": (sig.CompareAcrossParams, PydCompareAcrossParams),
            "export": (sig.ExportParams, PydExportParams),
        }
        if tool_name not in mapping:
            raise ValueError(f"No signature mapping for tool {tool_name}")
        return mapping[tool_name]

    def _session_context(self) -> str:
        focus = self.session.get_focus() if hasattr(self.session, "get_focus") else None
        sets = self.session.list_sets() if hasattr(self.session, "list_sets") else []
        set_names = [s.name for s in sets if hasattr(s, "name")]
        return f"focus={focus.entity_id if focus else 'none'}; sets={set_names}"

    def _focus_hint(self) -> str:
        """Provide minimal focus hint without concrete IDs to reduce bias."""
        focus = self.session.get_focus() if hasattr(self.session, "get_focus") else None
        if not focus:
            return "no current focus"
        return f"current entity type={focus.entity_type}"

    def _synthesize_response(self, query: str, tool_results: List[ToolResult]) -> str:
        lines = [f"Query: {query}"]
        for tr in tool_results:
            lines.append(f"- {tr.summary}")
            if tr.warnings:
                lines.append(f"  warnings: {', '.join(tr.warnings)}")
        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    # Post-processing helpers
    # ------------------------------------------------------------------ #
    def _postprocess_params(self, tool_name: str, params, query: Optional[str] = None):
        """Fill missing fields from session context where sensible."""
        focus_id = None
        focus_type = None
        working_sets = []
        if hasattr(self.session, "get_focus"):
            focus = self.session.get_focus()
            if focus:
                focus_id = focus.entity_id
                focus_type = getattr(focus, "entity_type", None)
        if hasattr(self.session, "list_sets"):
            working_sets = self.session.list_sets()
        first_set_name = working_sets[0].name if working_sets else None

        if tool_name in {"manage_sets"}:
            # Default set target from focus or first available set
            if getattr(params, "set_name", None) in (None, "",) and focus_type == "set" and focus_id:
                setattr(params, "set_name", focus_id)
            if getattr(params, "set_name", None) in (None, "",) and first_set_name:
                setattr(params, "set_name", first_set_name)
            # If a protein is focused and no proteins provided, seed with focus protein
            if getattr(params, "protein_ids", None) == [] and focus_type == "protein" and focus_id:
                setattr(params, "protein_ids", [focus_id])

        if tool_name in {"export"}:
            if getattr(params, "source", None) in (None, ""):
                setattr(params, "source", "set")
            if getattr(params, "source_id", None) in (None, "",) and focus_type == "set" and focus_id:
                setattr(params, "source_id", focus_id)
            if getattr(params, "source_id", None) in (None, "",) and first_set_name:
                setattr(params, "source_id", first_set_name)

        if tool_name in {"get_context"}:
            if getattr(params, "protein_id", None) in (None, "",):
                setattr(params, "protein_id", focus_id or "")
            # Prefer gene window if neither set
            if getattr(params, "window_genes", None) in (None, 0):
                setattr(params, "window_genes", 10)

        if tool_name in {"find_similar"}:
            if getattr(params, "query_id", None) in (None, "",):
                setattr(params, "query_id", focus_id or "")
            if getattr(params, "n", None) in (None, 0):
                setattr(params, "n", 20)

        if tool_name in {"find_proteins"}:
            if hasattr(params, "similar_to") and getattr(params, "similar_to", None) in (None, "",) and focus_id:
                setattr(params, "similar_to", focus_id)
            if getattr(params, "in_set", None) in (None, "",) and first_set_name:
                setattr(params, "in_set", first_set_name)

            # Clean function keywords (drop stopwords/noise)
            stopwords = {
                "dataset",
                "data",
                "protein",
                "proteins",
                "find",
                "me",
                "this",
                "that",
                "these",
                "those",
                "in",
                "on",
                "for",
                "all",
            }
            funcs = []
            for f in getattr(params, "functions", []) or []:
                clean = f.strip().lower()
                if clean and clean not in stopwords:
                    funcs.append(clean)
            # de-duplicate preserving order
            seen = set()
            funcs = [x for x in funcs if not (x in seen or seen.add(x))]
            setattr(params, "functions", funcs)

            # Heuristic mapping: integrase keyword -> Pfam domains if none given
            if "integrase" in funcs and getattr(params, "domains", []) == []:
                setattr(params, "domains", ["PF00589", "PF09003"])

        if tool_name in {"find_anomalies"}:
            if getattr(params, "scope", None) in (None, ""):
                setattr(params, "scope", "genome")
            if params.scope == "genome" and focus_type in {"bin", "contig"}:
                setattr(params, "scope", focus_type)
                if focus_id:
                    setattr(params, "scope_id", focus_id)
            if getattr(params, "scope_id", None) in (None, "",) and first_set_name and params.scope == "set":
                setattr(params, "scope_id", first_set_name)
            if getattr(params, "scope_id", None) in (None, "",) and focus_id and params.scope in {"bin", "contig"}:
                setattr(params, "scope_id", focus_id)

        if tool_name in {"detect_loci"}:
            if getattr(params, "scope", None) in (None, ""):
                setattr(params, "scope", "genome")
            if params.scope == "genome" and focus_type in {"bin", "contig"}:
                setattr(params, "scope", focus_type)
                if focus_id:
                    setattr(params, "scope_id", focus_id)
            if getattr(params, "scope_id", None) in (None, "",) and focus_id and params.scope in {"bin", "contig"}:
                setattr(params, "scope_id", focus_id)
            if getattr(params, "scope_id", None) in (None, "",) and params.scope == "set" and first_set_name:
                setattr(params, "scope_id", first_set_name)

        return params

    def _coerce_tools_needed(self, value):
        """
        Ensure tools_needed is a list of registered tool names.
        Accepts:
          - list[str]
          - JSON array string
          - JSON object with "value": [...]
        Raises ValueError if parsing fails or tools are unknown.
        """
        allowed = {t["name"] for t in self.registry.list_tools()}

        tools_list = None
        if isinstance(value, list):
            tools_list = value
        elif isinstance(value, str):
            import json

            try:
                parsed = json.loads(value)
                if isinstance(parsed, dict) and "value" in parsed and isinstance(parsed["value"], list):
                    tools_list = parsed["value"]
                elif isinstance(parsed, list):
                    tools_list = parsed
            except Exception:
                tools_list = None

        if not tools_list or not isinstance(tools_list, list):
            raise ValueError(f"Invalid tools_needed from router: {value!r}")

        # validate tool names
        tools_list = [t for t in tools_list if isinstance(t, str)]
        unknown = [t for t in tools_list if t not in allowed]
        if unknown:
            raise ValueError(f"Router returned unknown tools: {unknown}")
        return tools_list

    # ------------------------------------------------------------------ #
    # Heuristics
    # ------------------------------------------------------------------ #
    def _heuristic_params(self, tool_name: str, query: str, pyd_cls):
        """Very simple keyword-based extraction for no-LLM mode."""
        if tool_name != "find_proteins":
            return pyd_cls()

        q = query.lower()
        params = pyd_cls()

        tokens = [t.strip(" .!?") for t in q.replace(",", " ").split() if t.strip(" .!?")]

        # Pfam accessions like PFxxxxx
        pfams = []
        for token in tokens:
            up = token.upper()
            if up.startswith("PF") and up[2:].isdigit():
                pfams.append(up)
        if pfams:
            setattr(params, "domains", pfams)

        # function keywords: naive stemming (drop trailing s/es)
        keywords = []
        for token in tokens:
            if len(token) < 5 or not token.isalpha():
                continue
            kw = token.lower()
            if kw.endswith("ses") or kw.endswith("ases"):
                kw = kw[:-1]  # drop single trailing s (integrases -> integrase)
            elif kw.endswith("es"):
                kw = kw[:-2]
            elif kw.endswith("s"):
                kw = kw[:-1]
            keywords.append(kw)
        if keywords:
            setattr(params, "functions", list(dict.fromkeys(keywords)))

        return params

    # ------------------------------------------------------------------ #
    # Debug logging
    # ------------------------------------------------------------------ #
    def _debug(self, msg: str):
        typer.echo(f"[debug] {msg}", err=True)


__all__ = ["BennuAgent"]
