"""Tool registry and ToolResult (Part 4.1)."""

from __future__ import annotations

from typing import Any, Dict, Literal, Protocol, runtime_checkable

from pydantic import BaseModel


class ToolResult(BaseModel):
    """Standard result format for all tools."""

    success: bool
    data: Any
    summary: str
    count: int = 0
    truncated: bool = False
    warnings: list[str] = []
    suggestions: list[str] = []


@runtime_checkable
class Tool(Protocol):
    """Protocol for Bennu tools."""

    name: str
    description: str
    tier: Literal[1, 2, 3]  # 1=Semantic, 2=Query, 3=Sandbox

    def execute(self, params: BaseModel, session: Any) -> ToolResult:
        ...


class ToolRegistry:
    """Simple registry for tool instances."""

    def __init__(self) -> None:
        self._tools: Dict[str, Tool] = {}

    def register(self, tool: Tool) -> None:
        self._tools[tool.name] = tool

    def get(self, name: str) -> Tool:
        return self._tools[name]

    def list_tools(self) -> list[dict]:
        return [{"name": t.name, "description": t.description, "tier": t.tier} for t in self._tools.values()]


__all__ = ["ToolResult", "Tool", "ToolRegistry"]
