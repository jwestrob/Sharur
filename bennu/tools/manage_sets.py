"""Working set management tool (Part 4.8)."""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field

from bennu.core.types import WorkingSet
from bennu.core.session import ExplorationSession
from bennu.tools.registry import Tool, ToolResult


class ManageSetsParams(BaseModel):
    action: str = Field(
        default="list",
        description="create|add|remove|delete|list|info|union|intersect|difference",
    )
    set_name: Optional[str] = None
    protein_ids: list[str] = Field(default_factory=list)
    description: Optional[str] = None
    other_set: Optional[str] = None


class ManageSetsTool:
    name = "manage_sets"
    description = "Create, modify, and list working sets of proteins."
    tier: int = 1

    def execute(self, params: ManageSetsParams, session: ExplorationSession) -> ToolResult:
        action = params.action.lower()
        summary = ""
        data = None
        warnings: list[str] = []

        try:
            if action == "create":
                ws = session.create_set(
                    name=params.set_name or "unnamed_set",
                    protein_ids=params.protein_ids,
                    description=params.description,
                )
                data = ws
                summary = f"Created set '{ws.name}' with {len(ws.protein_ids)} proteins."

            elif action == "add":
                session.add_to_set(params.set_name, params.protein_ids)
                ws = session.get_set(params.set_name)
                data = ws
                summary = f"Added {len(params.protein_ids)} proteins to '{params.set_name}'. Now {len(ws.protein_ids)}."

            elif action == "remove":
                session.remove_from_set(params.set_name, params.protein_ids)
                ws = session.get_set(params.set_name)
                data = ws
                summary = f"Removed {len(params.protein_ids)} proteins from '{params.set_name}'. Now {len(ws.protein_ids)}."

            elif action == "delete":
                session.delete_set(params.set_name)
                data = None
                summary = f"Deleted set '{params.set_name}'."

            elif action == "list":
                sets = session.list_sets()
                data = sets
                summary = f"{len(sets)} working sets."

            elif action == "info":
                ws = session.get_set(params.set_name)
                if not ws:
                    raise KeyError(f"Set '{params.set_name}' not found")
                data = ws
                summary = f"Set '{ws.name}': {len(ws.protein_ids)} proteins."

            elif action in {"union", "intersect", "difference"}:
                ws1 = session.get_set(params.set_name)
                ws2 = session.get_set(params.other_set or "")
                if not ws1 or not ws2:
                    raise KeyError("Both sets must exist for set operations")
                set1 = set(ws1.protein_ids)
                set2 = set(ws2.protein_ids)
                if action == "union":
                    result_ids = sorted(set1 | set2)
                elif action == "intersect":
                    result_ids = sorted(set1 & set2)
                else:
                    result_ids = sorted(set1 - set2)
                new_name = params.description or f"{params.set_name}_{action}_{params.other_set}"
                ws = session.create_set(new_name, result_ids, description=f"{action} of {ws1.name} and {ws2.name}")
                data = ws
                summary = f"Created set '{ws.name}' with {len(ws.protein_ids)} proteins ({action})."

            else:
                raise ValueError(f"Unknown action '{params.action}'")

            count = len(data.protein_ids) if isinstance(data, WorkingSet) else (len(data) if data else 0)
            return ToolResult(success=True, data=data, summary=summary, count=count, warnings=warnings, suggestions=[])

        except Exception as exc:
            return ToolResult(
                success=False,
                data=None,
                summary=f"manage_sets failed: {exc}",
                count=0,
                warnings=[str(exc)],
                suggestions=[],
            )


__all__ = ["ManageSetsTool", "ManageSetsParams"]
