"""Agent preamble bundler.

Subagents launched via the Task tool receive ONLY the `prompt` parameter. They
do NOT inherit CLAUDE.md, the project's MEMORY.md, the skill spec, or any
feedback memories the parent learned during the session.

This is the leak that causes Claude(child) to reinvent mistakes that
Claude(parent) was already corrected on — per-chunk mmseqs2, raw DefenseFinder
HMM hits, KOFAM on the login node, the iconic-member trap, etc. The parent has
the lesson; the child does not.

`bundle(skill_name)` builds a single string preamble that captures:

  1. Skill spec content from `.claude/skills/<skill_name>.md`
  2. Shared validation protocols (`.claude/skills/_validation_protocols.md`)
  3. Core rules and SLURM rules from `CLAUDE.md`
  4. Feedback memories — every `feedback_*.md` in the project's memory dir
  5. (optional) Dataset context — counts, columns, schema sketch

Prepend the result to your subagent's `prompt` argument before dispatch.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

# Default project root (works when imported from within the Sharur install).
# Override via SHARUR_REPO env var if you call this from elsewhere.
DEFAULT_REPO = Path(os.environ.get(
    "SHARUR_REPO",
    "/groups/banfield/users/jwestrob/bin/Sharur",
))

# Default memory dir for the Sharur project. Override via SHARUR_MEMORY env var.
DEFAULT_MEMORY_DIR = Path(os.environ.get(
    "SHARUR_MEMORY",
    "/home/jwestrob/jwestrob/.claude/projects/"
    "-groups-banfield-users-jwestrob-bin-Sharur/memory",
))

# User-level CLAUDE.md (research-track context, vault pointers, hardware notes,
# folding/ELSA conventions). Auto-injected for the main Claude Code session but
# NOT inherited by Task() subagents — bundle pulls it in.
DEFAULT_USER_CLAUDE_MD = Path(os.environ.get(
    "SHARUR_USER_CLAUDE_MD",
    str(Path.home() / ".claude" / "CLAUDE.md"),
))


# Sections of CLAUDE.md worth including verbatim in every subagent preamble.
# Headings to match (level-2 markdown headings under "Core Rules" and the SLURM
# block). Listed in the order they appear in CLAUDE.md.
CORE_HEADINGS = (
    "Disposition",
    "Data Integrity",
    "Reproducible Findings",
    "Database Queries",
    "DuckDB Performance",
    "Scientific Rigor",
    "Genome Completeness",
    "Check for Functional Detail",
)


@dataclass
class BundleParts:
    """Structured view of what went into a bundle — useful for debugging."""

    skill: Optional[str]
    skill_chars: int
    validation_chars: int
    core_rules_chars: int
    feedback_files: list[str]
    feedback_chars: int
    dataset_context_chars: int
    total_chars: int


def _read(path: Path) -> str:
    try:
        return path.read_text()
    except (FileNotFoundError, PermissionError):
        return ""


def _extract_core_rules(claude_md: str) -> str:
    """Pull the level-3 headings listed in CORE_HEADINGS, plus the SLURM block.

    CLAUDE.md uses ### for sub-rules under '## Core Rules (always apply)'.
    The SLURM block is at '## SLURM Rules — READ EVERY TIME'.
    """
    out: list[str] = []

    # Core rules sub-sections (### headings under ## Core Rules)
    in_core = False
    current_heading: Optional[str] = None
    buffer: list[str] = []
    for line in claude_md.splitlines():
        if line.startswith("## "):
            if in_core and buffer:
                out.append("\n".join(buffer).rstrip())
                buffer = []
            in_core = "Core Rules" in line
            current_heading = None
            continue
        if not in_core:
            continue
        if line.startswith("### "):
            # Flush previous heading if it was wanted
            if current_heading and buffer:
                out.append("\n".join(buffer).rstrip())
                buffer = []
            heading_text = line.removeprefix("### ").strip()
            if heading_text in CORE_HEADINGS:
                current_heading = heading_text
                buffer = [line]
            else:
                current_heading = None
                buffer = []
            continue
        if current_heading:
            buffer.append(line)
    if current_heading and buffer:
        out.append("\n".join(buffer).rstrip())

    # SLURM rules block (entire ## section)
    slurm_match = re.search(
        r"^## SLURM Rules.*?(?=^## )",
        claude_md,
        flags=re.MULTILINE | re.DOTALL,
    )
    if slurm_match:
        out.append(slurm_match.group(0).rstrip())

    return "\n\n".join(out).strip()


def _list_feedback_memories(memory_dir: Path) -> list[Path]:
    if not memory_dir.exists():
        return []
    return sorted(memory_dir.glob("feedback_*.md"))


def _format_feedback(paths: Iterable[Path]) -> str:
    """Pack feedback memories into a compact section.

    Each memory is short. We strip frontmatter and prefix with the slug so
    Claude can recall the source if needed.
    """
    chunks: list[str] = []
    for p in paths:
        body = _read(p)
        # Strip YAML frontmatter
        body = re.sub(r"^---\s*\n.*?\n---\s*\n", "", body, count=1, flags=re.DOTALL)
        slug = p.stem.removeprefix("feedback_")
        chunks.append(f"### {slug}\n{body.strip()}")
    return "\n\n".join(chunks).strip()


def _dataset_context(dataset_path: Path) -> str:
    """Generate a small context block describing the dataset.

    Currently lazy: just reports presence of key files. A future version can
    inspect the duckdb schema and write a one-line summary per table.
    """
    if not dataset_path.exists():
        return f"(Dataset path {dataset_path} does not exist.)"

    interesting = [
        "sharur.duckdb", "manifest.json", "sharur_ops.db",
        "survey", "exploration", "annotations", "embeddings",
        "structures", "synteny", "stage00_prepared", "stage03_prodigal",
        "stage04_astra", "stage05c_crispr",
    ]
    present, absent = [], []
    for name in interesting:
        if (dataset_path / name).exists():
            present.append(name)
        else:
            absent.append(name)
    lines = [
        f"**Dataset path:** `{dataset_path}`",
        f"**Present:** {', '.join(present) if present else '(none of the expected files)'}",
    ]
    if absent:
        lines.append(f"**Not yet generated:** {', '.join(absent)}")
    return "\n".join(lines)


def bundle(
    skill_name: Optional[str] = None,
    *,
    dataset_path: Optional[Path | str] = None,
    repo_root: Optional[Path | str] = None,
    memory_dir: Optional[Path | str] = None,
    user_claude_md: Optional[Path | str] = None,
    include_memories: bool = True,
    include_skill: bool = True,
    include_core: bool = True,
    include_validation: bool = True,
    include_user_claude_md: bool = True,
) -> str:
    """Build a subagent preamble blob.

    Args:
        skill_name: short name like "defense", "hydrogenase", "explore". The
            corresponding `.claude/skills/{skill_name}.md` is included verbatim.
            If None, no skill spec is included.
        dataset_path: path to the dataset directory (`data/<name>`). If given,
            a small "what's present in this dataset" block is appended.
        repo_root: Sharur repo path. Defaults to SHARUR_REPO env or hardcoded.
        memory_dir: project memory dir. Defaults to SHARUR_MEMORY env or hardcoded.
        include_memories: include feedback_*.md memories. Default True.
        include_skill: include the skill spec. Default True.
        include_core: include CLAUDE.md core rules + SLURM block. Default True.
        include_validation: include `_validation_protocols.md`. Default True.

    Returns:
        A single string ready to prepend to a Task() prompt. Sections are
        separated by horizontal rules so Claude reads them as distinct context.
    """
    repo = Path(repo_root) if repo_root else DEFAULT_REPO
    mem = Path(memory_dir) if memory_dir else DEFAULT_MEMORY_DIR
    user_md = Path(user_claude_md) if user_claude_md else DEFAULT_USER_CLAUDE_MD

    parts: list[str] = [
        "# Sharur Subagent Preamble",
        "",
        "You are a subagent dispatched within the Sharur metagenomic exploration "
        "system. The blocks below are mandatory context — read them BEFORE acting "
        "on your task. They were learned by the parent agent and would otherwise "
        "not reach you.",
    ]

    if include_user_claude_md:
        user_md_text = _read(user_md)
        if user_md_text:
            parts.append("\n---\n")
            parts.append(
                "## User-Level CLAUDE.md (`~/.claude/CLAUDE.md`)\n"
                "Auto-injected for the parent session, NOT inherited by Task() "
                "subagents. Includes: user profile, vault layout, hardware/SLURM "
                "conventions, related software (fold, ELSA, etc.). When the task "
                "touches Omnitrophota / SR-VP / DPANN / giant proteins / "
                "lanthanide binding / structural biology, the vault project "
                "directories pointed at here are likely the deepest context.\n"
            )
            parts.append(user_md_text)

    if include_core:
        claude_md = _read(repo / "CLAUDE.md")
        core = _extract_core_rules(claude_md) if claude_md else ""
        if core:
            parts.append("\n---\n")
            parts.append("## Project Core Rules (from `<repo>/CLAUDE.md`)\n")
            parts.append(core)

    if include_validation:
        vproto = _read(repo / ".claude/skills/_validation_protocols.md")
        if vproto:
            parts.append("\n---\n")
            parts.append("## Shared Validation Protocols\n")
            parts.append(vproto)

        # V2 is now the normal predicate backend (commit ebc141a). Subagents
        # need the V2 data-model overview so they can query the right tables
        # and read atom relations correctly instead of falling back to V1's
        # flat boolean predicates.
        v2_doc = _read(repo / "docs/predicates_v2.md")
        if v2_doc:
            parts.append("\n---\n")
            parts.append("## Predicate V2 (the normal backend)\n")
            parts.append(v2_doc)

    if include_skill and skill_name:
        spec = _read(repo / f".claude/skills/{skill_name}.md")
        if spec:
            parts.append("\n---\n")
            parts.append(f"## Skill Spec: {skill_name}\n")
            parts.append(spec)

    if include_memories:
        feedback_paths = _list_feedback_memories(mem)
        if feedback_paths:
            parts.append("\n---\n")
            parts.append(
                "## Feedback Memories (lessons the parent already learned)\n"
                "These are session-spanning corrections the user gave the parent. "
                "Apply them; do not relearn the lesson the hard way."
            )
            parts.append(_format_feedback(feedback_paths))

    if dataset_path is not None:
        parts.append("\n---\n")
        parts.append("## Dataset Context\n")
        parts.append(_dataset_context(Path(dataset_path)))

    parts.append("\n---\n")
    parts.append(
        "## Task\n"
        "(The specific task follows below this preamble.)"
    )

    return "\n".join(parts)


def bundle_parts(*args, **kwargs) -> BundleParts:
    """Like `bundle`, but also returns a sizing report.

    Useful for figuring out token budget before dispatching.
    """
    # Re-implementation that records sizes per section
    repo = Path(kwargs.get("repo_root") or DEFAULT_REPO)
    mem = Path(kwargs.get("memory_dir") or DEFAULT_MEMORY_DIR)
    skill_name = kwargs.get("skill_name") or (args[0] if args else None)

    claude_md = _read(repo / "CLAUDE.md")
    core_text = _extract_core_rules(claude_md) if claude_md else ""
    vproto = _read(repo / ".claude/skills/_validation_protocols.md")
    skill_text = _read(repo / f".claude/skills/{skill_name}.md") if skill_name else ""
    feedback_paths = _list_feedback_memories(mem)
    feedback_text = _format_feedback(feedback_paths)
    dataset_context_text = (
        _dataset_context(Path(kwargs["dataset_path"]))
        if kwargs.get("dataset_path") else ""
    )

    blob = bundle(*args, **kwargs)
    return BundleParts(
        skill=skill_name,
        skill_chars=len(skill_text),
        validation_chars=len(vproto),
        core_rules_chars=len(core_text),
        feedback_files=[p.name for p in feedback_paths],
        feedback_chars=len(feedback_text),
        dataset_context_chars=len(dataset_context_text),
        total_chars=len(blob),
    )


__all__ = ["bundle", "bundle_parts", "BundleParts", "DEFAULT_REPO", "DEFAULT_MEMORY_DIR"]
