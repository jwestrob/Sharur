"""Headless model-CLI drivers for Sharur worker executors.

One driver per provider in the review policy. Each takes a frozen, bounded,
sequence-free prompt and returns a `ModelRun` carrying the parsed structured
record plus the provenance the task contract requires: stderr, exit status,
model metadata, prompt hash, and result hash.

Design note — the packet is fetched by the worker and piped in on **stdin**,
never fetched by the model as a tool call. Two reasons: tool output is subject
to truncation, and a tool-calling turn decouples the provider request count
from the Atlas packet count (one packet should be ~one inference turn).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


class ModelError(RuntimeError):
    """Model invocation failed in a way that is not a rate limit."""


class ModelRateLimited(RuntimeError):
    """Provider refused the call for quota/rate reasons. Retry after backoff."""


@dataclass
class ModelRun:
    """One model invocation and everything needed to audit it later."""

    provider: str
    model: str
    reasoning_effort: str
    record: dict[str, Any]
    prompt_sha256: str
    result_sha256: str
    exit_status: int
    stderr: str
    usage: dict[str, Any] = field(default_factory=dict)
    raw_stdout_bytes: int = 0

    def provenance(self) -> dict[str, Any]:
        return {
            "provider": self.provider,
            "model": self.model,
            "reasoning_effort": self.reasoning_effort,
            "prompt_sha256": self.prompt_sha256,
            "result_sha256": self.result_sha256,
            "exit_status": self.exit_status,
            "stderr_tail": self.stderr[-2000:] if self.stderr else "",
            "usage": self.usage,
        }


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


_RATE_LIMIT_MARKERS = (
    "rate limit",
    "rate_limit",
    "429",
    "quota",
    "usage limit",
    "too many requests",
    "resource_exhausted",
)


def _looks_rate_limited(stderr: str, stdout: str) -> bool:
    blob = f"{stderr}\n{stdout}".lower()
    return any(marker in blob for marker in _RATE_LIMIT_MARKERS)


def _extract_json(text: str) -> dict[str, Any]:
    """Pull the structured record out of a CLI's stdout.

    Tries whole-text JSON first, then the last balanced ``{...}`` span. Models
    occasionally wrap JSON in prose even under a schema, so this stays lenient
    rather than failing a whole genome on a stray sentence.
    """
    text = text.strip()
    if not text:
        raise ModelError("model returned empty stdout")
    try:
        parsed = json.loads(text)
        if isinstance(parsed, dict):
            return parsed
    except json.JSONDecodeError:
        pass

    depth = 0
    start = -1
    best: str | None = None
    for idx, ch in enumerate(text):
        if ch == "{":
            if depth == 0:
                start = idx
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0 and start >= 0:
                best = text[start : idx + 1]
    if best is None:
        raise ModelError(f"no JSON object in model stdout (first 400 chars): {text[:400]}")
    try:
        parsed = json.loads(best)
    except json.JSONDecodeError as exc:
        raise ModelError(f"model stdout was not valid JSON: {exc}") from exc
    if not isinstance(parsed, dict):
        raise ModelError("model returned a non-object JSON value")
    return parsed


def _run(argv: list[str], stdin_text: str, timeout: int) -> tuple[int, str, str]:
    proc = subprocess.run(
        argv,
        input=stdin_text,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=os.environ.copy(),
    )
    return proc.returncode, proc.stdout, proc.stderr


def run_openai(
    *,
    model: str,
    reasoning_effort: str,
    system_prompt: str,
    payload_text: str,
    output_schema: dict[str, Any],
    timeout: int = 1800,
) -> ModelRun:
    """Drive `codex exec` headlessly, stdin-fed, with a forced output schema."""
    codex = shutil.which("codex")
    if codex is None:
        raise ModelError("codex CLI not found on PATH")

    prompt = f"{system_prompt}\n\n=== GENOME PACKET (JSON) ===\n{payload_text}\n"

    with tempfile.TemporaryDirectory(prefix="sharur-scan-") as tmp:
        schema_path = Path(tmp) / "schema.json"
        schema_path.write_text(json.dumps(output_schema), encoding="utf-8")
        argv = [
            codex,
            "exec",
            "--skip-git-repo-check",
            "-s",
            "read-only",
            "--json",
            "-m",
            model,
            "-c",
            f"model_reasoning_effort={reasoning_effort}",
            "--output-schema",
            str(schema_path),
            "-",
        ]
        code, stdout, stderr = _run(argv, prompt, timeout)

    if code != 0 and _looks_rate_limited(stderr, stdout):
        raise ModelRateLimited(stderr.strip()[-500:] or "codex reported a rate limit")

    # `--json` emits JSONL events; the record is in the final agent_message and
    # usage is on turn.completed. Fall back to raw stdout if events are absent.
    usage: dict[str, Any] = {}
    message_text: str | None = None
    for line in stdout.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("type") == "turn.completed":
            usage = event.get("usage") or {}
        item = event.get("item") or {}
        if item.get("type") == "agent_message" and item.get("text"):
            message_text = item["text"]

    if code != 0 and message_text is None:
        raise ModelError(f"codex exited {code}: {stderr.strip()[-800:]}")

    record = _extract_json(message_text if message_text is not None else stdout)
    result_text = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return ModelRun(
        provider="openai",
        model=model,
        reasoning_effort=reasoning_effort,
        record=record,
        prompt_sha256=_sha256(prompt),
        result_sha256=_sha256(result_text),
        exit_status=code,
        stderr=stderr,
        usage=usage,
        raw_stdout_bytes=len(stdout.encode("utf-8")),
    )


def run_anthropic(
    *,
    model: str,
    reasoning_effort: str,
    system_prompt: str,
    payload_text: str,
    output_schema: dict[str, Any],
    timeout: int = 1800,
) -> ModelRun:
    """Drive the `claude` CLI headlessly (`-p`, JSON output).

    The CLI has no `--output-schema`, so the schema is stated in the prompt and
    the response is parsed leniently. Reasoning effort is not a CLI flag either;
    it is recorded for provenance so the task record still reflects the policy.
    """
    claude = shutil.which("claude")
    if claude is None:
        raise ModelError("claude CLI not found on PATH")

    schema_text = json.dumps(output_schema, indent=2)
    prompt = (
        f"{system_prompt}\n\n"
        "Respond with a single JSON object and no other text. It must validate "
        f"against this JSON Schema:\n{schema_text}\n\n"
        f"=== GENOME PACKET (JSON) ===\n{payload_text}\n"
    )

    argv = [
        claude,
        "-p",
        "--model",
        model,
        "--output-format",
        "json",
        "--permission-mode",
        "plan",
    ]
    code, stdout, stderr = _run(argv, prompt, timeout)

    if code != 0 and _looks_rate_limited(stderr, stdout):
        raise ModelRateLimited(stderr.strip()[-500:] or "claude reported a rate limit")
    if code != 0:
        raise ModelError(f"claude exited {code}: {stderr.strip()[-800:]}")

    usage: dict[str, Any] = {}
    body = stdout
    try:
        envelope = json.loads(stdout)
        if isinstance(envelope, dict):
            usage = envelope.get("usage") or {}
            if isinstance(envelope.get("result"), str):
                body = envelope["result"]
    except json.JSONDecodeError:
        pass

    record = _extract_json(body)
    result_text = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return ModelRun(
        provider="anthropic",
        model=model,
        reasoning_effort=reasoning_effort,
        record=record,
        prompt_sha256=_sha256(prompt),
        result_sha256=_sha256(result_text),
        exit_status=code,
        stderr=stderr,
        usage=usage,
        raw_stdout_bytes=len(stdout.encode("utf-8")),
    )


PROVIDERS = {"openai": run_openai, "anthropic": run_anthropic}


def run_profile(
    *,
    provider: str,
    model: str,
    reasoning_effort: str,
    system_prompt: str,
    payload_text: str,
    output_schema: dict[str, Any],
    timeout: int = 1800,
) -> ModelRun:
    driver = PROVIDERS.get(provider)
    if driver is None:
        raise ModelError(f"no driver for provider {provider!r}")
    return driver(
        model=model,
        reasoning_effort=reasoning_effort,
        system_prompt=system_prompt,
        payload_text=payload_text,
        output_schema=output_schema,
        timeout=timeout,
    )
