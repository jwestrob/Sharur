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

import contextlib
import hashlib
import json
import re
import os
import shutil
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


class ModelError(RuntimeError):
    """Model invocation failed in a way that is not a rate limit."""


class ModelRateLimited(RuntimeError):
    """Provider refused the call for quota/rate reasons. Retry after backoff."""


class ModelQuotaExhausted(ModelRateLimited):
    """The subscription window is spent and names a far-future reset.

    Distinct from a throttle because backing off cannot help. The Dormibacteria
    pilot hit "You've hit your usage limit ... try again at Aug 1st" six days
    out; the backoff caps at 1800s, so eight workers would have spun for six
    days, each waking every half hour to be refused again, holding leases and
    burning attempts for nothing.

    A subclass so existing `except ModelRateLimited` handlers still catch it if
    a caller does not care about the distinction.
    """

    def __init__(self, message: str, reset_at: str | None = None) -> None:
        super().__init__(message)
        self.reset_at = reset_at


# "try again at <when>" is the provider telling us the window is spent rather
# than that we are being throttled. Capture it so the worker can say when to
# come back instead of rediscovering the wall every 30 minutes.
_QUOTA_RESET = re.compile(
    r"(?:try again (?:at|after)|resets? (?:at|on))\s+([^\".]{4,60})", re.IGNORECASE
)
_QUOTA_MARKERS = ("usage limit", "purchase more credits", "quota exceeded")


def _quota_message(blob: str) -> str | None:
    """The provider's own sentence, when stderr is empty and it is in stdout."""
    for line in blob.splitlines():
        if any(marker in line.lower() for marker in _QUOTA_MARKERS):
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                return line.strip()[:300]
            if isinstance(event, dict):
                msg = event.get("message")
                if isinstance(msg, str):
                    return msg[:300]
                err = event.get("error")
                if isinstance(err, dict) and isinstance(err.get("message"), str):
                    return err["message"][:300]
    return None


def _quota_exhausted(blob: str) -> tuple[bool, str | None]:
    lowered = blob.lower()
    if not any(marker in lowered for marker in _QUOTA_MARKERS):
        return False, None
    match = _QUOTA_RESET.search(blob)
    return True, (match.group(1).strip() if match else None)


class ModelTransient(RuntimeError):
    """Transient transport failure (DNS, socket, timeout).

    Distinct from :class:`ModelError` because it must NOT consume a task
    attempt. Biotite's login node resolves against external nameservers with
    `options timeout:2 attempts:3`, so lookups fail intermittently under load
    and the CLI gives up immediately. Treating that as a task failure would
    let five DNS blips permanently kill a genome at the default max_attempts.
    """


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
    "quota",
    "usage limit",
    "too many requests",
    "resource_exhausted",
)

# `429` needs context. As a bare substring it matched any occurrence anywhere in
# a verbose CLI stderr -- request ids, token counts, byte offsets, URLs -- and
# rate-limit classification runs BEFORE transient, so one incidental "429" turned
# a DNS blip into a rate-limit backoff. That backoff doubles and persists across
# tasks (60s -> 1800s), so on a host where DNS failures are endemic a worker
# could idle for half an hour over a network hiccup it should have retried in
# seconds. Observed once in the Dormibacteria pilot: a websocket DNS failure
# logged as "rate limited".
_HTTP_429 = re.compile(r"(?:status|code|http|error)\D{0,12}\b429\b|\b429\b\s*(?:too many|rate)")


def _looks_rate_limited(stderr: str, stdout: str) -> bool:
    blob = f"{stderr}\n{stdout}".lower()
    if any(marker in blob for marker in _RATE_LIMIT_MARKERS):
        return True
    return _HTTP_429.search(blob) is not None


_TRANSIENT_MARKERS = (
    "failed to lookup address information",
    "temporary failure in name resolution",
    "failed to connect to websocket",
    "connection reset",
    "connection refused",
    "network is unreachable",
    "no route to host",
    "broken pipe",
    "timed out",
    "timeout",
    "eof while parsing",
    "stream closed",
    "no output; treating as a stalled connection",
)


def _looks_transient(stderr: str, stdout: str) -> bool:
    """Transport-level failure worth an immediate in-process retry."""
    blob = f"{stderr}\n{stdout}".lower()
    return any(marker in blob for marker in _TRANSIENT_MARKERS)


def _classify_failure(code: int, stderr: str, stdout: str) -> None:
    """Raise the right exception class for a non-zero CLI exit.

    Order matters: a rate limit can also mention a timeout, and the caller
    handles the two very differently (release-and-back-off vs retry in place).
    """
    blob = f"{stderr}\n{stdout}"
    exhausted, reset_at = _quota_exhausted(blob)
    if exhausted:
        detail = stderr.strip()[-500:] or _quota_message(blob) or "usage limit reached"
        raise ModelQuotaExhausted(detail, reset_at=reset_at)
    if _looks_rate_limited(stderr, stdout):
        raise ModelRateLimited(stderr.strip()[-500:] or "provider reported a rate limit")
    if _looks_transient(stderr, stdout):
        raise ModelTransient(stderr.strip()[-500:] or "transient transport failure")
    raise ModelError(f"CLI exited {code}: {stderr.strip()[-800:]}")


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
    """Run a model CLI, bounded by SILENCE rather than by total duration.

    A hard wall-clock cap is the wrong instrument here: a hard genome at high
    effort can legitimately think for a long time, and killing it discards
    everything it has done. What actually indicates a dead call is the absence
    of progress, and `codex exec --json` streams JSONL events throughout, so
    output activity is a direct liveness signal.

    `timeout` is therefore interpreted as a STALL timeout — the maximum silence
    tolerated between writes — not a total runtime limit. A call producing
    output runs as long as it likes.
    """
    proc = subprocess.Popen(
        argv,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        env=os.environ.copy(),
    )

    out_chunks: list[str] = []
    err_chunks: list[str] = []
    last_activity = [time.monotonic()]
    lock = threading.Lock()

    def pump(stream, sink: list[str]) -> None:
        try:
            for line in iter(stream.readline, ""):
                with lock:
                    sink.append(line)
                    last_activity[0] = time.monotonic()
        finally:
            with contextlib.suppress(Exception):
                stream.close()

    threads = [
        threading.Thread(target=pump, args=(proc.stdout, out_chunks), daemon=True),
        threading.Thread(target=pump, args=(proc.stderr, err_chunks), daemon=True),
    ]
    for t in threads:
        t.start()

    try:
        proc.stdin.write(stdin_text)
        proc.stdin.close()
    except BrokenPipeError:
        pass

    stalled = False
    while proc.poll() is None:
        with lock:
            silent_for = time.monotonic() - last_activity[0]
        if silent_for > timeout:
            stalled = True
            proc.terminate()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
            break
        time.sleep(2.0)

    for t in threads:
        t.join(timeout=10.0)
    code = proc.poll()
    stdout, stderr = "".join(out_chunks), "".join(err_chunks)
    if stalled:
        stderr += (
            f"\n[sharur] terminated after {timeout}s with no output; "
            "treating as a stalled connection"
        )
        return (code if code is not None else 1), stdout, stderr
    return (code if code is not None else 0), stdout, stderr


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

    if code != 0:
        _classify_failure(code, stderr, stdout)

    # `--json` emits JSONL events; the record is in the final agent_message and
    # usage is on turn.completed. Fall back to raw stdout if events are absent.
    usage: dict[str, Any] = {}
    message_text: str | None = None
    turn_error: str | None = None
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
        if event.get("type") in ("turn.failed", "error"):
            # codex exits 0 on a failed turn, so the exit code alone is not a
            # reliable success signal — the event stream is authoritative.
            detail = event.get("message") or json.dumps(event.get("error") or {})
            turn_error = str(detail)[:1500]
        item = event.get("item") or {}
        if item.get("type") == "agent_message" and item.get("text"):
            message_text = item["text"]

    if turn_error and message_text is None:
        _classify_failure(code or 1, f"{stderr}\n{turn_error}", stdout)
    if code != 0 and message_text is None:
        _classify_failure(code, stderr, stdout)

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


def _parse_anthropic_stream(stdout: str) -> tuple[dict[str, Any], str]:
    """Extract (usage, result text) from a `--output-format stream-json` run.

    The stream is JSONL; the terminal `result` event carries both. Falls back to
    treating the whole of stdout as the body, so a format change degrades to
    lenient parsing rather than losing the answer outright.
    """
    usage: dict[str, Any] = {}
    body: str | None = None
    for line in stdout.splitlines():
        line = line.strip()
        if not line or not line.startswith("{"):
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(event, dict):
            continue
        if event.get("type") == "result":
            if isinstance(event.get("usage"), dict):
                usage = event["usage"]
            if isinstance(event.get("result"), str):
                body = event["result"]
            # An error result must not be mistaken for an empty answer.
            if event.get("is_error") and event.get("subtype") not in (None, "success"):
                raise ModelError(
                    f"claude turn failed: {event.get('subtype')} "
                    f"{str(event.get('result'))[:300]}"
                )
    return usage, body if body is not None else stdout


def run_anthropic(
    *,
    model: str,
    reasoning_effort: str,
    system_prompt: str,
    payload_text: str,
    output_schema: dict[str, Any],
    timeout: int = 1800,
) -> ModelRun:
    """Drive the `claude` CLI headlessly (`-p`, streaming JSON output).

    The CLI has no `--output-schema`, so the schema is stated in the prompt and
    the response is parsed leniently.

    Three things this driver previously got wrong, all of which made the
    provenance record disagree with what actually ran:

    * **Effort was recorded but never applied.** `--effort` does exist and takes
      exactly the policy vocabulary (low/medium/high/xhigh/max). Recording the
      requested effort while not passing it meant the task record asserted a
      setting the run never used -- worse than not recording it at all.
    * **The silence timeout was really a wall-clock timeout.** `--output-format
      json` emits a single blob at the end, so a healthy long call looks
      identical to a hung one and gets killed at the stall deadline. Streaming
      makes silence mean silence.
    * **The packet-only boundary leaked.** Project settings, MCP servers and
      tools were all live, so the model could read the filesystem instead of
      only the frozen packet it was handed. That silently violates the input
      contract every other part of this pipeline enforces.
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
        "--effort",
        reasoning_effort,
        # Stream so the stall detector measures silence rather than duration.
        "--output-format",
        "stream-json",
        "--include-partial-messages",
        "--verbose",
        # Packet-only boundary: no tools, no MCP servers, no user/project/local
        # settings. The model must reason from the packet it was handed and
        # nothing else, which is the same contract the OpenAI path enforces.
        "--tools",
        "",
        "--strict-mcp-config",
        "--setting-sources",
        "",
        "--permission-mode",
        "plan",
    ]
    code, stdout, stderr = _run(argv, prompt, timeout)

    if code != 0:
        _classify_failure(code, stderr, stdout)

    usage, body = _parse_anthropic_stream(stdout)
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
