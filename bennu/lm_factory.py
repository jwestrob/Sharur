"""LM factory aligned with microbial_claude_matter/src/llm/lm_factory.py.

Key goals:
- Use chat endpoint for GPT-5 models.
- Strip legacy params (max_tokens, response_format, drop_params) before hitting OpenAI.
- Default temperature 1.0 (GPT-5 rejects 0.0) and max_completion_tokens 16000.
- Preserve reasoning_effort for gpt-5-* aliases.

We cannot rely on dspy.LM in this environment (not shipped in dspy-ai 2.4.0),
so we subclass dsp.modules.gpt3.GPT3 to get the LM interface that DSPy expects.
"""

from __future__ import annotations

import os
from typing import Any, Optional

import dsp
from dsp.modules.gpt3 import GPT3
import json

# -----------------------
# Alias + effort helpers
# -----------------------


def _resolve_alias(model_id: str) -> Optional[str]:
    s = (model_id or "").strip().lower()
    if not s:
        return None
    name = s.split("/", 1)[1] if "/" in s else s
    alias = {
        "gpt-5": "openai/gpt-5-2025-08-07",
        "gpt-5-high": "openai/gpt-5-2025-08-07",
        "gpt-5-premium": "openai/gpt-5-2025-08-07",
        "gpt-5-minimal": "openai/gpt-5-2025-08-07",
        "gpt-5-mini": "openai/gpt-5-2025-08-07",
        "gpt-5-medium": "openai/gpt-5-2025-08-07",
        "gpt-5-low": "openai/gpt-5-2025-08-07",
        "gpt-4.1-mini": "openai/gpt-4.1-mini",
        "4.1-mini": "openai/gpt-4.1-mini",
    }
    return alias.get(name)


def _normalize_model(model_id: str) -> str:
    mid = (model_id or "").strip()
    if not mid:
        return "gpt-5-mini-2025-08-07"
    resolved = _resolve_alias(mid)
    if resolved:
        mid = resolved
    if "/" not in mid:
        return mid
    if mid.lower().startswith("openai/"):
        return mid.split("/", 1)[1]
    if "/" not in mid:
        return mid
    return mid


def _extract_gpt5_effort(model_id: str) -> Optional[str]:
    if not model_id:
        return None
    alias = model_id.split("/", 1)[1] if "/" in model_id else model_id
    a = alias.strip().lower()
    if a.startswith("gpt-5-"):
        if a.endswith("-minimal"):
            return "minimal"
        if a.endswith("-low"):
            return "low"
        if a.endswith("-medium"):
            return "medium"
        if a.endswith("-high"):
            return "high"
    return None


# -----------------------
# LM implementation
# -----------------------

class BennuChatLM(GPT3):
    """Chat LM that cleans unsupported params for GPT-5 style chat models."""

    def __init__(
        self,
        model: str,
        temperature: float = 1.0,
        max_completion_tokens: int = 16000,
        reasoning_effort: Optional[str] = None,
    ):
        # Force chat model_type to hit chat endpoint in GPT3.basic_request
        super().__init__(model=model, model_type="chat")

        # Reset kwargs to avoid legacy defaults
        self.kwargs = {
            "model": model,
            "temperature": temperature,
            "top_p": 1,
            "frequency_penalty": 0,
            "presence_penalty": 0,
            "n": 1,
            # Keep max_tokens for DSPy internal retry logic; we will strip it before API call
            "max_tokens": max_completion_tokens,
            "max_completion_tokens": max_completion_tokens,
        }
        if reasoning_effort:
            self.kwargs["reasoning_effort"] = reasoning_effort

    def basic_request(self, prompt: str, **kwargs):
        # Drop unsupported/legacy params that may get forwarded
        for k in (
            "max_tokens",
            "max_output_tokens",
            "response_format",
            "drop_params",
            "additional_drop_params",
        ):
            kwargs.pop(k, None)

        payload = {**self.kwargs, **kwargs}
        # Remove legacy max_tokens before hitting chat API
        payload.pop("max_tokens", None)
        # GPT-5 chat only accepts default temperature; force 1.0
        payload["temperature"] = 1.0
        payload["messages"] = [{"role": "user", "content": prompt}]
        # GPT3.chat_request expects 'stringify_request' wrapper for caching
        payload = {"stringify_request": json.dumps(payload)}
        response = dsp.modules.gpt3.chat_request(**payload)
        self.history.append({"prompt": prompt, "response": response, "kwargs": payload})
        return response


# -----------------------
# Factory
# -----------------------

def make_lm(model_id: str, temperature: float = 1.0, max_completion_tokens: int = 16000) -> Any:
    model = _normalize_model(model_id)
    effort = _extract_gpt5_effort(model_id)
    lm = BennuChatLM(
        model=model,
        temperature=temperature,
        max_completion_tokens=max_completion_tokens,
        reasoning_effort=effort,
    )
    return lm
