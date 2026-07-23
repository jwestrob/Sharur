"""Authentication bridge from the query data plane to Sharur Ops."""

from __future__ import annotations

import hashlib
import secrets
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any

import requests


ALLOWED_ROLES = frozenset({"reader", "worker", "coordinator", "operator"})


class InvalidQueryCredentialError(RuntimeError):
    """Raised when a bearer credential cannot identify an allowed principal."""


class AuthServiceUnavailableError(RuntimeError):
    """Raised when configured Ops token introspection is unavailable."""


@dataclass(frozen=True)
class QueryPrincipal:
    agent_id: str
    role: str
    bootstrap: bool = False


@dataclass(frozen=True)
class _CachedPrincipal:
    principal: QueryPrincipal
    expires_at: float


def is_loopback(host: str | None) -> bool:
    if host is None:
        return False
    return host.strip().lower() in {
        "127.0.0.1",
        "::1",
        "localhost",
        "testclient",
    }


class OpsTokenIntrospector:
    """Resolve per-agent bearer tokens through the sole Ops HTTP owner."""

    def __init__(
        self,
        ops_url: str,
        *,
        timeout_seconds: float = 3.0,
        cache_ttl_seconds: float = 30.0,
        max_cache_entries: int = 4_096,
        session: requests.Session | None = None,
    ):
        if (
            timeout_seconds <= 0
            or cache_ttl_seconds < 0
            or max_cache_entries < 1
        ):
            raise ValueError(
                "Auth timeout/cache size must be positive and cache TTL non-negative"
            )
        self.url = f"{ops_url.rstrip('/')}/auth/whoami"
        self.timeout_seconds = timeout_seconds
        self.cache_ttl_seconds = cache_ttl_seconds
        self.max_cache_entries = max_cache_entries
        self._owns_session = session is None
        self._session = session or requests.Session()
        self._cache: OrderedDict[str, _CachedPrincipal] = OrderedDict()
        self._lock = threading.Lock()

    def close(self) -> None:
        if self._owns_session:
            self._session.close()

    @staticmethod
    def _cache_key(token: str) -> str:
        return hashlib.sha256(token.encode("utf-8")).hexdigest()

    def introspect(self, token: str) -> QueryPrincipal | None:
        key = self._cache_key(token)
        now = time.monotonic()
        with self._lock:
            cached = self._cache.get(key)
            if cached is not None and cached.expires_at >= now:
                self._cache.move_to_end(key)
                return cached.principal
            self._cache.pop(key, None)

        try:
            response = self._session.get(
                self.url,
                headers={"Authorization": f"Bearer {token}"},
                timeout=self.timeout_seconds,
            )
        except requests.RequestException as exc:
            raise AuthServiceUnavailableError(
                f"Sharur Ops token introspection failed: {type(exc).__name__}"
            ) from exc
        if response.status_code in {401, 403}:
            return None
        try:
            response.raise_for_status()
            payload: Any = response.json()
        except (requests.RequestException, ValueError) as exc:
            raise AuthServiceUnavailableError(
                f"Sharur Ops token introspection returned {response.status_code}"
            ) from exc
        if not isinstance(payload, dict):
            raise AuthServiceUnavailableError(
                "Sharur Ops returned a malformed principal"
            )
        agent_id = payload.get("agent_id")
        role = payload.get("role")
        if not isinstance(agent_id, str) or role not in ALLOWED_ROLES:
            raise AuthServiceUnavailableError(
                "Sharur Ops returned an unsupported principal"
            )
        principal = QueryPrincipal(
            agent_id=agent_id,
            role=str(role),
            bootstrap=bool(payload.get("bootstrap", False)),
        )
        with self._lock:
            self._cache[key] = _CachedPrincipal(
                principal=principal,
                expires_at=now + self.cache_ttl_seconds,
            )
            self._cache.move_to_end(key)
            while len(self._cache) > self.max_cache_entries:
                self._cache.popitem(last=False)
        return principal


class QueryAuthenticator:
    """Authenticate a query request with a service token or Ops agent token."""

    def __init__(
        self,
        *,
        api_token: str | None = None,
        ops_url: str | None = None,
        introspector: OpsTokenIntrospector | None = None,
    ):
        self.api_token = api_token.strip() if api_token else None
        self.introspector = introspector or (
            OpsTokenIntrospector(ops_url) if ops_url else None
        )
        self._owns_introspector = introspector is None

    @property
    def auth_required(self) -> bool:
        return self.api_token is not None or self.introspector is not None

    def close(self) -> None:
        if self._owns_introspector and self.introspector is not None:
            self.introspector.close()

    def authenticate(
        self,
        bearer_token: str | None,
        *,
        client_host: str | None,
    ) -> QueryPrincipal:
        if (
            self.api_token is not None
            and bearer_token is not None
            and secrets.compare_digest(self.api_token, bearer_token)
        ):
            return QueryPrincipal("query-bootstrap-operator", "operator", bootstrap=True)
        if bearer_token is not None and self.introspector is not None:
            principal = self.introspector.introspect(bearer_token)
            if principal is not None:
                return principal
        if self.auth_required:
            raise InvalidQueryCredentialError(
                "Missing or invalid Sharur Query bearer token"
            )
        if is_loopback(client_host):
            return QueryPrincipal("local-query-operator", "operator", bootstrap=True)
        raise InvalidQueryCredentialError(
            "Unauthenticated Sharur Query access is restricted to loopback"
        )


__all__ = [
    "ALLOWED_ROLES",
    "AuthServiceUnavailableError",
    "InvalidQueryCredentialError",
    "OpsTokenIntrospector",
    "QueryAuthenticator",
    "QueryPrincipal",
    "is_loopback",
]
