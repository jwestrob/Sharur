"""Bounded HTTP client for the Sharur read-only query service."""

from __future__ import annotations

import uuid
from typing import Any
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


DEFAULT_TIMEOUT = (3.05, 330.0)


class SharurQuery:
    """Typed client for the query data plane; arbitrary SQL is intentionally absent."""

    def __init__(
        self,
        base_url: str = "http://localhost:8812",
        *,
        api_token: str | None = None,
        timeout: float | tuple[float, float] = DEFAULT_TIMEOUT,
        retries: int = 2,
        session: requests.Session | None = None,
        verify_connection: bool = True,
    ):
        if retries < 0:
            raise ValueError("retries must be non-negative")
        if isinstance(timeout, tuple):
            if len(timeout) != 2 or any(value <= 0 for value in timeout):
                raise ValueError(
                    "timeout tuple must contain positive connect/read values"
                )
        elif timeout <= 0:
            raise ValueError("timeout must be positive")
        self.base = base_url.rstrip("/")
        self.timeout = timeout
        self._owns_session = session is None
        self._session = session or requests.Session()
        if api_token:
            self._session.headers.update({"Authorization": f"Bearer {api_token}"})
        if self._owns_session:
            retry_policy = Retry(
                total=retries,
                connect=retries,
                read=retries,
                status=retries,
                allowed_methods=frozenset({"GET", "HEAD", "OPTIONS"}),
                status_forcelist=(429, 502, 503, 504),
                backoff_factor=0.25,
                respect_retry_after_header=True,
            )
            adapter = HTTPAdapter(max_retries=retry_policy)
            self._session.mount("http://", adapter)
            self._session.mount("https://", adapter)
        if verify_connection:
            self.health()

    def close(self) -> None:
        if self._owns_session:
            self._session.close()

    def __enter__(self) -> SharurQuery:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def _request(self, method: str, path: str, **kwargs: Any) -> requests.Response:
        headers = dict(kwargs.pop("headers", {}))
        headers.setdefault("X-Sharur-Query-ID", str(uuid.uuid4()))
        response = self._session.request(
            method,
            f"{self.base}{path}",
            headers=headers,
            timeout=kwargs.pop("timeout", self.timeout),
            **kwargs,
        )
        response.raise_for_status()
        return response

    def health(self) -> dict[str, Any]:
        return self._request("GET", "/health").json()

    def overview(self) -> dict[str, Any]:
        return self._request("GET", "/v1/overview").json()

    def list_genomes(self, **filters: Any) -> dict[str, Any]:
        return self._request("POST", "/v1/genomes/search", json=filters).json()

    def list_proteins(self, **filters: Any) -> dict[str, Any]:
        return self._request("POST", "/v1/proteins/search", json=filters).json()

    def get_genome(self, genome_id: str, *, verbosity: int = 1) -> dict[str, Any]:
        return self._request(
            "GET",
            f"/v1/genomes/{quote(genome_id, safe='')}",
            params={"verbosity": verbosity},
        ).json()

    def get_protein(self, protein_id: str, *, verbosity: int = 1) -> dict[str, Any]:
        return self._request(
            "GET",
            f"/v1/proteins/{quote(protein_id, safe='')}",
            params={"verbosity": verbosity},
        ).json()

    def neighborhood(
        self,
        protein_id: str,
        *,
        window: int = 10,
        verbosity: int = 1,
        all_annotations: bool = False,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/v1/neighborhood",
            json={
                "protein_id": protein_id,
                "window": window,
                "verbosity": verbosity,
                "all_annotations": all_annotations,
            },
        ).json()

    def search_predicates(
        self,
        *,
        has: list[str] | None = None,
        lacks: list[str] | None = None,
        limit: int = 50,
        offset: int = 0,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/v1/predicates/search",
            json={
                "has": has or [],
                "lacks": lacks or [],
                "limit": limit,
                "offset": offset,
            },
        ).json()

    def search_atoms(self, **filters: Any) -> dict[str, Any]:
        return self._request("POST", "/v1/atoms/search", json=filters).json()

    def search_by_atoms(
        self,
        *,
        has: list[str] | None = None,
        lacks: list[str] | None = None,
        limit: int = 50,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/v1/atoms/proteins",
            json={"has": has or [], "lacks": lacks or [], "limit": limit},
        ).json()

    def cancel(self, query_id: str) -> dict[str, Any]:
        return self._request(
            "POST",
            f"/v1/queries/{quote(query_id, safe='')}/cancel",
        ).json()


__all__ = ["DEFAULT_TIMEOUT", "SharurQuery"]
