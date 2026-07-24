"""Bounded read-only HTTP data plane for coordinated Sharur agents."""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import json
import os
import re
import tempfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Any, Literal

from fastapi import (
    APIRouter,
    Depends,
    FastAPI,
    HTTPException,
    Query,
    Request,
)
from fastapi import (
    Path as FastAPIPath,
)
from fastapi.responses import JSONResponse, PlainTextResponse, Response
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel, ConfigDict, Field, field_validator

from sharur import __version__
from sharur.dataset_seal import DEFAULT_SEAL_NAME, DatasetSealError, verify_dataset_seal
from sharur.operators.base import SharurResult
from sharur.operators.contigs import (
    DEFAULT_GENOME_PACKET_BYTES,
    DEFAULT_GENOME_PACKET_CONTIGS,
    DEFAULT_GENOME_PACKET_PROTEINS,
    MAX_GENOME_PACKET_BYTES,
    MAX_GENOME_PACKET_CONTIGS,
    MAX_GENOME_PACKET_PROTEINS,
    get_contig,
    get_contig_packet,
    get_genome_packet,
    list_contigs,
)
from sharur.operators.introspection import overview
from sharur.operators.navigation import (
    get_genome,
    get_neighborhood,
    get_protein,
    list_genomes,
    list_proteins,
)
from sharur.operators.predicates_v2 import search_atoms, search_by_atoms
from sharur.operators.search import search_by_predicates
from sharur.query.auth import (
    AuthServiceUnavailableError,
    InvalidQueryCredentialError,
    QueryAuthenticator,
    QueryPrincipal,
    is_loopback,
)
from sharur.query.runtime import (
    AdmissionRejectedError,
    AdmissionTimeoutError,
    DuplicateQueryIDError,
    QueryCancelledError,
    QueryTimeoutError,
    ReadOnlyDuckDBRuntime,
    WeightedAdmissionController,
)
from sharur.query.staging import DEFAULT_RESERVE_BYTES, StagedDatabase, stage_database


if TYPE_CHECKING:
    from collections.abc import Callable


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8812
DEFAULT_THREADS = 8
DEFAULT_MEMORY_LIMIT = "16GB"
DEFAULT_MAX_TEMP_SIZE = "128GB"
DEFAULT_CAPACITY_UNITS = 12
DEFAULT_HEAVY_WEIGHT = 3
DEFAULT_MAX_QUEUE = 128
DEFAULT_QUEUE_TIMEOUT_SECONDS = 30.0
DEFAULT_LIGHT_TIMEOUT_SECONDS = 30.0
DEFAULT_HEAVY_TIMEOUT_SECONDS = 300.0
DEFAULT_MAX_REQUEST_BYTES = 1024 * 1024
DEFAULT_MAX_RESULT_BYTES = 2 * 1024 * 1024
MAX_ROWS = 500
MAX_TERMS = 64

_QUERY_ID = re.compile(r"[A-Za-z0-9_.:-]{1,128}")


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class GenomeSearch(StrictModel):
    taxonomy_filter: str | None = Field(default=None, max_length=1_024)
    min_completeness: float | None = Field(default=None, ge=0, le=100)
    max_contamination: float | None = Field(default=None, ge=0, le=100)
    limit: int = Field(default=20, ge=1, le=MAX_ROWS)
    offset: int = Field(default=0, ge=0, le=2**63 - 1)


class ProteinSearch(StrictModel):
    genome_id: str | None = Field(default=None, max_length=1_024)
    contig_id: str | None = Field(default=None, max_length=2_048)
    min_length: int | None = Field(default=None, ge=0)
    max_length: int | None = Field(default=None, ge=0)
    has_annotation: bool | None = None
    limit: int = Field(default=50, ge=1, le=MAX_ROWS)
    offset: int = Field(default=0, ge=0, le=2**63 - 1)


class NeighborhoodQuery(StrictModel):
    protein_id: str = Field(min_length=1, max_length=2_048)
    window: int = Field(default=10, ge=0, le=100)
    verbosity: int = Field(default=1, ge=0, le=1)
    all_annotations: bool = False


class ContigPacketQuery(StrictModel):
    genome_id: str = Field(min_length=1, max_length=1_024)
    contig_id: str = Field(min_length=1, max_length=2_048)
    cursor: str | None = Field(default=None, max_length=4_096)
    limit: int = Field(default=100, ge=1, le=250)
    all_annotations: bool = True


class GenomePacketQuery(StrictModel):
    genome_id: str = Field(min_length=1, max_length=1_024)
    cursor: str | None = Field(default=None, max_length=4_096)
    max_contigs: int = Field(
        default=DEFAULT_GENOME_PACKET_CONTIGS,
        ge=1,
        le=MAX_GENOME_PACKET_CONTIGS,
    )
    max_proteins: int = Field(
        default=DEFAULT_GENOME_PACKET_PROTEINS,
        ge=1,
        le=MAX_GENOME_PACKET_PROTEINS,
    )
    max_model_payload_bytes: int = Field(
        default=DEFAULT_GENOME_PACKET_BYTES,
        ge=1_024,
        le=MAX_GENOME_PACKET_BYTES,
    )
    all_annotations: bool = True


class PredicateSearch(StrictModel):
    has: list[str] = Field(default_factory=list, max_length=MAX_TERMS)
    lacks: list[str] = Field(default_factory=list, max_length=MAX_TERMS)
    limit: int = Field(default=50, ge=1, le=MAX_ROWS)
    offset: int = Field(default=0, ge=0, le=2**63 - 1)

    @field_validator("has", "lacks")
    @classmethod
    def bounded_terms(cls, terms: list[str]) -> list[str]:
        if any(not term or len(term) > 256 for term in terms):
            raise ValueError("Predicate terms must contain 1-256 characters")
        return terms


class AtomSearch(StrictModel):
    atom_id: str | list[str] | None = None
    facet: str | None = Field(default=None, max_length=256)
    relation: str | None = Field(default=None, max_length=256)
    source_db: str | None = Field(default=None, max_length=256)
    source_accession: str | None = Field(default=None, max_length=512)
    limit: int = Field(default=50, ge=1, le=MAX_ROWS)

    @field_validator("atom_id")
    @classmethod
    def bounded_atoms(cls, value: str | list[str] | None):
        terms = [value] if isinstance(value, str) else (value or [])
        if len(terms) > MAX_TERMS or any(not term or len(term) > 256 for term in terms):
            raise ValueError("Atom IDs require 1-256 characters and at most 64 values")
        return value


class AtomProteinSearch(StrictModel):
    has: list[str] = Field(default_factory=list, max_length=MAX_TERMS)
    lacks: list[str] = Field(default_factory=list, max_length=MAX_TERMS)
    limit: int = Field(default=50, ge=1, le=MAX_ROWS)

    @field_validator("has", "lacks")
    @classmethod
    def bounded_terms(cls, terms: list[str]) -> list[str]:
        if any(not term or len(term) > 256 for term in terms):
            raise ValueError("Atom terms must contain 1-256 characters")
        return terms


class RequestSizeLimitMiddleware:
    """Apply the same bound to declared and streamed request bodies."""

    def __init__(self, app: Any, *, max_bytes: int):
        self.app = app
        self.max_bytes = max_bytes

    async def __call__(self, scope: dict, receive: Any, send: Any) -> None:
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        headers = dict(scope.get("headers") or [])
        raw_length = headers.get(b"content-length")
        if raw_length is not None:
            with contextlib.suppress(ValueError):
                if int(raw_length) > self.max_bytes:
                    await JSONResponse(
                        {"detail": f"Request body exceeds {self.max_bytes} bytes"},
                        status_code=413,
                    )(scope, receive, send)
                    return
        messages: list[dict] = []
        received = 0
        while True:
            message = await receive()
            messages.append(message)
            if message["type"] != "http.request":
                break
            received += len(message.get("body", b""))
            if received > self.max_bytes:
                await JSONResponse(
                    {"detail": f"Request body exceeds {self.max_bytes} bytes"},
                    status_code=413,
                )(scope, receive, send)
                return
            if not message.get("more_body", False):
                break

        async def replay_receive():
            if messages:
                return messages.pop(0)
            return await receive()

        await self.app(scope, replay_receive, send)


def _operator_payload(value: Any, *, operator: str) -> dict[str, Any]:
    if isinstance(value, SharurResult):
        return json.loads(value.to_json())
    rows = len(value) if isinstance(value, list) else int(value is not None)
    raw_bytes = len(
        json.dumps(value, default=str, separators=(",", ":")).encode("utf-8")
    )
    return {
        "status": "empty" if rows == 0 else "ok",
        "warnings": [],
        "meta": {
            "rows": rows,
            "total_rows": None,
            "bytes": raw_bytes,
            "time_ms": 0,
            "truncated": False,
            "index_used": None,
        },
        "ref": None,
        "trace": {
            "operator": operator,
            "params": {},
            "dataset_version": "query-service",
            "schema_version": "live",
            "trace_hash": None,
        },
        "raw": value,
    }


def _prometheus(
    runtime: ReadOnlyDuckDBRuntime,
    admission: WeightedAdmissionController,
) -> str:
    snapshot = runtime.snapshot()
    admitted = admission.snapshot()
    lines = [
        f"sharur_query_active_queries {snapshot['active_queries']}",
        f"sharur_query_cursors {snapshot['cursor_count']}",
        f"sharur_query_admission_active_units {admitted['active_units']}",
        f"sharur_query_admission_queued {admitted['queued_queries']}",
        f"sharur_query_admission_admitted_total {admitted['admitted_total']}",
        f"sharur_query_admission_rejected_total {admitted['rejected_total']}",
        (
            "sharur_query_admission_timeouts_total "
            f"{admitted['queue_timeouts_total']}"
        ),
        (
            "sharur_query_admission_wait_seconds_total "
            f"{admitted['wait_seconds_total']}"
        ),
    ]
    for operator, metrics in snapshot["operators"].items():
        label = operator.replace("\\", "\\\\").replace('"', '\\"')
        lines.extend(
            [
                f'sharur_query_requests_total{{operator="{label}"}} {metrics["queries"]}',
                f'sharur_query_errors_total{{operator="{label}"}} {metrics["errors"]}',
                f'sharur_query_timeouts_total{{operator="{label}"}} {metrics["timeouts"]}',
                (
                    f'sharur_query_cancellations_total{{operator="{label}"}} '
                    f'{metrics["cancellations"]}'
                ),
                (
                    f'sharur_query_duration_seconds_sum{{operator="{label}"}} '
                    f'{metrics["duration_seconds"]}'
                ),
            ]
        )
    return "\n".join(lines) + "\n"


def create_app(
    *,
    db_path: str | Path,
    api_token: str | None = None,
    ops_url: str | None = None,
    authenticator: QueryAuthenticator | None = None,
    staged_database: StagedDatabase | None = None,
    allow_unsealed: bool = False,
    threads: int = DEFAULT_THREADS,
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
    temp_directory: str | Path | None = None,
    max_temp_directory_size: str = DEFAULT_MAX_TEMP_SIZE,
    capacity_units: int = DEFAULT_CAPACITY_UNITS,
    heavy_weight: int = DEFAULT_HEAVY_WEIGHT,
    max_queue: int = DEFAULT_MAX_QUEUE,
    queue_timeout_seconds: float = DEFAULT_QUEUE_TIMEOUT_SECONDS,
    light_timeout_seconds: float = DEFAULT_LIGHT_TIMEOUT_SECONDS,
    heavy_timeout_seconds: float = DEFAULT_HEAVY_TIMEOUT_SECONDS,
    max_request_bytes: int = DEFAULT_MAX_REQUEST_BYTES,
    max_result_bytes: int = DEFAULT_MAX_RESULT_BYTES,
) -> FastAPI:
    """Build a lazy, single-owner query application."""

    resolved_db = Path(db_path).expanduser().resolve()
    if staged_database is None and not allow_unsealed:
        raise ValueError(
            "create_app requires a verified StagedDatabase; "
            "allow_unsealed=True is reserved for tests and controlled development"
        )
    if (
        staged_database is not None
        and staged_database.path.expanduser().resolve() != resolved_db
    ):
        raise ValueError("staged_database.path must match db_path")
    if heavy_weight < 1 or heavy_weight > capacity_units:
        raise ValueError("heavy_weight must fit within capacity_units")
    for value, name in (
        (queue_timeout_seconds, "queue_timeout_seconds"),
        (light_timeout_seconds, "light_timeout_seconds"),
        (heavy_timeout_seconds, "heavy_timeout_seconds"),
    ):
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    if max_request_bytes < 1 or max_result_bytes < 1:
        raise ValueError("HTTP request and result bounds must be positive")
    configured_temp = (
        Path(temp_directory).expanduser().resolve()
        if temp_directory is not None
        else Path(tempfile.gettempdir())
        / "sharur-query-spill"
        / uuid.uuid5(uuid.NAMESPACE_URL, str(resolved_db)).hex[:16]
    )
    runtime = ReadOnlyDuckDBRuntime(
        resolved_db,
        threads=threads,
        memory_limit=memory_limit,
        temp_directory=configured_temp,
        max_temp_directory_size=max_temp_directory_size,
        dataset_id=staged_database.dataset_id if staged_database else None,
    )
    admission = WeightedAdmissionController(
        capacity=capacity_units,
        max_queue=max_queue,
    )
    query_auth = authenticator or QueryAuthenticator(
        api_token=api_token,
        ops_url=ops_url,
    )
    owns_authenticator = authenticator is None

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        executor = ThreadPoolExecutor(
            max_workers=capacity_units,
            thread_name_prefix="sharur-query",
        )
        application.state.query_executor = executor
        loop = asyncio.get_running_loop()
        try:
            await loop.run_in_executor(executor, runtime.open)
            yield
        finally:
            await loop.run_in_executor(executor, runtime.close)
            executor.shutdown(wait=True, cancel_futures=True)
            application.state.query_executor = None
            if owns_authenticator:
                query_auth.close()

    application = FastAPI(
        title="Sharur Query",
        version=__version__,
        lifespan=lifespan,
    )
    application.state.query_runtime = runtime
    application.state.query_admission = admission
    application.state.query_authenticator = query_auth
    application.state.query_executor = None
    application.state.staged_database = staged_database
    application.add_middleware(
        RequestSizeLimitMiddleware,
        max_bytes=max_request_bytes,
    )

    bearer = HTTPBearer(auto_error=False)

    def authorize(
        request: Request,
        credentials: HTTPAuthorizationCredentials | None = Depends(bearer),
    ) -> QueryPrincipal:
        token = (
            credentials.credentials
            if credentials is not None and credentials.scheme.lower() == "bearer"
            else None
        )
        try:
            principal = query_auth.authenticate(
                token,
                client_host=request.client.host if request.client else None,
            )
        except AuthServiceUnavailableError as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc
        except InvalidQueryCredentialError as exc:
            status = 401 if query_auth.auth_required else 403
            raise HTTPException(
                status_code=status,
                detail=str(exc),
                headers={"WWW-Authenticate": "Bearer"} if status == 401 else None,
            ) from exc
        request.state.query_principal = principal
        return principal

    router = APIRouter(dependencies=[Depends(authorize)])

    def principal(request: Request) -> QueryPrincipal:
        return request.state.query_principal

    def query_id(request: Request) -> str:
        supplied = request.headers.get("X-Sharur-Query-ID")
        candidate = supplied or str(uuid.uuid4())
        if _QUERY_ID.fullmatch(candidate) is None:
            raise HTTPException(
                status_code=422,
                detail="X-Sharur-Query-ID must contain 1-128 safe identifier characters",
            )
        return candidate

    async def run_operator(
        request: Request,
        *,
        operator: str,
        query_class: Literal["light", "heavy"],
        operation: Callable[[Any], Any],
    ) -> Response:
        identifier = query_id(request)
        actor = principal(request)
        weight = 1 if query_class == "light" else heavy_weight
        timeout = (
            light_timeout_seconds
            if query_class == "light"
            else heavy_timeout_seconds
        )
        try:
            wait_seconds = await admission.acquire(
                weight,
                timeout_seconds=queue_timeout_seconds,
            )
        except AdmissionRejectedError as exc:
            raise HTTPException(
                status_code=429,
                detail=str(exc),
                headers={"Retry-After": "1"},
            ) from exc
        except AdmissionTimeoutError as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc

        executor = request.app.state.query_executor
        if executor is None:
            raise HTTPException(status_code=503, detail="Query executor is unavailable")
        task = asyncio.get_running_loop().run_in_executor(
            executor,
            partial(
                runtime.execute,
                query_id=identifier,
                principal=actor,
                operator=operator,
                timeout_seconds=timeout,
                operation=operation,
            ),
        )
        try:
            while True:
                done, _pending = await asyncio.wait({task}, timeout=0.25)
                if done:
                    execution = task.result()
                    break
                if await request.is_disconnected():
                    runtime.cancel(identifier, principal=actor, reason="disconnect")
                    with contextlib.suppress(Exception):
                        await task
                    raise HTTPException(
                        status_code=499,
                        detail=f"Client disconnected; query {identifier} was cancelled",
                    )
        except asyncio.CancelledError:
            runtime.cancel(identifier, principal=actor, reason="disconnect")
            with contextlib.suppress(Exception):
                await asyncio.shield(task)
            raise
        except DuplicateQueryIDError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        except QueryTimeoutError as exc:
            raise HTTPException(status_code=504, detail=str(exc)) from exc
        except QueryCancelledError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(
                status_code=422,
                detail=f"{type(exc).__name__}: {exc}",
            ) from exc
        finally:
            await admission.release(weight)

        payload = _operator_payload(execution.value, operator=operator)
        payload["service"] = {
            "query_id": identifier,
            "operator": operator,
            "class": query_class,
            "dataset_id": runtime.dataset_id,
            "queue_wait_ms": round(wait_seconds * 1000, 3),
            "execution_ms": round(execution.execution_seconds * 1000, 3),
        }
        encoded = json.dumps(
            payload,
            default=str,
            separators=(",", ":"),
        ).encode("utf-8")
        if len(encoded) > max_result_bytes:
            raise HTTPException(
                status_code=413,
                detail=(
                    f"Serialized result is {len(encoded):,} bytes; "
                    f"service limit is {max_result_bytes:,} bytes"
                ),
            )
        return Response(
            content=encoded,
            media_type="application/json",
            headers={"X-Sharur-Query-ID": identifier},
        )

    @router.get("/health")
    async def health(request: Request):
        snapshot = runtime.snapshot()
        staged = request.app.state.staged_database
        return {
            "status": "ok",
            "database": str(resolved_db),
            "dataset_id": staged.dataset_id if staged else None,
            "seal_strength": staged.seal_strength if staged else None,
            "read_only": True,
            "duckdb_owner": "single-http-server",
            "workers": 1,
            "auth_required": query_auth.auth_required,
            "resource_budget": runtime.resource_budget,
            "settings": snapshot["settings"],
            "admission": admission.snapshot(),
            "ts": time.time(),
        }

    @router.get("/auth/whoami")
    async def whoami(request: Request):
        actor = principal(request)
        return {
            "agent_id": actor.agent_id,
            "role": actor.role,
            "bootstrap": actor.bootstrap,
        }

    @router.get("/metrics")
    async def metrics():
        return PlainTextResponse(
            _prometheus(runtime, admission),
            media_type="text/plain; version=0.0.4",
        )

    @router.get("/v1/overview")
    async def dataset_overview(request: Request):
        return await run_operator(
            request,
            operator="overview",
            query_class="heavy",
            operation=overview,
        )

    @router.post("/v1/genomes/search")
    async def genomes_search(body: GenomeSearch, request: Request):
        return await run_operator(
            request,
            operator="list_genomes",
            query_class="light",
            operation=lambda store: list_genomes(store, **body.model_dump()),
        )

    @router.post("/v1/proteins/search")
    async def proteins_search(body: ProteinSearch, request: Request):
        if (
            body.min_length is not None
            and body.max_length is not None
            and body.min_length > body.max_length
        ):
            raise HTTPException(
                status_code=422,
                detail="min_length must be less than or equal to max_length",
            )
        return await run_operator(
            request,
            operator="list_proteins",
            query_class="light",
            operation=lambda store: list_proteins(store, **body.model_dump()),
        )

    @router.get("/v1/genomes/{genome_id}")
    async def genome_detail(
        request: Request,
        genome_id: Annotated[str, FastAPIPath(min_length=1, max_length=1_024)],
        verbosity: Annotated[int, Query(ge=0, le=1)] = 1,
    ):
        return await run_operator(
            request,
            operator="get_genome",
            query_class="light",
            operation=lambda store: get_genome(store, genome_id, verbosity),
        )

    @router.get("/v1/genomes/{genome_id}/contigs")
    async def contigs_list(
        request: Request,
        genome_id: Annotated[str, FastAPIPath(min_length=1, max_length=1_024)],
        limit: Annotated[int, Query(ge=1, le=1_000)] = 100,
        cursor: Annotated[str | None, Query(max_length=4_096)] = None,
    ):
        return await run_operator(
            request,
            operator="list_contigs",
            query_class="light",
            operation=lambda store: list_contigs(
                store,
                genome_id,
                limit=limit,
                cursor=cursor,
            ),
        )

    @router.get("/v1/genomes/{genome_id}/contigs/{contig_id}")
    async def contig_detail(
        request: Request,
        genome_id: Annotated[str, FastAPIPath(min_length=1, max_length=1_024)],
        contig_id: Annotated[str, FastAPIPath(min_length=1, max_length=2_048)],
        verbosity: Annotated[int, Query(ge=0, le=1)] = 1,
    ):
        return await run_operator(
            request,
            operator="get_contig",
            query_class="light",
            operation=lambda store: get_contig(
                store,
                genome_id,
                contig_id,
                verbosity=verbosity,
            ),
        )

    @router.post("/v1/contigs/packet")
    async def contig_packet(body: ContigPacketQuery, request: Request):
        return await run_operator(
            request,
            operator="get_contig_packet",
            query_class="light",
            operation=lambda store: get_contig_packet(
                store,
                **body.model_dump(),
            ),
        )

    @router.post("/v1/genomes/packet")
    async def genome_packet(body: GenomePacketQuery, request: Request):
        return await run_operator(
            request,
            operator="get_genome_packet",
            query_class="light",
            operation=lambda store: get_genome_packet(
                store,
                **body.model_dump(),
            ),
        )

    @router.get("/v1/proteins/{protein_id}")
    async def protein_detail(
        request: Request,
        protein_id: Annotated[str, FastAPIPath(min_length=1, max_length=2_048)],
        verbosity: Annotated[int, Query(ge=0, le=1)] = 1,
    ):
        return await run_operator(
            request,
            operator="get_protein",
            query_class="light",
            operation=lambda store: get_protein(store, protein_id, verbosity),
        )

    @router.post("/v1/neighborhood")
    async def neighborhood(body: NeighborhoodQuery, request: Request):
        return await run_operator(
            request,
            operator="get_neighborhood",
            query_class="light",
            operation=lambda store: get_neighborhood(
                store,
                entity_id=body.protein_id,
                window=body.window,
                verbosity=body.verbosity,
                all_annotations=body.all_annotations,
            ),
        )

    @router.post("/v1/predicates/search")
    async def predicates_search(body: PredicateSearch, request: Request):
        return await run_operator(
            request,
            operator="search_by_predicates",
            query_class="heavy",
            operation=lambda store: search_by_predicates(
                store,
                has=body.has,
                lacks=body.lacks,
                limit=body.limit,
                offset=body.offset,
            ),
        )

    @router.post("/v1/atoms/search")
    async def atoms_search(body: AtomSearch, request: Request):
        return await run_operator(
            request,
            operator="search_atoms",
            query_class="heavy",
            operation=lambda store: search_atoms(store, **body.model_dump()),
        )

    @router.post("/v1/atoms/proteins")
    async def atom_protein_search(body: AtomProteinSearch, request: Request):
        return await run_operator(
            request,
            operator="search_by_atoms",
            query_class="heavy",
            operation=lambda store: search_by_atoms(
                store,
                has=body.has,
                lacks=body.lacks,
                limit=body.limit,
            ),
        )

    @router.post("/v1/queries/{query_id_value}/cancel", status_code=202)
    async def cancel_query(
        request: Request,
        query_id_value: Annotated[
            str,
            FastAPIPath(min_length=1, max_length=128),
        ],
    ):
        if _QUERY_ID.fullmatch(query_id_value) is None:
            raise HTTPException(status_code=422, detail="Invalid query identifier")
        try:
            cancelled = runtime.cancel(
                query_id_value,
                principal=principal(request),
            )
        except PermissionError as exc:
            raise HTTPException(status_code=403, detail=str(exc)) from exc
        if not cancelled:
            raise HTTPException(status_code=404, detail="Active query was not found")
        return {"query_id": query_id_value, "status": "cancelling"}

    application.include_router(router)
    return application


def _direct_sealed_database(
    database: Path,
    seal_path: Path,
) -> tuple[str, str]:
    verification = verify_dataset_seal(seal_path, db_path=database)
    if not verification.valid:
        changed = ", ".join(verification.changed_sections) or "dataset identity"
        raise DatasetSealError(f"Dataset seal verification failed; changed: {changed}")
    payload = json.loads(seal_path.read_text())
    return str(payload["dataset_id"]), str(payload.get("seal_strength", "unknown"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument(
        "--stage-dir",
        type=Path,
        default=(
            Path(os.environ["SHARUR_QUERY_STAGE_DIR"])
            if os.environ.get("SHARUR_QUERY_STAGE_DIR")
            else None
        ),
        help="Campaign-local storage for the atomic immutable replica",
    )
    parser.add_argument(
        "--direct",
        action="store_true",
        help="Open the verified immutable source database in place",
    )
    parser.add_argument("--seal", type=Path)
    parser.add_argument(
        "--reserve-gb",
        type=float,
        default=DEFAULT_RESERVE_BYTES / 1024**3,
    )
    parser.add_argument("--host", default=os.environ.get("SHARUR_QUERY_HOST", DEFAULT_HOST))
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("SHARUR_QUERY_PORT", DEFAULT_PORT)),
    )
    parser.add_argument(
        "--ops-url",
        default=os.environ.get("SHARUR_OPS_URL"),
        help="Sharur Ops base URL for per-agent token introspection",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=int(os.environ.get("SHARUR_QUERY_THREADS", DEFAULT_THREADS)),
    )
    parser.add_argument(
        "--memory-limit",
        default=os.environ.get("SHARUR_QUERY_MEMORY_LIMIT", DEFAULT_MEMORY_LIMIT),
    )
    parser.add_argument("--temp-dir", type=Path)
    parser.add_argument(
        "--max-temp-size",
        default=os.environ.get("SHARUR_QUERY_MAX_TEMP_SIZE", DEFAULT_MAX_TEMP_SIZE),
    )
    parser.add_argument(
        "--capacity-units",
        type=int,
        default=int(
            os.environ.get("SHARUR_QUERY_CAPACITY_UNITS", DEFAULT_CAPACITY_UNITS)
        ),
    )
    parser.add_argument(
        "--heavy-weight",
        type=int,
        default=int(os.environ.get("SHARUR_QUERY_HEAVY_WEIGHT", DEFAULT_HEAVY_WEIGHT)),
    )
    parser.add_argument(
        "--max-queue",
        type=int,
        default=int(os.environ.get("SHARUR_QUERY_MAX_QUEUE", DEFAULT_MAX_QUEUE)),
    )
    parser.add_argument(
        "--light-timeout",
        type=float,
        default=float(
            os.environ.get(
                "SHARUR_QUERY_LIGHT_TIMEOUT_SECONDS",
                DEFAULT_LIGHT_TIMEOUT_SECONDS,
            )
        ),
    )
    parser.add_argument(
        "--heavy-timeout",
        type=float,
        default=float(
            os.environ.get(
                "SHARUR_QUERY_HEAVY_TIMEOUT_SECONDS",
                DEFAULT_HEAVY_TIMEOUT_SECONDS,
            )
        ),
    )
    args = parser.parse_args()

    database = args.db.expanduser().resolve()
    seal_path = (
        args.seal.expanduser().resolve()
        if args.seal is not None
        else database.parent / DEFAULT_SEAL_NAME
    )
    if args.direct and args.stage_dir is not None:
        parser.error("--direct and --stage-dir are mutually exclusive")
    if not args.direct and args.stage_dir is None:
        parser.error(
            "--stage-dir (or SHARUR_QUERY_STAGE_DIR) is required; "
            "use --direct to serve the verified immutable source in place"
        )
    token = os.environ.get("SHARUR_QUERY_TOKEN")
    if not is_loopback(args.host) and not token and not args.ops_url:
        parser.error(
            "remote binding requires SHARUR_QUERY_TOKEN or --ops-url token introspection"
        )

    staged: StagedDatabase | None
    if args.direct:
        dataset_id, seal_strength = _direct_sealed_database(database, seal_path)
        staged = StagedDatabase(
            path=database,
            source_path=database,
            seal_path=seal_path,
            dataset_id=dataset_id,
            seal_strength=seal_strength,
            artifact_digest={},
            reused=True,
            staged_at="direct",
        )
    else:
        staged = stage_database(
            database,
            args.stage_dir,
            seal_path=seal_path,
            reserve_bytes=int(args.reserve_gb * 1024**3),
        )
        database = staged.path

    temp_directory = (
        args.temp_dir.expanduser().resolve()
        if args.temp_dir is not None
        else database.parent / "spill" / staged.dataset_id[:12]
    )

    import uvicorn  # noqa: PLC0415

    server_app = create_app(
        db_path=database,
        api_token=token,
        ops_url=args.ops_url,
        staged_database=staged,
        threads=args.threads,
        memory_limit=args.memory_limit,
        temp_directory=temp_directory,
        max_temp_directory_size=args.max_temp_size,
        capacity_units=args.capacity_units,
        heavy_weight=args.heavy_weight,
        max_queue=args.max_queue,
        light_timeout_seconds=args.light_timeout,
        heavy_timeout_seconds=args.heavy_timeout,
    )
    uvicorn.run(server_app, host=args.host, port=args.port, workers=1)


__all__ = ["create_app", "main"]


if __name__ == "__main__":
    main()
