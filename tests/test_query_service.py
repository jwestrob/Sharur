"""Correctness, concurrency, security, and bounds for Sharur Query."""

from __future__ import annotations

import asyncio
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb
import pytest
import requests
from fastapi.testclient import TestClient

from sharur.query.auth import (
    AuthServiceUnavailableError,
    OpsTokenIntrospector,
    QueryPrincipal,
)
from sharur.query.client import SharurQuery
from sharur.query.runtime import (
    AdmissionRejectedError,
    QueryCancelledError,
    QueryTimeoutError,
    ReadOnlyDuckDBRuntime,
    WeightedAdmissionController,
)
from sharur.query.server import create_app
from sharur.query.staging import StagedDatabase
from sharur.storage.duckdb_store import DuckDBStore


TOKEN = "query-test-secret"


def _query_database(path):
    store = DuckDBStore(path)
    store.execute(
        """
        INSERT INTO bins(
            bin_id, completeness, contamination, taxonomy, n_contigs, total_length
        ) VALUES ('bin1', 95, 1, 'd__Bacteria;p__Testota', 1, 900)
        """
    )
    store.execute(
        """
        INSERT INTO contigs(contig_id, bin_id, length, gc_content)
        VALUES ('contig1', 'bin1', 900, 0.5)
        """
    )
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence, sequence_length
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 1, 'MAAA', 100)
        """
    )
    store.execute(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name, description,
            evalue, score
        ) VALUES (1, 'p1', 'pfam', 'PFTEST', 'Observed_domain',
                  'Fixture observation', 1e-20, 80)
        """
    )
    store.execute(
        """
        INSERT INTO protein_predicates(protein_id, predicates)
        VALUES ('p1', ['giant'])
        """
    )
    store.execute(
        """
        INSERT INTO semantic_atoms(
            protein_id, atom_id, facet, relation, source_accession, source_db,
            evidence_evalue, evidence_score
        ) VALUES ('p1', 'pfam:PFTEST', 'architecture', 'observed',
                  'PFTEST', 'pfam', 1e-20, 80)
        """
    )
    store.close()
    return path


def _app(database, tmp_path, **kwargs):
    return create_app(
        db_path=database,
        temp_directory=tmp_path / "spill",
        threads=2,
        memory_limit="256MB",
        max_temp_directory_size="512MB",
        allow_unsealed=True,
        **kwargs,
    )


def test_query_service_is_lazy_read_only_hardened_and_typed(tmp_path):
    database = _query_database(tmp_path / "sharur.duckdb")
    spill = tmp_path / "spill"
    with pytest.raises(ValueError, match="verified StagedDatabase"):
        create_app(db_path=database)
    app = _app(database, tmp_path)
    assert not spill.exists()

    with TestClient(app) as client:
        health = client.get("/health")
        assert health.status_code == 200
        body = health.json()
        assert body["read_only"] is True
        assert body["duckdb_owner"] == "single-http-server"
        assert body["workers"] == 1
        assert body["settings"]["enable_external_access"] == "false"
        assert body["settings"]["autoload_known_extensions"] == "false"
        assert body["settings"]["autoinstall_known_extensions"] == "false"
        assert body["settings"]["lock_configuration"] == "true"

        genomes = client.post("/v1/genomes/search", json={"limit": 10})
        assert genomes.status_code == 200
        assert genomes.json()["raw"][0]["bin_id"] == "bin1"
        assert genomes.headers["X-Sharur-Query-ID"]

        predicates = client.post(
            "/v1/predicates/search",
            json={"has": ["giant"], "limit": 10},
        )
        assert predicates.status_code == 200
        assert predicates.json()["raw"][0]["protein_id"] == "p1"

        atoms = client.post(
            "/v1/atoms/search",
            json={"atom_id": "pfam:PFTEST"},
        )
        assert atoms.status_code == 200
        assert atoms.json()["raw"][0]["protein_id"] == "p1"

        protein = client.get("/v1/proteins/p1")
        assert protein.status_code == 200
        assert protein.json()["raw"]["protein_id"] == "p1"
        assert "sequence" not in protein.json()["raw"]

        contigs = client.get("/v1/genomes/bin1/contigs")
        assert contigs.status_code == 200
        assert contigs.json()["raw"][0]["contig_id"] == "contig1"

        contig = client.get("/v1/genomes/bin1/contigs/contig1")
        assert contig.status_code == 200
        assert contig.json()["raw"]["protein_count"] == 1

        packet = client.post(
            "/v1/contigs/packet",
            json={
                "genome_id": "bin1",
                "contig_id": "contig1",
                "limit": 1,
            },
        )
        assert packet.status_code == 200
        assert packet.json()["raw"]["proteins"][0]["protein_id"] == "p1"
        assert "sequence" not in packet.text.lower()

        genome_packet = client.post(
            "/v1/genomes/packet",
            json={
                "genome_id": "bin1",
                "max_contigs": 1_000,
                "max_proteins": 1_000,
                "max_model_payload_bytes": 524_288,
            },
        )
        assert genome_packet.status_code == 200
        assert genome_packet.json()["raw"]["bin_id"] == "bin1"
        assert {
            contig["bin_id"]
            for contig in genome_packet.json()["raw"]["model_payload"]["contigs"]
        } == {"bin1"}
        assert "sequence" not in genome_packet.text.lower()

        assert client.get("/v1/proteins/p1", params={"verbosity": 2}).status_code == 422
        assert client.post("/v1/sql", json={"sql": "SELECT 1"}).status_code == 404
        assert client.post(
            "/v1/proteins/search",
            json={"limit": 501},
        ).status_code == 422

        metrics = client.get("/metrics")
        assert metrics.status_code == 200
        assert 'sharur_query_requests_total{operator="list_genomes"} 1' in metrics.text
    assert spill.is_dir()


def test_query_service_propagates_sealed_dataset_identity(tmp_path):
    database = _query_database(tmp_path / "sharur.duckdb").resolve()
    staged = StagedDatabase(
        path=database,
        source_path=database,
        seal_path=Path(tmp_path / "dataset.seal.json"),
        dataset_id="dataset-identity-1234567890",
        seal_strength="full",
        artifact_digest={"algorithm": "sha256", "value": "fixture"},
        reused=True,
        staged_at="2026-07-23T00:00:00+00:00",
    )
    app = _app(database, tmp_path, staged_database=staged)

    with TestClient(app) as client:
        response = client.get("/v1/genomes/bin1")
        body = response.json()

    assert response.status_code == 200
    assert body["trace"]["dataset_version"] == staged.dataset_id
    assert body["service"]["dataset_id"] == staged.dataset_id


def test_query_service_auth_remote_guard_and_payload_bounds(tmp_path):
    database = _query_database(tmp_path / "sharur.duckdb")
    secured = _app(database, tmp_path, api_token=TOKEN, max_request_bytes=256)
    with TestClient(secured) as client:
        assert client.get("/health").status_code == 401
        headers = {"Authorization": f"Bearer {TOKEN}"}
        who = client.get("/auth/whoami", headers=headers)
        assert who.status_code == 200
        assert who.json()["role"] == "operator"
        oversized = client.post(
            "/v1/genomes/search",
            headers=headers,
            content=b'{"taxonomy_filter":"' + b"x" * 1_000 + b'"}',
        )
        assert oversized.status_code == 413

    remote_app = _app(database, tmp_path)
    with TestClient(
        remote_app,
        client=("198.51.100.20", 50000),
    ) as remote:
        assert remote.get("/health").status_code == 403

    tiny_results = _app(database, tmp_path, max_result_bytes=100)
    with TestClient(tiny_results) as client:
        response = client.post("/v1/genomes/search", json={})
        assert response.status_code == 413
        assert "Serialized result" in response.json()["detail"]


@pytest.mark.asyncio
async def test_weighted_admission_is_fifo_bounded_and_releases_capacity():
    admission = WeightedAdmissionController(capacity=3, max_queue=1)
    await admission.acquire(3, timeout_seconds=1)
    waiting = asyncio.create_task(admission.acquire(1, timeout_seconds=1))
    for _ in range(100):
        if admission.snapshot()["queued_queries"] == 1:
            break
        await asyncio.sleep(0)
    with pytest.raises(AdmissionRejectedError):
        await admission.acquire(1, timeout_seconds=1)
    await admission.release(3)
    assert await waiting >= 0
    assert admission.snapshot()["active_units"] == 1
    await admission.release(1)
    assert admission.snapshot()["active_units"] == 0
    assert admission.snapshot()["rejected_total"] == 1


def test_runtime_shares_one_database_instance_across_thread_local_cursors(tmp_path):
    database = _query_database(tmp_path / "sharur.duckdb")
    runtime = ReadOnlyDuckDBRuntime(
        database,
        threads=2,
        memory_limit="256MB",
        temp_directory=tmp_path / "spill",
        max_temp_directory_size="512MB",
    )
    runtime.open()
    duplicate_owner = ReadOnlyDuckDBRuntime(
        database,
        threads=1,
        memory_limit="256MB",
        temp_directory=tmp_path / "other-spill",
        max_temp_directory_size="512MB",
    )
    with pytest.raises(RuntimeError, match="already owns"):
        duplicate_owner.open()
    barrier = threading.Barrier(2)
    principal = QueryPrincipal("worker", "worker")

    def query(identifier):
        return runtime.execute(
            query_id=identifier,
            principal=principal,
            operator="parallel_count",
            timeout_seconds=2,
            operation=lambda store: (
                barrier.wait(timeout=1),
                store.execute("SELECT COUNT(*) FROM proteins")[0][0],
            )[1],
        ).value

    with ThreadPoolExecutor(max_workers=2) as executor:
        assert list(executor.map(query, ["q1", "q2"])) == [1, 1]
    snapshot = runtime.snapshot()
    assert snapshot["cursor_count"] == 2
    assert snapshot["operators"]["parallel_count"]["queries"] == 2
    runtime.close()


def test_runtime_enforces_timeout_cancellation_and_read_only_mode(tmp_path):
    database = _query_database(tmp_path / "sharur.duckdb")
    runtime = ReadOnlyDuckDBRuntime(
        database,
        threads=1,
        memory_limit="256MB",
        temp_directory=tmp_path / "spill",
        max_temp_directory_size="512MB",
    )
    runtime.open()
    owner = QueryPrincipal("owner", "worker")
    other = QueryPrincipal("other", "worker")

    with pytest.raises(duckdb.Error, match="read-only"):
        runtime.execute(
            query_id="write-attempt",
            principal=owner,
            operator="forbidden_write",
            timeout_seconds=1,
            operation=lambda store: store.execute("CREATE TABLE forbidden(i INT)"),
        )

    with pytest.raises(QueryTimeoutError):
        runtime.execute(
            query_id="timeout",
            principal=owner,
            operator="slow_fixture",
            timeout_seconds=0.01,
            operation=lambda _store: time.sleep(0.03),
        )
    recovered = runtime.execute(
        query_id="after-timeout",
        principal=owner,
        operator="recovery_fixture",
        timeout_seconds=1,
        operation=lambda store: store.execute("SELECT 1")[0][0],
    )
    assert recovered.value == 1

    release = threading.Event()
    outcome = {}

    def cancellable():
        try:
            runtime.execute(
                query_id="cancel-me",
                principal=owner,
                operator="cancel_fixture",
                timeout_seconds=2,
                operation=lambda _store: release.wait(1),
            )
        except Exception as exc:  # captured for assertion in the parent thread
            outcome["error"] = exc

    thread = threading.Thread(target=cancellable)
    thread.start()
    for _ in range(100):
        if runtime.snapshot()["active_queries"] == 1:
            break
        time.sleep(0.005)
    with pytest.raises(PermissionError):
        runtime.cancel("cancel-me", principal=other)
    assert runtime.cancel("cancel-me", principal=owner)
    release.set()
    thread.join(timeout=2)
    assert isinstance(outcome["error"], QueryCancelledError)
    runtime.close()


class _Response:
    status_code = 200

    def raise_for_status(self):
        return None

    def json(self):
        return {"agent_id": "agent-a", "role": "reader", "bootstrap": False}


class _IntrospectionSession:
    def __init__(self):
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        return _Response()


def test_ops_token_introspection_caches_hash_keyed_principals():
    session = _IntrospectionSession()
    introspector = OpsTokenIntrospector(
        "http://ops:8811",
        session=session,
        cache_ttl_seconds=30,
    )
    first = introspector.introspect("agent-secret")
    second = introspector.introspect("agent-secret")
    assert first == second == QueryPrincipal("agent-a", "reader")
    assert len(session.calls) == 1
    assert session.calls[0][1]["headers"]["Authorization"] == "Bearer agent-secret"
    assert "agent-secret" not in repr(introspector._cache)

    bounded = OpsTokenIntrospector(
        "http://ops:8811",
        session=session,
        cache_ttl_seconds=30,
        max_cache_entries=1,
    )
    bounded.introspect("first-secret")
    bounded.introspect("second-secret")
    assert len(bounded._cache) == 1


def test_ops_token_introspection_fails_closed_when_ops_is_unavailable():
    class FailingSession:
        def get(self, *_args, **_kwargs):
            raise requests.ConnectionError("offline")

    introspector = OpsTokenIntrospector(
        "http://ops:8811",
        session=FailingSession(),
    )
    with pytest.raises(AuthServiceUnavailableError, match="ConnectionError"):
        introspector.introspect("agent-secret")


def test_query_client_generates_query_ids_and_escapes_entity_paths():
    class ClientResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return {"status": "ok"}

    class ClientSession:
        def __init__(self):
            self.headers = {}
            self.calls = []

        def request(self, method, url, **kwargs):
            self.calls.append((method, url, kwargs))
            return ClientResponse()

    session = ClientSession()
    query = SharurQuery(
        "http://query:8812/",
        api_token="agent-token",
        session=session,
        verify_connection=False,
    )
    assert query.get_protein("protein/with/slash") == {"status": "ok"}
    method, url, kwargs = session.calls[0]
    assert method == "GET"
    assert url == "http://query:8812/v1/proteins/protein%2Fwith%2Fslash"
    assert kwargs["headers"]["X-Sharur-Query-ID"]
    assert session.headers["Authorization"] == "Bearer agent-token"

    assert query.list_contigs("genome/with/slash", limit=5) == {"status": "ok"}
    assert (
        session.calls[1][1]
        == "http://query:8812/v1/genomes/genome%2Fwith%2Fslash/contigs"
    )
    assert query.get_contig("genome", "contig/with/slash") == {"status": "ok"}
    assert (
        session.calls[2][1]
        == "http://query:8812/v1/genomes/genome/contigs/contig%2Fwith%2Fslash"
    )
    assert query.contig_packet("genome", "contig", limit=25) == {"status": "ok"}
    assert session.calls[3][1] == "http://query:8812/v1/contigs/packet"
    assert session.calls[3][2]["json"]["limit"] == 25
    assert query.genome_packet("genome", max_proteins=500) == {"status": "ok"}
    assert session.calls[4][1] == "http://query:8812/v1/genomes/packet"
    assert session.calls[4][2]["json"]["genome_id"] == "genome"
    assert session.calls[4][2]["json"]["max_proteins"] == 500
    assert "max_contigs" not in session.calls[4][2]["json"]

    with pytest.raises(ValueError, match="timeout"):
        SharurQuery(session=session, timeout=0, verify_connection=False)
