"""Semantic-term providers: one contract over ``semantic_terms`` (SQL) or a compact generation.

Both providers answer the same questions with the same results:

- ``search_compatible(has, lacks, limit)`` — the live ``search_by_atoms`` contract
  (unordered limited protein IDs; empty request and repeated ``has`` values select
  nothing; lacks-only starts from the active universe; unknown positives select
  nothing; unknown exclusions change nothing).
- ``search(request, offset, limit)`` — the browser extension: AND/OR/NOT over distinct
  terms, an optional genome scope, an exact total and a page of protein IDs in
  canonical UTF-8 order. ``genomes=None`` skips the scope; an empty tuple selects none.
- ``rich_rows(protein_id)`` — every stored row as an exact six-field tuple
  ``(term_id, term_kind, facet, relation, source_db, source_accession)`` in the API order
  ``term_kind, term_id, source_db, source_accession`` with ``facet, relation`` breaking ties,
  NULLs last. NULL stays None and blank strings stay blank.
- ``term_catalog()`` — every active term with its distinct active protein count.

Active membership is exactly ``(term_kind != 'atom' OR relation != 'excludes')`` with SQL
NULL semantics. Attach a provider to a store with :func:`attach`; operators find it with
:func:`attached`.
"""

from __future__ import annotations

import contextlib
import re
import threading
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from sharur.semantic_index.generation import Generation, GenerationError, resolve, verify


if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


ACTIVE = "(term_kind != 'atom' OR relation != 'excludes')"
RICH_ORDER = "term_kind, term_id, source_db, source_accession, facet, relation"
MAX_TERMS = 64


def _distinct(values: Iterable[str] | None) -> tuple[str, ...]:
    return tuple(dict.fromkeys(v for v in (values or ()) if v))


@dataclass(frozen=True)
class SearchRequest:
    """A normalized semantic-term search; operands are distinct and order-preserving."""

    has: tuple[str, ...] = ()
    any_of: tuple[str, ...] = ()
    lacks: tuple[str, ...] = ()
    genomes: tuple[str, ...] | None = None

    @classmethod
    def build(cls, has=None, any_of=None, lacks=None, genomes=None) -> SearchRequest:
        if genomes is not None:
            if isinstance(genomes, str) or any(not isinstance(g, str) for g in genomes):
                raise ValueError("genomes takes a list of genome ID strings; NULL ownership is not a genome")
            genomes = tuple(dict.fromkeys(genomes))
        request = cls(_distinct(has), _distinct(any_of), _distinct(lacks), genomes)
        if max(len(request.has), len(request.any_of), len(request.lacks)) > MAX_TERMS:
            raise ValueError(f"At most {MAX_TERMS} terms per operand")
        return request

    @property
    def empty(self) -> bool:
        return not (self.has or self.any_of or self.lacks)

    def terms(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(self.has + self.any_of + self.lacks))


@dataclass
class SearchPage:
    total: int
    protein_ids: list[str]
    elapsed_ms: float


def rank_terms(catalog: list[tuple[str, int]], q: str, limit: int = 20) -> list[tuple[str, int]]:
    """Terms whose ID contains ``q`` (case-insensitive): exact, prefix, token start, substring; then by count."""
    needle = q.strip().lower()
    if not needle:
        return []
    token = re.compile(r"[:_.\-]" + re.escape(needle))
    scored = []
    for term, count in catalog:
        low = term.lower()
        if needle not in low:
            continue
        if low == needle:
            score = 0
        elif low.startswith(needle):
            score = 1
        elif token.search(low):
            score = 2
        else:
            score = 3
        scored.append((score, -count, term, count))
    scored.sort()
    return [(term, count) for _, _, term, count in scored[:limit]]


class _Counters:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self.calls: dict[str, int] = {}

    def hit(self, name: str) -> None:
        with self._lock:
            self.calls[name] = self.calls.get(name, 0) + 1

    def snapshot(self) -> dict[str, int]:
        with self._lock:
            return dict(self.calls)


class SqlSemanticProvider:
    """``semantic_terms`` in the dataset's DuckDB; the default and the comparison baseline."""

    backend = "sql"
    needs_store_lock = True
    has_rich_rows = True
    # The term catalog is a GROUP BY over every row; the main search box skips it.
    cheap_catalog = False

    def __init__(self, store) -> None:
        self.store = store
        self.counters = _Counters()
        self._catalog: list[tuple[str, int]] | None = None
        self._catalog_lock = threading.Lock()
        self._available: bool | None = None

    def identity(self) -> dict[str, Any]:
        return {"backend": self.backend, "table": "semantic_terms"}

    def available(self) -> bool:
        if self._available is None:
            self._available = bool(self.store.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_name = 'semantic_terms'"))
        return self._available

    def search_compatible(self, has, lacks, limit: int) -> list[str]:
        from sharur.operators import predicates_v2  # noqa: PLC0415  (operators import this module)

        self.counters.hit("search_compatible")
        has, lacks = list(has or []), list(lacks or [])
        if not has and not lacks:
            return []
        return predicates_v2._search_by_atoms_from_semantic_terms(self.store, has, lacks, limit)

    def _result_sql(self, request: SearchRequest) -> tuple[str, list]:
        ctes, params, filters = [], [], []

        def placeholders(values):
            params.extend(values)
            return ",".join("?" * len(values))

        if request.has:
            ctes.append(f"has_set AS (SELECT protein_id FROM semantic_terms WHERE {ACTIVE} "
                        f"AND term_id IN ({placeholders(request.has)}) GROUP BY protein_id "
                        "HAVING COUNT(DISTINCT term_id) = ?)")
            params.append(len(request.has))
            base = "has_set"
        if request.any_of:
            ctes.append(f"any_set AS (SELECT DISTINCT protein_id FROM semantic_terms WHERE {ACTIVE} "
                        f"AND term_id IN ({placeholders(request.any_of)}))")
            if request.has:
                filters.append("protein_id IN (SELECT protein_id FROM any_set)")
            else:
                base = "any_set"
        if not request.has and not request.any_of:
            ctes.append(f"universe AS (SELECT DISTINCT protein_id FROM semantic_terms WHERE {ACTIVE})")
            base = "universe"
        if request.lacks:
            ctes.append(f"lacks_set AS (SELECT DISTINCT protein_id FROM semantic_terms WHERE {ACTIVE} "
                        f"AND term_id IN ({placeholders(request.lacks)}) AND protein_id IS NOT NULL)")
            filters.append("protein_id NOT IN (SELECT protein_id FROM lacks_set)")
        if request.genomes is not None:
            ctes.append("scope AS (SELECT protein_id FROM proteins WHERE bin_id IN (SELECT UNNEST(?::VARCHAR[])))")
            params.append(list(request.genomes))
            filters.append("protein_id IN (SELECT protein_id FROM scope)")
        where = f" WHERE {' AND '.join(filters)}" if filters else ""
        return f"WITH {', '.join(ctes)}, result AS (SELECT protein_id FROM {base}{where})", params

    def search(self, request: SearchRequest, *, offset: int = 0, limit: int = 50) -> SearchPage:
        self.counters.hit("search")
        started = time.perf_counter()
        if request.empty or request.genomes == ():
            return SearchPage(0, [], (time.perf_counter() - started) * 1000)
        prefix, params = self._result_sql(request)
        total = self.store.execute(f"{prefix} SELECT COUNT(*) FROM result", params)[0][0]
        rows = self.store.execute(f"{prefix} SELECT protein_id FROM result ORDER BY protein_id LIMIT ? OFFSET ?",
                                  [*params, limit, offset]) if total > offset else []
        return SearchPage(int(total), [r[0] for r in rows], (time.perf_counter() - started) * 1000)

    def rich_rows(self, protein_id: str) -> list[tuple]:
        self.counters.hit("rich_rows")
        return [tuple(r) for r in self.store.execute(
            "SELECT term_id, term_kind, facet, relation, source_db, source_accession FROM semantic_terms "
            f"WHERE protein_id = ? ORDER BY {RICH_ORDER}", [protein_id])]

    def rich_row_count(self, protein_id: str) -> int:
        return int(self.store.execute("SELECT COUNT(*) FROM semantic_terms WHERE protein_id = ?", [protein_id])[0][0])

    def rich_rows_many(self, protein_ids: list[str]) -> dict[str, list[tuple]]:
        self.counters.hit("rich_rows_many")
        out: dict[str, list[tuple]] = {pid: [] for pid in protein_ids}
        if protein_ids:
            for row in self.store.execute(
                    "SELECT protein_id, term_id, term_kind, facet, relation, source_db, source_accession "
                    "FROM semantic_terms WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) "
                    f"ORDER BY protein_id, {RICH_ORDER}", [list(protein_ids)]):
                out[row[0]].append(tuple(row[1:]))
        return out

    def term_catalog(self) -> list[tuple[str, int]]:
        with self._catalog_lock:
            if self._catalog is None:
                self.counters.hit("term_catalog")
                self._catalog = [(t, int(n)) for t, n in self.store.execute(
                    f"SELECT term_id, COUNT(DISTINCT protein_id) FROM semantic_terms WHERE {ACTIVE} "
                    "GROUP BY term_id ORDER BY encode(term_id)")]
            return self._catalog

    def term_counts(self, terms: Iterable[str]) -> dict[str, int]:
        terms = list(dict.fromkeys(terms))
        if self._catalog is not None:
            counts = dict(self._catalog)
        elif terms:
            counts = dict(self.store.execute(
                f"SELECT term_id, COUNT(DISTINCT protein_id) FROM semantic_terms WHERE {ACTIVE} "
                "AND term_id IN (SELECT UNNEST(?::VARCHAR[])) GROUP BY term_id", [terms]))
        else:
            counts = {}
        return {t: int(counts.get(t, 0)) for t in terms}

    def suggest_terms(self, q: str, limit: int = 20) -> list[tuple[str, int]]:
        self.counters.hit("suggest_terms")
        return rank_terms(self.term_catalog(), q, limit)

    def stats(self) -> dict[str, Any]:
        rows = self.store.execute("SELECT COUNT(*) FROM semantic_terms")[0][0]
        return {"rich_rows": int(rows), "terms": len(self.term_catalog())}

    def close(self) -> None:
        self._catalog = None


class CompactSemanticProvider:
    """A verified compact generation: Roaring/list membership, rich CSR rows and genome postings."""

    backend = "compact"
    needs_store_lock = False
    cheap_catalog = True

    def __init__(self, generation: Generation, membership, forward, scope) -> None:
        self.generation = generation
        self.membership, self.forward, self.scope = membership, forward, scope
        # A generation without the rows component serves membership only; callers read rows from SQL.
        self.has_rich_rows = forward is not None
        self.counters = _Counters()
        self._catalog: list[tuple[str, int]] | None = None
        self._closed = False

    @classmethod
    def open(cls, index_dir: str | Path, db_path: str | Path, *, verify_payloads: bool = True,
             store=None, content_match: bool = False, seal_path: str | Path | None = None,
             ) -> CompactSemanticProvider:
        """Resolve CURRENT, verify the source binding and payloads, then open the readers once."""
        from sharur.semantic_index.readers import (  # noqa: PLC0415
            ForwardIndex,
            GenomeScope,
            PostingIndex,
        )

        generation = resolve(index_dir)
        verify(generation, db_path, payloads=verify_payloads, content_match=content_match, seal_path=seal_path)
        membership = PostingIndex(generation.component("membership"))
        try:
            forward = (ForwardIndex(generation.component("forward"), generation.component("membership"))
                       if "forward" in generation.record["components"] else None)
            try:
                scope = GenomeScope(generation.component("genome_scope"), generation.component("membership"))
            except BaseException:
                if forward is not None:
                    forward.close()
                raise
        except BaseException:
            membership.close()
            raise
        provider = cls(generation, membership, forward, scope)
        if store is not None:
            try:
                provider.check_store(store)
            except BaseException:
                provider.close()
                raise
        return provider

    def check_store(self, store) -> None:
        """Counts the open dataset must agree with: proteins always, rich rows while the SQL table exists."""
        counts = self.generation.record["counts"]
        proteins = store.execute("SELECT COUNT(*) FROM proteins")[0][0]
        if proteins != counts["proteins"]:
            raise GenerationError(f"Dataset has {proteins} proteins; generation has {counts['proteins']}")
        tables = {r[0] for r in store.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'")}
        if "semantic_terms" in tables and counts["rich_rows"] is not None:
            rows = store.execute("SELECT COUNT(*) FROM semantic_terms")[0][0]
            if rows != counts["rich_rows"]:
                raise GenerationError(f"semantic_terms has {rows} rows; generation has {counts['rich_rows']}")

    def _open(self) -> None:
        if self._closed:
            raise GenerationError("Compact semantic provider is closed")

    def _rows(self) -> None:
        self._open()
        if self.forward is None:
            raise GenerationError("This generation serves membership only; read stored rows from SQL")

    def identity(self) -> dict[str, Any]:
        return {"backend": self.backend, **self.generation.identity(), "root": str(self.generation.root),
                "rich_rows": "compact" if self.has_rich_rows else "sql"}

    def available(self) -> bool:
        return not self._closed

    def search_compatible(self, has, lacks, limit: int) -> list[str]:
        self._open()
        self.counters.hit("search_compatible")
        result = self.membership.search_compatible(has=list(has or []), lacks=list(lacks or []))
        return self.membership.decode(result, limit=limit)

    def search(self, request: SearchRequest, *, offset: int = 0, limit: int = 50) -> SearchPage:
        self._open()
        self.counters.hit("search")
        started = time.perf_counter()
        if request.empty or request.genomes == ():
            return SearchPage(0, [], (time.perf_counter() - started) * 1000)
        result = self.membership.search(has=request.has, any_of=request.any_of or None, lacks=request.lacks)
        if request.genomes is not None and result:
            result &= self.scope.scope(request.genomes)
        ids = self.membership.decode(result, limit=limit, offset=offset) if len(result) > offset else []
        return SearchPage(len(result), ids, (time.perf_counter() - started) * 1000)

    def rich_rows(self, protein_id: str) -> list[tuple]:
        self._rows()
        self.counters.hit("rich_rows")
        return self.forward.rows_for_protein(protein_id)

    def rich_row_count(self, protein_id: str) -> int:
        self._rows()
        pid = self.forward.proteins.find(protein_id)
        return 0 if pid is None else self.forward.row_count(pid)

    def rich_rows_many(self, protein_ids: list[str]) -> dict[str, list[tuple]]:
        self._rows()
        self.counters.hit("rich_rows_many")
        return {pid: self.forward.rows_for_protein(pid) for pid in protein_ids}

    def term_catalog(self) -> list[tuple[str, int]]:
        self._open()
        if self._catalog is None:
            terms, cardinality = self.membership.terms, self.membership.catalog["cardinality"].tolist()
            self._catalog = [(terms.get(i), int(cardinality[i])) for i in range(len(terms))]
        return self._catalog

    def term_counts(self, terms: Iterable[str]) -> dict[str, int]:
        self._open()
        return {t: self.membership.count(t) for t in terms}

    def suggest_terms(self, q: str, limit: int = 20) -> list[tuple[str, int]]:
        self.counters.hit("suggest_terms")
        return rank_terms(self.term_catalog(), q, limit)

    def stats(self) -> dict[str, Any]:
        self._open()
        rows = int(self.forward.manifest["rows"]) if self.forward is not None else None
        return {"rich_rows": rows, "terms": len(self.membership.terms)}

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._catalog = None
        for reader in (self.scope, self.forward, self.membership):
            if reader is not None:
                reader.close()


def attach(store, provider) -> None:
    """Route this store's semantic-term operators through ``provider``."""
    store.semantic_provider = provider


def attached(store):
    """The provider attached to ``store`` (through any proxy), or None."""
    return getattr(store, "semantic_provider", None)


def provider_for(store):
    """The attached provider, else an SQL provider over ``semantic_terms`` cached on the store."""
    provider = attached(store)
    if provider is None:
        provider = getattr(store, "_sql_semantic_provider", None)
        if provider is None:
            provider = SqlSemanticProvider(store)
            with contextlib.suppress(AttributeError):
                store._sql_semantic_provider = provider
    return provider
