"""Readers for the compact semantic-term artifacts (membership, rich rows, genome scope).

Ported from the certified isolated prototype (``search_index_prototype_2026-10-04``)
with its file formats unchanged, so certified generations open as-is:

- membership (``sharur-isolated-membership-v1``): per-term sorted uint32 lists or
  portable Roaring bitmaps of numeric protein IDs, plus the active universe.
- rich forward rows (``sharur-isolated-rich-forward-v1``): every original
  ``semantic_terms`` row as a dictionary-coded six-field record, CSR by protein.
- genome scope (``sharur-isolated-genome-scope-v1``): per-genome protein postings.

Protein and term IDs are dense uint32 codes in UTF-8 byte order of the canonical
strings, so numeric order is canonical ID order. All payloads are memory-mapped
read-only; ``verify_hashes=True`` rehashes every payload against its manifest.
"""

from __future__ import annotations

import json
import mmap
import os
import sys
from pathlib import Path

import numpy as np
from pyroaring import BitMap

from sharur.semantic_index.formats import (
    CONTROL_FILES,
    FIELDS,
    FORWARD_FORMAT,
    MEMBERSHIP_FORMAT,
    SCOPE_FORMAT,
    sha256_file,
)


CATALOG = np.dtype([("offset", "<u8"), ("bytes", "<u8"), ("cardinality", "<u4"),
                    ("encoding", "u1"), ("padding", "u1", (3,))])


class ArtifactError(ValueError):
    """An artifact is missing, malformed, or differs from its manifest."""


def _check_files(entries, verify_hashes: bool) -> None:
    for root, name, expected in entries:
        path = Path(root) / name
        if not path.is_file():
            raise ArtifactError(f"Artifact file missing: {path}")
        if path.stat().st_size != expected["bytes"]:
            raise ArtifactError(f"Artifact size differs from its manifest: {path}")
        if verify_hashes and sha256_file(path) != expected["sha256"]:
            raise ArtifactError(f"Artifact checksum differs from its manifest: {path}")


def _map(path: Path):
    handle = path.open("rb")
    if os.fstat(handle.fileno()).st_size == 0:
        return handle, b""
    return handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)


class PackedDictionary:
    """Dense zero-based string dictionary: UTF-8 blob plus uint64 offsets with a final sentinel."""

    def __init__(self, root: str | Path, prefix: str) -> None:
        root = Path(root)
        self.offsets = np.memmap(root / f"{prefix}.offsets.u64", dtype="<u8", mode="r")
        self.handle, self.blob = _map(root / f"{prefix}.utf8")

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def raw(self, i: int) -> bytes:
        if i < 0 or i >= len(self):
            raise IndexError(i)
        return self.blob[int(self.offsets[i]):int(self.offsets[i + 1])]

    def get(self, i: int) -> str:
        return self.raw(i).decode("utf-8")

    def find(self, value: str) -> int | None:
        needle, low, high = value.encode("utf-8"), 0, len(self)
        while low < high:
            mid = (low + high) // 2
            if self.raw(mid) < needle:
                low = mid + 1
            else:
                high = mid
        return low if low < len(self) and self.raw(low) == needle else None

    def close(self) -> None:
        self.offsets._mmap.close()
        if hasattr(self.blob, "close"):
            self.blob.close()
        self.handle.close()


class PostingIndex:
    """Active protein/term membership: per-term numeric protein sets."""

    def __init__(self, root: str | Path, *, verify_hashes: bool = False) -> None:
        self.root = Path(root)
        self.manifest = json.loads((self.root / "manifest.json").read_text())
        if self.manifest.get("format") != MEMBERSHIP_FORMAT:
            raise ArtifactError("Unknown membership artifact format")
        _check_files(((self.root, n, e) for n, e in self.manifest["files"].items() if n not in CONTROL_FILES),
                     verify_hashes)
        self.proteins = PackedDictionary(self.root, "proteins")
        self.terms = PackedDictionary(self.root, "terms")
        n = self.manifest["term_count"]
        self.catalog = (np.memmap(self.root / "catalog.bin", dtype=CATALOG, mode="r") if n
                        else np.empty(0, dtype=CATALOG))
        if len(self.proteins) != self.manifest["protein_count"] or len(self.terms) != n or len(self.catalog) != n:
            raise ArtifactError("Membership dictionary/catalog count mismatch")
        self.handle, self.postings = _map(self.root / "postings.bin")
        self._universe: BitMap | None = None

    def term_id(self, term: str) -> int | None:
        return self.terms.find(term)

    def count(self, term: str) -> int:
        tid = self.terms.find(term)
        return 0 if tid is None else int(self.catalog[tid]["cardinality"])

    def bitmap_for_tid(self, tid: int) -> BitMap:
        row = self.catalog[tid]
        offset, cardinality = int(row["offset"]), int(row["cardinality"])
        if int(row["encoding"]) == 0:
            result = BitMap(np.frombuffer(self.postings, dtype="<u4", count=cardinality, offset=offset))
        else:
            result = BitMap.deserialize(memoryview(self.postings)[offset:offset + int(row["bytes"])])
        if len(result) != cardinality:
            raise ArtifactError("Posting cardinality differs from its catalog entry")
        return result

    def bitmap(self, term: str) -> BitMap:
        tid = self.terms.find(term)
        return BitMap() if tid is None else self.bitmap_for_tid(tid)

    def universe(self) -> BitMap:
        if self._universe is None:
            universe = BitMap.deserialize((self.root / "active_universe.roaring").read_bytes())
            if len(universe) != self.manifest["active_protein_count"]:
                raise ArtifactError("Active universe size differs from its manifest")
            self._universe = universe
        return self._universe.copy()

    def union(self, terms) -> BitMap:
        result = BitMap()
        for term in dict.fromkeys(terms):
            result |= self.bitmap(term)
        return result

    def intersection(self, terms) -> BitMap:
        names = sorted(set(terms), key=self.count)
        if not names:
            return self.universe()
        result = self.bitmap(names[0])
        for name in names[1:]:
            if not result:
                break
            result &= self.bitmap(name)
        return result

    def search(self, *, has=(), any_of=None, lacks=()) -> BitMap:
        """Set algebra with idempotent operands: AND over ``has``, OR over ``any_of``, minus ``lacks``."""
        result = self.intersection(has)
        if any_of is not None:
            result &= self.union(any_of)
        if lacks and result:
            result -= self.union(lacks)
        return result

    def search_compatible(self, *, has=(), lacks=()) -> BitMap:
        """The live ``search_by_atoms`` contract: empty request and repeated ``has`` values select nothing."""
        has, lacks = list(has), list(lacks)
        if (not has and not lacks) or len(has) != len(set(has)):
            return BitMap()
        return self.search(has=has, lacks=lacks)

    def decode(self, result: BitMap, *, limit: int | None = None, offset: int = 0) -> list[str]:
        stop = None if limit is None else offset + limit
        return [self.proteins.get(pid) for pid in result[offset:stop]] if len(result) else []

    def close(self) -> None:
        if hasattr(self.catalog, "_mmap"):
            self.catalog._mmap.close()
        self.proteins.close()
        self.terms.close()
        if hasattr(self.postings, "close"):
            self.postings.close()
        self.handle.close()


class ForwardIndex:
    """Every original rich term row, contiguous per numeric protein ID (CSR)."""

    def __init__(self, root: str | Path, matching: str | Path, *, verify_hashes: bool = False) -> None:
        self.root, self.matching = Path(root), Path(matching)
        self.manifest = json.loads((self.root / "manifest.json").read_text())
        if self.manifest.get("format") != FORWARD_FORMAT:
            raise ArtifactError("Unknown rich forward artifact format")
        entries = [(self.root, n, e) for n, e in self.manifest["files"].items()]
        entries.extend((self.matching, n, e) for n, e in self.manifest["shared_protein_dictionary"]["files"].items())
        _check_files(entries, verify_hashes)
        self.proteins = PackedDictionary(self.matching, "proteins")
        self.dictionaries = {f: PackedDictionary(self.root, f) for f in FIELDS}
        self.dtype = np.dtype([tuple(v) for v in self.manifest["row_dtype"]])
        if self.dtype.names != FIELDS:
            raise ArtifactError("Rich row fields differ from the expected six-field layout")
        self.offsets = np.memmap(self.root / "protein_offsets.u64", dtype="<u8", mode="r")
        count = self.manifest["rows"]
        self.rows = (np.memmap(self.root / "rows.bin", dtype=self.dtype, mode="r") if count
                     else np.empty(0, dtype=self.dtype))
        if len(self.offsets) != len(self.proteins) + 1 or int(self.offsets[-1]) != count or len(self.rows) != count:
            raise ArtifactError("Rich forward shape/count mismatch")
        # Metadata dictionaries (38,578 labels for DPANN) decode once; protein IDs stay mapped.
        caches = []
        for field in FIELDS:
            dictionary = self.dictionaries[field]
            declared = self.manifest["dictionaries"][field]["non_null_values"]
            if len(dictionary) != declared:
                raise ArtifactError("Rich metadata dictionary differs from its manifest")
            caches.append((None, *(sys.intern(dictionary.get(i)) for i in range(declared))))
        self._labels = tuple(caches)

    def raw_rows_for_pid(self, pid: int):
        if pid < 0 or pid >= len(self.proteins):
            raise IndexError(pid)
        return self.rows[int(self.offsets[pid]):int(self.offsets[pid + 1])]

    def row_count(self, pid: int) -> int:
        if pid < 0 or pid >= len(self.proteins):
            raise IndexError(pid)
        return int(self.offsets[pid + 1] - self.offsets[pid])

    def rows_for_pid(self, pid: int) -> list[tuple]:
        """Original six-field tuples; code 0 decodes to None (SQL NULL), blank strings stay blank."""
        labels = self._labels
        return [tuple(values[code] for values, code in zip(labels, row, strict=True))
                for row in self.raw_rows_for_pid(pid).tolist()]

    def rows_for_protein(self, protein_id: str) -> list[tuple]:
        pid = self.proteins.find(protein_id)
        return [] if pid is None else self.rows_for_pid(pid)

    def close(self) -> None:
        self._labels = ()
        self.offsets._mmap.close()
        if hasattr(self.rows, "_mmap"):
            self.rows._mmap.close()
        for dictionary in self.dictionaries.values():
            dictionary.close()
        self.proteins.close()


class GenomeScope:
    """Protein ownership by genome (non-NULL ``proteins.bin_id``), sharing the membership protein IDs."""

    def __init__(self, root: str | Path, matching: str | Path, *, verify_hashes: bool = False) -> None:
        self.root, self.matching = Path(root), Path(matching)
        self.manifest = json.loads((self.root / "manifest.json").read_text())
        if self.manifest.get("format") != SCOPE_FORMAT:
            raise ArtifactError("Unknown genome scope artifact format")
        entries = [(self.root, n, e) for n, e in self.manifest["files"].items()]
        entries.extend((self.matching, n, e) for n, e in self.manifest["shared_protein_dictionary"]["files"].items())
        _check_files(entries, verify_hashes)
        self.genomes = PackedDictionary(self.root, "genomes")
        count = self.manifest["genome_count"]
        self.catalog = (np.memmap(self.root / "catalog.bin", dtype=CATALOG, mode="r") if count
                        else np.empty(0, dtype=CATALOG))
        if len(self.catalog) != count or len(self.genomes) != count:
            raise ArtifactError("Genome scope dictionary/catalog count mismatch")
        self.handle, self.postings = _map(self.root / "postings.bin")

    def bitmap(self, genome: str) -> BitMap:
        gid = self.genomes.find(genome)
        if gid is None:
            return BitMap()
        row = self.catalog[gid]
        offset, cardinality = int(row["offset"]), int(row["cardinality"])
        if int(row["encoding"]) == 0:
            result = BitMap(np.frombuffer(self.postings, dtype="<u4", offset=offset, count=cardinality))
        else:
            result = BitMap.deserialize(memoryview(self.postings)[offset:offset + int(row["bytes"])])
        if len(result) != cardinality:
            raise ArtifactError("Genome posting cardinality differs from its catalog entry")
        return result

    def scope(self, genomes) -> BitMap:
        """Union of the named genomes' proteins. Callers skip filtering for an absent scope (``None``)."""
        if genomes is None:
            return BitMap()
        if isinstance(genomes, str):
            raise TypeError("scope expects an iterable of genome IDs")
        result = BitMap()
        for genome in dict.fromkeys(genomes):
            if genome is None:
                raise ValueError("NULL ownership is explicit via unassigned(), not a genome ID")
            result |= self.bitmap(genome)
        return result

    def unassigned(self) -> BitMap:
        result = BitMap.deserialize((self.root / "null_ownership.roaring").read_bytes())
        if len(result) != self.manifest["null_ownership_count"]:
            raise ArtifactError("NULL ownership bitmap count differs from its manifest")
        return result

    def close(self) -> None:
        if hasattr(self.catalog, "_mmap"):
            self.catalog._mmap.close()
        self.genomes.close()
        if hasattr(self.postings, "close"):
            self.postings.close()
        self.handle.close()
