"""Immutable compact-index generations bound to one exact source database.

Layout of an index directory::

    INDEX/
      CURRENT                         {"generation": "<id>"}, replaced atomically
      generations/<id>/
        generation.json               source binding + component manifest hashes
        membership/  forward/  genome_scope/

A generation is complete and immutable once its directory is published. Opening one
verifies, before any request is served: the generation and component formats; each
component manifest hash and the shared protein-dictionary binding; the source's size,
mtime and inode (and its seal's dataset ID, when recorded) against the adopted state; and,
by default, every payload checksum. A mismatch raises :class:`StaleGenerationError`
or :class:`GenerationError`; there is no silent fallback.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from sharur.semantic_index.formats import (
    CONTROL_FILES,
    FORWARD_FORMAT,
    MEMBERSHIP_FORMAT,
    SCOPE_FORMAT,
    sha256_file,
)


GENERATION_FORMAT = "sharur-semantic-index-generation-v1"
COMPONENTS = {"membership": MEMBERSHIP_FORMAT, "forward": FORWARD_FORMAT, "genome_scope": SCOPE_FORMAT}
POINTER = "CURRENT"


class GenerationError(RuntimeError):
    """A generation is missing, malformed, or internally inconsistent."""


class StaleGenerationError(GenerationError):
    """A generation was adopted for a different state of the source database."""


@dataclass(frozen=True)
class Generation:
    root: Path
    generation_id: str
    record: dict[str, Any]

    def component(self, name: str) -> Path:
        return self.root / self.record["components"][name]["dir"]

    def identity(self) -> dict[str, Any]:
        source = self.record["source"]
        return {"generation_id": self.generation_id, "format": self.record["format"],
                "source_sha256": source["sha256"], "seal_dataset_id": (source.get("seal") or {}).get("dataset_id"),
                "components": {k: v["manifest_sha256"] for k, v in self.record["components"].items()}}


def source_state(db_path: str | Path) -> dict[str, int]:
    p = Path(db_path).resolve()
    s = p.stat()
    return {"bytes": s.st_size, "mtime_ns": s.st_mtime_ns, "inode": s.st_ino}


def _write_json_atomic(path: Path, value: Any) -> None:
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(json.dumps(value, indent=2) + "\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _copy(src: Path, dst: Path) -> None:
    """Copy one immutable payload; APFS clones share blocks until either side changes."""
    if sys.platform == "darwin":
        if subprocess.run(["cp", "-c", str(src), str(dst)], check=False, capture_output=True).returncode == 0:
            return
    shutil.copy2(src, dst)


def _manifest(directory: Path, expected_format: str) -> tuple[dict, str]:
    path = directory / "manifest.json"
    if not path.is_file():
        raise GenerationError(f"Component manifest missing: {path}")
    manifest = json.loads(path.read_text())
    if manifest.get("format") != expected_format:
        raise GenerationError(f"Component format differs from {expected_format}: {directory}")
    return manifest, sha256_file(path)


def _check_bindings(manifests: dict[str, tuple[dict, str]]) -> None:
    membership, membership_sha = manifests["membership"]
    shared = {name: membership["files"][name] for name in ("proteins.utf8", "proteins.offsets.u64")}
    for name in ("forward", "genome_scope"):
        manifest = manifests[name][0]
        if manifest["matching_manifest_sha256"] != membership_sha:
            raise GenerationError(f"{name} was built against a different membership manifest")
        if manifest["shared_protein_dictionary"]["files"] != shared:
            raise GenerationError(f"{name} binds a different shared protein dictionary")
        if manifest["source"]["sha256"] != membership["source"]["sha256"]:
            raise GenerationError(f"{name} was built from a different source")


def _verify_payloads(directory: Path, manifest: dict) -> None:
    for name, expected in manifest["files"].items():
        if name in CONTROL_FILES:
            continue
        path = directory / name
        if not path.is_file() or path.stat().st_size != expected["bytes"] or sha256_file(path) != expected["sha256"]:
            raise GenerationError(f"Payload differs from its manifest: {path}")


def assemble(index_dir: str | Path, db_path: str | Path, *, membership: str | Path, forward: str | Path,
             genome_scope: str | Path, seal_path: str | Path | None = None,
             provenance: dict | None = None) -> Generation:
    """Adopt built components for ``db_path`` as a new immutable generation (CURRENT is left unchanged).

    The full source SHA-256 must equal the one recorded by the components; every deployed
    payload is checksum-verified while it is copied into the generation.
    """
    index_dir = Path(index_dir).resolve()
    db_path = Path(db_path).resolve()
    sources = {"membership": Path(membership), "forward": Path(forward), "genome_scope": Path(genome_scope)}
    manifests = {name: _manifest(path, COMPONENTS[name]) for name, path in sources.items()}
    _check_bindings(manifests)
    membership_manifest = manifests["membership"][0]
    if Path(str(db_path) + ".wal").exists():
        raise StaleGenerationError(f"Source has a pending WAL: {db_path}")
    state = source_state(db_path)
    source_sha = sha256_file(db_path)
    if source_sha != membership_manifest["source"]["sha256"]:
        raise StaleGenerationError("Source SHA-256 differs from the one the components were built from")
    if source_state(db_path) != state:
        raise StaleGenerationError("Source changed while it was being hashed")
    seal = None
    if seal_path is not None:
        seal_path = Path(seal_path).resolve()
        seal = {"name": seal_path.name, "sha256": sha256_file(seal_path),
                "dataset_id": json.loads(seal_path.read_text()).get("dataset_id")}
        recorded = membership_manifest["source"].get("seal")
        if recorded and recorded.get("dataset_id") != seal["dataset_id"]:
            raise StaleGenerationError("Seal dataset ID differs from the one recorded at membership build time")
    digest = hashlib.sha256("\n".join(manifests[n][1] for n in COMPONENTS).encode()).hexdigest()
    generation_id = f"g1-{digest[:16]}"
    generations = index_dir / "generations"
    generations.mkdir(parents=True, exist_ok=True)
    final = generations / generation_id
    if final.exists():
        raise FileExistsError(f"Generation already exists: {final}")
    staging = Path(tempfile.mkdtemp(prefix=f".{generation_id}.", dir=generations))
    try:
        components = {}
        for name, src in sources.items():
            manifest, manifest_sha = manifests[name]
            dst = staging / name
            dst.mkdir()
            for fname in manifest["files"]:
                if fname not in CONTROL_FILES:
                    _copy(src / fname, dst / fname)
            _copy(src / "manifest.json", dst / "manifest.json")
            _verify_payloads(dst, manifest)
            if sha256_file(dst / "manifest.json") != manifest_sha:
                raise GenerationError(f"Copied {name} manifest differs from its source")
            components[name] = {"dir": name, "format": COMPONENTS[name], "manifest_sha256": manifest_sha,
                                "payload_bytes": sum(e["bytes"] for f, e in manifest["files"].items()
                                                     if f not in CONTROL_FILES)}
        forward_manifest = manifests["forward"][0]
        scope_manifest = manifests["genome_scope"][0]
        record = {"format": GENERATION_FORMAT, "generation_id": generation_id, "components": components,
                  "source": {"name": db_path.name, "sha256": source_sha, "state": state, "seal": seal},
                  "counts": {"proteins": membership_manifest["protein_count"],
                             "active_proteins": membership_manifest["active_protein_count"],
                             "terms": membership_manifest["term_count"],
                             "memberships": membership_manifest["memberships"],
                             "rich_rows": forward_manifest["rows"],
                             "rich_term_ids": forward_manifest["dictionaries"]["term_id"]["non_null_values"],
                             "genomes": scope_manifest["genome_count"],
                             "null_owned_proteins": scope_manifest["null_ownership_count"]},
                  "adopted_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                  "provenance": provenance or {}}
        _write_json_atomic(staging / "generation.json", record)
        os.chmod(staging, 0o755)
        os.rename(staging, final)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return Generation(final, generation_id, record)


def select(index_dir: str | Path, generation_id: str) -> None:
    """Point CURRENT at an existing generation (atomic replace)."""
    index_dir = Path(index_dir).resolve()
    if not (index_dir / "generations" / generation_id / "generation.json").is_file():
        raise GenerationError(f"No such generation: {generation_id}")
    _write_json_atomic(index_dir / POINTER, {"generation": generation_id})


def resolve(index_dir: str | Path) -> Generation:
    """The generation CURRENT names (or ``index_dir`` itself when it is a generation directory)."""
    root = Path(index_dir).resolve()
    if (root / "generation.json").is_file():
        gen_dir = root
    else:
        pointer = root / POINTER
        if not pointer.is_file():
            raise GenerationError(f"No CURRENT pointer or generation.json in {root}")
        try:
            generation_id = json.loads(pointer.read_text())["generation"]
        except (ValueError, KeyError) as exc:
            raise GenerationError(f"Malformed CURRENT pointer: {pointer}") from exc
        gen_dir = root / "generations" / generation_id
        if not (gen_dir / "generation.json").is_file():
            raise GenerationError(f"CURRENT names a missing generation: {generation_id}")
    record = json.loads((gen_dir / "generation.json").read_text())
    if record.get("format") != GENERATION_FORMAT:
        raise GenerationError(f"Unknown generation format in {gen_dir}")
    if set(record.get("components", {})) != set(COMPONENTS):
        raise GenerationError(f"Generation lacks a component: {gen_dir}")
    return Generation(gen_dir, record["generation_id"], record)


def verify(generation: Generation, db_path: str | Path, *, payloads: bool = True,
           seal_path: str | Path | None = None, content_match: bool = False) -> None:
    """Refuse a generation whose source binding, manifests or (optionally) payloads differ.

    The source must keep the size, mtime and inode it had at adoption. With
    ``content_match`` (for verified copies such as a staged query replica), a file whose
    state differs is accepted when its full SHA-256 equals the adopted one.
    """
    db_path = Path(db_path).resolve()
    source = generation.record["source"]
    if Path(str(db_path) + ".wal").exists():
        raise StaleGenerationError(f"Source has a pending WAL: {db_path}")
    if source_state(db_path) != source["state"] and (
            not content_match or sha256_file(db_path) != source["sha256"]):
        raise StaleGenerationError(
            f"{db_path.name} differs from the state generation {generation.generation_id} was adopted for "
            "(size, mtime or inode changed); rebuild or re-adopt the index")
    if source.get("seal"):
        # Bind the dataset identity; a reseal that only refreshes software provenance keeps it.
        seal = Path(seal_path) if seal_path else db_path.parent / source["seal"]["name"]
        if not seal.is_file() or json.loads(seal.read_text()).get("dataset_id") != source["seal"]["dataset_id"]:
            raise StaleGenerationError(f"Dataset seal identity differs from the one bound to {generation.generation_id}")
    manifests = {}
    for name, info in generation.record["components"].items():
        manifest, manifest_sha = _manifest(generation.component(name), COMPONENTS[name])
        if manifest_sha != info["manifest_sha256"]:
            raise GenerationError(f"{name} manifest differs from the generation record")
        if manifest["source"]["sha256"] != source["sha256"]:
            raise StaleGenerationError(f"{name} was built from a different source")
        manifests[name] = (manifest, manifest_sha)
    _check_bindings(manifests)
    if payloads:
        for name, (manifest, _) in manifests.items():
            _verify_payloads(generation.component(name), manifest)
