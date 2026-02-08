"""
Foldseek structure search operators.

Searches predicted structures against AlphaFold DB, PDB, and ESM Atlas
using local Foldseek binary (preferred) or the web API (fallback).
"""

from __future__ import annotations

import csv
import io
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import requests

from bennu.operators.base import BennuResult, OperatorContext

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


# Foldseek API endpoints
FOLDSEEK_BASE = "https://search.foldseek.com/api"
FOLDSEEK_DATABASES = f"{FOLDSEEK_BASE}/databases"
FOLDSEEK_SUBMIT = f"{FOLDSEEK_BASE}/ticket"

# Reliable databases (some like esm30_folddisco can be flaky)
DEFAULT_DATABASES = ["afdb50", "afdb-swissprot", "pdb100"]

# Local Foldseek database search path
LOCAL_DB_DIR = Path.home() / ".foldseek"

# Output format for local foldseek easy-search
LOCAL_OUTPUT_COLUMNS = (
    "query,target,evalue,bits,prob,qstart,qend,tstart,tend,taxname,theader"
)


def _find_foldseek_binary() -> Optional[str]:
    """Find the local foldseek binary. Returns path or None."""
    return shutil.which("foldseek")


def _find_local_database(db_name: str) -> Optional[Path]:
    """Check if a local Foldseek database exists. Returns path or None."""
    db_path = LOCAL_DB_DIR / db_name / db_name
    if db_path.exists():
        return db_path
    # Also check without the nested directory
    alt_path = LOCAL_DB_DIR / db_name
    if alt_path.exists() and not alt_path.is_dir():
        return alt_path
    return None


def search_foldseek_local(
    pdb_path: str,
    database: str = "pdb100",
    top_k: int = 10,
) -> Optional[list[dict]]:
    """Run local foldseek easy-search.

    Returns list of hit dicts, or None if local foldseek is not available
    (caller should fall back to web API).
    """
    binary = _find_foldseek_binary()
    if binary is None:
        return None

    db_path = _find_local_database(database)
    if db_path is None:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        result_path = Path(tmpdir) / "result.tsv"
        tmp_path = Path(tmpdir) / "tmp"
        tmp_path.mkdir()

        cmd = [
            binary,
            "easy-search",
            pdb_path,
            str(db_path),
            str(result_path),
            str(tmp_path),
            "--format-output",
            LOCAL_OUTPUT_COLUMNS,
        ]

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return None

        if proc.returncode != 0:
            return None

        if not result_path.exists():
            return []

        content = result_path.read_text().strip()
        if not content:
            return []

        hits = []
        reader = csv.reader(io.StringIO(content), delimiter="\t")
        for row in reader:
            if len(row) < 11:
                continue
            try:
                evalue = float(row[2]) if row[2] else None
            except ValueError:
                evalue = None
            try:
                score = float(row[3]) if row[3] else None
            except ValueError:
                score = None

            hits.append({
                "database": database,
                "target": row[1],
                "evalue": evalue,
                "score": score,
                "seq_identity": None,  # Not in this output format
                "description": row[10] if len(row) > 10 else "",
                "taxon": row[9] if len(row) > 9 else "",
                "qstart": int(row[5]) if row[5].isdigit() else None,
                "qend": int(row[6]) if row[6].isdigit() else None,
                "tstart": int(row[7]) if row[7].isdigit() else None,
                "tend": int(row[8]) if row[8].isdigit() else None,
            })

        hits.sort(
            key=lambda x: x["evalue"] if x["evalue"] is not None else float("inf")
        )
        return hits[:top_k]


def _ensure_pdb_format(pdb_content: str) -> str:
    """Ensure PDB has TER and END lines (required by Foldseek)."""
    pdb_content = pdb_content.rstrip()
    if "TER" not in pdb_content:
        pdb_content += "\nTER"
    if not pdb_content.endswith("END"):
        pdb_content += "\nEND"
    return pdb_content + "\n"


def _submit_foldseek_job(pdb_path: str, databases: list[str]) -> str:
    """Submit a Foldseek search job and return ticket ID."""
    # Read and fix PDB
    pdb_content = Path(pdb_path).read_text()
    pdb_content = _ensure_pdb_format(pdb_content)

    # Submit job
    files = {"q": ("protein.pdb", pdb_content, "chemical/x-pdb")}
    data = {"mode": "3diaa", "database[]": databases}

    response = requests.post(FOLDSEEK_SUBMIT, files=files, data=data, timeout=30)

    if response.status_code != 200:
        raise Exception(f"Foldseek submit failed: {response.text}")

    result = response.json()
    return result["id"]


def _poll_foldseek_job(ticket_id: str, max_wait: int = 300) -> dict:
    """Poll for job completion and return results."""
    status_url = f"{FOLDSEEK_BASE}/ticket/{ticket_id}"
    result_url = f"{FOLDSEEK_BASE}/result/{ticket_id}/0"  # Note: /0 suffix required!

    start_time = time.time()
    while time.time() - start_time < max_wait:
        status_response = requests.get(status_url, timeout=10)
        status = status_response.json().get("status")

        if status == "COMPLETE":
            # Get results
            results_response = requests.get(result_url, timeout=30)
            return results_response.json()
        elif status == "ERROR":
            raise Exception("Foldseek server error")

        time.sleep(3)  # Poll every 3 seconds

    raise Exception(f"Foldseek job timed out after {max_wait}s")


def _parse_foldseek_results(data: dict, top_k: int = 10) -> list[dict]:
    """Parse Foldseek results into a clean list of hits."""
    hits = []

    for db_result in data.get("results", []):
        db_name = db_result.get("db", "unknown")
        alignments = db_result.get("alignments", [])

        if alignments and alignments[0]:
            for hit in alignments[0][:top_k]:
                # seqId from Foldseek is already a fraction (0-1), not percentage
                seq_id = hit.get("seqId")
                if seq_id is not None and seq_id > 1:
                    seq_id = seq_id / 100  # Convert if it's a percentage

                hits.append({
                    "database": db_name,
                    "target": hit.get("target", ""),
                    "evalue": hit.get("eval"),
                    "score": hit.get("score"),
                    "seq_identity": seq_id,
                    "description": hit.get("tDescription", ""),
                    "taxon": hit.get("taxName", ""),
                    "qstart": hit.get("qStartPos"),
                    "qend": hit.get("qEndPos"),
                    "tstart": hit.get("tStartPos"),
                    "tend": hit.get("tEndPos"),
                })

    # Sort by e-value
    hits.sort(key=lambda x: x["evalue"] if x["evalue"] is not None else float("inf"))
    return hits[:top_k]


def search_foldseek(
    pdb_path: str,
    databases: Optional[list[str]] = None,
    top_k: int = 10,
    prefer_local: bool = True,
) -> BennuResult:
    """
    Search a PDB structure against Foldseek databases.

    Uses local foldseek binary when available (faster, no rate limits).
    Falls back to web API for databases not available locally.

    Args:
        pdb_path: Path to PDB file
        databases: List of databases to search (default: afdb50, afdb-swissprot, pdb100)
        top_k: Number of top hits to return
        prefer_local: Try local foldseek first, fall back to web API (default True)

    Returns:
        BennuResult with structural homologs
    """
    if databases is None:
        databases = DEFAULT_DATABASES

    params = {"pdb_path": pdb_path, "databases": databases, "top_k": top_k}

    with OperatorContext("search_foldseek", params) as ctx:
        if not Path(pdb_path).exists():
            return ctx.make_result(
                data=f"PDB file not found: {pdb_path}",
                rows=0,
            )

        all_hits: list[dict] = []
        web_databases: list[str] = []

        # Try local search first for each database
        if prefer_local:
            for db in databases:
                local_hits = search_foldseek_local(pdb_path, database=db, top_k=top_k)
                if local_hits is not None:
                    all_hits.extend(local_hits)
                else:
                    web_databases.append(db)
        else:
            web_databases = list(databases)

        # Fall back to web API for databases not available locally
        if web_databases:
            try:
                ticket_id = _submit_foldseek_job(pdb_path, web_databases)
                results = _poll_foldseek_job(ticket_id)
                web_hits = _parse_foldseek_results(results, top_k)
                all_hits.extend(web_hits)
            except Exception as e:
                if not all_hits:
                    return ctx.make_result(
                        data=f"Foldseek search failed: {e}",
                        rows=0,
                    )
                # If we have local hits, continue with those

        if not all_hits:
            return ctx.make_result(
                data=[],
                rows=0,
            )

        # Sort combined results by e-value
        all_hits.sort(
            key=lambda x: x["evalue"] if x["evalue"] is not None else float("inf")
        )
        all_hits = all_hits[:top_k]

        return ctx.make_result(
            data=all_hits,
            rows=len(all_hits),
        )


def format_foldseek_hits(hits: list[dict], query_name: str = "query") -> str:
    """
    Format Foldseek hits as a readable markdown table.

    Args:
        hits: List of hit dictionaries from search_foldseek
        query_name: Name of query for header

    Returns:
        Formatted markdown string
    """
    if not hits:
        return "No structural homologs found."

    lines = [
        f"# Foldseek Structure Search",
        f"**Query:** {query_name}",
        f"**Hits found:** {len(hits)}",
        "",
        "| Rank | Database | E-value | Identity | Region | Target |",
        "|------|----------|---------|----------|--------|--------|",
    ]

    for i, hit in enumerate(hits, 1):
        db = hit.get("database", "?")[:12]
        evalue = f"{hit['evalue']:.1e}" if hit.get("evalue") else "N/A"
        identity = f"{hit['seq_identity']:.1%}" if hit.get("seq_identity") else "N/A"
        region = f"{hit.get('qstart', '?')}-{hit.get('qend', '?')}"
        target = (hit.get("target") or "")[:45]
        lines.append(f"| {i} | {db} | {evalue} | {identity} | {region} | {target} |")

    # Add top hit details
    top = hits[0]
    lines.extend([
        "",
        "## Top Hit",
        f"- **Target:** {top.get('target', 'N/A')}",
        f"- **E-value:** {top.get('evalue', 'N/A')}",
        f"- **Identity:** {top.get('seq_identity', 0):.1%}" if top.get('seq_identity') else "",
        f"- **Query region:** {top.get('qstart', '?')}-{top.get('qend', '?')}",
        f"- **Taxon:** {top.get('taxon', 'N/A')}" if top.get('taxon') else "",
    ])

    return "\n".join(line for line in lines if line)  # Filter empty lines


def search_foldseek_for_protein(
    store: "DuckDBStore",
    protein_id: str,
    pdb_path: Optional[str] = None,
    databases: Optional[list[str]] = None,
    top_k: int = 10,
) -> BennuResult:
    """
    Search Foldseek for a protein (requires existing PDB or will predict structure).

    Args:
        store: DuckDB store
        protein_id: Protein ID
        pdb_path: Optional path to existing PDB file
        databases: Foldseek databases to search
        top_k: Number of hits to return

    Returns:
        BennuResult with structural homologs (data is list of hit dicts)
    """
    params = {"protein_id": protein_id, "databases": databases}

    with OperatorContext("search_foldseek_for_protein", params) as ctx:
        # If no PDB provided, check for existing or predict
        if pdb_path is None:
            # Check common locations using protein_id
            safe_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in protein_id[:40])
            possible_paths = [
                f"/tmp/bennu_structures/{safe_name}.pdb",
                f"/tmp/bennu_structure_{safe_name}.pdb",
            ]
            for path in possible_paths:
                if Path(path).exists():
                    pdb_path = path
                    break

        if pdb_path is None or not Path(pdb_path).exists():
            return ctx.make_result(
                data=[],
                rows=0,
            )

        # Run Foldseek search - returns structured data
        result = search_foldseek(pdb_path, databases, top_k)

        # Add protein_id to each hit for traceability
        if result.data and isinstance(result.data, list):
            for hit in result.data:
                hit["query_protein_id"] = protein_id

        return result


def list_foldseek_databases() -> BennuResult:
    """List available Foldseek databases. Returns list of database info dicts."""
    with OperatorContext("list_foldseek_databases", {}) as ctx:
        try:
            response = requests.get(FOLDSEEK_DATABASES, timeout=10)
            databases = response.json()

            # Return structured data
            return ctx.make_result(
                data=databases,  # List of {"name": ..., "description": ...}
                rows=len(databases),
            )

        except Exception as e:
            return ctx.make_result(
                data=[],
                rows=0,
            )


__all__ = [
    "search_foldseek",
    "search_foldseek_local",
    "search_foldseek_for_protein",
    "list_foldseek_databases",
    "format_foldseek_hits",
]
