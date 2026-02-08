#!/usr/bin/env python3
"""
Rate-limited Foldseek batch search for ESM3 predicted structures.
Run with: python scripts/foldseek_batch.py

Respects Foldseek API rate limits with 60-second delays between queries.
"""

import requests
import time
import tarfile
import json
import io
from pathlib import Path
import sys

STRUCTURES_DIR = Path("data/thorarchaeota_production/structures")
FOLDSEEK_DIR = Path("data/thorarchaeota_production/foldseek_results")
API = "https://search.foldseek.com/api"
DELAY_SECONDS = 60  # Delay between queries to avoid rate limit

def search_structure(pdb_path: Path) -> dict | None:
    """Submit PDB to Foldseek and return parsed results."""
    pid = pdb_path.stem
    result_file = FOLDSEEK_DIR / f"{pid}.json"

    if result_file.exists():
        print(f"[{pid}] SKIP - already done")
        return None

    print(f"[{pid}] Submitting...", flush=True)

    try:
        # Submit
        with open(pdb_path, 'rb') as f:
            resp = requests.post(
                f"{API}/ticket",
                files={'q': (pdb_path.name, f)},
                data={'database[]': ['afdb50', 'pdb100'], 'mode': '3diaa'},
                timeout=60
            )

        if resp.status_code == 429:
            print(f"[{pid}] RATE LIMITED - waiting 5 minutes")
            time.sleep(300)
            return None

        if resp.status_code != 200:
            print(f"[{pid}] Submit failed: {resp.status_code}")
            return None

        ticket_id = resp.json().get('id')
        print(f"[{pid}] Ticket: {ticket_id}")

        # Poll for completion
        for attempt in range(90):  # 7.5 min max
            time.sleep(5)
            r = requests.get(f"{API}/ticket/{ticket_id}", timeout=30)
            if r.status_code == 200:
                status = r.json().get('status')
                if status == 'COMPLETE':
                    print(f"[{pid}] Complete after {(attempt+1)*5}s")
                    break
                elif status == 'ERROR':
                    print(f"[{pid}] Job error")
                    return None
        else:
            print(f"[{pid}] Timeout")
            return None

        # Download tar.gz
        r = requests.get(f"{API}/result/download/{ticket_id}", timeout=120)
        content = r.content

        if len(content) < 100:
            print(f"[{pid}] Empty result")
            return None

        # Parse tar.gz
        results = {'query': pid, 'hits': []}
        tar_buffer = io.BytesIO(content)

        with tarfile.open(fileobj=tar_buffer, mode='r:gz') as tar:
            for member in tar.getmembers():
                if 'alis_' in member.name:
                    f = tar.extractfile(member)
                    if f:
                        data = f.read().decode('utf-8')
                        for line in data.strip().split('\n'):
                            if line:
                                fields = line.split('\t')
                                if len(fields) >= 11:
                                    results['hits'].append({
                                        'database': member.name.split('_')[1] if '_' in member.name else 'unknown',
                                        'target': fields[1],
                                        'identity': float(fields[2]) if fields[2] else 0,
                                        'evalue': float(fields[10]) if fields[10] else 1,
                                        'qcov': float(fields[3]) if fields[3] else 0,
                                        'tcov': float(fields[4]) if fields[4] else 0,
                                    })

        # Save
        with open(result_file, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"[{pid}] Saved {len(results['hits'])} hits")
        return results

    except Exception as e:
        print(f"[{pid}] Error: {type(e).__name__}: {e}")
        return None


def main():
    FOLDSEEK_DIR.mkdir(parents=True, exist_ok=True)

    # Get all PDBs sorted by name
    pdbs = sorted(STRUCTURES_DIR.glob("*.pdb"))
    print(f"Found {len(pdbs)} structures to search")
    print(f"Results will be saved to {FOLDSEEK_DIR}/")
    print(f"Using {DELAY_SECONDS}s delay between queries\n")

    for i, pdb in enumerate(pdbs, 1):
        print(f"\n=== [{i}/{len(pdbs)}] ===")
        result = search_structure(pdb)

        if result is not None:
            # Successful query - wait before next
            print(f"Waiting {DELAY_SECONDS}s before next query...")
            time.sleep(DELAY_SECONDS)

    # Summary
    results = list(FOLDSEEK_DIR.glob("*.json"))
    print(f"\n=== COMPLETE ===")
    print(f"Searched: {len(results)} structures")

    # Show top hits
    print("\nTop hits by structure:")
    for rf in sorted(results)[:10]:
        with open(rf) as f:
            data = json.load(f)
        hits = data.get('hits', [])
        if hits:
            best = min(hits, key=lambda x: x.get('evalue', 1))
            print(f"  {rf.stem}: {len(hits)} hits, best={best['target'][:40]}")


if __name__ == '__main__':
    main()
