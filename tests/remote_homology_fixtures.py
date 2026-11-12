"""Synthetic HH-suite output for tests/test_remote_homology.py (standard library only).

The stand-in hhsearch/hhblits scripts of those tests import this module, so it imports
nothing from Sharur.
"""

from __future__ import annotations


def render_hhr(query: str, match_columns: int, hits: list[dict], *, alignments: int | None = None) -> str:
    """HH-suite 3 result text for ``hits`` (dicts with name, prob, evalue, score, cols, q, t, tlen)."""
    lines = [f"Query         {query}", f"Match_columns {match_columns}", "No_of_seqs    5 out of 9",
             "Neff          2.1", "Searched_HMMs 42", "Date          Thu Jan  1 00:00:00 2026",
             "Command       hhsearch -i x -d y", "",
             " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM"]
    for i, h in enumerate(hits, 1):
        qs, qe = h["q"]
        ts, te = h["t"]
        lines.append(f"{i:3d} {h['name'][:30]:<30} {h['prob']:5.1f} {h['evalue']:7.2G} {h['evalue'] / 1e5:7.2G} "
                     f"{h['score']:6.1f} {0.0:5.1f} {h['cols']:4d} {qs:4d}-{qe:<4d} {ts:4d}-{te:<4d}({h['tlen']:4d})")
    lines.append("")
    for i, h in enumerate(hits[:alignments] if alignments is not None else hits, 1):
        lines += [f"No {i}", f">{h['name']}",
                  f"Probab={h['prob']:.2f}  E-value={h['evalue']:.2g}  Score={h['score']:.2f}  "
                  f"Aligned_cols={h['cols']}  Identities={h.get('ident', 25)}%  Similarity=0.300  "
                  f"Sum_probs={h['cols'] * 0.8:.1f}  Template_Neff=3.000", "",
                  f"Q {query:14s} {h['q'][0]:>4d} ACDEFGHIKL{'':10s} {h['q'][1]:>4d} ({match_columns})",
                  f"T {h['name'].split()[0]:14s} {h['t'][0]:>4d} ACDEFGHIKL{'':10s} {h['t'][1]:>4d} ({h['tlen']})", ""]
    lines.append("Done!")
    return "\n".join(lines) + "\n"


PFAM_HITS = [
    {"name": "PF90001.3 SynthA; synthetic family A", "prob": 99.6, "evalue": 1e-25, "score": 120.0, "cols": 100,
     "q": (5, 110), "t": (1, 105), "tlen": 108},
    {"name": "PF90003.1 DUF9999; Domain of unknown function (DUF9999)", "prob": 95.2, "evalue": 1e-6,
     "score": 50.0, "cols": 40, "q": (60, 100), "t": (3, 44), "tlen": 50},
    {"name": "PF90004.2 WeakDom; weakly similar domain", "prob": 61.0, "evalue": 2.0, "score": 20.0, "cols": 25,
     "q": (10, 35), "t": (1, 25), "tlen": 1162},
    {"name": "PF90005.1 NoiseDom; noise", "prob": 30.0, "evalue": 50.0, "score": 10.0, "cols": 12,
     "q": (1, 12), "t": (1, 12), "tlen": 80},
]
PDB_HITS = [
    {"name": "9zzz_A mol:protein length:205  SYNTHETIC ENZYME", "prob": 98.1, "evalue": 1e-12, "score": 90.0,
     "cols": 95, "q": (3, 104), "t": (10, 110), "tlen": 205},
    {"name": "9yyy_B", "prob": 91.0, "evalue": 1e-4, "score": 40.0, "cols": 31, "q": (20, 52), "t": (1, 33),
     "tlen": 60},
    {"name": "9xxx_C 1-233", "prob": 88.0, "evalue": 1e-3, "score": 35.0, "cols": 50, "q": (20, 72), "t": (1, 50),
     "tlen": 233},
]
UNIREF_HITS = [
    {"name": f"UniRef100_SYN{i:03d} Uncharacterized protein n=1 Tax=Synthetica TaxID=1 RepID=SYN{i:03d}",
     "prob": 100.0 - i, "evalue": 1e-30, "score": 150.0, "cols": 110, "q": (1, 110), "t": (1, 110), "tlen": 115}
    for i in range(30)
]
TITLES = "9yyy\tCrystal structure of a protein of unknown function\n9xxx\tStructure of a synthetic transporter\n"
