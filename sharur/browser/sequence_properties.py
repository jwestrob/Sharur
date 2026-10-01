"""Sequence properties computed from a protein's own sequence, for the protein page.

All of these are rules applied to the sequence, stated as such on the page:

- ``hydropathy``: Kyte-Doolittle values averaged over a 19-residue window.
- ``hydrophobic_segments``: runs of windows averaging at least 1.6 (the
  Kyte-Doolittle membrane-spanning rule); each segment covers its windows.
- ``low_complexity``: SEG-like regions — 12-residue windows with Shannon
  entropy at most 2.2 bits seed a region, which extends while windows stay at
  or below 2.5 bits.
- ``n_terminal``: a hydrophobic stretch (8 residues averaging at least 2.0)
  within the first 35 residues, after a positive residue in the first 10 —
  the shape of a signal peptide or N-terminal anchor; a heuristic.
- ``composition``: length, mass, net charge (K + R − D − E), cysteines, and the
  residues most enriched over the average Swiss-Prot composition.
- ``self_similarity``: k-mer self-matches in a 10-letter reduced alphabet
  (k from 4 to 8, the shortest keeping chance matches at or below 0.05 per
  cell given the sequence's composition) on a grid of at most 600 × 600 cells,
  keeping cells with a neighbour along the diagonal (internal repeats show as
  off-diagonal lines), and the strongest offset between matching k-mers as a
  repeat period estimate.
"""

from __future__ import annotations

import base64
import math
import struct
import zlib
from typing import Any

import numpy as np

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# average residue masses (Da), for a protein mass estimate
MASS = {"A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "Q": 128.13, "E": 129.12, "G": 57.05,
        "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08,
        "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13}
# UniProtKB/Swiss-Prot average amino-acid composition (%)
SWISSPROT = {"A": 8.25, "R": 5.53, "N": 4.06, "D": 5.45, "C": 1.37, "Q": 3.93, "E": 6.75, "G": 7.07, "H": 2.27,
             "I": 5.96, "L": 9.66, "K": 5.84, "M": 2.42, "F": 3.86, "P": 4.70, "S": 6.56, "T": 5.34, "W": 1.08,
             "Y": 2.92, "V": 6.87}
AMINO = "ACDEFGHIKLMNPQRSTVWY"

TM_WINDOW, TM_THRESHOLD = 19, 1.6
LC_WINDOW, LC_TRIGGER, LC_EXTEND = 12, 2.2, 2.5
GRID = 600
KMER_CAP = 400          # k-mers occurring more often than this (low complexity) are left out of the dot plot
OFFSET_SAMPLE = 120     # positions per k-mer used for the offset histogram


def clean(sequence: str) -> str:
    return "".join(c for c in (sequence or "").upper() if c.isalpha() and c != "*")


def window_mean(values: np.ndarray, window: int) -> np.ndarray:
    """Mean of every full window (len(values) - window + 1 entries)."""
    if len(values) < window:
        return np.zeros(0)
    c = np.concatenate([[0.0], np.cumsum(values)])
    return (c[window:] - c[:-window]) / window


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """(start, end) index pairs (end exclusive) of True runs."""
    if not len(mask):
        return []
    edges = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return list(zip(np.flatnonzero(edges == 1).tolist(), np.flatnonzero(edges == -1).tolist()))


def hydropathy(seq: str) -> np.ndarray:
    values = np.array([KD.get(c, 0.0) for c in seq])
    return window_mean(values, TM_WINDOW)


def hydrophobic_segments(profile: np.ndarray) -> list[dict[str, Any]]:
    """Runs of 19-residue windows averaging >= 1.6, as 1-based residue spans."""
    out = []
    for s, e in runs(profile >= TM_THRESHOLD):
        out.append({"start": s + 1, "end": e - 1 + TM_WINDOW, "peak": round(float(profile[s:e].max()), 2)})
    return out


def entropy_windows(seq: str, window: int = LC_WINDOW) -> np.ndarray:
    n = len(seq)
    if n < window:
        return np.zeros(0)
    codes = np.array([AMINO.find(c) for c in seq])
    onehot = np.zeros((n, 21), dtype=np.int32)
    onehot[np.arange(n), np.where(codes < 0, 20, codes)] = 1
    cum = np.vstack([np.zeros((1, 21), dtype=np.int32), np.cumsum(onehot, axis=0)])
    counts = (cum[window:] - cum[:-window]).astype(float)
    p = counts / window
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.nansum(np.where(p > 0, p * np.log2(p), 0.0), axis=1)
    return h


def low_complexity(seq: str) -> list[dict[str, Any]]:
    """SEG-like regions: windows at <= 2.2 bits seed, extended while <= 2.5 bits."""
    h = entropy_windows(seq)
    if not len(h):
        return []
    out = []
    for s, e in runs(h <= LC_EXTEND):
        if (h[s:e] <= LC_TRIGGER).any():
            start, end = s + 1, e - 1 + LC_WINDOW
            region = seq[start - 1:end]
            top = max(set(region), key=region.count)
            out.append({"start": start, "end": end, "top": top, "top_share": round(region.count(top) / len(region), 2)})
    return out


def n_terminal(seq: str) -> dict[str, Any] | None:
    """A signal-peptide-like N-terminus: K/R in the first 10, then 8 hydrophobic residues within the first 35."""
    head = seq[:35]
    if len(head) < 15 or not any(c in "KR" for c in head[:10]):
        return None
    values = np.array([KD.get(c, 0.0) for c in head])
    means = window_mean(values, 8)
    best = int(np.argmax(means)) if len(means) else 0
    if len(means) and means[best] >= 2.0 and best >= 1:
        return {"start": best + 1, "end": best + 8, "mean": round(float(means[best]), 2)}
    return None


def composition(seq: str) -> dict[str, Any]:
    n = len(seq) or 1
    counts = {a: seq.count(a) for a in AMINO}
    share = {a: 100.0 * counts[a] / n for a in AMINO}
    enriched = sorted(((a, share[a] / SWISSPROT[a]) for a in AMINO if counts[a]), key=lambda x: -x[1])[:4]
    return {
        "length": len(seq),
        "mass_kda": round(sum(MASS.get(c, 110.0) for c in seq) / 1000 + 0.018, 1),
        "net_charge": counts["K"] + counts["R"] - counts["D"] - counts["E"],
        "charged": round(100.0 * sum(counts[a] for a in "DEKR") / n, 1),
        "hydrophobic": round(100.0 * sum(counts[a] for a in "AILMFVW") / n, 1),
        "cysteines": counts["C"],
        "enriched": [{"aa": a, "share": round(share[a], 1), "ratio": round(r, 1)} for a, r in enriched if r >= 1.5],
    }


# Murphy et al. (2000) 10-letter reduced alphabet: diverged repeat copies keep sharing words
REDUCED = {a: i for i, group in enumerate(("LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H")) for a in group}


def _kmer_codes(seq: str, k: int) -> np.ndarray:
    codes = np.array([REDUCED.get(c, -1) for c in seq], dtype=np.int64)
    valid = codes >= 0
    codes = np.where(valid, codes, 0)
    n = len(seq) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.int64)
    out = np.zeros(n, dtype=np.int64)
    ok = np.ones(n, dtype=bool)
    for j in range(k):
        out = out * 10 + codes[j:j + n]
        ok &= valid[j:j + n]
    return np.where(ok, out, -1)


def diagonal_runs(cells: np.ndarray) -> np.ndarray:
    """Cells with a non-empty neighbour along the diagonal; isolated chance matches drop out."""
    filled = cells > 0
    keep = np.zeros_like(filled)
    keep[1:, 1:] |= filled[:-1, :-1]
    keep[:-1, :-1] |= filled[1:, 1:]
    return np.where(filled & keep, cells, 0.0)


CHANCE_PER_CELL = 0.05


def word_length(seq: str, bins: int) -> int:
    """The shortest k (4-8) keeping chance matches per plot cell at or below CHANCE_PER_CELL.

    q is the chance two residues match in the reduced alphabet, from this
    sequence's own composition; a cell spans (n / bins)² residue pairs.
    """
    groups = [REDUCED[c] for c in seq if c in REDUCED]
    if not groups:
        return 4
    freq = np.bincount(groups, minlength=10) / len(groups)
    q = float((freq ** 2).sum())
    pairs = max(1.0, (len(seq) / bins) ** 2)
    k = math.ceil(math.log(CHANCE_PER_CELL / pairs) / math.log(q)) if 0 < q < 1 else 4
    return int(min(max(k, 4), 8))


def self_similarity(seq: str, *, grid: int = GRID) -> dict[str, Any] | None:
    """Binned k-mer self-match grid (reduced alphabet) and the strongest repeat offset."""
    n = len(seq)
    if n < 30:
        return None
    k = word_length(seq, min(grid, n))
    codes = _kmer_codes(seq, k)
    pos = np.flatnonzero(codes >= 0)
    codes = codes[pos]
    order = np.argsort(codes, kind="stable")
    codes, pos = codes[order], pos[order]
    bounds = np.flatnonzero(np.diff(codes)) + 1
    starts = np.concatenate([[0], bounds])
    ends = np.concatenate([bounds, [len(codes)]])
    sizes = ends - starts
    bins = min(grid, n)
    scale = bins / n
    cells = np.zeros((bins, bins), dtype=np.float64)
    offsets = np.zeros(n, dtype=np.float64)
    skipped = 0
    for s, e in zip(starts[sizes >= 2], ends[sizes >= 2]):
        p = pos[s:e]
        if len(p) > KMER_CAP:
            skipped += 1
            continue
        b, c = np.unique((p * scale).astype(np.int64), return_counts=True)
        block = np.outer(c, c).astype(np.float64)
        np.fill_diagonal(block, np.diag(block) - c)        # drop each position's match with itself
        cells[np.ix_(b, b)] += block
        sample = p if len(p) <= OFFSET_SAMPLE else p[np.linspace(0, len(p) - 1, OFFSET_SAMPLE).astype(int)]
        d = (sample[None, :] - sample[:, None])
        d = d[d > 0]
        np.add.at(offsets, d, len(p) / len(sample) if len(p) > len(sample) else 1.0)
    period = _period(offsets, n, k, pos, codes, starts, ends, sizes)
    return {"k": k, "bins": bins, "cells": diagonal_runs(cells), "skipped_kmers": skipped, "period": period}


def _period(offsets: np.ndarray, n: int, k: int, pos, codes, starts, ends, sizes) -> dict[str, Any] | None:
    """The smallest strong offset between matching k-mers, and the stretch where it repeats.

    An offset qualifies when its matches (smoothed over ±1) reach 25 and eight
    times the mean over all offsets; the smallest offset with at least half the
    strongest signal is the repeat unit. The repeat stretch is the densest
    cluster of positions whose k-mer recurs at that offset (± 2).
    """
    hi = n // 2
    if hi < 8:
        return None
    o = offsets[:hi + 1].copy()
    o[:max(k, 4)] = 0                                 # overlapping self-matches
    smooth = np.convolve(o, np.ones(3), mode="same")
    # mean over every offset, zeros included: random exact k-mer matches at one offset are rare
    background = float(smooth[max(k, 4):].mean())
    candidates = np.flatnonzero(smooth >= max(25.0, 8 * background))
    if not len(candidates):
        return None
    strongest = float(smooth.max())
    d = int(candidates[smooth[candidates] >= 0.5 * strongest][0])
    d = int(np.argmax(o[max(0, d - 2):d + 3]) + max(0, d - 2))     # the raw peak; smoothing ties neighbours
    hits = []
    for s, e in zip(starts[sizes >= 2], ends[sizes >= 2]):
        p = pos[s:e]
        if len(p) > KMER_CAP:
            continue
        diff = p[None, :] - p[:, None]
        rows = np.flatnonzero(((diff >= d - 2) & (diff <= d + 2)).any(axis=1))
        hits.extend(p[rows].tolist())
    if len(hits) < 6:
        return None
    hits = np.sort(np.array(hits))
    gap = max(8 * d, 200)              # diverged copies share exact k-mers only now and then
    breaks = np.flatnonzero(np.diff(hits) > gap) + 1
    clusters = np.split(hits, breaks)
    best = max(clusters, key=len)
    if len(best) < 6:
        return None
    lo, hi_pos = int(best[0]), min(n, int(best[-1]) + d + k)
    return {"period": d, "start": lo + 1, "end": hi_pos, "copies": round((hi_pos - lo) / d, 1),
            "matches": int(smooth[d]), "background": round(background, 1), "hits": int(len(best))}


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #

def _png_rgba(alpha: np.ndarray, rgb: tuple[int, int, int]) -> bytes:
    """A minimal RGBA PNG: one colour, per-pixel alpha (0-255)."""
    h, w = alpha.shape
    raw = bytearray()
    r, g, b = rgb
    for y in range(h):
        raw.append(0)
        row = np.empty((w, 4), dtype=np.uint8)
        row[:, 0], row[:, 1], row[:, 2], row[:, 3] = r, g, b, alpha[y]
        raw.extend(row.tobytes())

    def chunk(tag: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(bytes(raw), 6)) + chunk(b"IEND", b""))


def dotplot_png(cells: np.ndarray, rgb: tuple[int, int, int] = (47, 168, 148)) -> str:
    """Base64 PNG of the self-match grid: log-scaled opacity, symmetric."""
    grid = cells + cells.T
    top = grid.max()
    if top <= 0:
        alpha = np.zeros(grid.shape, dtype=np.uint8)
    else:
        alpha = (255 * np.sqrt(np.log1p(grid) / math.log1p(top))).clip(0, 255).astype(np.uint8)
    return base64.b64encode(_png_rgba(alpha, rgb)).decode()


def analyse(sequence: str) -> dict[str, Any] | None:
    seq = clean(sequence)
    if not seq:
        return None
    profile = hydropathy(seq)
    return {
        "length": len(seq),
        "profile": profile,
        "segments": hydrophobic_segments(profile),
        "low_complexity": low_complexity(seq),
        "n_terminal": n_terminal(seq),
        "composition": composition(seq),
        "self": self_similarity(seq),
    }
