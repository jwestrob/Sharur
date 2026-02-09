#!/usr/bin/env python3
"""Visualize ESM2 protein embeddings using UMAP dimensionality reduction.

Usage:
    python scripts/visualize_embeddings.py \
        --db data/DATASET/sharur.duckdb \
        --output figures/embedding_umap.html \
        --color-by genome \
        --limit 10000

    python scripts/visualize_embeddings.py \
        --db data/DATASET/sharur.duckdb \
        --output figures/embedding_umap.png \
        --color-by predicate \
        --predicate hydrogenase \
        --limit 5000

Requirements:
    pip install umap-learn plotly  # for interactive HTML
    pip install umap-learn matplotlib  # for static PNG
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def load_embeddings(db_path: str, limit: int = 10000) -> tuple[np.ndarray, list[str]]:
    """Load ESM2 embeddings from LanceDB alongside the Sharur database.

    Returns (embeddings_array, protein_ids).
    """
    import lancedb

    db_dir = Path(db_path).parent
    lance_path = db_dir / "embeddings"

    if not lance_path.exists():
        # Try alternative locations
        for alt in ["lancedb", "embeddings.lance"]:
            alt_path = db_dir / alt
            if alt_path.exists():
                lance_path = alt_path
                break
        else:
            print(f"Error: No LanceDB embeddings found near {db_path}")
            print(f"  Checked: {db_dir / 'embeddings'}, {db_dir / 'lancedb'}")
            sys.exit(1)

    db = lancedb.connect(str(lance_path))
    table = db.open_table("embeddings")

    # Read embeddings
    df = table.to_pandas()
    if limit and len(df) > limit:
        df = df.sample(n=limit, random_state=42)

    protein_ids = df["protein_id"].tolist()
    vectors = np.stack(df["vector"].values)

    print(f"Loaded {len(protein_ids)} embeddings ({vectors.shape[1]}-dim)")
    return vectors, protein_ids


def get_color_data(
    db_path: str,
    protein_ids: list[str],
    color_by: str,
    predicate: str | None = None,
    annotation_source: str | None = None,
) -> tuple[list[str], str]:
    """Get color labels for each protein.

    Returns (labels, legend_title).
    """
    import duckdb

    conn = duckdb.connect(db_path, read_only=True)

    if color_by == "genome":
        # Color by bin_id
        placeholders = ",".join(["?"] * len(protein_ids))
        rows = conn.execute(
            f"SELECT protein_id, bin_id FROM proteins WHERE protein_id IN ({placeholders})",
            protein_ids,
        ).fetchall()
        pid_to_label = {r[0]: r[1] or "unknown" for r in rows}
        labels = [pid_to_label.get(pid, "unknown") for pid in protein_ids]
        return labels, "Genome"

    elif color_by == "predicate":
        if not predicate:
            # Use annotation status
            predicate = "unannotated"

        placeholders = ",".join(["?"] * len(protein_ids))
        rows = conn.execute(
            f"""
            SELECT protein_id, list_contains(predicates, ?) as has_pred
            FROM protein_predicates
            WHERE protein_id IN ({placeholders})
            """,
            [predicate] + protein_ids,
        ).fetchall()
        pid_to_label = {r[0]: predicate if r[1] else "other" for r in rows}
        labels = [pid_to_label.get(pid, "other") for pid in protein_ids]
        return labels, f"Predicate: {predicate}"

    elif color_by == "annotation":
        source = annotation_source or "pfam"
        placeholders = ",".join(["?"] * len(protein_ids))
        rows = conn.execute(
            f"""
            SELECT DISTINCT a.protein_id, a.name
            FROM annotations a
            WHERE a.protein_id IN ({placeholders})
              AND a.source = ?
            """,
            protein_ids + [source],
        ).fetchall()

        # Use the first annotation per protein
        pid_to_label: dict[str, str] = {}
        for pid, name in rows:
            if pid not in pid_to_label:
                pid_to_label[pid] = name or "unknown"

        labels = [pid_to_label.get(pid, "unannotated") for pid in protein_ids]
        return labels, f"Annotation ({source})"

    else:
        raise ValueError(f"Unknown color-by: {color_by}")


def run_umap(vectors: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    """Run UMAP dimensionality reduction."""
    import umap

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        random_state=42,
        metric="cosine",
    )
    embedding_2d = reducer.fit_transform(vectors)
    print(f"UMAP complete: {embedding_2d.shape}")
    return embedding_2d


def plot_interactive(
    embedding_2d: np.ndarray,
    labels: list[str],
    protein_ids: list[str],
    legend_title: str,
    output_path: str,
) -> None:
    """Create interactive HTML plot with plotly."""
    import plotly.express as px
    import pandas as pd

    df = pd.DataFrame({
        "UMAP1": embedding_2d[:, 0],
        "UMAP2": embedding_2d[:, 1],
        "label": labels,
        "protein_id": protein_ids,
    })

    # Limit legend to top 20 categories
    top_labels = df["label"].value_counts().head(20).index.tolist()
    df["label_display"] = df["label"].apply(
        lambda x: x if x in top_labels else "other"
    )

    fig = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color="label_display",
        hover_data=["protein_id", "label"],
        title=f"ESM2 Embedding UMAP ({len(df)} proteins)",
        labels={"label_display": legend_title},
        opacity=0.6,
        width=1000,
        height=800,
    )
    fig.update_traces(marker=dict(size=3))
    fig.write_html(output_path)
    print(f"Interactive plot saved: {output_path}")


def plot_static(
    embedding_2d: np.ndarray,
    labels: list[str],
    legend_title: str,
    output_path: str,
) -> None:
    """Create static PNG plot with matplotlib."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    unique_labels = sorted(set(labels))

    # Limit to top 15 categories for readability
    from collections import Counter

    label_counts = Counter(labels)
    top_labels = [l for l, _ in label_counts.most_common(15)]

    fig, ax = plt.subplots(figsize=(12, 10))

    for label in top_labels:
        mask = np.array([l == label for l in labels])
        ax.scatter(
            embedding_2d[mask, 0],
            embedding_2d[mask, 1],
            s=2,
            alpha=0.5,
            label=f"{label} ({mask.sum()})",
        )

    # Plot "other" category
    other_mask = np.array([l not in top_labels for l in labels])
    if other_mask.any():
        ax.scatter(
            embedding_2d[other_mask, 0],
            embedding_2d[other_mask, 1],
            s=1,
            alpha=0.2,
            color="gray",
            label=f"other ({other_mask.sum()})",
        )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"ESM2 Embedding UMAP ({len(labels)} proteins)")
    ax.legend(
        title=legend_title,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=7,
        markerscale=4,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Static plot saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize ESM2 protein embeddings using UMAP"
    )
    parser.add_argument("--db", required=True, help="Path to sharur.duckdb")
    parser.add_argument("--output", required=True, help="Output path (.html or .png)")
    parser.add_argument(
        "--color-by",
        choices=["genome", "predicate", "annotation"],
        default="genome",
        help="How to color points (default: genome)",
    )
    parser.add_argument(
        "--predicate",
        help="Predicate to highlight (with --color-by predicate)",
    )
    parser.add_argument(
        "--annotation-source",
        help="Annotation source to color by (with --color-by annotation, default: pfam)",
    )
    parser.add_argument(
        "--limit", type=int, default=10000, help="Max proteins to plot (default: 10000)"
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="UMAP n_neighbors parameter (default: 15)",
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.1,
        help="UMAP min_dist parameter (default: 0.1)",
    )

    args = parser.parse_args()

    # Load embeddings
    vectors, protein_ids = load_embeddings(args.db, limit=args.limit)

    # Get color labels
    labels, legend_title = get_color_data(
        args.db,
        protein_ids,
        args.color_by,
        predicate=args.predicate,
        annotation_source=args.annotation_source,
    )

    # Run UMAP
    embedding_2d = run_umap(vectors, n_neighbors=args.n_neighbors, min_dist=args.min_dist)

    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Plot
    if output_path.suffix == ".html":
        plot_interactive(embedding_2d, labels, protein_ids, legend_title, str(output_path))
    else:
        plot_static(embedding_2d, labels, legend_title, str(output_path))


if __name__ == "__main__":
    main()
