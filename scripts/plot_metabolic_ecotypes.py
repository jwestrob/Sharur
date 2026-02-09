#!/usr/bin/env python3
"""
Regenerate metabolic ecotypes bar chart for Omnitrophota manuscript.

Phase 4 review corrected numbers (non-overlapping decomposition, n=1,831 genomes).
"""

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

OUTPUT = "data/omni_production/figures/metabolic_ecotypes.png"

# Phase 4 corrected ecotype counts (non-overlapping, sorted descending)
ecotypes = [
    ("Core fermenter\n(G4+ Rnf\u2212 Cox\u2212 WL\u2212)",   764),
    ("Enhanced fermenter\n(G4+ Rnf+ Cox\u2212 WL\u2212)",     325),
    ("No-hydrogenase\n(G4\u2212 Rnf\u2212)",                   300),
    ("Microaerobe\n(G4+ Cox+)",                                 223),
    ("Rnf-without-G4\n(G4\u2212 Rnf+)",                        131),
    ("Other\ncombinations",                                      58),
    ("Autotroph\n(G4+ Rnf+ WL+)",                               30),
]

labels = [e[0] for e in ecotypes]
counts = [e[1] for e in ecotypes]

# Distinct colors matching the spirit of the original figure
colors = [
    "#1E88E5",  # blue  - core fermenter
    "#43A047",  # green - enhanced fermenter
    "#9E9E9E",  # grey  - no-hydrogenase
    "#FB8C00",  # orange - microaerobe
    "#8E24AA",  # purple - Rnf-without-G4
    "#6D4C41",  # brown - other
    "#E53935",  # red   - autotroph
]

fig, ax = plt.subplots(figsize=(12, 7))

bars = ax.bar(range(len(counts)), counts, color=colors, edgecolor="black",
              linewidth=0.6, width=0.72)

# Count labels above each bar
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 12,
            f"{count}", ha="center", va="bottom", fontsize=13, fontweight="bold")

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=10, ha="center")
ax.set_ylabel("Number of Genomes", fontsize=13)
ax.set_title("Metabolic Ecotypes in Omnitrophota (n=1,831 genomes)",
             fontsize=16, fontweight="bold", pad=14)

# Y-axis range with headroom for labels
ax.set_ylim(0, max(counts) * 1.15)
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))

# Legend with abbreviation key
legend_text = (
    "G4 = Group 4 NiFe hydrogenase\n"
    "Rnf = Rnf complex\n"
    "Cox = Cytochrome oxidase\n"
    "WL = Wood\u2013Ljungdahl pathway"
)
props = dict(boxstyle="round,pad=0.5", facecolor="wheat", alpha=0.85,
             edgecolor="gray")
ax.text(0.98, 0.97, legend_text, transform=ax.transAxes, fontsize=9.5,
        verticalalignment="top", horizontalalignment="right", bbox=props)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT, dpi=150, bbox_inches="tight")
print(f"Saved: {OUTPUT}")
