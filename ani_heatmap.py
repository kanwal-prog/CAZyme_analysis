#!/usr/bin/env python3
"""
Draw an ANI heatmap (genome vs genome) as a PNG for a paper or poster.

Usage (from ~/Desktop/CAZyme_analysis):
    python ani_heatmap.py ecoli
    python ani_heatmap.py bfragilis
Reads  <sp>/anvio/ANI/ANIb_percentage_identity.txt
Writes <sp>_ANI_heatmap.png
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sp = sys.argv[1]
m = pd.read_csv(f"{sp}/anvio/ANI/ANIb_percentage_identity.txt", sep="\t", index_col=0) * 100
m = (m + m.T) / 2          # ANIb is directional, so average the two directions
n = len(m)

off = m.values[~np.eye(n, dtype=bool)]
vmin = np.floor(off.min()) - 0.5

fig, ax = plt.subplots(figsize=(1.2 * n + 3, 1.2 * n + 2))
im = ax.imshow(m.values, cmap="viridis", vmin=vmin, vmax=100)
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(m.columns, rotation=45, ha="right")
ax.set_yticklabels(m.index)
for i in range(n):
    for j in range(n):
        v = m.values[i, j]
        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                color="white" if v < (vmin + 100) / 2 else "black", fontsize=9)
fig.colorbar(im, ax=ax, label="ANI (%)", shrink=0.8)
ax.set_title(f"{sp}: average nucleotide identity")
fig.tight_layout()
fig.savefig(f"{sp}_ANI_heatmap.png", dpi=300)
print("saved", f"{sp}_ANI_heatmap.png")
