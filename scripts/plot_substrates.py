#!/usr/bin/env python3
"""Substrate distribution of CAZyme gene clusters (CGCs) from dbCAN easy_substrate output.

One genome  -> horizontal bar chart of CGC counts per predicted substrate.
Several     -> grouped bars, one colour per genome, for comparison.
--style heatmap -> substrate x genome grid with counts (easier to read for 4+ genomes).
--style stacked -> one stacked bar per genome, coloured by substrate (teacher's figure style).

  conda activate dbcan
  python scripts/plot_substrates.py -i dbcan_teacher/Bfrag_NCTC9343/4_sub_ncbi_gff -l "B. fragilis NCTC9343" --fmt png
  python scripts/plot_substrates.py --style stacked \
    -i dbcan_teacher/Bfrag_NCTC9343/4_sub_ncbi_gff dbcan_teacher/Ecoli_k12_MG1655/4_sub_ncbi_gff \
    -l "B. fragilis NCTC9343" "E. coli K12" --fmt png

Substrate per CGC: the dbCAN-PUL substrate, or the dbCAN-sub substrate when no PUL matched.
Output: <out>.<fmt> and <out>.tsv (the counts behind the figure).
"""
import argparse
import csv
import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# categorical slots, fixed order (one per genome)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"


def read_substrates(d, min_score=0):
    path = os.path.join(d, "substrate_prediction.tsv")
    if not os.path.exists(path):
        raise SystemExit(f"not found: {path}")
    counts = Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pul = row.get("dbCAN-PUL substrate", "").strip()
            score = float(row.get("bitscore") or 0)
            sub = pul if pul and score >= min_score else row.get("dbCAN-sub substrate", "").strip()
            if not sub or sub == "-":
                continue
            for s in {x.strip().lower() for x in sub.replace(";", ",").split(",")}:
                if s and s != "-":
                    counts[s] += 1
    return counts


def stacked(subs, data, labels, path):
    cmap = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    fig, ax = plt.subplots(figsize=(1.6 * len(labels) + 6.5, 7))
    bottoms = [0] * len(data)
    for k, s in enumerate(subs):
        vals = [c[s] for c in data]
        name = s.title() if len(s) <= 45 else s.title()[:44] + "..."
        ax.bar(range(len(data)), vals, bottom=bottoms, width=0.7, color=cmap[k % len(cmap)],
               edgecolor="white", linewidth=0.5, label=name)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
    ax.set_ylabel("Number of predicted PUL substrates", fontsize=12)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def heatmap(subs, data, labels, path):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "blue", ["#f7f7f5", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    m = [[c[s] for c in data] for s in subs]
    vmax = max(max(r) for r in m) or 1
    fig, ax = plt.subplots(figsize=(1.1 * len(labels) + 3, 0.3 * len(subs) + 1.8))
    ax.imshow(m, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i, row in enumerate(m):
        for j, v in enumerate(row):
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=8,
                        color="white" if v > vmax * 0.5 else INK)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8, color=INK)
    ax.set_yticks(range(len(subs)))
    ax.set_yticklabels(subs, fontsize=8, color=INK)
    ax.set_xticks([x - 0.5 for x in range(1, len(labels))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(subs))], minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Predicted CGC substrates (number of CGCs)", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True, help="dbCAN substrate output folder(s)")
    ap.add_argument("-l", "--label", nargs="+", help="genome label(s), same order as -i")
    ap.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"])
    ap.add_argument("-o", "--out", help="output file name without extension")
    ap.add_argument("--top", type=int, default=0, help="show only the N most frequent substrates")
    ap.add_argument("--style", default="bar", choices=["bar", "heatmap", "stacked"])
    ap.add_argument("--min-score", type=float, default=0, help="ignore PUL matches below this bitscore")
    a = ap.parse_args()

    labels = a.label or [os.path.basename(os.path.normpath(d)) for d in a.input]
    if len(labels) != len(a.input):
        raise SystemExit("give one -l label per -i folder")
    if a.style == "bar" and len(labels) > len(COLORS):
        raise SystemExit(f"at most {len(COLORS)} genomes per figure; split into several figures")

    data = [read_substrates(d, a.min_score) for d in a.input]
    total = Counter()
    for c in data:
        total.update(c)
    subs = [s for s, _ in sorted(total.items(), key=lambda kv: (-kv[1], kv[0]))]
    if a.top:
        subs = subs[:a.top]
    if not subs:
        raise SystemExit("no substrate predictions found")

    single = len(data) == 1
    out = a.out or ("substrate_" + labels[0].replace(" ", "_").replace(".", "") if single
                    else "substrate_comparison")

    with open(out + ".tsv", "w") as fh:
        fh.write("substrate\t" + "\t".join(labels) + "\n")
        for s in subs:
            fh.write(s + "\t" + "\t".join(str(c[s]) for c in data) + "\n")

    if a.style == "stacked":
        stacked(sorted(total), data, labels, f"{out}.{a.fmt}")
        print(f"wrote {out}.{a.fmt} and {out}.tsv")
        return

    if a.style == "heatmap":
        heatmap(subs, data, labels, f"{out}.{a.fmt}")
        print(f"wrote {out}.{a.fmt} and {out}.tsv")
        return

    n, k = len(subs), len(data)
    bar_h = 0.8 / k
    fig, ax = plt.subplots(figsize=(7.5, max(3, 0.28 * n * max(1, k ** 0.6) + 1.2)))
    for j, (c, lab) in enumerate(zip(data, labels)):
        ys = [i + (j - (k - 1) / 2) * bar_h for i in range(n)]
        ax.barh(ys, [c[s] for s in subs], height=bar_h * 0.85, color=COLORS[j],
                label=lab, edgecolor="white", linewidth=0.6)
        if single:
            for y, s in zip(ys, subs):
                ax.text(c[s] + 0.1, y, str(c[s]), va="center", fontsize=8, color=MUTED)

    ax.set_yticks(range(n))
    ax.set_yticklabels(subs, fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel("Number of CGCs", fontsize=9, color=INK)
    ax.xaxis.get_major_locator().set_params(integer=True)
    ax.tick_params(axis="x", labelsize=8, colors=MUTED)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    title = (f"Predicted CGC substrates - {labels[0]}" if single
             else "Predicted CGC substrates")
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    if not single:
        ax.legend(fontsize=8, frameon=False, loc="lower right")

    fig.tight_layout()
    fig.savefig(f"{out}.{a.fmt}", dpi=300)
    print(f"wrote {out}.{a.fmt} and {out}.tsv")


if __name__ == "__main__":
    main()
