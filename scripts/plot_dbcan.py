#!/usr/bin/env python3
"""
Auto-generate dbCAN plots from output directory.

Usage:
  Copy/place this script into a dbCAN output directory (or multiple) and run:

    # From inside a single output directory:
    cd b_xylanisolvens_output_fna_sub
    python ../plot_dbcan.py

    # Or specify one or more output directories:
    python plot_dbcan.py -i b_xylanisolvens_output_fna_sub
    python plot_dbcan.py -i dir1 dir2 dir3 -l "B. xylanisolvens" "E. coli" "B. fragilis"

  Generates:
    1) cazyme_abundance.pdf    - Barplot of CAZyme counts per genome (Fig. 9 style)
    2) pul_substrates.pdf      - Stacked barplot of PUL substrates per genome (Fig. 10 style)
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from pathlib import Path

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"


# ---------------------------------------------------------------------------
# Taxonomy color map (extend as needed for your species)
# ---------------------------------------------------------------------------
TAXONOMY_COLORS = {
    "B. caecimuris":    "#6aae4a",
    "B. fragilis":      "#999999",
    "B. muris":         "#4472c4",
    "B. xylanisolvens":  "#c0504d",
    "P. ovatus":        "#d4c36a",
    "P. sartorii":      "#e8943a",
    "B. ovatus":        "#d4c36a",
    "B. thetaiotaomicron": "#8c64a8",
    "E. coli":          "#f4a582",
}
DEFAULT_BAR_COLOR = "#5b9bd5"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def detect_output_dirs(paths):
    """Return list of valid dbCAN output directories."""
    dirs = []
    for p in paths:
        p = Path(p)
        if p.is_dir() and (p / "overview.tsv").exists():
            dirs.append(p)
        elif p.is_file() and p.name == "overview.tsv":
            dirs.append(p.parent)
    return dirs


def count_cazymes(output_dir, min_tools=2):
    """Count CAZymes from overview.tsv (genes found by >= min_tools)."""
    overview = output_dir / "overview.tsv"
    df = pd.read_csv(overview, sep="\t")
    return int((df["#ofTools"] >= min_tools).sum())


def load_substrate_counts(output_dir):
    """Load substrate_prediction.tsv and return substrate value counts."""
    sp = output_dir / "substrate_prediction.tsv"
    if not sp.exists() or sp.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(sp, sep="\t", comment="#")
    except Exception:
        return None
    if df.empty:
        return None
    for col in ["dbCAN-PUL substrate", "substrate"]:
        if col in df.columns:
            return df[col].dropna().str.strip().value_counts()
    # fallback: 3rd column
    if df.shape[1] >= 3:
        return df.iloc[:, 2].dropna().str.strip().value_counts()
    return None


def get_color_for_label(label):
    """Match label to taxonomy color map."""
    for taxon, color in TAXONOMY_COLORS.items():
        if taxon.lower() in label.lower():
            return color
    return DEFAULT_BAR_COLOR


# ---------------------------------------------------------------------------
# Plot 1: CAZyme abundance barplot (Fig. 1 style)
# ---------------------------------------------------------------------------

def plot_cazyme_abundance(labels, counts, output_file, width=None, height=8, dpi=300):
    n = len(labels)
    colors = [get_color_for_label(l) for l in labels]

    fig_w = width if width else max(4, n * 0.9 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, height))

    x = np.arange(n)
    ax.bar(x, counts, color=colors, edgecolor="white", linewidth=0.5, width=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Number of CAZymes", fontsize=12)
    ax.set_xlim(-0.5, n - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Build legend from unique taxonomy colors used
    seen = {}
    for label, color in zip(labels, colors):
        for taxon in TAXONOMY_COLORS:
            if taxon.lower() in label.lower() and taxon not in seen:
                seen[taxon] = color
                break
    if seen:
        handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in seen.values()]
        legend_labels = [f"$\\it{{{t}}}$" for t in seen.keys()]
        ax.legend(handles, legend_labels, loc="upper right", fontsize=9,
                  frameon=True, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
    print(f"  -> {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 2: PUL substrate stacked barplot (Fig. 2 style)
# ---------------------------------------------------------------------------

def plot_pul_substrates(labels, substrate_data, output_file, width=None, height=8, dpi=300):
    # Build DataFrame: rows=substrates, columns=genomes
    df = pd.DataFrame(substrate_data).fillna(0).astype(int)
    df.index = df.index.str.replace("_", " ").str.title()
    df.index = df.index.map(lambda x: x[:50] + "..." if len(x) > 50 else x)
    df = df.sort_index()

    n_substrates = len(df.index)
    n_genomes = len(df.columns)

    if n_substrates <= 20:
        cmap = matplotlib.colormaps.get_cmap("tab20")
        colors = [cmap(i) for i in range(n_substrates)]
    else:
        cmap = matplotlib.colormaps.get_cmap("tab20")
        cmap2 = matplotlib.colormaps.get_cmap("tab20b")
        colors = [cmap(i) for i in range(20)] + \
                 [cmap2(i) for i in range(n_substrates - 20)]

    fig_w = width if width else max(4, n_genomes * 1.5 + 4)
    fig, ax = plt.subplots(figsize=(fig_w, height))

    bottom = np.zeros(n_genomes)
    for i, substrate in enumerate(df.index):
        ax.bar(range(n_genomes), df.loc[substrate].values,
               bottom=bottom, color=colors[i], edgecolor="white",
               linewidth=0.5, label=substrate, width=0.7)
        bottom += df.loc[substrate].values

    ax.set_xticks(range(n_genomes))
    ax.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Number of predicted PUL substrates", fontsize=12)
    ax.set_xlim(-0.5, n_genomes - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
              frameon=True, framealpha=0.9, ncol=1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
    print(f"  -> {output_file}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Auto-generate dbCAN plots from output directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "-i", "--input", nargs="*", default=None,
        help="dbCAN output directory/directories. "
             "If omitted, uses current directory.",
    )
    parser.add_argument(
        "-l", "--labels", nargs="*", default=None,
        help="Labels for each genome (same order as -i). "
             "If omitted, directory names are used.",
    )
    parser.add_argument(
        "-o", "--outdir", default=None,
        help="Directory for output figures (default: current directory).",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--width", type=float, default=None)
    parser.add_argument("--height", type=float, default=8)
    parser.add_argument(
        "--min_tools", type=int, default=2,
        help="Minimum number of tools for a CAZyme to be counted (default: 2).",
    )
    parser.add_argument(
        "--fmt", default="pdf", choices=["pdf", "png", "svg"],
        help="Output figure format (default: pdf).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Resolve input directories
    input_paths = args.input if args.input else ["."]
    dirs = detect_output_dirs(input_paths)

    if not dirs:
        print("ERROR: No valid dbCAN output directory found.", file=sys.stderr)
        print("  Either run from inside an output dir or use: -i <dir1> [dir2 ...]",
              file=sys.stderr)
        sys.exit(1)

    # Labels
    if args.labels and len(args.labels) == len(dirs):
        labels = args.labels
    else:
        labels = [d.resolve().name for d in dirs]

    outdir = Path(args.outdir) if args.outdir else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Found {len(dirs)} output directory(s): {[str(d) for d in dirs]}")

    # ---- Plot 1: CAZyme abundance ----
    cazyme_counts = [count_cazymes(d, args.min_tools) for d in dirs]
    cazyme_file = outdir / f"fig1_cazyme_abundance.{args.fmt}"
    print(f"\n[Fig 1] CAZyme abundance (>={args.min_tools} tools):")
    for lab, cnt in zip(labels, cazyme_counts):
        print(f"  {lab}: {cnt} CAZymes")
    plot_cazyme_abundance(labels, cazyme_counts, str(cazyme_file),
                          args.width, args.height, args.dpi)

    # ---- Plot 2: PUL substrates ----
    substrate_data = {}
    has_substrates = False
    for lab, d in zip(labels, dirs):
        counts = load_substrate_counts(d)
        if counts is not None and len(counts) > 0:
            substrate_data[lab] = counts
            has_substrates = True
        else:
            print(f"  [skip] No substrate_prediction.tsv in {d}")

    if has_substrates:
        sub_file = outdir / f"fig2_pul_substrates.{args.fmt}"
        print(f"\n[Fig 2] PUL substrate prediction:")
        plot_pul_substrates(list(substrate_data.keys()), substrate_data,
                            str(sub_file), args.width, args.height, args.dpi)
    else:
        print("\n[Fig 2] Skipped — no substrate prediction data found.")

    print("\nDone!")


if __name__ == "__main__":
    main()
