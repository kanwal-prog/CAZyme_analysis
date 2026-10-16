#!/usr/bin/env python3
"""PUL figure: gene maps of selected CAZyme gene clusters (CGCs) across genomes.

For each PUL (a dbCAN-PUL reference ID), every genome's best-matching CGC is drawn
as gene arrows coloured by type and labelled (GH/GT/PL/CE family, SusC, SusD, ...).
Each locus is oriented so most of its CAZymes point right, to make strains comparable.

  conda activate dbcan
  python scripts/pul_figure.py                                  # default PUL panel
  python scripts/pul_figure.py --pul PUL0380 PUL0479 --fmt pdf  # your own selection

Find PUL IDs in dbcan_teacher/<genome>/4_sub_ncbi_gff/substrate_prediction.tsv (column PULID).

Output: figures/pul_figure.<fmt> and figures/pul_figure.tsv (genes drawn).
"""
import argparse
import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon

RESULTS = "dbcan_teacher"
def load_genomes(path="genomes.tsv"):
    """genome name -> figure label, from genomes.tsv"""
    with open(path) as fh:
        return {r["name"]: r["label"] for r in csv.DictReader(fh, delimiter="	")}


GENOMES = load_genomes()
DEFAULT_PULS = ["PUL0380", "PUL0556", "PUL0479", "PUL0166", "PUL0230"]

# fixed categorical order
TYPES = {
    "CAZyme": "#2a78d6", "SusC": "#eb6834", "SusD": "#eda100", "Transporter": "#1baf7a",
    "Regulator (TF/STP)": "#e87ba4", "Sulfatase": "#6250d6", "Peptidase": "#e34948",
    "Other": "#d6d5cf",
}
INK, MUTED = "#1f1f1e", "#6b6a64"


def products(genome):
    """protein_id -> NCBI product name."""
    out = {}
    path = f"{RESULTS}/{genome}/{genome}.gff"
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if "\tCDS\t" not in line:
                continue
            pid = re.search(r"protein_id=([^;\n]+)", line)
            prod = re.search(r"product=([^;\n]+)", line)
            if pid and prod:
                out[pid.group(1)] = prod.group(1).replace("%2C", ",")
    return out


def classify(gtype, annot, product):
    """-> (type, label)"""
    p = product.lower()
    if gtype == "CAZyme":
        fams = re.findall(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+)", annot)
        fams = list(dict.fromkeys(fams))
        return "CAZyme", "+".join(fams[:2])
    if "susd" in p or "ragb" in p:
        return "SusD", "SusD"
    if gtype == "TC" and "1.B.14" in annot:
        return "SusC", "SusC"
    if gtype == "TC":
        return "Transporter", "TC"
    if gtype in ("TF", "STP"):
        return "Regulator (TF/STP)", gtype
    if gtype == "Sulfatase":
        return "Sulfatase", "Sulf"
    if gtype == "Peptidase":
        return "Peptidase", "Pep"
    return "Other", ""


def best_cgc(genome, pul):
    """Highest-bitscore CGC in this genome matched to the given PUL."""
    path = f"{RESULTS}/{genome}/4_sub_ncbi_gff/substrate_prediction.tsv"
    best = None
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["PULID"] == pul:
                score = float(row["bitscore"] or 0)
                if best is None or score > best[1]:
                    best = (row["#cgcid"], score, row["dbCAN-PUL substrate"])
    return best


def cgc_genes(genome, cgcid, prod):
    contig, cgc = cgcid.split("|")
    genes = []
    with open(f"{RESULTS}/{genome}/4_sub_ncbi_gff/cgc_standard_out.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Contig ID"] == contig and row["CGC#"] == cgc:
                t, lab = classify(row["Gene Type"], row["Gene Annotation"],
                                  prod.get(row["Protein ID"], ""))
                genes.append(dict(start=int(row["Gene Start"]), end=int(row["Gene Stop"]),
                                  strand=row["Gene Strand"], type=t, label=lab,
                                  pid=row["Protein ID"], product=prod.get(row["Protein ID"], "")))
    genes.sort(key=lambda g: g["start"])
    # orient so most CAZymes point right
    caz = [g["strand"] for g in genes if g["type"] == "CAZyme"] or [g["strand"] for g in genes]
    if caz.count("-") > caz.count("+"):
        hi = max(g["end"] for g in genes)
        for g in genes:
            g["start"], g["end"] = hi - g["end"], hi - g["start"]
            g["strand"] = "+" if g["strand"] == "-" else "-"
        genes.sort(key=lambda g: g["start"])
    lo = genes[0]["start"]
    for g in genes:
        g["x0"], g["x1"] = (g["start"] - lo) / 1000, (g["end"] - lo) / 1000
    return genes


def arrow(ax, x0, x1, y, strand, color, h=0.32):
    head = min(0.35, (x1 - x0) * 0.4)
    if strand == "+":
        pts = [(x0, y - h / 2), (x1 - head, y - h / 2), (x1, y), (x1 - head, y + h / 2), (x0, y + h / 2)]
    else:
        pts = [(x1, y - h / 2), (x0 + head, y - h / 2), (x0, y), (x0 + head, y + h / 2), (x1, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, fc=color, ec="white", lw=0.8))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pul", nargs="+", default=DEFAULT_PULS, help="dbCAN-PUL IDs to show")
    ap.add_argument("--genomes", nargs="+", default=list(GENOMES), help="genome folder names")
    ap.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"])
    ap.add_argument("-o", "--out", default="figures/pul_figure")
    a = ap.parse_args()

    prods = {g: products(g) for g in a.genomes}
    rows = []  # (pul, substrate, genome, cgcid, score, genes)
    for pul in a.pul:
        for g in a.genomes:
            hit = best_cgc(g, pul)
            if hit:
                rows.append((pul, hit[2], g, hit[0], hit[1], cgc_genes(g, hit[0], prods[g])))
    if not rows:
        raise SystemExit("none of these PULs were found")

    width = max(max(x["x1"] for x in r[5]) for r in rows)
    n_pul = len(dict.fromkeys(r[0] for r in rows))
    fig_h = 0.62 * len(rows) + 0.5 * n_pul + 1.2
    fig, ax = plt.subplots(figsize=(13, fig_h))

    y, last_pul, used = 0, None, set()
    for pul, substrate, g, cgcid, score, genes in rows:
        if pul != last_pul:
            y += 0.55
            ax.text(-0.2, y, f"{substrate}  ({pul})", fontsize=10, fontweight="bold",
                    color=INK, ha="left", va="center")
            y += 0.7
            last_pul = pul
        ax.text(-0.4, y, GENOMES.get(g, g), fontsize=8.5, style="italic", color=INK,
                ha="right", va="center")
        ax.text(width + 0.4, y, f"{cgcid.split('|')[1]}\nbitscore {score:.0f}",
                fontsize=7, color=MUTED, ha="left", va="center")
        ax.plot([0, genes[-1]["x1"]], [y, y], color="#b9b8b1", lw=0.8, zorder=0)
        for k, gene in enumerate(genes):
            arrow(ax, gene["x0"], gene["x1"], y, gene["strand"], TYPES[gene["type"]])
            used.add(gene["type"])
            if gene["label"]:
                ax.text((gene["x0"] + gene["x1"]) / 2, y + (0.3 if k % 2 == 0 else -0.3),
                        gene["label"], fontsize=6.5, color=INK, ha="center",
                        va="bottom" if k % 2 == 0 else "top")
        y += 0.95

    ax.set_xlim(-0.3, width + 2.2)
    ax.set_ylim(y, 0)
    ax.set_yticks([])
    ax.set_xlabel("kb", fontsize=8, color=MUTED)
    ax.tick_params(axis="x", labelsize=7, colors=MUTED)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#d6d5cf")
    handles = [Patch(fc=c, label=t) for t, c in TYPES.items() if t in used]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06 * 12 / fig_h),
              ncol=len(handles), fontsize=8, frameon=False)
    ax.set_title("Polysaccharide utilization loci (PULs) across genomes",
                 fontsize=12, color=INK, loc="left")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(f"{a.out}.{a.fmt}", dpi=300, bbox_inches="tight")
    with open(f"{a.out}.tsv", "w") as fh:
        fh.write("PUL\tsubstrate\tgenome\tCGC\tbitscore\tprotein_id\tstart\tend\tstrand\ttype\tlabel\tNCBI_product\n")
        for pul, sub, g, cgcid, score, genes in rows:
            for x in genes:
                fh.write(f"{pul}\t{sub}\t{g}\t{cgcid}\t{score:.0f}\t{x['pid']}\t{x['start']}\t"
                         f"{x['end']}\t{x['strand']}\t{x['type']}\t{x['label']}\t{x['product']}\n")
    print(f"wrote {a.out}.{a.fmt} and {a.out}.tsv")


if __name__ == "__main__":
    main()
