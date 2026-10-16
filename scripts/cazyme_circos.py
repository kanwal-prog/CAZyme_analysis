#!/usr/bin/env python3
"""CAZyme circos plot for one genome, from dbcan_teacher results.

Tracks (outside -> inside):
  1. contigs with Mb ticks
  2. CAZymes on the + strand, coloured by class (GH, GT, PL, CE, AA, CBM)
  3. CAZymes on the - strand
  4. CGCs coloured by predicted substrate group (grey = no prediction)
  5. GC skew

  conda activate dbcan
  python scripts/cazyme_circos.py all                                              # every genome in genomes.tsv
  python scripts/cazyme_circos.py Bfrag_NCTC9343 bfragilis/Bfrag_NCTC9343.fna      # one genome

Output: figures/circos_<genome>.png and .svg
"""
import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos

RESULTS = "dbcan_teacher"
OUTDIR = "figures"

def load_genomes(path="genomes.tsv"):
    """genome name -> .fna path, from genomes.tsv"""
    with open(path) as fh:
        return {r["name"]: f"{r['folder']}/{r['name']}.fna" for r in csv.DictReader(fh, delimiter="	")}

CLASS_COLORS = {
    "GH": "#2a78c4", "GT": "#e08a1e", "PL": "#c2378a",
    "CE": "#2e9c5a", "AA": "#8a5cc7", "CBM": "#8c8c8c",
}
CLASS_ORDER = ["GH", "PL", "CE", "AA", "GT", "CBM"]  # which class wins for multi-domain genes

SUBSTRATE_GROUPS = {
    "Host glycan / mucin": ["host glycan", "mucin", "glycosaminoglycan",
                            "human milk oligosaccharide", "sialic acid"],
    "Plant fibre": ["xylan", "pectin", "cellulose", "hemicellulose", "arabinan",
                    "beta-mannan", "galactomannan", "alpha-mannan", "beta-glucan",
                    "xyloglucan", "arabinoxylan", "arabinogalactan", "fructan",
                    "inulin", "levan", "rhamnogalacturonan", "glucomannan"],
    "Storage glucan": ["starch", "glycogen", "trehalose", "pullulan", "dextran"],
    "Oligo/disaccharide": ["lactose", "melibiose", "cellobiose", "raffinose",
                           "galactooligosaccharide", "beta-galactooligosaccharide",
                           "beta-glucoside", "sucrose", "maltose"],
    "Chitin / peptidoglycan": ["chitin", "peptidoglycan"],
}
GROUP_COLORS = {
    "Host glycan / mucin": "#d1495b", "Plant fibre": "#3a9d5d",
    "Storage glucan": "#e6a23c", "Oligo/disaccharide": "#4f8fd1",
    "Chitin / peptidoglycan": "#8a5cc7", "Other substrate": "#5c5c5c",
    "No prediction": "#cfcfcf",
}


def substrate_group(sub):
    s = sub.strip().lower()
    for group, members in SUBSTRATE_GROUPS.items():
        if s in members:
            return group
    return "Other substrate"


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def gc_skew(seq, window):
    xs, ys = [], []
    for i in range(0, len(seq), window):
        w = seq[i:i + window]
        g, c = w.count("G"), w.count("C")
        xs.append(min(i + window / 2, len(seq)))
        ys.append((g - c) / (g + c) if g + c else 0)
    return xs, ys


def read_cazymes(d):
    """CAZyme genes (>=2 tools) with coordinates: list of (contig, start, end, strand, class)."""
    keep = {}
    with open(f"{d}/overview.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["#ofTools"]) < 2:
                continue
            fams = re.findall(r"(GH|GT|PL|CE|AA|CBM)\d", row["Recommend Results"])
            if fams:
                keep[row["Gene ID"]] = next(c for c in CLASS_ORDER if c in fams)
    genes = []
    with open(f"{d}/cgc.gff") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pid = re.search(r"protein_id=([^;]+)", f[8])
            if pid and pid.group(1) in keep:
                genes.append((f[0], int(f[3]), int(f[4]), f[6], keep[pid.group(1)]))
    return genes


def read_cgcs(d):
    """CGC spans with substrate group: list of (contig, start, end, group)."""
    subs = {}
    if os.path.exists(f"{d}/substrate_prediction.tsv"):
        with open(f"{d}/substrate_prediction.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                s = row["dbCAN-PUL substrate"] or row.get("dbCAN-sub substrate", "")
                if s:
                    subs[row["#cgcid"]] = substrate_group(s.split(",")[0])
    spans = {}
    with open(f"{d}/cgc_standard_out.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["Contig ID"], row["CGC#"])
            s, e = int(row["Gene Start"]), int(row["Gene Stop"])
            lo, hi = spans.get(key, (s, e))
            spans[key] = (min(lo, s), max(hi, e))
    return [(c, s, e, subs.get(f"{c}|{n}", "No prediction")) for (c, n), (s, e) in spans.items()]


def plot(genome, fna):
    d = f"{RESULTS}/{genome}/4_sub_ncbi_gff"
    if not os.path.exists(f"{d}/cgc_standard_out.tsv"):
        print(f"skip {genome}: results not finished")
        return
    seqs = read_fasta(fna)
    genes = read_cazymes(d)
    cgcs = read_cgcs(d)

    sectors = {name: len(s) for name, s in seqs.items()}
    circos = Circos(sectors, space=3 if len(sectors) > 1 else 0)
    total = sum(sectors.values())
    win = max(5000, total // 1000)

    for sector in circos.sectors:
        big = sector.size > total * 0.05
        outer = sector.add_track((96, 100))
        outer.axis(fc="#404040")
        if big:
            outer.xticks_by_interval(500_000, outer=True, label_size=7,
                                     label_formatter=lambda v: f"{v / 1e6:.1f} Mb")
        sector.text(sector.name if big else "", r=110, size=8)

        plus, minus = sector.add_track((86, 94)), sector.add_track((77, 85))
        for c, s, e, strand, cls in genes:
            if c == sector.name:
                (plus if strand == "+" else minus).rect(s, max(e, s + total // 2000),
                                                        fc=CLASS_COLORS[cls], lw=0)

        cgc_track = sector.add_track((66, 75))
        cgc_track.axis(fc="#f4f4f4", ec="none")
        for c, s, e, group in cgcs:
            if c == sector.name:
                cgc_track.rect(s, max(e, s + total // 1000), fc=GROUP_COLORS[group], lw=0)

        skew = sector.add_track((45, 63))
        xs, ys = gc_skew(seqs[sector.name], win)
        if len(xs) > 1:
            pos = [max(y, 0) for y in ys]
            neg = [min(y, 0) for y in ys]
            lim = max(abs(y) for y in ys) or 1
            skew.fill_between(xs, pos, y2=0, vmin=-lim, vmax=lim, fc="#4d4d4d", lw=0)
            skew.fill_between(xs, neg, y2=0, vmin=-lim, vmax=lim, fc="#a8a8a8", lw=0)

    n_cgc_pred = sum(g != "No prediction" for *_, g in cgcs)
    circos.text(f"{genome.replace('_', ' ')}\n{len(genes)} CAZymes\n"
                f"{len(cgcs)} CGCs ({n_cgc_pred} with substrate)", size=10)

    fig = circos.plotfig(figsize=(8, 8))
    used_cls = [c for c in CLASS_COLORS if any(g[4] == c for g in genes)]
    used_grp = [g for g in GROUP_COLORS if any(x[3] == g for x in cgcs)]
    leg1 = fig.legend([Patch(fc=CLASS_COLORS[c]) for c in used_cls], used_cls,
                      title="CAZyme class (rings 2-3: +/- strand)",
                      loc="upper left", fontsize=8, title_fontsize=8, frameon=False)
    fig.legend([Patch(fc=GROUP_COLORS[g]) for g in used_grp], used_grp,
               title="CGC substrate (ring 4)", loc="upper right",
               fontsize=8, title_fontsize=8, frameon=False)
    fig.add_artist(leg1)
    fig.legend([Patch(fc="#4d4d4d"), Patch(fc="#a8a8a8")], ["GC skew +", "GC skew -"],
               loc="lower right", fontsize=8, frameon=False)

    os.makedirs(OUTDIR, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(f"{OUTDIR}/circos_{genome}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR}/circos_{genome}.png/.svg")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "all":
        for g, f in load_genomes().items():
            plot(g, f)
    elif len(sys.argv) == 3:
        plot(sys.argv[1], sys.argv[2])
    else:
        sys.exit(__doc__)
