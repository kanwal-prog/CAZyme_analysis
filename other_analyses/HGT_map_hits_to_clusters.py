#!/usr/bin/env python3
"""
Put AMR / virulence / plasmid hits (from abricate) onto the pangenome:
which gene cluster is each hit in, and is that cluster core, accessory or singleton?
Also counts clusters with mobile-element-like annotations (transposase, integrase, ...)
in each bin.

Usage (run from ~/Desktop/CAZyme_analysis, inside the anvio-9 environment):
    python map_hits_to_clusters.py ecoli ECOLI
    python map_hits_to_clusters.py bfragilis BFRAGILIS

Needs:
    <sp>/hits/<genome>_<db>.tsv                  abricate output
    <sp>/anvio/gene_calls/<genome>.txt           from anvi-export-gene-calls
    <sp>/anvio/<SPU>_PAN_SUMMARY_final/<SPU>_PAN_gene_clusters_summary.txt.gz
"""
import glob
import os
import sys

import pandas as pd

sp, spu = sys.argv[1], sys.argv[2]
base = f"{sp}/anvio"
summary_file = f"{base}/{spu}_PAN_SUMMARY_final/{spu}_PAN_gene_clusters_summary.txt.gz"

summ = pd.read_csv(summary_file, sep="\t")
genes = summ[["genome_name", "gene_callers_id", "gene_cluster_id", "bin_name"]].copy()
genes["gene_callers_id"] = genes["gene_callers_id"].astype(float)

# ---------------------------------------------------------------- 1. abricate hits
rows = []
for f in sorted(glob.glob(f"{sp}/hits/*.tsv")):
    genome, db = os.path.basename(f)[:-4].rsplit("_", 1)
    try:
        h = pd.read_csv(f, sep="\t")
    except pd.errors.EmptyDataError:
        continue
    if h.empty:
        continue
    h["genome_name"] = genome
    h["db"] = db
    rows.append(h)

if not rows:
    print("No hits found in", f"{sp}/hits/")
else:
    hits = pd.concat(rows, ignore_index=True).rename(
        columns={"%IDENTITY": "pct_identity", "%COVERAGE": "pct_coverage"}
    )

    # ------------------------------------------- 2. which gene call does a hit overlap?
    calls = {}

    def load_calls(g):
        if g not in calls:
            calls[g] = pd.read_csv(f"{base}/gene_calls/{g}.txt", sep="\t")
        return calls[g]

    ids = []
    for _, r in hits.iterrows():
        gc = load_calls(r["genome_name"])
        c = gc[gc["contig"] == r["SEQUENCE"]]
        # abricate coordinates are 1-based, anvi'o gene calls are 0-based
        ov = c["stop"].clip(upper=r["END"]) - c["start"].clip(lower=r["START"] - 1)
        ov = ov[ov > 0]
        ids.append(float(c.loc[ov.idxmax(), "gene_callers_id"]) if len(ov) else float("nan"))
    hits["gene_callers_id"] = ids

    # ----------------------------------------------- 3. gene cluster and bin of each hit
    m = hits.merge(genes, on=["genome_name", "gene_callers_id"], how="left")
    keep = [c for c in ["genome_name", "db", "GENE", "PRODUCT", "RESISTANCE",
                        "pct_identity", "pct_coverage", "SEQUENCE", "START", "END",
                        "gene_callers_id", "gene_cluster_id", "bin_name"] if c in m.columns]
    out = m[keep]
    out.to_csv(f"{sp}/hits_mapped_to_clusters.tsv", sep="\t", index=False)

    print("\n=== Hits per database and pangenome bin (a hit = one gene found in one genome) ===")
    print(out.groupby(["db", "bin_name"], dropna=False).size().unstack(fill_value=0))

    print("\n=== Which genomes carry each gene (1 = present) ===")
    pres = (out.assign(v=1)
               .pivot_table(index=["db", "GENE"], columns="genome_name",
                            values="v", aggfunc="max", fill_value=0))
    # show genomes without any hit as 0 columns too
    pres = pres.reindex(columns=sorted(genes["genome_name"].unique()), fill_value=0)
    print(pres.to_string())
    pres.to_csv(f"{sp}/hits_presence_absence.tsv", sep="\t")

    unmapped = out["gene_cluster_id"].isna().sum()
    if unmapped:
        print(f"\nNote: {unmapped} hit(s) did not overlap a called gene (often partial "
              "genes or RNA genes), so they have no cluster.")

# ----------------------------------------------- 4. mobile-element-like annotations
text_cols = [c for c in ["COG24_FUNCTION", "KOfam", "KEGG_BRITE"] if c in summ.columns]
summ["_t"] = summ[text_cols].fillna("").astype(str).agg(" ".join, axis=1)
pattern = (r"transposase|integrase|recombinase|resolvase|phage|conjug|relaxase|"
           r"mobiliz|insertion sequence")
summ["mobile"] = summ["_t"].str.contains(pattern, case=False, regex=True)
cl = summ.groupby("gene_cluster_id").agg(bin=("bin_name", "first"), mobile=("mobile", "any"))
tab = cl.groupby("bin")["mobile"].agg(["sum", "count"])
tab["percent"] = (100 * tab["sum"] / tab["count"]).round(1)
print("\n=== Gene clusters with a transposase/integrase/phage-like annotation, per bin ===")
print(tab)
print("(This is a rough proxy for mobile DNA, not proof of horizontal transfer.)")
