"""Estimate which cell types are present in each spot.

Since a spot holds several cells, we do not give each spot one label.
Instead we score every spot for a set of marker genes per cell type. The
score tells us how much that cell type's program is switched on in the spot
compared with the rest of the tissue. This is simpler than full deconvolution
tools such as cell2location, needs no single cell reference, and gives a clear
first picture of the tissue.

We then ask which cell types tend to show up in the same spots, by correlating
the scores across all spots.
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import squidpy as sq

warnings.filterwarnings("ignore")
sys.stdout = sys.stderr = open(snakemake.log[0], "w")
seed = snakemake.config["seed"]

adata = sc.read_h5ad(snakemake.input[0])

MARKERS = {
    "Tumour epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "ERBB2"],
    "Luminal hormone sensing": ["ESR1", "PGR", "AR", "FOXA1", "TFF1", "TFF3"],
    "Myoepithelial": ["KRT5", "KRT14", "KRT17", "TP63", "ACTA2", "MYLK"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5"],
    "T cells": ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "IL7R"],
    "B and plasma cells": ["MS4A1", "CD79A", "IGKC", "JCHAIN", "MZB1", "IGHG1"],
    "Macrophages": ["CD68", "CD163", "C1QA", "C1QB", "C1QC", "AIF1"],
    "Adipocytes": ["ADIPOQ", "PLIN1", "FABP4", "LEP", "CIDEA"],
}
MARKERS = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
for k, v in MARKERS.items():
    print(f"{k}: {len(v)} genes found ({', '.join(v)})")

score_cols = []
for cell_type, genes in MARKERS.items():
    col = f"score_{cell_type}"
    sc.tl.score_genes(adata, genes, score_name=col, random_state=seed)
    score_cols.append(col)

titles = list(MARKERS)
sq.pl.spatial_scatter(adata, color=score_cols, title=titles, ncols=3, figsize=(5, 5),
                      cmap="viridis")
plt.savefig(snakemake.output.tissue, dpi=150, bbox_inches="tight")
plt.close("all")

by_cluster = adata.obs.groupby("cluster", observed=True)[score_cols].mean()
by_cluster.columns = titles
by_cluster.to_csv(snakemake.output.table, sep="\t")
z = (by_cluster - by_cluster.mean()) / by_cluster.std()
fig, ax = plt.subplots(figsize=(9, 0.45 * len(z) + 2))
sns.heatmap(z, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "Scaled mean score"})
ax.set_ylabel("Cluster")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_title("Cell type scores in each cluster")
fig.tight_layout()
fig.savefig(snakemake.output.by_cluster, dpi=150)
plt.close(fig)

corr = adata.obs[score_cols].corr(method="spearman")
corr.index = corr.columns = titles
g = sns.clustermap(corr, cmap="RdBu_r", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f",
                   figsize=(8, 8), cbar_kws={"label": "Spearman correlation"})
g.figure.suptitle("Which cell types share spots", y=1.02)
g.savefig(snakemake.output.colocal, dpi=150, bbox_inches="tight")
plt.close("all")

adata.write_h5ad(snakemake.output.adata, compression="gzip")
