"""Cluster spots by their expression and draw the clusters on the tissue.

This is the same recipe used for single cell data: normalize, log transform,
pick variable genes, PCA, neighbour graph and Leiden clustering. The new part
is that we can now see where each cluster sits in the tissue. Clusters that
form clear regions are often called spatial domains.
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

warnings.filterwarnings("ignore")
sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config
seed = cfg["seed"]

adata = sc.read_h5ad(snakemake.input[0])

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=cfg["n_top_genes"], flavor="seurat")
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable", random_state=seed)
sc.pp.neighbors(adata, n_pcs=cfg["n_pcs"], random_state=seed)
sc.tl.umap(adata, random_state=seed)
sc.tl.leiden(adata, resolution=cfg["leiden_resolution"], flavor="igraph",
             n_iterations=2, random_state=seed, key_added="cluster")
print(adata.obs["cluster"].value_counts().sort_index())

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1, 1.2]})
sc.pl.umap(adata, color="cluster", ax=axes[0], show=False, frameon=False,
           legend_loc="on data", title="Clusters on UMAP")
sq.pl.spatial_scatter(adata, color="cluster", ax=axes[1], title="Clusters on the tissue")
fig.savefig(snakemake.output.tissue, dpi=150, bbox_inches="tight")
plt.close(fig)

sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers = markers[markers["pvals_adj"] < 0.05].groupby("group", observed=True).head(30)
markers.to_csv(snakemake.output.markers, sep="\t", index=False)

dp = sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, standard_scale="var",
                                     show=False, return_fig=True)
dp.savefig(snakemake.output.dotplot, dpi=150, bbox_inches="tight")
plt.close("all")

adata.write_h5ad(snakemake.output.adata, compression="gzip")
