"""Spatial statistics with Squidpy.

1. Spatially variable genes. Moran's I asks whether a gene's level in a spot
   looks like its level in the spots around it. A high value means the gene
   forms patches or gradients in the tissue rather than random speckles.
2. Neighbourhood enrichment. For every pair of clusters we count how often
   their spots touch, and compare that to random shuffles of the labels.
   Positive scores mean two clusters sit next to each other more than chance.
3. Co occurrence. How likely we are to find a cluster at growing distances
   from another cluster.
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq

warnings.filterwarnings("ignore")
sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config
seed = cfg["seed"]

adata = sc.read_h5ad(snakemake.input[0])

# Build the spot graph from the Visium hexagonal grid
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=cfg["n_rings"])

# Moran's I on the variable genes
genes = adata.var_names[adata.var["highly_variable"]]
sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=100, seed=seed,
                       n_jobs=snakemake.threads, show_progress_bar=False)
moran = adata.uns["moranI"].sort_values("I", ascending=False)
moran.index.name = "gene"
moran.to_csv(snakemake.output.svg_table, sep="\t")
print(moran.head(20))

top = moran.index[: cfg["n_svg_to_plot"]].tolist()
sq.pl.spatial_scatter(adata, color=top, ncols=4, figsize=(4, 4), cmap="magma",
                      title=[f"{g} (I = {moran.loc[g, 'I']:.2f})" for g in top])
plt.savefig(snakemake.output.svg_plot, dpi=150, bbox_inches="tight")
plt.close("all")

sq.gr.nhood_enrichment(adata, cluster_key="cluster", n_perms=1000, seed=seed,
                       n_jobs=snakemake.threads, show_progress_bar=False)
# A cluster always touches itself the most, so we hide the diagonal.
# That lets the colour scale show the more interesting pairs.
z = pd.DataFrame(adata.uns["cluster_nhood_enrichment"]["zscore"],
                 index=adata.obs["cluster"].cat.categories,
                 columns=adata.obs["cluster"].cat.categories)
mask = np.eye(len(z), dtype=bool)
lim = np.nanmax(np.abs(z.values[~mask]))
fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(z, mask=mask, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".0f",
            ax=ax, cbar_kws={"label": "z score (positive = touch more than chance)"})
ax.set_xlabel("Cluster")
ax.set_ylabel("Cluster")
ax.set_title("Neighbourhood enrichment between clusters")
fig.tight_layout()
fig.savefig(snakemake.output.nhood, dpi=150)
plt.close(fig)

sq.gr.co_occurrence(adata, cluster_key="cluster", show_progress_bar=False, n_jobs=snakemake.threads)
occ = adata.uns["cluster_co_occurrence"]["occ"]
interval = adata.uns["cluster_co_occurrence"]["interval"][1:]
# Distances are in pixels of the full resolution image. The scale factors file
# tells us how many pixels one 55 micrometre spot covers.
lib = list(adata.uns["spatial"].keys())[0]
um_per_px = 55 / adata.uns["spatial"][lib]["scalefactors"]["spot_diameter_fullres"]
dist_um = interval * um_per_px

cats = adata.obs["cluster"].cat.categories
colours = dict(zip(cats, adata.uns["cluster_colors"]))
ncols = 3
nrows = int(np.ceil(len(cats) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.6 * nrows), sharex=True)
for ax, (i, c) in zip(axes.flat, enumerate(cats)):
    for j, other in enumerate(cats):
        ax.plot(dist_um, occ[i, j, :], color=colours[other], lw=2 if other == c else 1, label=other)
    ax.axhline(1, ls="--", c="grey", lw=0.8)
    ax.set_title(f"Around cluster {c}")
    ax.set_ylabel("Observed / expected")
for ax in axes.flat[len(cats):]:
    ax.axis("off")
for ax in axes[-1]:
    ax.set_xlabel("Distance (micrometres)")
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, title="Cluster", loc="center right")
fig.suptitle("How likely each cluster is to be found at a given distance from another")
fig.tight_layout(rect=(0, 0, 0.93, 0.97))
fig.savefig(snakemake.output.cooc, dpi=150)
plt.close(fig)

adata.write_h5ad(snakemake.output.adata, compression="gzip")
