"""Load the Visium data and remove poor quality spots.

Each Visium spot is 55 micrometres wide and usually holds somewhere between
1 and 10 cells. So a spot is a small mix of cells and not a single cell.
We check the number of reads, the number of genes and the share of
mitochondrial reads per spot, and we look at where on the tissue the poor
spots are. Spots with few reads often sit at tissue edges or in torn areas.
"""

import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq

warnings.filterwarnings("ignore")
sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config

adata = sq.read.visium("data/visium", counts_file="filtered_feature_bc_matrix.h5",
                       library_id=cfg["sample"])
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
print(adata)

# The H&E image on its own
fig, ax = plt.subplots(figsize=(7, 7))
sq.pl.spatial_scatter(adata, color=None, ax=ax, title="H&E image of the section")
fig.savefig(snakemake.output.tissue, dpi=150, bbox_inches="tight")
plt.close(fig)

# QC values drawn on the tissue
sq.pl.spatial_scatter(adata, color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                      title=["UMIs per spot", "Genes per spot", "Mitochondrial reads (%)"],
                      img_alpha=0.5, ncols=3, figsize=(6, 6))
plt.savefig(snakemake.output.qc_spatial, dpi=150, bbox_inches="tight")
plt.close("all")

# Histograms with the cutoffs we use
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, col, cut, label in [(axes[0], "total_counts", cfg["min_counts"], "UMIs per spot"),
                            (axes[1], "n_genes_by_counts", cfg["min_genes"], "Genes per spot"),
                            (axes[2], "pct_counts_mt", cfg["max_pct_mito"], "Mitochondrial reads (%)")]:
    ax.hist(adata.obs[col], bins=60, color="#4C72B0")
    ax.axvline(cut, color="#C44E52", ls="--")
    ax.set_title(label)
fig.suptitle("Spot quality. Red lines show the cutoffs.")
fig.tight_layout()
fig.savefig(snakemake.output.qc_hist, dpi=150)
plt.close(fig)

steps = [("Spots under tissue", adata.n_obs)]
keep = (adata.obs["total_counts"] >= cfg["min_counts"]) & \
       (adata.obs["n_genes_by_counts"] >= cfg["min_genes"]) & \
       (adata.obs["pct_counts_mt"] <= cfg["max_pct_mito"])
adata = adata[keep].copy()
steps.append(("Spots after filtering", adata.n_obs))
sc.pp.filter_genes(adata, min_cells=cfg["min_spots_per_gene"])
steps.append(("Genes kept", adata.n_vars))

summary = pd.DataFrame(steps, columns=["step", "count"])
summary.to_csv(snakemake.output.table, sep="\t", index=False)
print(summary)

adata.layers["counts"] = adata.X.copy()
adata.write_h5ad(snakemake.output.adata, compression="gzip")
