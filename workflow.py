#!/usr/bin/env python3
"""
Single-Cell RNA-Seq Workflow: Zebrafish Retina Regeneration (Project 1)
Framework: Scanpy / AnnData
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

def main():
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)

    data_dir = Path("data")
    sample_dirs = sorted([
        p for p in data_dir.iterdir() 
        if p.is_dir() and (p / "filtered_feature_bc_matrix").exists()
    ])

    print(f"[INFO] Found {len(sample_dirs)} samples.")
    adatas = {}
    for sample_path in sample_dirs:
        sample_name = sample_path.name
        matrix_dir = sample_path / "filtered_feature_bc_matrix"
        adata_sample = sc.read_10x_mtx(path=matrix_dir, var_names="gene_symbols", cache=False)
        adata_sample.obs["sample"] = sample_name

        if sample_name.startswith("ctrl"):
            adata_sample.obs["condition"] = "ctrl"
            adata_sample.obs["replicate"] = sample_name.replace("ctrl", "")
        elif "dp" in sample_name:
            cond, rep = sample_name.split("dp")
            adata_sample.obs["condition"] = f"{cond}dp"
            adata_sample.obs["replicate"] = rep
        else:
            adata_sample.obs["condition"] = sample_name
            adata_sample.obs["replicate"] = "1"

        adatas[sample_name] = adata_sample

    adata = ad.concat(adatas, label="sample_batch", index_unique="-", join="outer")
    adata.var_names_make_unique()
    print(f"[SUCCESS] Combined AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

    # QC Metrics
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl"))
    adata.var["hb"] = adata.var_names.str.startswith(("hba", "hbb", "hb-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)

    # Filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 15, :].copy()

    # Normalization & Scaling
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # HVG & PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True)

    # Neighbors, UMAP & Clustering
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.6, key_added="leiden_res_0.6")

    print("[SUCCESS] Project 1 workflow completed.")
    return adata

if __name__ == "__main__":
    main()
