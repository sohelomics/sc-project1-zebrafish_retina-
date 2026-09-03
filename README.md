# sc-project1-zebrafish_retina-
# Single-Cell RNA-Seq Processing and Concatenation Pipeline

A Python-based workflow for loading, annotating, and merging multi-sample single-cell RNA sequencing (scRNA-seq) datasets into a unified `AnnData` object using **Scanpy** and **AnnData**.

---

## 📌 Project Overview

This repository provides an end-to-end data ingestion pipeline for scRNA-seq analyses. It automates:
- Automated sample directory discovery across multiple experimental conditions.
- Loading raw or filtered 10x Genomics matrices (`filtered_feature_bc_matrix`) or flat `.tsv` files.
- Sample metadata extraction (e.g., conditions like `ctrl`, `3dp`, `7dp` and replicate numbers) from directory identifiers.
- Resolution of duplicate gene symbols and cell barcodes across individual samples.
- Seamless concatenation into a single unified `AnnData` object (`adata`).

---

## 📁 Repository & Data Structure

Place your raw input matrices inside a `data/` directory. The automated loader expects folder structures matching the following format:

```text
.
├── data/
│   ├── ctrl1/
│   │   └── filtered_feature_bc_matrix/
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   ├── 3dp1/
│   │   └── filtered_feature_bc_matrix/
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   └── 7dp1/
│       └── filtered_feature_bc_matrix/
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── load_and_concat.py
├── README.md
└── .gitignore
