# Single-Cell RNA-seq Analysis Pipeline

A complete scRNA-seq analysis workflow implemented in R using [Seurat](https://satijalab.org/seurat/), applied to a human brain organoid dataset (DS1). The pipeline covers everything from raw count matrix loading through quality control, dimensionality reduction, clustering, marker gene visualization, and cell type annotation. The original tutorial is present at [QuadBio](https://github.com/quadbio/scRNAseq_analysis_vignette)


## Overview

This pipeline was used to analyze single-cell transcriptomic data and identify major brain cell populations including dorsal and ventral telencephalon progenitors, intermediate progenitors, and various neuronal subtypes spanning telencephalic and diencephalic/midbrain regions.


## Pipeline Steps

| Step | Description |
|------|-------------|
| **Data loading** | Read 10X Genomics sparse matrix (`.mtx`, barcodes, features) |
| **Quality control** | Filter by gene count (500–5000) and mitochondrial content (<5%) |
| **Normalization** | Log-normalization via `NormalizeData` |
| **Feature selection** | Top 3,000 highly variable genes |
| **Scaling & regression** | Regress out `nFeature_RNA` and `percent.mt` |
| **PCA** | 50 PCs; significant dimensions selected via ElbowPlot |
| **Clustering** | KNN graph + Louvain clustering (resolution = 1) → 15 clusters |
| **Embedding** | tSNE and UMAP using top 20 PCs |
| **Marker visualization** | FeaturePlot and VlnPlot for canonical and cluster-specific markers |
| **Cell type annotation** | Manual annotation based on known marker gene expression |
| **Trajectory subset** | Dorsal telencephalon subset (clusters 0, 2, 5, 6, 10) |
| **Cell cycle scoring** | S/G2M scoring and regression for trajectory analysis |


## Key Results

### Cell Type Annotation (UMAP)

Fifteen clusters were identified and annotated into the following cell types:

- **Dorsal telencephalon:** NPC, G2M NPC, IP, Neuron
- **Diencephalon & midbrain:** NPC, G2M NPC, IP, Excitatory neuron, Inhibitory neuron, Early inhibitory neuron
- **Other:** MGE-like neuron, Ventral telen. neuron, Midbrain-hindbrain boundary neuron
- **Unresolved:** Unknown 1, Unknown 2

### Canonical Marker Genes Used

| Marker | Cell identity |
|--------|--------------|
| `MKI67` | Cycling cells (G2M) |
| `NES` | Neural progenitor cells |
| `DCX` | Neurons |
| `FOXG1` | Telencephalon |
| `EMX1`, `PAX6`, `GLI3`, `EOMES`, `NEUROD6` | Dorsal telencephalon |
| `DLX2`, `DLX5`, `NKX2.1`, `SOX6`, `NR2F2` | Ventral telencephalon / MGE |
| `OTX2`, `LHX9`, `RSPO3`, `TFAP2A` | Non-telencephalic regions |
| `NEUROD2`, `NEUROD6` | Dorsal excitatory neurons |
| `LHX9` | Diencephalon/midbrain (clusters 4–5) |


## Requirements

```r
library(Seurat)   # >= 4.x (v3 assay format)
library(Matrix)
library(patchwork)
library(dplyr)
```


## Usage

### 1. Prepare your data

Place 10X Genomics output files in a `data/` directory:

```
data/
├── matrix.mtx.gz
├── barcodes.tsv.gz
└── features.tsv.gz
```

### 2. Run the pipeline

```r
source("scrna_pipeline.R")
```

Outputs are saved to `results/`:

```
results/
├── seurat_obj_all.rds       # Full annotated Seurat object
└── seurat_obj_dorsal.rds    # Dorsal telencephalon subset
```

### 3. Load saved objects

```r
seurat        <- readRDS("results/seurat_obj_all.rds")
seurat_dorsal <- readRDS("results/seurat_obj_dorsal.rds")
```


## Output Figures

| Figure | Description |
|--------|-------------|
| UMAP (annotated) | Cell type labels across all 15 clusters |
| tSNE / UMAP (unlabeled) | Initial embedding colored by sample |
| Cluster tSNE / UMAP | Numeric cluster identity plots |
| Canonical marker FeaturePlots | 9-gene panel across tSNE and UMAP |
| Cell-type marker heatmap | 22-gene panel across all clusters |
| NEUROD2 / NEUROD6 plots | FeaturePlot + VlnPlot |
| LHX9 plots | FeaturePlot + VlnPlot |


## Notes

- Seurat v3 assay format is enforced via `options(Seurat.object.assay.version = "v3")` for compatibility.
- Cell cycle genes (`cc.genes`) are removed from HVGs prior to trajectory analysis to prevent cell cycle effects from dominating the embedding.
- Clustering resolution of `1.0` yielded 15 biologically interpretable clusters for this dataset; adjust as needed for other datasets.
