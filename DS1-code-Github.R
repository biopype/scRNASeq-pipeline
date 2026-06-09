# ============================================================
# Single-Cell RNA-seq Analysis Pipeline using Seurat
# ============================================================

# Load Seurat package
library(Seurat)

set.seed(42)

data_dir <- "data"
results_dir <- "results"

dir.create(results_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# Seurat Object Setup
# ------------------------------------------------------------

# Tell Seurat to use Assay v3 format instead of Assay5 (v5)
options(Seurat.object.assay.version = "v3")

# Read 10X Genomics formatted dataset
counts <- Read10X(data.dir = data_dir)

# Create Seurat object to store expression data and analysis results
seurat <- CreateSeuratObject(counts, project = "DS1")


# ------------------------------------------------------------
# Manual Loading of Matrix Files (Alternative Method)
# ------------------------------------------------------------

library(Matrix)

# Read sparse count matrix
counts <- readMM(file.path(data_dir, "matrix.mtx.gz"))

# Read cell barcodes
barcodes <- read.table(
  file.path(data_dir, "barcodes.tsv.gz"),
  stringsAsFactors = FALSE
)[,1]

# Read gene/feature annotations
features <- read.csv(
  file.path(data_dir, "features.tsv.gz"),
  stringsAsFactors = FALSE,
  sep = "\t",
  header = FALSE
)

# Assign gene names and cell barcodes
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes

# Recreate Seurat object using manually loaded data
seurat <- CreateSeuratObject(counts, project = "DS1")


# ============================================================
# Quality Control (QC)
# ============================================================

# ------------------------------------------------------------
# Calculate Percentage of Mitochondrial Gene Expression
# High mitochondrial content may indicate stressed or dying cells
# ------------------------------------------------------------

seurat[["percent.mt"]] <- PercentageFeatureSet(
  seurat,
  pattern = "^MT[-\\.]"
)

# ------------------------------------------------------------
# Visualize QC Metrics
# ------------------------------------------------------------
# nFeature_RNA = number of detected genes per cell
# nCount_RNA   = total RNA molecules (UMIs) per cell
# percent.mt   = mitochondrial transcript percentage

VlnPlot(
  seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# Scatter plots for QC relationships
library(patchwork)

plot1 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

plot2 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

plot1 + plot2


# ------------------------------------------------------------
# Filter Low-Quality Cells
# ------------------------------------------------------------
# Keep cells with:
# - More than 500 genes
# - Less than 5000 genes
# - Less than 5% mitochondrial expression

seurat <- subset(
  seurat,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)


# ============================================================
# Normalization and Feature Selection
# ============================================================

# ------------------------------------------------------------
# Normalize Data
# Makes expression values comparable across cells
# ------------------------------------------------------------

seurat <- NormalizeData(seurat)

# ------------------------------------------------------------
# Identify Highly Variable Genes (HVGs)
# These genes are most informative for downstream analysis
# ------------------------------------------------------------

seurat <- FindVariableFeatures(
  seurat,
  nfeatures = 3000
)

# ------------------------------------------------------------
# Visualize Highly Variable Features
# ------------------------------------------------------------

top_features <- head(VariableFeatures(seurat), 20)

plot1 <- VariableFeaturePlot(seurat)

plot2 <- LabelPoints(
  plot = plot1,
  points = top_features,
  repel = TRUE
)

plot1 + plot2


# ============================================================
# Data Scaling and Regression
# ============================================================

# ------------------------------------------------------------
# Scale Data
# Centers and scales expression values
# ------------------------------------------------------------

seurat <- ScaleData(seurat)

# ------------------------------------------------------------
# Regress Out Unwanted Sources of Variation
# Removes effects of:
# - total detected genes
# - mitochondrial transcript percentage
# ------------------------------------------------------------

seurat <- ScaleData(
  seurat,
  vars.to.regress = c("nFeature_RNA", "percent.mt")
)


# ============================================================
# Dimensionality Reduction
# ============================================================

# ------------------------------------------------------------
# Principal Component Analysis (PCA)
# Linear dimensionality reduction
# ------------------------------------------------------------

seurat <- RunPCA(
  seurat,
  npcs = 50
)

# ------------------------------------------------------------
# Determine Number of Significant PCs
# ------------------------------------------------------------

ElbowPlot(
  seurat,
  ndims = ncol(Embeddings(seurat, "pca"))
)

# ------------------------------------------------------------
# Visualize Genes Contributing to Top PCs
# ------------------------------------------------------------

PCHeatmap(
  seurat,
  dims = 1:20,
  cells = 500,
  balanced = TRUE,
  ncol = 4
)


# ============================================================
# Non-Linear Dimensionality Reduction
# ============================================================

# ------------------------------------------------------------
# Run tSNE and UMAP using top 20 PCs
# ------------------------------------------------------------

seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

# ------------------------------------------------------------
# Visualize tSNE and UMAP Embeddings
# ------------------------------------------------------------

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)

plot1 + plot2


# ============================================================
# Marker Gene Visualization
# ============================================================

# ------------------------------------------------------------
# Visualize Canonical Marker Genes
# Helps identify cell types or cell states
# ------------------------------------------------------------

plot1 <- FeaturePlot(
  seurat,
  features = c(
    "MKI67", "NES", "DCX", "FOXG1",
    "DLX2", "EMX1", "OTX2", "LHX9", "TFAP2A"
  ),
  ncol = 3,
  reduction = "tsne"
)

plot2 <- FeaturePlot(
  seurat,
  features = c(
    "MKI67", "NES", "DCX", "FOXG1",
    "DLX2", "EMX1", "OTX2", "LHX9", "TFAP2A"
  ),
  ncol = 3,
  reduction = "umap"
)

plot1 / plot2


# ============================================================
# Clustering
# ============================================================

# ------------------------------------------------------------
# Construct KNN Graph
# ------------------------------------------------------------

seurat <- FindNeighbors(
  seurat,
  dims = 1:20
)

# ------------------------------------------------------------
# Cluster Cells
# Resolution controls cluster granularity
# ------------------------------------------------------------

seurat <- FindClusters(
  seurat,
  resolution = 1
)

# ------------------------------------------------------------
# Visualize Clusters
# ------------------------------------------------------------

plot1 <- DimPlot(
  seurat,
  reduction = "tsne",
  label = TRUE
)

plot2 <- DimPlot(
  seurat,
  reduction = "umap",
  label = TRUE
)

plot1 + plot2


# ============================================================
# Heatmap of Known Cell-Type Marker Genes
# ============================================================

ct_markers <- c(
  "MKI67", "NES", "DCX", "FOXG1",        # G2M, NPC, neuron, telencephalon
  "DLX2", "DLX5", "ISL1", "SIX3",
  "NKX2.1", "SOX6", "NR2F2",             # ventral telencephalon markers
  "EMX1", "PAX6", "GLI3", "EOMES",
  "NEUROD6",                             # dorsal telencephalon markers
  "RSPO3", "OTX2", "LHX9", "TFAP2A",
  "RELN", "HOXB2", "HOXB5"              # non-telencephalon markers
)

DoHeatmap(
  seurat,
  features = ct_markers
) + NoLegend()


# ============================================================
# Differential Gene Analysis (DGA)
# ============================================================

# ------------------------------------------------------------
# Identify Marker Genes for Each Cluster
# ------------------------------------------------------------

cl_markers <- FindAllMarkers(
  seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = log(1.2)
)

library(dplyr)

# ------------------------------------------------------------
# Extract Top 2 Markers per Cluster
# ------------------------------------------------------------

top_markers <- cl_markers %>%
  group_by(cluster) %>%
  slice_max(
    n = 2,
    order_by = avg_log2FC
  )


# ============================================================
# Heatmap of Top Cluster Markers
# ============================================================

library(dplyr)
library(Seurat)

# ------------------------------------------------------------
# Extract Top 10 Markers per Cluster
# ------------------------------------------------------------

top10_cl_markers <- cl_markers %>%
  group_by(cluster) %>%
  slice_max(
    n = 10,
    order_by = avg_log2FC
  )

# Extract unique gene names
top10_genes <- unique(top10_cl_markers$gene)

# Plot heatmap
DoHeatmap(
  seurat,
  features = top10_genes
) + NoLegend()


# ============================================================
# Detailed Marker Visualization
# ============================================================

# ------------------------------------------------------------
# NEUROD2 and NEUROD6 Marker Visualization
# ------------------------------------------------------------

plot1 <- FeaturePlot(
  seurat,
  c("NEUROD2", "NEUROD6"),
  ncol = 1
)

plot2 <- VlnPlot(
  seurat,
  features = c("NEUROD2", "NEUROD6"),
  pt.size = 0
)

plot1 + plot2 + plot_layout(widths = c(1, 2))


# ------------------------------------------------------------
# LHX9 Marker Visualization
# ------------------------------------------------------------

plot1 <- FeaturePlot(
  seurat,
  c("LHX9"),
  ncol = 1
)

plot2 <- VlnPlot(
  seurat,
  features = c("LHX9"),
  pt.size = 0
)

plot1 + plot2 + plot_layout(widths = c(1, 2))


# ============================================================
# Optional Cell Type Annotation
# ============================================================

library(Seurat)

# Rename cluster identities with biological annotations
new_ident <- setNames(
  c(
    "Dorsal telen. NPC",
    "Midbrain-hindbrain boundary neuron",
    "Dorsal telen. neuron",
    "Dien. and midbrain excitatory neuron",
    "MGE-like neuron",
    "G2M dorsal telen. NPC",
    "Dorsal telen. IP",
    "Dien. and midbrain NPC",
    "Dien. and midbrain IP and excitatory early neuron",
    "G2M Dien. and midbrain NPC",
    "G2M dorsal telen. NPC",
    "Dien. and midbrain inhibitory neuron",
    "Dien. and midbrain IP and early inhibitory neuron",
    "Ventral telen. neuron",
    "Unknown 1",
    "Unknown 2"
  ),
  levels(seurat)
)

# Apply annotations
seurat <- RenameIdents(seurat, new_ident)

# Visualize annotated clusters
DimPlot(
  seurat,
  reduction = "umap",
  label = TRUE
) + NoLegend()


# ============================================================
# Pseudotime / Trajectory Analysis
# ============================================================

# ------------------------------------------------------------
# Subset Dorsal Cell Populations
# ------------------------------------------------------------

seurat_dorsal <- subset(
  seurat,
  subset = RNA_snn_res.1 %in% c(0, 2, 5, 6, 10)
)

# ------------------------------------------------------------
# Identify Highly Variable Genes in Subset
# ------------------------------------------------------------

seurat_dorsal <- FindVariableFeatures(
  seurat_dorsal,
  nfeatures = 2000
)

# ------------------------------------------------------------
# Remove Cell-Cycle Related Genes from HVGs
# ------------------------------------------------------------

VariableFeatures(seurat) <- setdiff(
  VariableFeatures(seurat),
  unlist(cc.genes)
)

# ------------------------------------------------------------
# PCA + UMAP on Subset
# ------------------------------------------------------------

seurat_dorsal <- RunPCA(seurat_dorsal) %>%
  RunUMAP(dims = 1:20)

# Visualize developmental markers
FeaturePlot(
  seurat_dorsal,
  c("MKI67", "GLI3", "EOMES", "NEUROD6"),
  ncol = 4
)


# ============================================================
# Cell Cycle Scoring
# ============================================================

# ------------------------------------------------------------
# Assign Cell Cycle Scores
# ------------------------------------------------------------

seurat_dorsal <- CellCycleScoring(
  seurat_dorsal,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE
)

# ------------------------------------------------------------
# Regress Out Cell Cycle Effects
# ------------------------------------------------------------

seurat_dorsal <- ScaleData(
  seurat_dorsal,
  vars.to.regress = c("S.Score", "G2M.Score")
)

# Recompute PCA and UMAP after regression
seurat_dorsal <- RunPCA(seurat_dorsal) %>%
  RunUMAP(dims = 1:20)

# Visualize markers again
FeaturePlot(
  seurat_dorsal,
  c("MKI67", "GLI3", "EOMES", "NEUROD6"),
  ncol = 4
)

#Save
saveRDS(
  seurat,
  file = file.path(results_dir, "seurat_obj_all.rds")
)

saveRDS(
  seurat_dorsal,
  file = file.path(results_dir, "seurat_obj_dorsal.rds")
)

seurat <- readRDS(
  file.path(results_dir, "seurat_obj_all.rds")
)

seurat_dorsal <- readRDS(
  file.path(results_dir, "seurat_obj_dorsal.rds")
)