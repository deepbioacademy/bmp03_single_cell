library(tidyverse)
library(Seurat)

seu <- readRDS("results/01_filtered_data.rds")

seu <- NormalizeData(
    seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000
)

top10 <- head(VariableFeatures(seu), 10)
vf_plot <- VariableFeaturePlot(seu)
p_vf <- LabelPoints(vf_plot, points = top10, repel = TRUE)
ggsave("plots/02_variable_features.png", p_vf, width = 9, height = 6)

seu <- ScaleData(seu)
seu <- RunPCA(seu)

p_elbow <- ElbowPlot(seu, ndims = 40) +
    geom_vline(xintercept = 25, linetype = "dashed", color = "red", linewidth = 1)
ggsave("plots/02_pca_elbow.png", p_elbow, width = 8, height = 5)

seu <- RunUMAP(seu, dims = 1:25, reduction.name = "umap.unintegrated")

p_unint <- DimPlot(seu, reduction = "umap.unintegrated", group.by = "orig.ident")
ggsave("plots/02_umap_unintegrated.png", p_unint, width = 8, height = 6)

saveRDS(seu, "results/02_dimred_data.rds")
