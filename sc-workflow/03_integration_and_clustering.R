library(tidyverse)
library(Seurat)
library(clustree)

seu <- readRDS("results/02_dimred_data.rds")

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

seu <- IntegrateLayers(
    object = seu,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
)

seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:25)
seu <- FindClusters(seu, resolution = seq(0.1, 0.8, by = 0.1))
seu <- RunUMAP(
    seu,
    reduction = "integrated.cca",
    dims = 1:25,
    reduction.name = "umap.cca"
)

p_mix <- DimPlot(seu, reduction = "umap.cca", group.by = "orig.ident", shuffle = TRUE)
ggsave("plots/03_umap_integrated_mixing.png", p_mix, width = 8, height = 6)

p_markers <- FeaturePlot(
    seu,
    reduction = "umap.cca",
    features = c("CD3E", "MS4A1", "CD14", "HBB", "CD34", "NKG7"),
    ncol = 3
)
ggsave("plots/03_canonical_markers.png", p_markers, width = 12, height = 8)

p_tree <- clustree(
    seu@meta.data[, grep("RNA_snn_res", colnames(seu@meta.data))],
    prefix = "RNA_snn_res."
)
ggsave("plots/03_clustree.png", p_tree, width = 10, height = 10)

saveRDS(seu, "results/03_clustered_data.rds")
