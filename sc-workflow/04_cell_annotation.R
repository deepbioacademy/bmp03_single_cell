library(tidyverse)
library(Seurat)
library(celldex)
library(SingleR)
library(dittoSeq)

seu <- readRDS("results/03_clustered_data.rds")

Idents(seu) <- seu$RNA_snn_res.0.3
DefaultAssay(seu) <- "RNA"

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

seu <- AddModuleScore(
    seu,
    features = list(tcell_genes, monocyte_genes),
    name = c("tcell_score", "monocyte_score")
)

s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

seu <- CellCycleScoring(
    seu,
    s.features = s_genes,
    g2m.features = g2m_genes
)

p_phase <- DimPlot(seu, reduction = "umap.cca", group.by = "Phase")
ggsave("plots/04_cell_cycle_phase.png", p_phase, width = 8, height = 6)

ref <- celldex::NovershternHematopoieticData()
seu_SingleR <- SingleR(
    test = GetAssayData(seu, layer = "data"),
    ref = ref,
    labels = ref$label.main
)

singleR_labels <- seu_SingleR$labels
label_counts <- table(singleR_labels)
rare_labels <- names(label_counts)[label_counts < 10]
singleR_labels[singleR_labels %in% rare_labels] <- "Other"
seu$SingleR_annot <- singleR_labels

p_annot <- dittoDimPlot(seu, var = "SingleR_annot", reduction.use = "umap.cca", size = 0.7)
ggsave("plots/04_singler_annotations.png", p_annot, width = 9, height = 6)

saveRDS(seu, "results/04_annotated_data.rds")
