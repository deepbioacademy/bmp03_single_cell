library(tidyverse)
library(Seurat)
library(Matrix)

if (!dir.exists("data")) {
    download.file(
        url = "https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/course_data.tar.gz",
        destfile = "data.tar.gz"
    )
    untar("data.tar.gz")
    file.remove("data.tar.gz")
}

sample_info <- read.csv("data/sample_info.csv")
data_dirs <- file.path(
    "data",
    "count_matrices",
    sample_info$SampleName,
    "outs",
    "filtered_feature_bc_matrix"
)
names(data_dirs) <- gsub("_", "-", sample_info$SampleName)
data_dirs <- data_dirs[1:3]

sparse_matrix <- Read10X(data.dir = data_dirs)

seu <- CreateSeuratObject(
    counts = sparse_matrix,
    project = "pbmmc",
    min.cells = 3,
    min.features = 100
)

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, pattern = "^RP[SL]", col.name = "percent_ribo")
seu <- PercentageFeatureSet(seu, pattern = "^HB[^(P)]", col.name = "percent_globin")

plot_most_expressed <- function(object, ngenes = 20) {
    cts <- GetAssayData(object, assay = "RNA", layer = "counts")
    cts <- t(cts) / colSums(cts) * 100
    medians <- apply(cts, 2, median)

    most_expressed <- order(medians, decreasing = TRUE)[ngenes:1]
    most_exp_matrix <- as.matrix(cts[, most_expressed])

    most_exp_df <- stack(as.data.frame(most_exp_matrix))
    colnames(most_exp_df) <- c("perc_total", "gene")

    ggplot(most_exp_df, aes(x = gene, y = perc_total)) +
        geom_boxplot() +
        coord_flip() +
        theme_minimal()
}

p_qc_pre <- VlnPlot(
    seu,
    features = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_globin"),
    ncol = 3
)
ggsave("plots/01_qc_pre_filter.png", p_qc_pre, width = 12, height = 8)

p_top_genes <- plot_most_expressed(seu, 20)
ggsave("plots/01_top_expressed_genes.png", p_top_genes, width = 8, height = 6)

seu <- subset(
    seu,
    subset = nFeature_RNA > 200 &
        nFeature_RNA < 5000 &
        percent_mito < 8
)

p_qc_post <- VlnPlot(seu, features = c("nFeature_RNA", "percent_mito"))
ggsave("plots/01_qc_post_filter.png", p_qc_post, width = 8, height = 5)

saveRDS(seu, "results/01_filtered_data.rds")
