if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}

pak::pkg_install(c(
    "tidyverse",
    "Seurat",
    "patchwork",
    "Matrix",
    "celldex",
    "SingleR",
    "edgeR",
    "limma",
    "scuttle",
    "presto",
    "clusterProfiler",
    "ggiraph",
    "enrichplot",
    "clustree",
    "dittoSeq"
))
