library(scRNAseq)
library(scater)
library(iSEE)
library(iSEEindex)
library(BiocFileCache)

sce <- ReprocessedAllenData(assays="tophat_counts")
sce <- logNormCounts(sce, exprs_values="tophat_counts")
sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)

saveRDS(sce, file="allen_sce.rds")

# Test
bfc <- BiocFileCache(cache = tempdir())

dataset_fun <- function() {
  list(list(id="ReprocessedAllenData",
            title="ReprocessedAllenData",
            uri="https://raw.githubusercontent.com/csoneson/iSEE-example-data/refs/heads/main/allen/allen_sce.rds",
            description="Reprocessed Allen Data.\n"))
}
initial_fun <- function() {
  list(list(id="ReprocessedAllenData_Config1",
            title="InitialConfig1",
            dataset="ReprocessedAllenData",
            uri="https://raw.githubusercontent.com/csoneson/iSEE-example-data/refs/heads/main/allen/allen_initial_1.R",
            description="PCA plot + feature assay plot"),
       list(id="ReprocessedAllenData_Config2",
            title="InitialConfig2",
            dataset="ReprocessedAllenData",
            uri="https://raw.githubusercontent.com/csoneson/iSEE-example-data/refs/heads/main/allen/allen_initial_2.R",
            description="tSNE plot + column data table + feature assay plot"))
}
iSEEindex(bfc, FUN.datasets=dataset_fun, FUN.initial=initial_fun)
