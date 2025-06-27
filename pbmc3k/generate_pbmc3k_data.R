library(scater)
library(iSEE)
library(iSEEindex)
library(BiocFileCache)
library(TENxPBMCData)
library(scran)
library(BiocSingular)
library(SingleR)
library(celldex)

sce <- TENxPBMCData(dataset = "pbmc3k")
colnames(sce) <- paste0("Cell", seq_len(ncol(sce)))
rownames(sce) <- scater::uniquifyFeatureNames(
  ID = rowData(sce)$ENSEMBL_ID,
  names = rowData(sce)$Symbol_TENx
)
MT <- rownames(sce)[grep("^MT-", rownames(sce))]
sce <- scater::addPerCellQC(sce, subsets = list(MT = MT))
sce <- scater::addPerFeatureQC(sce)
sce$log10_total <- log10(sce$total)
rowData(sce)$n_cells <- as.integer(rowData(sce)$detected * ncol(sce))
rowData(sce)$log10_total <- log10(rowSums(assay(sce, "counts")) + 1)
sce <- sce[, sce$subsets_MT_percent < 5]
assay(sce, "counts") <- as(assay(sce, "counts"), "sparseMatrix") # dirty fix for DockerHub
sce <- scran::computeSumFactors(sce, min.mean = 0.1)
sce <- scater::logNormCounts(sce)
dec <- scran::modelGeneVar(sce)
top.dec <- dec[order(dec$bio, decreasing = TRUE), ]
set.seed(1000)
sce <- scran::denoisePCA(sce, technical = dec)
set.seed(1000)
sce <- scater::runTSNE(sce, dimred = "PCA", perplexity = 30)
sce <- scater::runUMAP(sce, dimred = "PCA")
snn.gr <- scran::buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
markers <- scran::findMarkers(sce, groups = sce$Cluster,
                              test.type = "t",
                              direction = "up", pval.type = "all")
for (i in names(markers)) {
  rowData(sce)[, paste0("FDR_cluster", i)] <-
    markers[[i]]$FDR[match(rownames(sce),
                           rownames(markers[[i]]))]
}
ref_monaco <- MonacoImmuneData()
pred_monaco_main <- SingleR(test = sce, ref = ref_monaco, labels = ref_monaco$label.main)
sce$labels_main <- pred_monaco_main$labels
pred_monaco_fine <- SingleR(test = sce, ref = ref_monaco, labels = ref_monaco$label.fine)
sce$labels_fine <- pred_monaco_fine$labels
pred_monaco_ont <- SingleR(test = sce, ref = ref_monaco, labels = ref_monaco$label.ont)
sce$labels_ont <- pred_monaco_ont$labels

saveRDS(sce, file="pbmc3k_sce.rds")

# Test
bfc <- BiocFileCache(cache = tempdir())

dataset_fun <- function() {
  list(list(id="TENxPBMC3k",
            title="TENxPBMC3k",
            uri="https://raw.githubusercontent.com/csoneson/iSEE-example-data/refs/heads/main/pbmc3k/pbmc3k_sce.rds",
            description="10x PBMC3k Data.\n"))
}
initial_fun <- function() {
  list(list(id="TENxPBMC3k_Config1",
            title="InitialConfig1",
            dataset="TENxPBMC3k",
            uri="https://raw.githubusercontent.com/csoneson/iSEE-example-data/refs/heads/main/pbmc3k/pbmc3k_initial_1.R",
            description="Config 1"))
}
iSEEindex(bfc, FUN.datasets=dataset_fun, FUN.initial=initial_fun)
