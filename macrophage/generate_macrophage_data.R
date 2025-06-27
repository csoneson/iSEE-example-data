library(SingleCellExperiment)
library(iSEE)
library(iSEEindex)
library(macrophage)
library(tximeta)
library(stringr)
library(org.Hs.eg.db)
library(scater)
library(DESeq2)
library(iSEEde)
library(iSEEpathways)
library(limma)
library(edgeR)
library(fgsea)

dir <- system.file("extdata", package = "macrophage")
coldata <- read.csv(file.path(dir, "coldata.csv"))[, c(1, 2, 3, 5)]
coldata$IFNg <- as.character(grepl("IFNg", coldata$condition_name))
coldata$SL1344 <- as.character(grepl("SL1344", coldata$condition_name))
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
se <- tximeta(coldata = coldata, type = "salmon", dropInfReps = TRUE)

seg <- summarizeToGene(se)
rownames(seg) <- str_replace(rownames(seg), "\\.\\d+$", "")
seg <- addIds(seg, "SYMBOL")
seg <- addIds(seg, "GOALL", multiVals = "list")
rownames(seg) <- scater::uniquifyFeatureNames(
  ID = rownames(seg), names = rowData(seg)$SYMBOL
)
dds <- DESeqDataSet(seg, design = ~ IFNg * SL1344)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
pca <- DESeq2::plotPCA(vst, intgroup = "condition_name", returnData = TRUE)
sce <- as(dds, "SingleCellExperiment")
assay(sce, "vst") <- assay(vst)
SingleCellExperiment::reducedDim(sce, "PCA") <- pca[, c("PC1", "PC2")]
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, name = "IFNgTRUE.SL1344TRUE",
                       lfcThreshold = 0)
sce <- embedContrastResults(res,
                            sce,
                            name = "IFNgTRUE.SL1344TRUE.DESeq2")
dge <- tximeta::makeDGEList(seg)
dge <- dge[rownames(dds), ]
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
design <- model.matrix(~ IFNg * SL1344, data = dge$samples)
fit <- limma::lmFit(logCPM, design = design)
fit <- eBayes(fit, trend = TRUE)
tt <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
sce <- embedContrastResults(tt,
                            sce,
                            name = "IFNgTRUE.SL1344TRUE.limma",
                            class = "limma")
pathways <- select(org.Hs.eg.db, keys(org.Hs.eg.db, "SYMBOL"), c("GOALL"),
                   keytype = "SYMBOL")
pathways <- subset(pathways, ONTOLOGYALL == "BP")
pathways <- unique(pathways[, c("SYMBOL", "GOALL")])
pathways <- split(pathways$SYMBOL, pathways$GOALL)
len_pathways <- lengths(pathways)
pathways <- pathways[len_pathways > 15 & len_pathways < 200]
feature_stats <- res$stat
names(feature_stats) <- rownames(res)
set.seed(42)
fgseaRes <- fgsea(pathways = pathways,
                  stats = feature_stats,
                  minSize = 15,
                  maxSize = 200)
fgseaRes <- fgseaRes[order(pval), ]
sce <- embedPathwaysResults(fgseaRes,
                            sce,
                            name = "IFNgTRUE.SL1344TRUE.DESeq2.fgsea",
                            class = "fgsea",
                            pathwayType = "GO",
                            pathwaysList = pathways,
                            featuresStats = feature_stats)
go_details <- function(x) {
  info <- select(GO.db, x, c("TERM", "ONTOLOGY", "DEFINITION"), "GOID")
  html <- list(p(strong(info$GOID), ":", info$TERM, paste0("(", info$ONTOLOGY, ")")))
  if (!is.na(info$DEFINITION)) {
    html <- append(html, list(p(info$DEFINITION)))
  }
  tagList(html)
}

## Define the mapping from GO terms to gene IDs
map_GO <- function(pathway_id, se) {
  pathway_symbol <- mapIds(org.Hs.eg.db, pathway_id, "SYMBOL",
                           keytype = "GOALL", multiVals = "CharacterList")[[pathway_id]]
  pathway_rownames <- rownames(se)[rowData(se)$SYMBOL %in% pathway_symbol]
  pathway_rownames
}
sce <- registerAppOptions(sce, Pathways.map.functions = list(GO = map_GO))
sce <- registerAppOptions(sce, PathwaysTable.select.details = go_details)

saveRDS(sce, "macrophage_sce.rds")
