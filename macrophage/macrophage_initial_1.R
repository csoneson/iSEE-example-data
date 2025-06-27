library(iSEE)
library(iSEEde)
library(iSEEpathways)

initial <- list(
  DETable(ContrastName = "IFNgTRUE.SL1344TRUE.DESeq2",
          HiddenColumns = c("baseMean", "lfcSE", "stat")),
  VolcanoPlot(ContrastName = "IFNgTRUE.SL1344TRUE.DESeq2"),
  MAPlot(ContrastName = "IFNgTRUE.SL1344TRUE.DESeq2"),
  PathwaysTable(ResultName = "IFNgTRUE.SL1344TRUE.limma.fgsea",
                Selected = "GO:0046324"),
  ComplexHeatmapPlot(RowSelectionSource = "PathwaysTable1",
                     CustomRows = FALSE, ColumnData = "condition_name",
                     ClusterRows = TRUE, Assay = "vst"),
  FgseaEnrichmentPlot(ResultName = "IFNgTRUE.SL1344TRUE.limma.fgsea",
                      PathwayId = "GO:0046324")
)
