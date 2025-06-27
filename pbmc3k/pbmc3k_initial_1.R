library(iSEE)
library(iSEEu)

initial <- list(
  ReducedDimensionPlot(PanelWidth=4L,
                       Type="TSNE",
                       ColorBy="Column data",
                       ColorByColumnData="labels_main"),
  DynamicMarkerTable(),
  FeatureAssayPlot(PanelWidth=4L,
                   XAxis="Column data",
                   XAxisColumnData="labels_main",
                   YAxisFeatureName="CD79A")
)
