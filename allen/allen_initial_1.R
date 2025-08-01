library(iSEE)

initial <- list(ReducedDimensionPlot(PanelWidth=6L,
                                     Type="PCA",
                                     ColorBy="Column data",
                                     ColorByColumnData="Core.Type"),
                FeatureAssayPlot(PanelWidth=6L,
                                 XAxis="Column data",
                                 XAxisColumnData="Core.Type",
                                 YAxisFeatureName="Cd80"))
