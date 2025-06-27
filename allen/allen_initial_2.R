library(iSEE)

initial <- list(ReducedDimensionPlot(PanelWidth=4L,
                                     Type="TSNE",
                                     ColorBy="Column data",
                                     ColorByColumnData="Core.Type"),
                ColumnDataTable(),
                FeatureAssayPlot(PanelWidth=4L,
                                 XAxis="Column data",
                                 XAxisColumnData="Core.Type",
                                 YAxisFeatureName="Cd80"))
