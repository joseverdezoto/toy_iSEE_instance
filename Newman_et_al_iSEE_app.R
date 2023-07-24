library(Seurat)
library(readr)
library(iSEE)
library(iSEEu)



# In this script will play around making an iSEE app for the Wirka data and Newman et al data  
newman_scrna = readr::read_rds("~/Desktop/SMC_modulation_project/mouse_scRNA_data/Nat_metabolism_scRNA_data/GFA_2020.01_VSKO-E002-1.rds")

# Convert to sce
newman_sce = Seurat::as.SingleCellExperiment(newman_scrna)

# Run default instance to configure panels
#iSEE::iSEE(newman_sce, appTitle = "Newman et al data")

# Define an empty instance
#iSEEu::modeEmpty(newman_sce)

################################################################################
# Configure panels
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, YAxis = 2L, 
                                          ColorByColumnData = "cell.type", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "orig.ident", 
                                          SizeByColumnData = "nCount_RNA", FacetByRow = "---", FacetByColumn = "---", 
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                          ColorByFeatureName = "Xkr4", ColorByFeatureSource = "---", 
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "AAACCCAAGCAATTCC_1", 
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent", 
                                          SelectionColor = "#5500FF", SelectionAlpha = 0.18, ZoomData = numeric(0), 
                                          BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = c("Color", 
                                                                                                       "Shape", "Size"), ContourAdd = FALSE, ContourColor = "#0000FF", 
                                          PointSize = 0.5, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                          CustomLabels = FALSE, CustomLabelsText = "AAACCCAAGCAATTCC_1", 
                                          FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                          HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                          LabelCentersColor = "#000000", PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                          RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                          ColumnSelectionType = "Saved", ColumnSelectionSaved = 0L, 
                                          SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "cell.type", XAxisFeatureName = "Xkr4", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "Myh11", YAxisFeatureSource = "---", YAxisFeatureDynamicSource = FALSE, 
                                      ColorByColumnData = "orig.ident", ColorByFeatureNameAssay = "logcounts", 
                                      ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "orig.ident", 
                                      SizeByColumnData = "nCount_RNA", FacetByRow = "---", FacetByColumn = "---", 
                                      ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "Xkr4", 
                                      ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                      ColorBySampleName = "AAACCCAAGCAATTCC_1", ColorBySampleSource = "---", 
                                      ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                      SelectionEffect = "Transparent", SelectionColor = "#FF0000", 
                                      SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                      VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
                                      ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                      Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                      CustomLabelsText = "AAACCCAAGCAATTCC_1", FontSize = 1, LegendPointSize = 1, 
                                      LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE, 
                                      LabelCentersBy = "orig.ident", LabelCentersColor = "#000000", 
                                      PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 500L, PanelWidth = 6L, 
                                      SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                      RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                      ColumnSelectionType = "Active", ColumnSelectionSaved = 0L, 
                                      SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "Fn1\nMyh11\nLtbp1\nIbsp", ClusterRows = TRUE, 
                                        ClusterRowsDistance = "euclidean", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("tissue", 
                                                                                                           "OriginalClusters"), RowData = character(0), CustomBounds = FALSE, 
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE, 
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow", 
                                        ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, SelectionEffect = "Color", SelectionColor = "#FF0000", 
                                        PanelId = c(ComplexHeatmapPlot = 1L), PanelHeight = 500L, 
                                        PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                        ColumnSelectionSource = "---", RowSelectionDynamicSource = FALSE, 
                                        RowSelectionType = "Active", RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                        ColumnSelectionType = "Active", ColumnSelectionSaved = 0L, 
                                        SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "SCT_snn_res.0.6", 
                                    XAxisColumnData = "orig.ident", ColorByColumnData = "cell.type", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    FacetByRow = "---", FacetByColumn = "---", ColorBy = "None", 
                                    ColorByDefaultColor = "#000000", ColorByFeatureName = "Xkr4", 
                                    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                    ColorBySampleName = "AAACCCAAGCAATTCC_1", ColorBySampleSource = "---", 
                                    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                    SelectionEffect = "Transparent", SelectionColor = "#FF0000", 
                                    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                    VisualBoxOpen = FALSE, VisualChoices = c("Color", "Facet"
                                    ), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "AAACCCAAGCAATTCC_1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", PanelId = c(ColumnDataPlot = 1L), 
                                    PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L, 
                                    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "Xkr4", Search = "", SearchColumns = "", 
                                  HiddenColumns = character(0), PanelId = c(RowDataTable = 1L), 
                                  PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                  RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                  DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                  RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                  ColumnSelectionType = "Active", ColumnSelectionSaved = 0L, 
                                  SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "AAACCCAAGCAATTCC_1", Search = "", 
                                     SearchColumns = c("", "", "", "", "", "", "", "", "", "", 
                                                       "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0), 
                                     PanelId = c(ColumnDataTable = 1L), PanelHeight = 500L, PanelWidth = 8L, 
                                     SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                     DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                     RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                     ColumnSelectionType = "Active", ColumnSelectionSaved = 0L, 
                                     SelectionHistory = list())

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "logcounts", XAxis = "None", XAxisRowData = "", 
                                     XAxisSampleName = "AAACCCAAGCAATTCC_1", XAxisSampleSource = "---", 
                                     XAxisSampleDynamicSource = FALSE, YAxisSampleName = "AAACCCAAGCAATTCC_1", 
                                     YAxisSampleSource = "---", YAxisSampleDynamicSource = FALSE, 
                                     ColorByRowData = "", ColorBySampleNameAssay = "logcounts", 
                                     ColorByFeatureNameColor = "#FF0000", ShapeByRowData = NA_character_, 
                                     SizeByRowData = "", FacetByRow = "---", FacetByColumn = "---", 
                                     ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "Xkr4", 
                                     ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                     ColorBySampleName = "AAACCCAAGCAATTCC_1", ColorBySampleSource = "---", 
                                     ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                     SelectionEffect = "Transparent", SelectionColor = "#FF0000", 
                                     SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                     VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
                                     ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                     Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                     CustomLabelsText = "Xkr4", FontSize = 1, LegendPointSize = 1, 
                                     LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE, 
                                     LabelCentersBy = "", LabelCentersColor = "#000000", PanelId = c(SampleAssayPlot = 1L), 
                                     PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                     RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                     DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active", 
                                     RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE, 
                                     ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
                                     SelectionHistory = list())

# Run the app with the newly configured panels
iSEE::iSEE(newman_sce, initial = initial, appTitle = "Newman et al scRNA data")





