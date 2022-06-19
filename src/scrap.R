# install.packages("devtools")
#
# if (!base::requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
#
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))
#
# install.packages("Seurat")
# install.packages("Se")
#
# devtools::install_github(repo = "kueckelj/confuns")
# devtools::install_github(repo = "theMILOlab/SPATA2")
#
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')

library(SPATA2)
library(patchwork)

data_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics"
data_304_path = file.path(data_dir, "304_T_SPATA_CNV_Pred.RDS")
data_334_path = file.path(data_dir, "334_T_SPATA_CNV_Pred.RDS")

# load objects
data_304 <- loadSpataObject(data_304_path)
data_334 <- loadSpataObject(data_334_path)

# example genes of interest
genes_a <- c("TUBA1B", "HOPX", "PLP1", "ACTB")
genes_b <- c("CARTPT", "OLIG1", "GFAP", "SYNPR", "HOPX", "CCK")

# plot a cluster feature
p1 <-
  plotSurface(
    object = data_304,
    color_by = "seurat_clusters",
    pt_clrp = "npg",
    pt_size = 1.2
  ) +
  ggplot2::labs(color = "Clusters")

# plot gene expression
p2 <-
  plotSurface(
    object = data_304,
    color_by = "CARTPT",
    pt_size = 1.2,
    pt_clrsp = "magma"
  )

p1 + legendTop() +
  p2 + legendTop()

# plot gene expression (spatially smoothed)
p1 <-
  plotSurface(
    object = data_334,
    color_by = "CARTPT",
    pt_size = 1.2,
    pt_clrsp = "magma"
  )

p2 <-
  plotSurface(
    object = data_334,
    color_by = "CARTPT",
    pt_size = 1.2,
    pt_clrsp = "magma",
    smooth = TRUE,
    smooth_span = 0.2
  )

# combine with patchwork
p1 + legendNone() +
  p2 + legendTop()

# surface plots interactive
plots <- plotSurfaceInteractive(object = data_334)
# save and return desired plots
names(plots)

plots$nCount_spatial +
  ggplot2::labs(title = "Number of spatial counts") +
  legendNone() +
  plots$pred_tumor +
  ggplot2::labs(title = "Predicted tumor probability")


##### Seurat pipeline demo
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# data_1512 is a Seurat object "brain", 3.3 GB
data_1512_path <-
  file.path(data_dir, "Sample-1512", "1512_rpca_integrated.rda")
load(data_1512_path)

# data_4774 is a Seurat object "brain.merge", 5.1 GB
data_4774_path <-
  file.path(data_dir, "Sample-4774", "4774_merged.rda")
load(data_4774_path)

plot1 <-
  VlnPlot(brain,
          idents = "0",
          features = "nCount_Spatial",
          pt.size = 0.1) + NoLegend()
plot2 <-
  SpatialFeaturePlot(brain, features = "nCount_Spatial", images = "brain1A") + theme(legend.position =  "right")

wrap_plots(plot1, plot2)

# normalize counts (variance is both technical and dependent on anatomy)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# gene expression visualization
SpatialFeaturePlot(brain,
                   features = c("CALB2", "CAMK2N1", "NRGN"),
                   images = "brain1A")

# dimensionality reduction, clustering and visualization
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <-
  SpatialDimPlot(brain,
                 label = TRUE,
                 label.size = 3,
                 images = "brain1A")
p1 + p2

SpatialDimPlot(
  brain,
  images = "brain1A",
  cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)),
  facet.highlight = TRUE,
  ncol = 3
)

# identify spatially variable features without prior annotation
brain <-
  FindSpatiallyVariableFeatures(
    brain,
    assay = "SCT",
    features = VariableFeatures(brain)[1:1000],
    selection.method = "markvariogram"
  )

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, images = "brain1A", features = top.features, ncol = 3, alpha = c(0.1, 1))

