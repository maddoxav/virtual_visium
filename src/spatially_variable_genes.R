library(Seurat)
library(SeuratData)
library(feather)
library(ggplot2)

get_save_spatial_vars <-
  function(var_genes_dir = "data/spatially_variable_genes",
           sample_dir = "data/Sample-4774",
           sample_name = "4774_merged",
           sample_id = "s4774",
           n = 100,
           assay = "SCT") {
    sample_path <- file.path(sample_dir,
                             paste(sample_name, ".rds", sep = ""))
    merged_seur_data <- readRDS(sample_path)
    var_feats <- VariableFeatures(merged_seur_data)[1:n]
    merged_seur_data <-
      FindSpatiallyVariableFeatures(
        merged_seur_data,
        asssay = assay,
        features = var_feats,
        selection.method = "markvariogram"
      )
    dest_path <- file.path(sample_dir,
                           paste(sample_name, "-spatvar", n, ".rds", sep = ""))
    var_genes_path <- file.path(var_genes_dir,
                                paste("t", n, "-spat_genes", sample_id, ".txt", sep = ""))
    spat_genes <-
      SpatiallyVariableFeatures(merged_seur_data, assay = "SCT", selection.method = "markvariogram")
    cat(spat_genes, file = var_genes_path, sep = "\n")
    return(merged_seur_data)
  }


plot_spatial_vars <-
  function(merged_seur_data,
           tissue_ids,
           sample_id,
           var_genes_dir) {
    spat_feats <-
      SpatiallyVariableFeatures(merged_seur_data, selection.method = "markvariogram")
    for (tissue_id in tissue_ids) {
      for (feat_idx in 1:10) {
        feat <- spat_feats[feat_idx]
        plt_file <-
          paste(paste(feat, tissue_id, sample_id, sep = "-"),
                ".pdf", sep = "")
        plt_path <- file.path(var_genes_dir, plt_file)
        plt <- SpatialFeaturePlot(
          merged_seur_data,
          images = tissue_id,
          features = feat,
          pt.size.factor = 1.6,
        )
        ggsave(
          plt_path,
          plot = plt,
          device = "pdf",
          width = 8,
          height = 4
        )
      }
    }
  }

## main
data_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics"
var_genes_dir <- file.path(data_dir, "spatially_variable_genes")
# sample 1512
# get/save top spatial vars
sample_dir <- file.path(data_dir, "Sample-1512")
sample_name <- "1512_rpca_integrated"
sample_id <- "s1512"
n <- 100
merged_seur_data <-
  get_save_spatial_vars(var_genes_dir, sample_dir, sample_name, sample_id, 100)
# plot top spatial vars
tissue_ids <- unique(merged_seur_data@meta.data$orig.ident)
plot_spatial_vars(merged_seur_data, tissue_ids, sample_id, var_genes_dir)
