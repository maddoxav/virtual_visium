library(Seurat)
library(SeuratData)
library(feather)

# get assay data from seurat object as data frame
assay_to_mat <-
  function(seur_data,
           assay = "Spatial",
           slot = "counts") {
    mat <-
      data.frame(GetAssayData(seur_data, assay = assay, slot = slot))
    mat$symbol <- rownames(mat)
    mat <- mat[, c(ncol(mat), 1:(ncol(mat) - 1))]
    return(mat)
  }

# get single tissue from merged seurat object
get_tissue <- function(merged_seur_data, tiss_id) {
  seur_data <-
    subset(merged_seur_data, subset = orig.ident == tiss_id)
  return(seur_data)
}

# write each sample's assay df
write_data <- function(dest_dir, data_path, data_id, tissue_ids) {
  merged_seur_data <- readRDS(data_path)
  # get data from each tissue in sample
  for (tissue_id in tissue_ids) {
    seur_data <- get_tissue(merged_seur_data, tissue_id)
    # get data from seurat object
    mat_data <- assay_to_mat(seur_data, assay = "SCT", slot = "counts")
    im_dims <- get_img_dims(seur_data, tissue_id)
    # output file name
    md_fname <-
      paste(tissue_id, data_id, "count_mat.feather", sep = "-")
    # output file path
    md_path <- file.path(dest_dir, md_fname)
    # write to file
    write_feather(mat_data, md_path)
  }
}

tissue_ids <- c("brain1A", "brain1B", "brain1C", "brain1D")
# sample 4774
sample_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics/Sample-4774"
dest_dir <- file.path(sample_dir, "assay_data")
data_path <-
  file.path(sample_dir, "4774_merged.rds")
data_id <- "s4774"
write_data(dest_dir, data_path, data_id, tissue_ids)
# sample 1512
sample_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics/Sample-1512"
dest_dir <- file.path(sample_dir, "assay_data")
data_path <-
  file.path(sample_dir, "1512_rpca_integrated.rds")
data_id <- "s1512"
write_data(dest_dir, data_path, data_id, tissue_ids)
