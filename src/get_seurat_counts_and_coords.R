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

# create data frame of tissue coordinates
get_tissue_coordinates <- function(merged_seur_data) {
  tc <- GetTissueCoordinates(merged_seur_data)
  tc$barcode <- rownames(tc)
  tc <- tc[, c(ncol(tc), 1:(ncol(tc) - 1))]
  return(tc)
}

# get tissue image dimensions
get_img_dims <- function(seur_data, img_name) {
  img <- seur_data@images[[img_name]]@image
  return(dim(img)[c(1, 2)])
}

# write each sample's assay df
# and barcode index to feather files, and img dims to txt
write_data <- function(dest_dir, data_path, data_id, tissue_ids) {
  merged_seur_data <- readRDS(data_path)
  # get tissue coordinates
  tc <- get_tissue_coordinates(merged_seur_data)
  # get data from each tissue in sample
  for (tissue_id in tissue_ids) {
    seur_data <- get_tissue(merged_seur_data, tissue_id)
    # get data from seurat object
    mat_data <- assay_to_mat(seur_data, assay = "SCT", slot = "counts")
    im_dims <- get_img_dims(seur_data, tissue_id)
    # output file names
    md_fname <-
      paste(tissue_id, data_id, "count_mat.feather", sep = "-")
    tc_fname <- paste(data_id, "tc.feather", sep = "-")
    im_dim_fname <-
      paste(tissue_id, data_id, "im_dims.txt", sep = "-")
    # output file paths
    md_path <- file.path(dest_dir, md_fname)
    tc_path <- file.path(dest_dir, tc_fname)
    im_dim_path <- file.path(dest_dir, im_dim_fname)
    # write to files
    write_feather(mat_data, md_path)
    write_feather(tc, tc_path)
    cat(im_dims, file=im_dim_path, sep=" ")
  }
}

tissue_ids <- c("brain1A", "brain1B", "brain1C", "brain1D")
sample_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics/Sample-4774"
dest_dir <- file.path(sample_dir, "assay_data")
# sample 4774
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
