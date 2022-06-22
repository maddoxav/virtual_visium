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
get_tissue_coordinates <- function(seur_data){
  tc <- GetTissueCoordinates(seur_data)
  tc$barcode <- rownames(tc)
  tc <- tc[,c(ncol(tc),1:(ncol(tc)-1))]
  return(tc)
}

# write each sample's assay df (just og spatial counts for now)
# and barcode index to feather files
tissue_ids <- c("brain1A", "brain1B", "brain1C", "brain1D")
data_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics"
# sample 4774
data_path <-
  file.path(data_dir, "Sample-4774/4774_merged.rds")
merged_seur_data <- readRDS(data_path)
for (tissue_id in tissue_ids) {
  seur_data <- get_tissue(merged_seur_data, tissue_id)
  mat_data <- assay_to_mat(seur_data)
  tc <- get_tissue_coordinates(seur_data)
  dest_dir <- file.path(data_dir, "assay_data")
  md_fname <- paste(tissue_id, "s4774", "count_mat.feather", sep="-")
  tc_fname <- paste(tissue_id, "s4774", "tc.feather", sep="-")
  md_path <- file.path(dest_dir, md_fname)
  tc_path <- file.path(dest_dir, tc_fname)
  write_feather(mat_data, md_path)
  write_feather(tc, tc_path)
}
# sample 1512
data_path <-
  file.path(data_dir, "Sample-1512/1512_rpca_integrated.rds")
merged_seur_data <- readRDS(data_path)
for (tissue_id in tissue_ids) {
  seur_data <- get_tissue(merged_seur_data, tissue_id)
  mat_data <- assay_to_mat(seur_data)
  tc <- get_tissue_coordinates(seur_data)
  dest_dir <- file.path(data_dir, "assay_data")
  md_fname <- paste(tissue_id, "s1512", "count_mat.feather", sep="-")
  tc_fname <- paste(tissue_id, "s1512", "tc.feather", sep="-")
  md_path <- file.path(dest_dir, md_fname)
  tc_path <- file.path(dest_dir, tc_fname)
  write_feather(mat_data, md_path)
  write_feather(tc, tc_path)
}
