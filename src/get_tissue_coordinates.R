library(Seurat)
library(SeuratData)
library(feather)

get_coordinates <- function(data, image){
  tc <- data@images[[image]]@coordinates
  tc <- cbind(rownames(tc), tc)
  names(tc)[1] <- "barcode"
  return(tc)
}

## Sample 1512
tissue_ids <- c("brain1A", "brain1B", "brain1C", "brain1D")
sample_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics/Sample-1512"
dest_dir <- file.path(sample_dir, "assay_data")

data_path <-
  file.path(sample_dir, "1512_rpca_integrated.rds")

merged_seur_data <- readRDS(data_path)

tcA <- get_coordinates(merged_seur_data, "brain1A")
tcB <- get_coordinates(merged_seur_data, "brain1B")
tcC <- get_coordinates(merged_seur_data, "brain1C")
tcD <- get_coordinates(merged_seur_data, "brain1D")

write_feather(tcA, file.path(dest_dir, "brain1A-s1512-tc.feather"))
write_feather(tcB, file.path(dest_dir, "brain1B-s1512-tc.feather"))
write_feather(tcC, file.path(dest_dir, "brain1C-s1512-tc.feather"))
write_feather(tcD, file.path(dest_dir, "brain1D-s1512-tc.feather"))

## Sample 4774
tissue_ids <- c("brain1A", "brain1B", "brain1C", "brain1D")
sample_dir <- "Rao_Lab/data/GBM_spatial_transcriptomics/Sample-4774"
dest_dir <- file.path(sample_dir, "assay_data")

data_path <-
  file.path(sample_dir, "4774_merged.rds")

merged_seur_data <- readRDS(data_path)

tcA <- get_coordinates(merged_seur_data, "brain1A")
tcB <- get_coordinates(merged_seur_data, "brain1B")
tcC <- get_coordinates(merged_seur_data, "brain1C")
tcD <- get_coordinates(merged_seur_data, "brain1D")

write_feather(tcA, file.path(dest_dir, "brain1A-s4774-tc.feather"))
write_feather(tcB, file.path(dest_dir, "brain1B-s4774-tc.feather"))
write_feather(tcC, file.path(dest_dir, "brain1C-s4774-tc.feather"))
write_feather(tcD, file.path(dest_dir, "brain1D-s4774-tc.feather"))

