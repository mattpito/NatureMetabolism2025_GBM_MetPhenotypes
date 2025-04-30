

############# The following chunk most likely will not be needed and is tailored to our HPC configuration ##################
library(curl)
# Load required libraries
lib_hdf5 <- "INSERT_PATH/libhdf5_hl.so.310"
lib_hdf5_main <- "INSERT_PATH/libhdf5.so.310"
dyn.load(lib_hdf5_main)
dyn.load(lib_hdf5)
#############################################################################################################################




# Core analysis libraries
library(hdf5r)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(BPCells)


# Utilities and visualization
library(ape)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scCustomize)
library(stringr)
library(tidyverse)
library(cowplot)
library(patchwork)

# --- 1) Define the correct base directory and discover cases ----
base_path <- "INSERT_PATH/Visium"
# find each immediate subdirectory under base_path (e.g. "ID-BRA5-FO-1", etc)   
case_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)

# prepare list to hold each slide’s Seurat object
spatial_obj_list <- list()

# --- 2) Loop through each case directory ----
for (case_dir in case_dirs) {
  case_name   <- basename(case_dir)                # e.g. "ID-BRA5-FO-1"
  case_prefix <- gsub("-", "_", case_name)         # e.g. "ID_BRA5_FO_1"
  
  # find spaceranger output folders within this case
  spatial_dirs <- Sys.glob(file.path(case_dir, "spaceranger*"))
  if (length(spatial_dirs) == 0) {
    warning("No spaceranger* folders in ", case_dir)
    next
  }
  
  # assign slice labels sequentially (A_, B_, C_, …)
  slice_letters <- LETTERS[seq_along(spatial_dirs)]
  
  for (i in seq_along(spatial_dirs)) {
    sd  <- spatial_dirs[i]                          # full path
    sl  <- paste0(slice_letters[i], "_")            # "A_", "B_", …
    spr <- basename(sd)                             # e.g. "spaceranger130_count_45132_…"
    sid <- paste0(case_prefix, "#", spr)            # e.g. "ID_BRA5_FO_1#spaceranger130_…"
    
    message("Loading: ", sid, " (slice=", sl, ")")
    
    # load + normalize
    so <- Load10X_Spatial(sd, slice = sl)
    so <- SCTransform(so, assay = "Spatial", verbose = FALSE)
    
    # rename barcodes: case#spaceranger#originalBarcode
    old_bc   <- colnames(so)
    new_bc   <- paste0(sid, "#", old_bc)
    colnames(so) <- new_bc
    
    # store slice label
    so$orig.sample <- sl
    
    # stash in list
    spatial_obj_list[[sid]] <- so
  }
}

# --- 3) Merge everything into one Seurat object ----
merged_visium <- merge(
  x       = spatial_obj_list[[1]],
  y       = spatial_obj_list[-1],
  project = "Merged_Visium"
)
merged_visium

DefaultAssay(merged_visium) <- "SCT"
VariableFeatures(merged_visium) <- c(VariableFeatures(spatial_obj_list[[1]]),
                                     VariableFeatures(spatial_obj_list[[2]]),
                                     VariableFeatures(spatial_obj_list[[3]]),
                                     VariableFeatures(spatial_obj_list[[4]]),
                                     VariableFeatures(spatial_obj_list[[5]]),
                                     VariableFeatures(spatial_obj_list[[6]]))

merged_visium <- RunPCA(merged_visium, verbose = FALSE)
merged_visium <- FindNeighbors(merged_visium, dims = 1:30)
merged_visium <- FindClusters(merged_visium, verbose = FALSE, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 2))
merged_visium <- RunUMAP(merged_visium, dims = 1:30)

#############################################
#############################################
#############################################

merged_visium$orig.ident
merged_visium$DonorID <- sub("#.*", "", rownames(merged_visium@meta.data))
merged_visium$DonorID
merged_visium$DonorID_rep <- sub("^[^#]*#([^#]*).*", "\\1", rownames(merged_visium@meta.data))
merged_visium$DonorID_rep
merged_visium$SCT_snn_res.0.05

# saveRDS(merged_visium, file = "INSERT_PATH/filename.rds")
