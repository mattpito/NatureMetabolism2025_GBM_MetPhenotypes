!!!! Deconvolution is tricky !!!! One should look at the bigger picture and not be fixated on spesific values.
!!!! This code will work reproducibly but keep in mind that there may be small changes in each run !!!


# Load processed spatial data
load(file = "PATH/spatial_data_ALL_filtered")
sessionInfo()

# Extract count and spatial coordinate data to proccess jointly
counts <- spatial_data_ALL_filtered[["Spatial"]]$counts
coords <- bind_rows(lapply(Images(spatial_data_ALL_filtered), function(x) GetTissueCoordinates(spatial_data_ALL_filtered, image = x)))
coords <- coords[!duplicated(paste0(coords$x, coords$y)), ]
colnames(coords) <- c("x", "y")
query <- SpatialRNA(coords, counts, colSums(counts))

# Load public dataset (metadata, raw counts, t-SNE) You can aquire these files form Darmanos et. al
gbm_metadata_public <- read.csv("PATH/GBM_metadata.csv", sep = " ")
gbm_metadata_public_rawcounts <- read.csv("PATH/GBM_raw_gene_counts.csv", sep = " ")
colnames(gbm_metadata_public_rawcounts) <- gsub("^[xX]", "", colnames(gbm_metadata_public_rawcounts))
gbm_metadata_public_TSNEcord <- as.matrix(read.csv("PATH/GBM_TSNE.csv", sep = " "))

# Create Seurat object from public data
seurat_object_Public <- CreateSeuratObject(counts = gbm_metadata_public_rawcounts, assay = "RNA")
seurat_object_Public <- AddMetaData(seurat_object_Public, metadata = gbm_metadata_public)
seurat_object_Public[["tsne"]] <- CreateDimReducObject(embeddings = gbm_metadata_public_TSNEcord, key = "tSNE_", assay = DefaultAssay(seurat_object_Public))

# Normalization, PCA, Clustering, and UMAP
seurat_object_Public <- seurat_object_Public %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.01, group.singletons = FALSE) %>%
  RunUMAP(dims = 1:30)

# Assign annotations
seurat_object_Public$Annotation <- case_when(
  seurat_object_Public$Cluster_2d == 1  ~ "Malignant",
  seurat_object_Public$Cluster_2d == 2  ~ "Oligodendrocytes",
  seurat_object_Public$Cluster_2d == 3  ~ "Vascular",
  seurat_object_Public$Cluster_2d == 4  ~ "Malignant",
  seurat_object_Public$Cluster_2d == 5  ~ "Neurons",
  seurat_object_Public$Cluster_2d == 6  ~ "Vascular",
  seurat_object_Public$Cluster_2d == 12 ~ "Vascular",
  seurat_object_Public$Cluster_2d %in% c(7, 8) ~ "Immune",
  seurat_object_Public$Cluster_2d == 11 ~ "Malignant",
  seurat_object_Public$Cluster_2d == 9  ~ "OPC",
  seurat_object_Public$Cluster_2d == 10  ~ "Astrocytes",
  TRUE ~ NA_character_  # Assign NA to any unexpected values
)

# Plot tSNE and UMAP
DimPlot(seurat_object_Public, reduction = "umap", group.by = "Annotation")
DimPlot(seurat_object_Public, reduction = "tsne", group.by = "Annotation", label = TRUE, label.box = TRUE, label.size = 2)


# Calculate transcriptomic similarity. !!Deconvolution methods are not perfect and may be prone to errors with high transcriptomic similarity.
Idents(seurat_object_Public) <- "Annotation"
avg_expression_matrix <- as.matrix(AverageExpression(seurat_object_Public)$RNA)
similarity_matrix <- cor(avg_expression_matrix, method = "spearman") #Try pearson or kendall to check differences.
pheatmap(similarity_matrix, main = "Transcriptomic Similarity Between Groups")
seurat_object_Public.cp <- seurat_object_Public

# ===============================
# RCTD Deconvolution Preparation
# ===============================
seurat_object_Public <- seurat_object_Public.cp
seurat_object_Public$Annotation2 <- ifelse(seurat_object_Public$Annotation == "Malignant", "Tumour", "TME")
seurat_object_Public$rownamesBarcode <- rownames(seurat_object_Public@meta.data)

# Balance cell sampling. We do not want to have 5000 Immune cells and 10 neurons. We set a threshold to 50. This is not exhaustive. For Tumour Cells, we use all.
ref_metadata <- seurat_object_Public@meta.data
balanced_data <- ref_metadata %>%
  mutate(Annotation2 = ifelse(Annotation == "Malignant", "Tumour", "TME")) %>%
  group_by(Annotation2, Annotation) %>%
  filter(Annotation2 == "Tumour" | (Annotation == "Immune" & row_number() <= 50) | (Annotation2 == "TME" & Annotation != "Immune" & row_number() <= 50)) %>%
  ungroup()

ref_metadata$in_sampled_data <- ifelse(ref_metadata$rownamesBarcode %in% balanced_data$rownamesBarcode, "Yes", "No")
seurat_object_Public$in_sampled_data <- ref_metadata$in_sampled_data

# Subset for RCTD
seurat_annodat_3500cellsEach_for_RCTD <- subset(seurat_object_Public, subset = in_sampled_data == "Yes")
counts <- seurat_annodat_3500cellsEach_for_RCTD[["RNA"]]$counts
cluster <- as.factor(seurat_annodat_3500cellsEach_for_RCTD$Annotation2)
names(cluster) <- colnames(seurat_annodat_3500cellsEach_for_RCTD)
nUMI <- seurat_annodat_3500cellsEach_for_RCTD$nCount_RNA
names(nUMI) <- colnames(seurat_annodat_3500cellsEach_for_RCTD)

# Create reference
reference <- Reference(counts, cluster, nUMI, min_UMI = 43)

# Run RCTD
RCTD_forfull <- create.RCTD(query, reference, max_cores = 10,CELL_MIN_INSTANCE = 20) # This is if using single cell reference. If you adjust the paramteres you will get slightly different results, depending on stringency
RCTD <- run.RCTD(RCTD_forfull, doublet_mode = "full")




# SAVE and LOAD
RCTD <- readRDS("PATH/RCTD_fromPUblic_darmanos.rds")

# Normalize and merge RCTD weights with Seurat spatial object
norm_weights_df <- as.data.frame(as.matrix(RCTD@results$norm_weights)) %>%
  tibble::rownames_to_column("rownames")

simple_weights_df <- as.data.frame(as.matrix(RCTD@results$weights)) %>%
  tibble::rownames_to_column("rownames")
colnames(simple_weights_df)[-1] <- paste0(colnames(simple_weights_df)[-1], "_weight")

spatial_norm <- spatial_data_ALL_filtered
spatial_norm@meta.data <- spatial_norm@meta.data %>%
  tibble::rownames_to_column("rownames") %>%
  left_join(norm_weights_df, by = "rownames") %>%
  left_join(simple_weights_df, by = "rownames") %>%
  tibble::column_to_rownames("rownames")

# Create MALIGNANT composite score. OPCs and Malignant in DARMANOS dataset have high transcriptmic similarity. We use them to block overly aggresive assigment of other TME to spots.
#We do not use the OPC population in inference.
spatial_norm$MALIGNANT <- spatial_norm$OPC + spatial_norm$Malignant













