# Set identities for plotting
Idents(spatial_data_RemoveMet) <- spatial_data_RemoveMet$CustomClusters   

# Initial spatial plot (with small legend)
SpatialDimPlot(spatial_data_RemoveMet, ncol = 3, label = TRUE, label.size = 1.5, pt.size.factor = 3.5, repel = TRUE) & 
  NoLegend() + 
  theme(legend.position = "right", legend.text = element_text(size = 8), legend.key.size = unit(0.2, 'cm'))

# Subset to remove unwanted cluster (e.g., 557728_R1_P2 should be excluded as is non-GBM)
spatial_data_RemoveMet <- subset(x = spatial_data_RemoveMet, subset = CustomClusters != "557728_R1_P2")

# Confirm subset worked
SpatialDimPlot(spatial_data_RemoveMet, ncol = 3, label = TRUE, label.size = 1.5, pt.size.factor = 3.5, repel = TRUE) & 
  NoLegend()

# Genes to keep — for kmeans (from Hallmark list, both Glycolysis and Oxphos)
genes2keep <- c("ABCB7", "ACAA1", "ACAA2", ..., "ZNF292")  # abreviated list, check script for UCELL score, or use your own sets.

# Extract expression matrix from Seurat object
counts <- GetAssayData(spatial_data_RemoveMet, assay = "SCT")
dim(counts)
counts <- counts[rownames(counts) %in% genes2keep, ]
counts_matrix <- as.matrix(counts)
counts_matrix <- t(counts_matrix)  # Transpose to have cells as rows

# K-means clustering
k <- kmeans(counts_matrix, centers = 3)

# Plot cluster centers
plot(k$centers[1,], type = 'l', lwd = 2, col = 'red2')
lines(k$centers[2,], type = 'l', lwd = 2, col = 'royalblue')
table(k$cluster)

# Save clustering result
df <- data.frame(cell_barcodes = names(k$cluster), kmeans_cluster = k$cluster)

# Assign back to full Seurat object
spatial_data_ALL_filtered$kmeans_clust3_updated <- NA
matching_indices <- match(rownames(spatial_data_ALL_filtered@meta.data), rownames(df))
spatial_data_ALL_filtered$kmeans_clust3_updated <- df$kmeans_cluster[matching_indices]

##### Plotting results #####
Idents(spatial_data_ALL_filtered) <- "kmeans_clust3_updated"
SpatialDimPlot(spatial_data_ALL_filtered, ncol = 3)

# Feature scatter plot with high confidence interval of the regression line "lm"
FeatureScatter(spatial_data_ALL_filtered, 
               feature1 = "Hallmark_Oxphos_UCell", 
               feature2 = "Hallmark_Glyco_UCell", 
               cells = rownames(df)) + 
  geom_smooth(method = "lm", se = TRUE, level = 0.999999999)

# Create interpretable group labels based on kmeans cluster
# 1, 2 , 3 MAY NOT BE THE SAME WITH YOUR RUN, SINCE THE PROCCESS IS STOCHASTIC WITHOUT A SEED!!!
spatial_data_ALL_filtered@meta.data <- spatial_data_ALL_filtered@meta.data %>%
  mutate(kmeans3_groups_updated = case_when(
    kmeans_clust3_updated == 2 ~ "LowBoth",
    kmeans_clust3_updated == 3 ~ "HighGlyco",
    kmeans_clust3_updated == 1 ~ "HighOxphos",
    TRUE ~ as.character(kmeans_clust3_updated)
  ))

!K-means results will differ in every run (unless you set up a seed). Small differences are OK, the big picture is unchanged. 



