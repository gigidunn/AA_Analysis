# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(Seurat)

################# Loading data & creating Seurat object #########################

# Define the data directory and subfolders for replicates
data_dir <- "C:/Users/DG1/Desktop/DALLAB/Experimenting/Data/Antibiotic resistance/BacDrop/RawData"
rep1_dir <- file.path(data_dir, "Replicate1", "og_names")
rep2_dir <- file.path(data_dir, "Replicate2")

# Collect all .tsv files from both replicate folders
rep1_files <- data.frame(file_path = list.files(path = rep1_dir, pattern = "\\.tsv$", full.names = TRUE)) %>%
  mutate(replicate = "Replicate1")

rep2_files <- data.frame(file_path = list.files(path = rep2_dir, pattern = "\\.tsv$", full.names = TRUE)) %>%
  mutate(replicate = "Replicate2")

# Combine into one data frame
tsv_files <- bind_rows(rep1_files, rep2_files)

# Add treatment and combined replicate_treatment labels
tsv_files <- tsv_files %>%
  mutate(
    treatment = case_when(
      grepl("P2", file_path) ~ "meropenem",
      grepl("P3", file_path) ~ "ciprofloxacin",
      grepl("P4", file_path) ~ "gentamicin",
      TRUE ~ NA_character_
    ),
    replicate_treatment = paste0(replicate, "_", treatment)
  )

# Read all files into a list of Seurat objects
data_list <- list()
for (file_num in seq_along(tsv_files$file_path)) {
  file <- tsv_files$file_path[file_num]
  message("Reading: ", file)
  
  data_table <- data.table::fread(file, sep = '\t', header = TRUE, nThread = 6, showProgress = TRUE)
  gene_names <- data_table[[1]]
  data_mat <- as.matrix(data_table[, -1])
  rownames(data_mat) <- gene_names
  
  sparse_mat <- Matrix::Matrix(data_mat, sparse = TRUE)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = sparse_mat, min.cells = 1, min.features = 1)
  
  # Add all 3 metadata columns
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(
      treatment = tsv_files$treatment[file_num],
      replicate = tsv_files$replicate[file_num],
      replicate_treatment = tsv_files$replicate_treatment[file_num]
    )
  
  data_list[[file_num]] <- seurat_obj
  remove(data_table, data_mat)
}

# Merge Seurat objects with unique cell IDs based on replicate_treatment
merged_seurat <- merge(
  x = data_list[[1]],
  y = data_list[-1],
  add.cell.ids = tsv_files$replicate_treatment
)

# Output summary
dim(merged_seurat)
table(merged_seurat@meta.data$treatment)
table(merged_seurat@meta.data$replicate)
table(merged_seurat@meta.data$replicate_treatment)


# Join layers (also changing variable name)
data_combined = merged_seurat
rm(merged_seurat)

data_combined[["RNA"]] = JoinLayers(data_combined[["RNA"]])
data_combined
table(data_combined@meta.data$project_id)

#save combined seurat object for next time
saveRDS(data_combined, "C:/Users/DG1/Desktop/DALLAB/Experimenting/Data/Antibiotic resistance/BacDrop/RawData/trial 2/mero_cip_gent_rep1_seurat_obj_250703_v0907.rds")


################# Visualize data ##########################

# Reload the object
data_combined <- readRDS("C:/Users/DG1/Desktop/DALLAB/Experimenting/Data/Antibiotic resistance/BacDrop/RawData/trial 2/mero_cip_gent_rep1_seurat_obj_250703_v0907.rds")
dim(data_combined)

# Visualize metrics using Seurat
Seurat::VlnPlot(data_combined, features = c("nFeature_RNA", "nCount_RNA"), group.by = "treatment", ncol = 2)
Seurat::FeatureScatter(data_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "treatment", raster = FALSE)

plot2 <- Seurat::FeatureScatter(data_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "treatment", raster = FALSE)
plot2 +
  ggplot2::xlim(0, 100) +
  ggplot2::ylim(0, 100)

# Check assay layers
head(data_combined@assays)

######### Subset to cells with ≥15 genes #########
data1_subset <- subset(data_combined, subset = nFeature_RNA >= 15)
dim(data1_subset)

# QC stats before and af"ter filtering
print ("Before filtering:")
mean(data_combined@meta.data$nFeature_RNA)
median(data_combined@meta.data$nFeature_RNA)
median(data_combined@meta.data$nCount_RNA)

print ("Removing Cells with <15 genes:")
mean(data1_subset@meta.data$nFeature_RNA)
median(data1_subset@meta.data$nFeature_RNA)
median(data1_subset@meta.data$nCount_RNA)

######## Removing cells with abnormally high numbers of mRNA detected ##############

# Calculate mean and standard deviation
mean_count <- mean(data1_subset$nCount_RNA)
sd_count <- sd(data1_subset$nCount_RNA)

# Set cutoff (e.g., 3x SD above the mean)
cutoff_high <- mean_count + 3 * sd_count  # You can also try 5 * sd_count

# Print cutoff for reference
cat("Mean:", mean_count, "\nSD:", sd_count, "\nCutoff (mean + 3*SD):", cutoff_high, "\n")

# Filter the dataset
data2_subset <- subset(data1_subset, subset = nCount_RNA <= cutoff_high)

# Show before/after
cat("Cells before filtering:", ncol(data1_subset), "\n")
cat("Cells after filtering:", ncol(data2_subset), "\n")


############ Normalize, Variable Features, PCA, Clustering, UMAP ############

# Global Normalisation - using data1_subset to not include the high gene removal
data <- Seurat::NormalizeData(data1_subset, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
top2000 <- head(VariableFeatures(data), 2000)

plot1 <- Seurat::VariableFeaturePlot(data)
plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale
all.genes <- rownames(data)
data <- Seurat::ScaleData(data, features = all.genes)

# PCA
data <- Seurat::RunPCA(data, features = top2000, npcs = 100)
Seurat::ElbowPlot(data, ndims = 50)
print(data[["pca"]], dims = 1:5, nfeatures = 5)
Seurat::VizDimLoadings(data, dims = 1:4, reduction = "pca")
Seurat::DimPlot(data, reduction = "pca", group.by = "treatment", raster = FALSE)

# Color by counts on PCA
# PCA colored by counts
pca_df <- Seurat::FetchData(data, vars = c("PC_1", "PC_2", "nCount_RNA", "nFeature_RNA", "treatment"))
ggplot(pca_df, aes(x = PC_1, y = PC_2, color = nCount_RNA)) +
  geom_point(size = 1, alpha = 0.3) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal()

# UMAP and Clustering
data <- Seurat::FindNeighbors(data, dims = 1:10, verbose = FALSE)
data <- Seurat::FindClusters(data, resolution = 0.5, verbose = FALSE)
data <- Seurat::RunUMAP(data, dims = 1:10, verbose = FALSE)

Seurat::DimPlot(data, reduction = "umap", label = TRUE)
Seurat::DimPlot(data, reduction = "umap", group.by = "treatment")

# Cluster stability (10 and 20 dims)
seurat_obj_10 <- Seurat::FindNeighbors(data, dims = 1:10, verbose = FALSE)
seurat_obj_10 <- Seurat::FindClusters(seurat_obj_10, resolution = 0.5, verbose = FALSE)
seurat_obj_10 <- Seurat::RunUMAP(seurat_obj_10, dims = 1:10)
Seurat::DimPlot(seurat_obj_10, reduction = "umap", label = TRUE)
Seurat::DimPlot(seurat_obj_10, reduction = "umap", group.by = "treatment") +
  ggtitle("UMAP coloured by treatment (10 dims)")


seurat_obj_20 <- Seurat::FindNeighbors(data, dims = 1:20, verbose = FALSE)
seurat_obj_20 <- Seurat::FindClusters(seurat_obj_20, resolution = 0.5, verbose = FALSE)
seurat_obj_20 <- Seurat::RunUMAP(seurat_obj_20, dims = 1:20)
Seurat::DimPlot(seurat_obj_20, reduction = "umap", group.by = "treatment")+
  ggtitle("UMAP coloured by treatment (20 dims)")


##########marker detection and annotation################
library(dplyr)

# Identify cluster markers with 25% expression minimum
markers <- Seurat::FindAllMarkers(
  seurat_obj_20,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = log2(4)  # corresponds to log2FC > 2
)

# Filter significant markers: adjusted p-value < 0.05 and log2FC > 2
sig_markers <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 2)

# View top markers per cluster
top_markers <- sig_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
print(top_markers)


# UMAP colored by SOS response marker (recA)
Seurat::FeaturePlot(
  seurat_obj_20,
  features = "cds-WP-002914769.1",
  cols = c("lightgrey", "blue"),
  pt.size = 0.4
) + ggplot2::ggtitle("recA (SOS response)")

# UMAP colored by heat shock marker (ibpB)
Seurat::FeaturePlot(
  seurat_obj_20,
  features = "cds-WP-000135058.1",
  cols = c("lightgrey", "blue"),
  pt.size = 0.4
) + ggplot2::ggtitle("ibpB (heat shock)")


###########converting the processed data into an AnnData object for MIDAA#################
#saving the data PRE PCA 

library(Matrix)
library(data.table)
library(Seurat)


#defining output path to save to 
output_dir <- "C:/Users/DG1/Desktop/DALLAB/Experimenting/Data/Antibiotic resistance/BacDrop/Seurat/Output exp1"


###Extracting expression matrix########

# Extract log-normalized matrix (pre-PCA)
norm_data <- GetAssayData(data, layer = "data")

# Convert to dense matrix (required for csv)
dense_norm <- as.matrix(GetAssayData(data, layer = "data"))

# write to CSV 
write.csv(dense_norm, file = "adata_X_lognorm_dense.csv")

# Save log-normalized matrix as CSV
write.csv(
  dense_norm,
  file = file.path(output_dir, "adata_X_lognorm_dense.csv")
)

######## Extracting Metadata ##############
# Add barcodes as a column to the metadata
metadata <- data@meta.data
metadata$barcode <- rownames(metadata)

# Save to CSV
write.csv(
  metadata,
  file = file.path(output_dir, "adata_metadata.csv"),
  row.names = FALSE
)


###extracting gene names for ParTIs ########
# Save gene names as plain text file (one per line)
writeLines(
  rownames(dense_norm),
  con = file.path(output_dir, "adata_genes.txt")
)
