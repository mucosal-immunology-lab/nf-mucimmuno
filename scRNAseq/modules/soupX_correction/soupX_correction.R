# Load necessary libraries
library(Seurat)
library(SoupX)
library(here)
library(DropletUtils)
library(glmGamPoi)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_dir <- args[1]
data_option <- args[2]
PCA_dims <- as.numeric(unlist(strsplit(args[3], ":")))
PCA_dims <- seq(PCA_dims[1], PCA_dims[2], 1)

# Define a function to run SoupX free mRNA contamination estimation calculations
soupX_correction <- function(sample_dir, data_option, PCA_dims) {
  # Construct paths
  filt_path <- file.path(sample_dir, data_option, 'filtered')
  raw_path <- file.path(sample_dir, data_option, 'raw')

  # Print paths for debugging
  cat("Filtered path:", filt_path, "\n")
  cat("Raw path:", raw_path, "\n")
  
  # Read in the filtered matrix (putative cells) and raw matrix (with empty droplets, a.k.a. 'soup')
  if (!dir.exists(filt_path)) stop("Filtered directory does not exist: ", filt_path)
  if (!dir.exists(raw_path)) stop("Raw directory does not exist: ", raw_path)

  filt_matrix <- Read10X(filt_path)
  raw_matrix <- Read10X(raw_path)

  # Handle feature names with underscores
  colnames(filt_matrix) <- gsub("_", "-", colnames(filt_matrix))
  rownames(filt_matrix) <- gsub("_", "-", rownames(filt_matrix))
  colnames(raw_matrix) <- gsub("_", "-", colnames(raw_matrix))
  rownames(raw_matrix) <- gsub("_", "-", rownames(raw_matrix))

  # Create a Seurat object from the sparse matrix
  options(Seurat.object.assay.version = "v3")
  srat <- CreateSeuratObject(counts = filt_matrix)
  
  # Print Seurat object summary for debugging
  print(srat)
  
  # Create the 'SoupChannel' needed by SoupX
  soup_channel <- SoupChannel(raw_matrix, filt_matrix)
  
  # Print SoupChannel summary for debugging
  print(soup_channel)
  
  # Create the clusters (in order to define marker genes) required by SoupX
  srat <- SCTransform(srat, verbose = TRUE)
  srat <- RunPCA(srat, verbose = FALSE)
  srat <- RunUMAP(srat, dims = PCA_dims, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = PCA_dims, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE)
  
  # Print Seurat object metadata for debugging
  print(head(srat@meta.data))
  
  # Add clusters to setClusters; setDR is useful for visualisation
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup_channel <- setClusters(soup_channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup_channel <- setDR(soup_channel, umap)
  
  # With defined clusters, run main SoupX function, calculating ambient RNA profile
  cat("Running autoEstCont...\n")
  soup_channel <- autoEstCont(soup_channel)
  
  # Use roundToInt to ensure you return an output integer matrix
  cat("Adjusting counts...\n")
  adj_matrix <- adjustCounts(soup_channel, roundToInt = TRUE)
  
  # Write the corrected read counts to the directory
  output_dir <- file.path(sample_dir, data_option, 'soupX_corrected')

  # Check if the output directory exists and remove it if it does
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
    message(paste("Existing directory", output_dir, "removed."))
  }

  # Write the corrected read counts to the directory
  DropletUtils::write10xCounts(output_dir, adj_matrix)
  
  # Clear unused RAM
  gc()
}

# Run the function with the parsed arguments
soupX_correction(sample_dir, data_option, PCA_dims)
