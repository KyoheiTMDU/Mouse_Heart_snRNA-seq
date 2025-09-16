# Clustering Analysis Script
# Perform clustering analysis on integrated data

# Load required libraries
library(Seurat)
library(tidyverse)
library(clustree)

set.seed(1234)
options(future.globals.maxSize = 4000 * 1024^2)
options(Seurat.object.assay.version = "v3")

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "integration", "data")
output_dir <- file.path(base_dir, "clustering")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)

# Step 1: Load integrated data
cat("Loading integrated Seurat object...\n")
seu <- readRDS(file.path(input_dir, "integrated_seurat.rds"))

cat("Loaded object summary:\n")
print(seu)

# Step 2: PCA analysis and visualization
cat("Performing PCA analysis...\n")

# Ensure PCA is computed
seu <- RunPCA(seu, verbose = FALSE)

# Create PC heatmap
pc_heatmap <- DimHeatmap(seu, 
                         dims = 1:9, 
                         cells = 500, 
                         balanced = TRUE)
ggsave(file.path(figures_output, "PC_heatmap.pdf"), 
       pc_heatmap, width = 12, height = 8)

# Create elbow plot
elbow_plot <- ElbowPlot(object = seu, ndims = 50)
ggsave(file.path(figures_output, "elbow_plot.pdf"), 
       elbow_plot, width = 8, height = 6)

# Step 3: Determine optimal number of PCs
cat("Determining optimal number of PCs...\n")

pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)

cat("PC analysis results:\n")
cat("90% variance + <5% PC variance threshold:", co1, "\n")
cat("Elbow method threshold:", co2, "\n")
cat("Selected PCs:", pcs, "\n")
cat("Cumulative variance at PC", pcs, ":", round(cumu[pcs], 2), "%\n")
cat("Cumulative variance at PC 20:", round(cumu[20], 2), "%\n")

# Step 4: Clustering analysis
cat("Performing clustering analysis...\n")

# Find neighbors using 20 PCs
seu <- FindNeighbors(object = seu, dims = 1:20, verbose = FALSE)

# Cluster at resolution 0.8 (as specified in the script)
resolution <- 0.8
seu <- FindClusters(object = seu, resolution = resolution, verbose = FALSE)

cat("Clustering completed at resolution:", resolution, "\n")

# Step 5: Create clustree visualization for available resolutions
cat("Creating clustree visualization...\n")

# Check what resolutions are available in metadata
available_resolutions <- grep("integrated_snn_res", colnames(seu@meta.data), value = TRUE)
cat("Available resolutions:", paste(available_resolutions, collapse = ", "), "\n")

if (length(available_resolutions) > 1) {
  clustree_plot <- clustree(seu, prefix = "integrated_snn_res.")
  ggsave(file.path(figures_output, "clustree.pdf"), 
         clustree_plot, width = 20, height = 12)
} else {
  cat("Only one resolution available, skipping clustree\n")
}

# Set optimal resolution
optimal_resolution <- resolution
cluster_column <- paste0("integrated_snn_res.", optimal_resolution)
Idents(object = seu) <- cluster_column

cat("Selected resolution:", optimal_resolution, "\n")
cat("Number of clusters:", length(levels(Idents(seu))), "\n")

# Step 6: UMAP visualization
cat("Creating UMAP visualizations...\n")

# Ensure UMAP is computed with consistent parameters
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20, verbose = FALSE)

# Basic cluster UMAP
umap_clusters <- DimPlot(seu, 
                         reduction = "umap", 
                         label = TRUE, 
                         label.size = 6) + 
  ggtitle(paste("Clusters (Resolution", optimal_resolution, ")"))

ggsave(file.path(figures_output, "UMAP_clusters_seurat_integrated.pdf"), 
       umap_clusters, width = 10, height = 8)

# UMAP colored by sample
umap_by_sample <- DimPlot(seu, 
                          reduction = "umap", 
                          label = FALSE, 
                          group.by = "sample") +
  ggtitle("Clusters Colored by Sample")

ggsave(file.path(figures_output, "UMAP_clusters_seurat_integrated_group_by_sample.pdf"), 
       umap_by_sample, width = 10, height = 8)

# Step 7: Generate cluster statistics
cat("Generating cluster statistics...\n")

# Extract identity and sample information
cluster_stats <- FetchData(seu, vars = c("ident", "sample")) %>% 
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n, fill = 0)

# Save cluster statistics
write.csv(cluster_stats, 
          file = file.path(data_output, "number_of_nuclei_per_cluster_per_sample.csv"),
          row.names = FALSE)

# Step 8: Additional visualizations
cat("Creating additional visualizations...\n")

# UMAP split by sample
umap_split_sample <- DimPlot(seu, 
                             label = TRUE, 
                             split.by = "sample") + 
  NoLegend() +
  ggtitle("Clusters Split by Sample")

ggsave(file.path(figures_output, "UMAP_cluster_split_by_sample.pdf"), 
       umap_split_sample, width = 20, height = 10)

# UMAP split by cell cycle phase (if available)
if ("Phase" %in% colnames(seu@meta.data)) {
  umap_cell_cycle <- DimPlot(seu,
                             label = TRUE, 
                             split.by = "Phase") + 
    NoLegend() +
    ggtitle("Clusters Split by Cell Cycle Phase")
  
  ggsave(file.path(figures_output, "UMAP_cell_cycle_phase.pdf"), 
         umap_cell_cycle, width = 15, height = 10)
} else {
  cat("Cell cycle phase information not available\n")
}

# Step 9: Quality control metrics visualization
cat("Creating quality control visualizations...\n")

# Define available metrics
all_metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
available_metrics <- all_metrics[all_metrics %in% colnames(seu@meta.data)]

cat("Available QC metrics:", paste(available_metrics, collapse = ", "), "\n")

if (length(available_metrics) > 0) {
  # Feature plot
  feature_plot <- FeaturePlot(seu, 
                              reduction = "umap", 
                              features = available_metrics,
                              pt.size = 0.4, 
                              order = TRUE,
                              min.cutoff = 'q10',
                              label = TRUE,
                              ncol = 2)
  
  ggsave(file.path(figures_output, "feature_plot.pdf"), 
         feature_plot, width = 12, height = ceiling(length(available_metrics)/2) * 4)
  
  # Violin plots for key metrics
  key_metrics <- intersect(c("nUMI", "nGene", "mitoRatio"), available_metrics)
  if (length(key_metrics) > 0) {
    violin_plot <- VlnPlot(seu, 
                           features = key_metrics, 
                           pt.size = 0, 
                           group.by = cluster_column, 
                           stack = TRUE) +
      ggtitle("QC Metrics by Cluster")
    
    ggsave(file.path(figures_output, "violin_plot_QC_metrics.pdf"), 
           violin_plot, width = 10, height = 8)
  }
}

# Step 10: PCA visualizations
cat("Creating PCA visualizations...\n")

# PCA grouped by sample
pca_group_sample <- PCAPlot(seu, group.by = "sample") +
  ggtitle("PCA Colored by Sample")
ggsave(file.path(figures_output, "PCA_group_by_sample.pdf"), 
       pca_group_sample, width = 8, height = 6)

# PCA split by sample
pca_split_sample <- PCAPlot(seu, split.by = "sample") +
  ggtitle("PCA Split by Sample")
ggsave(file.path(figures_output, "PCA_split_by_sample.pdf"), 
       pca_split_sample, width = 16, height = 8)

# Step 11: Save final object
cat("Saving final clustered object...\n")
saveRDS(seu, file = file.path(data_output, "clusters_seurat_integrated.rds"))

# Step 12: Create analysis summary
total_cells_per_cluster <- table(Idents(seu))
total_cells_per_sample <- table(seu$sample)

summary_stats <- list(
  total_cells = ncol(seu),
  total_genes = nrow(seu),
  resolution_used = optimal_resolution,
  number_of_clusters = length(levels(Idents(seu))),
  pcs_used = 20,
  cells_per_cluster = as.list(total_cells_per_cluster),
  cells_per_sample = as.list(total_cells_per_sample)
)

# Save summary
saveRDS(summary_stats, file.path(data_output, "clustering_summary.rds"))
write.csv(data.frame(
  Metric = names(unlist(summary_stats)),
  Value = unlist(summary_stats)
), file.path(data_output, "clustering_summary.csv"), row.names = FALSE)

cat("\n=== CLUSTERING SUMMARY ===\n")
cat("Total cells:", summary_stats$total_cells, "\n")
cat("Total genes:", summary_stats$total_genes, "\n")
cat("Resolution used:", summary_stats$resolution_used, "\n")
cat("Number of clusters:", summary_stats$number_of_clusters, "\n")
cat("PCs used:", summary_stats$pcs_used, "\n")
cat("Cluster distribution:\n")
for (cluster in names(total_cells_per_cluster)) {
  cat(sprintf("  Cluster %-2s: %4d cells\n", cluster, total_cells_per_cluster[cluster]))
}
cat("Output directory:", output_dir, "\n")
