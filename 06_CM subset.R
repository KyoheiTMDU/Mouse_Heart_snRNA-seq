# Cardiac Lineage Cells Subset Analysis Script
# Re-clustering and annotation of cardiac cell populations

# Load required libraries
library(Seurat)
library(tidyverse)
library(clustree)
library(sctransform)
library(patchwork)

set.seed(1234)
options(future.globals.maxSize = 4000 * 1024^2)  
options(Seurat.object.assay.version = "v3")

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "annotation", "data")
output_dir <- file.path(base_dir, "CM_subset")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_output, "clustering"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_output, "markers"), recursive = TRUE, showWarnings = FALSE)

# Step 1: Load annotated data
cat("Loading annotated Seurat object...\n")
Heart <- readRDS(file.path(input_dir, "annotation_seurat.rds"))

cat("Original object summary:\n")
print(Heart)

# Step 2: Visualize original cell types
cat("Creating original cell type visualization...\n")

DefaultAssay(Heart) <- "SCT"
original_plot <- DimPlot(Heart, group.by = "CellType", label = TRUE, repel = TRUE) +
  ggtitle("Original Cell Type Annotations")

ggsave(file.path(figures_output, "original_cell_types.pdf"), 
       original_plot, width = 12, height = 8)

# Step 3: Subset cardiac cells
cat("Subsetting cardiac cells...\n")

Idents(Heart) <- "CellType"
cardiac_cell_types <- c("Cardiac lineage cells_1", "Cardiac lineage cells_2")

CM <- subset(Heart, idents = cardiac_cell_types)

cat("Cardiac subset summary:\n")
print(CM)
cat("Number of cardiac cells:", ncol(CM), "\n")

# Initial visualization of subset
subset_plot <- DimPlot(CM, group.by = "CellType") +
  ggtitle("Cardiac Cells Subset")

ggsave(file.path(figures_output, "cardiac_cells_subset.pdf"), 
       subset_plot, width = 8, height = 6)

# Step 4: Re-normalize and process subset
cat("Re-normalizing cardiac cell subset...\n")

# SCTransform normalization
CM <- SCTransform(CM, variable.features.n = 3000, min_cells = 0, verbose = FALSE)

# Run PCA
CM <- RunPCA(object = CM, verbose = FALSE)

# Step 5: Determine optimal PCs
cat("Determining optimal number of PCs...\n")

pct <- CM[["pca"]]@stdev / sum(CM[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)

cat("PC analysis results:\n")
cat("90% variance + <5% PC variance threshold:", co1, "\n")
cat("Elbow method threshold:", co2, "\n")
cat("Selected PCs:", pcs, "\n")

# Create PC selection plot
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

pc_selection_plot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() +
  ggtitle("PC Selection for Cardiac Cells") +
  xlab("Cumulative Variance (%)") +
  ylab("Individual PC Variance (%)")

ggsave(file.path(figures_output, "clustering", "PC_selection_plot.pdf"), 
       pc_selection_plot, width = 8, height = 6)

# Step 6: UMAP and clustering
cat("Performing UMAP and clustering...\n")

# Use 10 PCs
n_pcs <- 10
CM <- RunUMAP(CM, dims = 1:n_pcs, reduction = "pca", verbose = FALSE)

# Initial UMAP plot
initial_umap <- DimPlot(CM) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +                           
  theme(legend.text = element_text(size = 10)) +
  ggtitle("Cardiac Cells UMAP")

ggsave(file.path(figures_output, "clustering", "initial_UMAP.pdf"), 
       initial_umap, width = 10, height = 8)

# Find neighbors and clusters
CM <- FindNeighbors(object = CM, dims = 1:n_pcs, verbose = FALSE)

# Test multiple resolutions
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
CM <- FindClusters(object = CM, resolution = resolutions, verbose = FALSE)

# Step 7: Resolution selection using clustree
cat("Creating clustree visualization...\n")

clustree_plot <- clustree(CM, prefix = "SCT_snn_res.")
ggsave(file.path(figures_output, "clustering", "clustree.pdf"), 
       clustree_plot, width = 20, height = 12)

# Step 8: Set optimal resolution and visualize
optimal_resolution <- 0.5
cluster_column <- paste0("SCT_snn_res.", optimal_resolution)
Idents(object = CM) <- cluster_column

cat("Selected resolution:", optimal_resolution, "\n")
cat("Number of clusters:", length(levels(Idents(CM))), "\n")

# Clustering visualizations
cluster_plots <- list(
  labeled = DimPlot(CM, reduction = "umap", label = TRUE, label.size = 6) +
    ggtitle(paste("Cardiac Clusters (Resolution", optimal_resolution, ")")),
  unlabeled = DimPlot(CM, reduction = "umap", label = FALSE, pt.size = 3) +
    ggtitle("Cardiac Clusters"),
  by_sample = DimPlot(CM, reduction = "umap", label = FALSE, group.by = "sample") +
    ggtitle("Cardiac Clusters by Sample")
)

if ("type" %in% colnames(CM@meta.data)) {
  cluster_plots$by_type <- DimPlot(CM, reduction = "umap", label = FALSE, group.by = "type") +
    ggtitle("Cardiac Clusters by Type")
}

# Save clustering plots
for (plot_name in names(cluster_plots)) {
  ggsave(file.path(figures_output, "clustering", paste0("clusters_", plot_name, ".pdf")), 
         cluster_plots[[plot_name]], width = 10, height = 8)
}

# Step 9: Prepare for marker analysis
cat("Preparing for marker gene analysis...\n")

# Re-run UMAP with adjusted parameters
Idents(object = CM) <- "CellType"
CM <- RunUMAP(CM, reduction = "pca", dims = 1:n_pcs, min.dist = 0.3, spread = 1, verbose = FALSE)

# Switch to RNA assay and scale data
DefaultAssay(CM) <- "RNA"
all.genes <- rownames(CM)
CM <- ScaleData(CM, features = all.genes, verbose = FALSE)

# Step 10: Annotate cardiac subclusters
cat("Annotating cardiac subclusters...\n")

Idents(object = CM) <- cluster_column

# Define subcluster annotations
subcluster_annotations <- c(
  "0" = "CMC_1",
  "1" = "CLCNM_2", 
  "2" = "CLCNM_3",
  "3" = "CMC_2",
  "4" = "CLCNM_4",
  "5" = "CLCNM_1",
  "6" = "CLCNM_5"
)

# Apply annotations
CM <- RenameIdents(object = CM, subcluster_annotations)

# Add to metadata
CM$SubCellType <- Idents(CM)
CM$SubCellType <- factor(CM$SubCellType,
                         levels = c("CMC_1", "CMC_2",
                                    "CLCNM_1", "CLCNM_2", "CLCNM_3", "CLCNM_4", "CLCNM_5"))
Idents(CM) <- CM$SubCellType

# Final annotated plot
final_plot <- DimPlot(CM, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Cardiac Cell Subtypes") +
  theme(legend.text = element_text(size = 10))

ggsave(file.path(figures_output, "cardiac_cell_subtypes.pdf"), 
       final_plot, width = 10, height = 8)

# Step 11: Save final object and generate statistics
cat("Saving final cardiac subset object...\n")
saveRDS(CM, file.path(data_output, "CM_subset.rds"))

# Generate statistics
subtype_counts <- table(CM$SubCellType)
subtype_by_sample <- FetchData(CM, vars = c("SubCellType", "sample")) %>%
  dplyr::count(SubCellType, sample) %>%
  tidyr::spread(sample, n, fill = 0)

# Save statistics
write.csv(subtype_by_sample, 
          file = file.path(data_output, "cardiac_subtypes_by_sample.csv"),
          row.names = FALSE)

subtype_prop_df <- data.frame(
  SubCellType = names(subtype_counts),
  Count = as.numeric(subtype_counts),
  Percentage = as.numeric(prop.table(subtype_counts) * 100)
)

write.csv(subtype_prop_df, 
          file = file.path(data_output, "cardiac_subtype_proportions.csv"),
          row.names = FALSE)

# Create summary
cm_summary <- list(
  total_cardiac_cells = ncol(CM),
  number_of_subtypes = length(unique(CM$SubCellType)),
  resolution_used = optimal_resolution,
  pcs_used = n_pcs,
  subtype_counts = as.list(subtype_counts)
)

saveRDS(cm_summary, file.path(data_output, "CM_subset_summary.rds"))

cat("\n=== CARDIAC SUBSET ANALYSIS SUMMARY ===\n")
cat("Total cardiac cells:", cm_summary$total_cardiac_cells, "\n")
cat("Number of subtypes:", cm_summary$number_of_subtypes, "\n")
cat("Resolution used:", cm_summary$resolution_used, "\n")
cat("PCs used:", cm_summary$pcs_used, "\n")
cat("\nCardiac subtype distribution:\n")
for (st in names(subtype_counts)) {
  cat(sprintf("  %-6s: %3d cells (%.1f%%)\n", 
              st, 
              subtype_counts[st], 
              subtype_prop_df$Percentage[subtype_prop_df$SubCellType == st]))
}
cat("Output directory:", output_dir, "\n")
