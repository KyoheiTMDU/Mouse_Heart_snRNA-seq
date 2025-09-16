# Cell Type Annotation Script
# Annotate clusters with cell type identities

# Load required libraries
library(Seurat)
library(tidyverse)

set.seed(1234)
options(future.globals.maxSize = 10000 * 1024^2)  
options(Seurat.object.assay.version = "v3")

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "clustering", "data")
output_dir <- file.path(base_dir, "annotation")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)

# Step 1: Load clustered data
cat("Loading clustered Seurat object...\n")
seu <- readRDS(file.path(input_dir, "clusters_seurat_integrated.rds"))

cat("Loaded object summary:\n")
print(seu)

# Step 2: Set up clustering resolution
cat("Setting up cluster identities...\n")

# Set cluster resolution (0.8 as specified in the script)
resolution <- 0.8
cluster_column <- paste0("integrated_snn_res.", resolution)
Idents(object = seu) <- cluster_column

# Get current cluster levels and order them numerically
current_clusters <- levels(seu[[cluster_column, drop = TRUE]])
ordered_clusters <- as.character(sort(as.numeric(current_clusters)))

# Set factor levels in numerical order
Idents(seu) <- factor(
  seu[[cluster_column, drop = TRUE]], 
  levels = ordered_clusters
)

cat("Number of clusters:", length(ordered_clusters), "\n")
cat("Cluster IDs:", paste(ordered_clusters, collapse = ", "), "\n")

# Step 3: Prepare for annotation
cat("Preparing data for annotation...\n")

# Switch to RNA assay for annotation
DefaultAssay(seu) <- "RNA"

# Scale all genes for marker gene analysis
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes, verbose = FALSE)

# Step 4: Create pre-annotation visualization
cat("Creating pre-annotation UMAP...\n")

pre_annotation_plot <- DimPlot(seu,
                               reduction = "umap",
                               label = TRUE,
                               label.size = 4) +
  ggtitle("Clusters Before Annotation")

ggsave(file.path(figures_output, "UMAP_clusters_before_annotation.pdf"), 
       pre_annotation_plot, width = 10, height = 8)

# Step 5: Define cell type annotations
cat("Applying cell type annotations...\n")

# Define cluster to cell type mapping (updated for resolution 0.8)
cluster_annotations <- c(
  "0" = "Cardiac lineage cells_1",
  "1" = "Endothelial_1",
  "2" = "Fibroblast_1", 
  "3" = "Fibroblast_2",
  "4" = "Endothelial_2",
  "5" = "Cardiac lineage cells_2",
  "6" = "Myeloid",
  "7" = "Fibroblast_3",
  "8" = "B_cell",
  "9" = "Pericyte_1",
  "10" = "Endothelial_3",
  "11" = "Pericyte_2",
  "12" = "Mesothelial",
  "13" = "Endothelial_4"
)

# Verify all clusters are covered
missing_clusters <- setdiff(ordered_clusters, names(cluster_annotations))
if (length(missing_clusters) > 0) {
  warning("Missing annotations for clusters: ", paste(missing_clusters, collapse = ", "))
}

# Apply annotations
Heart <- RenameIdents(object = seu, cluster_annotations)

# Step 6: Organize cell type factors
cat("Organizing cell type categories...\n")

# Add cell type to metadata
Heart$CellType <- Idents(Heart)

# Define logical cell type order
cell_type_order <- c(
  "Cardiac lineage cells_1", "Cardiac lineage cells_2",
  "Endothelial_1", "Endothelial_2", "Endothelial_3", "Endothelial_4",
  "Fibroblast_1", "Fibroblast_2", "Fibroblast_3", "Mesothelial",
  "Myeloid", "B_cell", "Pericyte_1", "Pericyte_2"
)

# Set factor levels
Heart$CellType <- factor(Heart$CellType, levels = cell_type_order)
Idents(Heart) <- Heart$CellType

# Verify annotation
cat("Cell type annotation summary:\n")
cell_type_counts <- table(Heart$CellType)
print(cell_type_counts)

# Step 7: Create visualizations
cat("Creating annotation visualizations...\n")

# Basic annotated UMAP
umap_annotated <- DimPlot(Heart, 
                          reduction = "umap", 
                          label = TRUE, 
                          repel = TRUE) + 
  NoLegend() +
  ggtitle("Cell Type Annotations") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave(file.path(figures_output, "UMAP_cell_type_annotations.pdf"), 
       umap_annotated, width = 10, height = 8)

# UMAP with legend
umap_with_legend <- DimPlot(Heart, 
                            reduction = "umap", 
                            label = FALSE) +
  ggtitle("Cell Type Annotations") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(figures_output, "UMAP_cell_type_annotations_with_legend.pdf"), 
       umap_with_legend, width = 12, height = 8)

# Split by sample
umap_split_sample <- DimPlot(Heart, 
                             reduction = "umap", 
                             split.by = "sample",
                             label = TRUE,
                             repel = TRUE) + 
  NoLegend() +
  ggtitle("Cell Types by Sample")

ggsave(file.path(figures_output, "UMAP_cell_types_split_by_sample.pdf"), 
       umap_split_sample, width = 20, height = 10)

# Split by type (CM vs nonCM) if available
if ("type" %in% colnames(Heart@meta.data)) {
  umap_split_type <- DimPlot(Heart, 
                             reduction = "umap", 
                             split.by = "type",
                             label = TRUE,
                             repel = TRUE) + 
    NoLegend() +
    ggtitle("Cell Types by Sample Type")
  
  ggsave(file.path(figures_output, "UMAP_cell_types_split_by_type.pdf"), 
         umap_split_type, width = 12, height = 8)
}

# Step 8: Generate annotation statistics
cat("Generating annotation statistics...\n")

# Cell type counts by sample
cell_type_by_sample <- FetchData(Heart, vars = c("CellType", "sample")) %>%
  dplyr::count(CellType, sample) %>%
  tidyr::spread(sample, n, fill = 0)

write.csv(cell_type_by_sample, 
          file = file.path(data_output, "cell_type_counts_by_sample.csv"),
          row.names = FALSE)

# Cell type proportions
cell_type_proportions <- prop.table(table(Heart$CellType)) * 100
cell_type_prop_df <- data.frame(
  CellType = names(cell_type_proportions),
  Count = as.numeric(table(Heart$CellType)),
  Percentage = as.numeric(cell_type_proportions)
)

write.csv(cell_type_prop_df, 
          file = file.path(data_output, "cell_type_proportions.csv"),
          row.names = FALSE)

# Cluster to cell type mapping table
mapping_table <- data.frame(
  Cluster = names(cluster_annotations),
  CellType = cluster_annotations,
  stringsAsFactors = FALSE
)

write.csv(mapping_table, 
          file = file.path(data_output, "cluster_to_celltype_mapping.csv"),
          row.names = FALSE)

# Step 9: Save annotated object
cat("Saving annotated Seurat object...\n")
saveRDS(Heart, file = file.path(data_output, "annotation_seurat.rds"))

# Step 10: Create summary report
cat("Creating annotation summary...\n")

annotation_summary <- list(
  total_cells = ncol(Heart),
  total_cell_types = length(unique(Heart$CellType)),
  cell_type_counts = as.list(cell_type_counts),
  resolution_used = resolution,
  clusters_annotated = length(cluster_annotations)
)

saveRDS(annotation_summary, file.path(data_output, "annotation_summary.rds"))

# Print summary
cat("\n=== ANNOTATION SUMMARY ===\n")
cat("Total cells:", annotation_summary$total_cells, "\n")
cat("Total cell types:", annotation_summary$total_cell_types, "\n")
cat("Resolution used:", annotation_summary$resolution_used, "\n")
cat("Clusters annotated:", annotation_summary$clusters_annotated, "\n")
cat("\nCell type distribution:\n")
for (ct in names(cell_type_counts)) {
  cat(sprintf("  %-18s: %4d cells (%.1f%%)\n", 
              ct, 
              cell_type_counts[ct], 
              cell_type_proportions[ct]))
}
cat("Output directory:", output_dir, "\n")
