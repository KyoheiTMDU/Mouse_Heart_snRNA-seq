# Integration Analysis Script
# Perform SCT normalization and integration across samples

# Load required libraries
library(Seurat)
library(tidyverse)
library(Matrix)

set.seed(1234)
options(future.globals.maxSize = 4000 * 1024^2)
options(Seurat.object.assay.version = "v3")

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "merge", "data")
output_dir <- file.path(base_dir, "integration")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_output, "PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_output, "UMAP"), recursive = TRUE, showWarnings = FALSE)

# Step 1: Load filtered data
cat("Loading filtered data...\n")
load(file.path(input_dir, "filtered_seurat_mt_genes_removed.RData"))
filtered_seurat <- subset_no_mt  # Using the corrected variable name from merge script

# Step 2: Cell cycle scoring
cat("Performing cell cycle scoring...\n")

# Normalize data for cell cycle scoring
seurat_phase <- NormalizeData(filtered_seurat)

# Load cell cycle markers
cycle_file <- file.path(base_dir, "cycle.rda")
if (file.exists(cycle_file)) {
  load(cycle_file)
  
  # Extract S phase genes
  s_genes <- cell_cycle_genes %>%
    dplyr::filter(phase == "S") %>%
    pull("gene_symbol")
  
  # Extract G2M phase genes
  g2m_genes <- cell_cycle_genes %>%
    dplyr::filter(phase == "G2/M") %>%
    pull("gene_symbol")
  
  # Score cells for cell cycle
  seurat_phase <- CellCycleScoring(seurat_phase, 
                                   g2m.features = g2m_genes, 
                                   s.features = s_genes)
} else {
  warning("Cell cycle markers file not found. Skipping cell cycle scoring.")
  seurat_phase$Phase <- "Unknown"
  seurat_phase$S.Score <- 0
  seurat_phase$G2M.Score <- 0
}

# Step 3: Initial processing for visualization
cat("Performing initial processing...\n")

# Find variable features
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale data
seurat_phase <- ScaleData(seurat_phase)

# Run PCA
seurat_phase <- RunPCA(seurat_phase)

# Create mitochondrial ratio categories for visualization
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks = c(-Inf, 0.003043, 0.109252, 0.287604, Inf), 
                                     labels = c("Low", "Medium", "Medium high", "High"))

# Step 4: Generate pre-integration visualizations
cat("Creating pre-integration plots...\n")

# PCA plots
pca_plots <- list(
  cell_cycle = DimPlot(seurat_phase, reduction = "pca", group.by = "Phase", split.by = "Phase"),
  sample = DimPlot(seurat_phase, reduction = "pca", group.by = "sample", split.by = "sample"),
  type = DimPlot(seurat_phase, reduction = "pca", group.by = "type", split.by = "type"),
  mito = DimPlot(seurat_phase, reduction = "pca", group.by = "mitoFr", split.by = "mitoFr")
)

# Save PCA plots
ggsave(file.path(figures_output, "PCA_colored_by_cell_cycle_phase.pdf"), 
       pca_plots$cell_cycle, width = 12, height = 8)
ggsave(file.path(figures_output, "PCA_colored_by_sample.pdf"), 
       pca_plots$sample, width = 12, height = 8)
ggsave(file.path(figures_output, "PCA_colored_by_type.pdf"), 
       pca_plots$type, width = 12, height = 8)
ggsave(file.path(figures_output, "PCA_colored_by_mito_fraction.pdf"), 
       pca_plots$mito, width = 12, height = 8)

# Step 5: Create unintegrated UMAP for comparison
cat("Creating unintegrated UMAP...\n")

umap_unintegrated <- RunUMAP(seurat_phase, dims = 1:40, reduction = "pca")

# Set factor levels for consistent ordering
Idents(umap_unintegrated) <- factor(Idents(umap_unintegrated), 
                                    levels = c("CM1", "CM2", "nonCM1", "nonCM2"))
umap_unintegrated$sample <- factor(umap_unintegrated$sample, 
                                   levels = c("CM1", "CM2", "nonCM1", "nonCM2"))

# UMAP plots
umap_plots <- list(
  combined = DimPlot(umap_unintegrated) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) + 
    theme(legend.text = element_text(size = 10)),
  split_sample = DimPlot(umap_unintegrated, split.by = "sample") + NoLegend(),
  split_type = DimPlot(umap_unintegrated, split.by = "type") + NoLegend()
)

# Save UMAP plots
ggsave(file.path(figures_output, "UMAP_unintegrated.pdf"), 
       umap_plots$combined, width = 10, height = 8)
ggsave(file.path(figures_output, "UMAP_unintegrated_split_by_sample.pdf"), 
       umap_plots$split_sample, width = 16, height = 8)
ggsave(file.path(figures_output, "UMAP_unintegrated_split_by_type.pdf"), 
       umap_plots$split_type, width = 12, height = 8)

# Individual sample UMAP plots
split_samples <- SplitObject(umap_unintegrated, split.by = "sample")
split_samples <- split_samples[c("CM1", "CM2", "nonCM1", "nonCM2")]

for (i in 1:length(split_samples)) {
  sample_name <- split_samples[[i]]@meta.data$sample[1]
  p <- DimPlot(split_samples[[i]]) + ggtitle(paste("Sample:", sample_name))
  ggsave(file.path(figures_output, paste0("UMAP_", sample_name, "_unintegrated.pdf")), 
         p, width = 8, height = 6)
}

# PCA split by sample
pca_unintegrated <- PCAPlot(seurat_phase, split.by = "sample") + NoLegend()
ggsave(file.path(figures_output, "PCA_unintegrated_split_by_sample.pdf"), 
       pca_unintegrated, width = 16, height = 8)

# Step 6: SCTransform and Integration
cat("Performing SCTransform and integration...\n")

# Split object by sample for SCTransform
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("CM1", "CM2", "nonCM1", "nonCM2")]

# Process each sample individually
for (i in 1:length(split_seurat)) {
  sample_name <- names(split_seurat)[i]
  cat("Processing sample:", sample_name, "\n")
  
  # Normalize data
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = FALSE)
  
  # Cell cycle scoring (if markers available)
  if (exists("s_genes") && exists("g2m_genes")) {
    split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], 
                                          g2m.features = g2m_genes, 
                                          s.features = s_genes)
  }
  
  # SCTransform normalization
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("mitoRatio"), 
                                   variable.features.n = 3000,
                                   verbose = FALSE)
}

# Save intermediate results
saveRDS(split_seurat, file.path(data_output, "split_seurat.rds"))

# Step 7: Integration
cat("Finding integration anchors...\n")

# Select integration features
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Prepare SCT integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find integration anchors
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features, 
                                        reference = 1,
                                        verbose = FALSE)

# Integrate data
cat("Integrating data...\n")
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT",
                                   verbose = FALSE)

# Step 8: Post-integration analysis
cat("Performing post-integration analysis...\n")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated, verbose = FALSE)

# Set identities
Idents(seurat_integrated) <- factor(seurat_integrated$sample, 
                                    levels = c("CM1", "CM2", "nonCM1", "nonCM2"))

# Determine optimal number of PCs
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

cat("Optimal number of PCs:", pcs, "\n")
cat("Cumulative variance at PC", pcs, ":", round(cumu[pcs], 2), "%\n")

# Run UMAP with determined PCs
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:20,
                             reduction = "pca",
                             verbose = FALSE)

# Step 9: Generate integrated visualizations
cat("Creating integrated plots...\n")

# PCA plot integrated
pca_integrated <- PCAPlot(seurat_integrated, split.by = "sample")
ggsave(file.path(figures_output, "PCA", "PCA_plot_integrated.pdf"), 
       pca_integrated, width = 16, height = 8)

# UMAP plots integrated
umap_integrated_plots <- list(
  combined = DimPlot(seurat_integrated) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) + 
    theme(legend.text = element_text(size = 10)),
  split_sample = DimPlot(seurat_integrated, split.by = "sample")
)

# Save integrated UMAP plots
ggsave(file.path(figures_output, "UMAP", "UMAP_integrated.pdf"), 
       umap_integrated_plots$combined, width = 10, height = 10)
ggsave(file.path(figures_output, "UMAP", "UMAP_integrated_split_by_sample.pdf"), 
       umap_integrated_plots$split_sample, width = 20, height = 10)

# Individual integrated sample plots
split_integrated <- SplitObject(seurat_integrated, split.by = "sample")
split_integrated <- split_integrated[c("CM1", "CM2", "nonCM1", "nonCM2")]

for (i in 1:length(split_integrated)) {
  sample_name <- split_integrated[[i]]@meta.data$sample[1]
  p <- DimPlot(split_integrated[[i]]) + ggtitle(paste("Integrated Sample:", sample_name))
  ggsave(file.path(figures_output, "UMAP", paste0("UMAP_", sample_name, "_integrated.pdf")), 
         p, width = 10, height = 10)
}

# Step 10: Save final results
cat("Saving final integrated object...\n")
saveRDS(seurat_integrated, file.path(data_output, "integrated_seurat.rds"))

# Print summary
cat("\n=== INTEGRATION SUMMARY ===\n")
cat("Input samples:", length(split_seurat), "\n")
cat("Integration features:", length(integ_features), "\n")
cat("Optimal PCs determined:", pcs, "\n")
cat("Final integrated object dimensions:", dim(seurat_integrated), "\n")
cat("Output directory:", output_dir, "\n")
