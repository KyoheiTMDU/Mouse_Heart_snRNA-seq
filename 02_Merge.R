# Merge Singlet Cells Script
# Process and merge samples after doublet removal

# Load required libraries
library(Seurat)
library(tidyverse)
library(Matrix)

set.seed(1234)
options(Seurat.object.assay.version = "v3")

# Set base directories
base_dir <- "/Heart"
output_dir <- file.path(base_dir, "merge")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories if they don't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)

# Define sample information
samples_info <- data.frame(
  sample_id = c("GSM5943175_CM_Sham_1", "GSM5943176_CM_Sham_2", 
                "GSM5355657_nonCM_Sham_1", "GSM5355658_nonCM_Sham_2"),
  sample_name = c("CM1", "CM2", "nonCM1", "nonCM2"),
  sample_type = c("CM", "CM", "nonCM", "nonCM"),
  stringsAsFactors = FALSE
)

# Function to create QC plots
create_qc_plots <- function(metadata, plot_name, output_path) {
  
  # Cell counts per sample
  p1 <- metadata %>% 
    ggplot(aes(x = sample, fill = sample)) + 
    geom_bar() +
    geom_text(aes(label = ..count..), stat = "count", vjust = -0.3, colour = "black", size = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Number of Cells")
  
  # UMI distribution
  p2 <- metadata %>% 
    ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # Gene count distribution
  p3 <- metadata %>% 
    ggplot(aes(color = sample, x = nGene, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  
  # Gene count boxplot
  p4 <- metadata %>% 
    ggplot(aes(x = sample, y = log10(nGene), fill = sample)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Cells vs Genes")
  
  # UMI vs gene correlation
  p5 <- metadata %>% 
    ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 300) +
    facet_wrap(~sample)
  
  # Mitochondrial gene distribution
  p6 <- metadata %>% 
    ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  
  # Gene complexity
  p7 <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  
  # Save individual plots
  ggsave(file.path(output_path, paste0(plot_name, "_cell_counts.pdf")), p1, width = 8, height = 6)
  ggsave(file.path(output_path, paste0(plot_name, "_UMI_distribution.pdf")), p2, width = 8, height = 6)
  ggsave(file.path(output_path, paste0(plot_name, "_gene_distribution.pdf")), p3, width = 8, height = 6)
  ggsave(file.path(output_path, paste0(plot_name, "_gene_boxplot.pdf")), p4, width = 8, height = 6)
  ggsave(file.path(output_path, paste0(plot_name, "_UMI_gene_correlation.pdf")), p5, width = 12, height = 8)
  ggsave(file.path(output_path, paste0(plot_name, "_mito_distribution.pdf")), p6, width = 8, height = 6)
  ggsave(file.path(output_path, paste0(plot_name, "_complexity.pdf")), p7, width = 8, height = 6)
}

# Step 1: Load 10X data and create Seurat objects
cat("Loading 10X data...\n")

seurat_objects <- list()
for (i in 1:nrow(samples_info)) {
  sample_info <- samples_info[i, ]
  data_path <- file.path(base_dir, sample_info$sample_id)
  
  cat("Loading", sample_info$sample_name, "...\n")
  data <- Read10X(data.dir = data_path)
  seurat_obj <- CreateSeuratObject(counts = data, min.features = 100)
  
  seurat_objects[[sample_info$sample_name]] <- seurat_obj
}

# Step 2: Create initial merged object for unfiltered analysis
cat("Creating initial merged object...\n")

merged_seurat <- merge(x = seurat_objects[[1]], 
                       y = seurat_objects[2:4],
                       add.cell.ids = names(seurat_objects))

# Add metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-") / 100

# Create comprehensive metadata
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)

# Assign sample and type information
metadata$sample <- case_when(
  str_detect(metadata$cells, "^CM1_") ~ "CM1",
  str_detect(metadata$cells, "^CM2_") ~ "CM2",
  str_detect(metadata$cells, "^nonCM1_") ~ "nonCM1",
  str_detect(metadata$cells, "^nonCM2_") ~ "nonCM2"
)

metadata$type <- case_when(
  str_detect(metadata$cells, "^CM1_|^CM2_") ~ "CM",
  str_detect(metadata$cells, "^nonCM1_|^nonCM2_") ~ "nonCM"
)

# Rename columns for consistency
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Update Seurat object metadata
merged_seurat@meta.data <- metadata
merged_seurat@meta.data$orig.ident <- merged_seurat@meta.data$cells

# Save unfiltered data
save(merged_seurat, file = file.path(data_output, "unfiltered_seurat_merged.RData"))

# Generate unfiltered statistics and plots
unfiltered_n_cells <- table(merged_seurat$sample)
write.csv(unfiltered_n_cells, file = file.path(data_output, "n_cells_unfiltered.csv"))

metadata$sample <- factor(metadata$sample, levels = c("CM1", "CM2", "nonCM1", "nonCM2"))
create_qc_plots(metadata, "unfiltered", figures_output)

# Step 3: Remove doublets from individual samples
cat("Removing doublets...\n")

clean_objects <- list()

for (i in 1:nrow(samples_info)) {
  sample_info <- samples_info[i, ]
  sample_name <- sample_info$sample_name
  
  cat("Processing doublets for", sample_name, "...\n")
  
  # Load doublet information
  doublet_file <- file.path(base_dir, sample_info$sample_id, "QC", 
                            paste0(sample_name, "_DoubletFinder_doublets.csv"))
  
  if (file.exists(doublet_file)) {
    doublets <- read.csv(doublet_file)
    
    # Add barcode information to original object
    seurat_objects[[sample_name]]$barcode_ident <- colnames(seurat_objects[[sample_name]])
    
    # Identify doublets
    predicted_doublet <- seurat_objects[[sample_name]]$barcode_ident %in% doublets$barcode
    seurat_objects[[sample_name]]$predicted_doublet <- predicted_doublet
    
    # Remove doublets
    clean_obj <- subset(seurat_objects[[sample_name]], subset = predicted_doublet == FALSE)
    clean_objects[[sample_name]] <- clean_obj
    
    cat("Removed", sum(predicted_doublet), "doublets from", sample_name, "\n")
  } else {
    warning("Doublet file not found for ", sample_name, ". Using original data.")
    clean_objects[[sample_name]] <- seurat_objects[[sample_name]]
  }
}

# Step 4: Create merged object from clean samples
cat("Creating merged object from singlet cells...\n")

merged_clean <- merge(x = clean_objects[[1]], 
                      y = clean_objects[2:4],
                      add.cell.ids = names(clean_objects))

# Add metadata to clean merged object
merged_clean$log10GenesPerUMI <- log10(merged_clean$nFeature_RNA) / log10(merged_clean$nCount_RNA)
merged_clean$mitoRatio <- PercentageFeatureSet(object = merged_clean, pattern = "^mt-") / 100

# Create metadata for clean object
metadata_clean <- merged_clean@meta.data
metadata_clean$cells <- rownames(metadata_clean)

# Assign sample and type information
metadata_clean$sample <- case_when(
  str_detect(metadata_clean$cells, "^CM1_") ~ "CM1",
  str_detect(metadata_clean$cells, "^CM2_") ~ "CM2",
  str_detect(metadata_clean$cells, "^nonCM1_") ~ "nonCM1",
  str_detect(metadata_clean$cells, "^nonCM2_") ~ "nonCM2"
)

metadata_clean$type <- case_when(
  str_detect(metadata_clean$cells, "^CM1_|^CM2_") ~ "CM",
  str_detect(metadata_clean$cells, "^nonCM1_|^nonCM2_") ~ "nonCM"
)

# Rename columns
metadata_clean <- metadata_clean %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Update merged object metadata
merged_clean@meta.data <- metadata_clean
merged_clean@meta.data$orig.ident <- merged_clean@meta.data$cells

# Step 5: Apply quality control filters
cat("Applying quality control filters...\n")

# Filter cells based on QC metrics
filtered_seurat <- subset(x = merged_clean, 
                          subset = (nGene >= 300) & (mitoRatio < 0.60))

# Filter genes (keep genes expressed in at least 10 cells)
counts <- GetAssayData(object = filtered_seurat, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

# Create final filtered object
filtered_seurat <- CreateSeuratObject(filtered_counts, 
                                      meta.data = filtered_seurat@meta.data)

# Update identities
Idents(filtered_seurat) <- filtered_seurat$orig.ident

# Generate filtered statistics and plots
filtered_n_cells <- table(filtered_seurat$sample)
write.csv(filtered_n_cells, file = file.path(data_output, "n_cells_filtered.csv"))

# Update metadata for plotting
metadata_filtered <- filtered_seurat@meta.data
metadata_filtered$sample <- factor(metadata_filtered$sample, levels = c("CM1", "CM2", "nonCM1", "nonCM2"))

create_qc_plots(metadata_filtered, "filtered", figures_output)

# Save filtered data
save(filtered_seurat, file = file.path(data_output, "filtered_seurat.RData"))

# Step 6: Remove mitochondrial genes
cat("Removing mitochondrial genes...\n")

# Identify mitochondrial genes
mt_genes <- grep("^mt-", rownames(filtered_seurat@assays$RNA@counts), value = TRUE)
cat("Found", length(mt_genes), "mitochondrial genes\n")

if (length(mt_genes) > 0) {
  # Get counts without mitochondrial genes
  counts_no_mt <- GetAssayData(filtered_seurat, assay = "RNA")
  counts_no_mt <- counts_no_mt[!rownames(counts_no_mt) %in% mt_genes, ]
  
  # Create subset without mitochondrial genes
  subset_no_mt <- CreateSeuratObject(counts_no_mt, meta.data = filtered_seurat@meta.data)
  
  # Save final object without mitochondrial genes
  save(subset_no_mt, file = file.path(data_output, "filtered_seurat_mt_genes_removed.RData"))
  
  cat("Mitochondrial genes removed. Final object has", nrow(subset_no_mt), "genes\n")
} else {
  cat("No mitochondrial genes found\n")
  subset_no_mt <- filtered_seurat
  save(subset_no_mt, file = file.path(data_output, "filtered_seurat_mt_genes_removed.RData"))
}

# Print summary statistics
cat("\n=== PROCESSING SUMMARY ===\n")
cat("Unfiltered cells:", sum(unfiltered_n_cells), "\n")
cat("Filtered cells:", sum(filtered_n_cells), "\n")
cat("Genes after filtering:", nrow(filtered_seurat), "\n")
cat("Genes after MT removal:", nrow(subset_no_mt), "\n")
cat("Output directory:", output_dir, "\n")