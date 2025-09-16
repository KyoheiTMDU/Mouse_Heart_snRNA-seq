# Integrated snRNA-seq QC Analysis Script
# Batch processing for multiple samples

options(future.globals.maxSize = 4000 * 1024^2)

# Load required libraries
suppressPackageStartupMessages({
  library("data.table")
  library("tidyverse")
  library("pheatmap")
  library("reldist")
  library("cowplot")
  library("ggrepel")
  library("Seurat")
  library("sctransform")
  library("patchwork")
  library("SingleCellExperiment")
  library("scales")
  library("RCurl")
  library("Matrix")
  library("DoubletFinder")
})

set.seed(1234)
options(Seurat.object.assay.version = "v3")

# Define sample information
samples_info <- data.frame(
  sample_id = c("GSM5355657_nonCM_Sham_1", "GSM5355658_nonCM_Sham_2", 
                "GSM5943175_CM_Sham_1", "GSM5943176_CM_Sham_2"),
  sample_name = c("nonCM_Sham_1", "nonCM_Sham_2", "CM_Sham_1", "CM_Sham_2"),
  cell_type = c("nonCM", "nonCM", "CM", "CM"),
  clustering_resolution = c(0.5, 0.5, 0.3, 0.3),
  stringsAsFactors = FALSE
)

# Set base directory (modify as needed)
base_dir <- "/Heart"

# Function to perform QC analysis
process_sample <- function(sample_id, sample_name, cell_type, clustering_resolution) {
  
  cat("Processing sample:", sample_name, "\n")
  
  # Set directory paths
  data_dir <- file.path(base_dir, sample_id)
  output_dir <- file.path(data_dir, "QC")
  
  # Load 10X data
  seurat_data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = 100)
  
  # Add cell type information (CM_Sham_2 only)
  if (sample_name == "CM_Sham_2") {
    seurat_obj$cells <- cell_type
  }
  
  # Add metadata
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  
  # Organize metadata
  metadata <- seurat_obj@meta.data
  metadata$cells <- rownames(metadata)
  metadata$sample <- sample_name
  metadata <- metadata %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  seurat_obj@meta.data <- metadata
  
  # Save initial data
  save(seurat_obj, file = file.path(output_dir, "filtered_seurat.RData"))
  
  # Create QC plots
  create_qc_plots(metadata, sample_name, output_dir)
  
  # Filter low quality cells
  filtered_seurat <- subset(x = seurat_obj, 
                            subset = (nGene >= 300) & (mitoRatio < 0.60))
  
  # Filter genes (keep genes expressed in at least 10 cells)
  counts <- GetAssayData(object = filtered_seurat, layer = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
  
  # Save filtered data
  save(filtered_seurat, file = file.path(output_dir, "seurat_filtered.RData"))
  
  # SCTransform normalization
  filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
  DefaultAssay(filtered_seurat) <- "SCT"
  
  # PCA analysis
  filtered_seurat <- RunPCA(filtered_seurat, assay = "SCT", npcs = 50)
  
  # Determine optimal number of PCs
  pct <- filtered_seurat[["pca"]]@stdev / sum(filtered_seurat[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  
  cat("Optimal PCs for", sample_name, ":", pcs, "\n")
  
  # UMAP, neighborhood graph, clustering
  filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", assay = "SCT", dims = 1:20)
  filtered_seurat <- FindNeighbors(object = filtered_seurat, reduction = "pca")
  filtered_seurat <- FindClusters(filtered_seurat, resolution = clustering_resolution)
  
  # Set clustering results
  resolution_col <- paste0("SCT_snn_res.", clustering_resolution)
  Idents(object = filtered_seurat) <- resolution_col
  
  # Save RDS file
  saveRDS(filtered_seurat, file = file.path(output_dir, paste0(sample_name, "_sct.rds")))
  
  # UMAP plot
  p <- DimPlot(filtered_seurat, reduction = "umap", label = TRUE, label.size = 6)
  ggsave(file.path(output_dir, paste0(sample_name, "_umap.pdf")), p, width = 8, height = 6)
  
  # DoubletFinder analysis
  doublet_results <- run_doubletfinder(filtered_seurat, sample_name, output_dir)
  
  return(filtered_seurat)
}

# Function to create QC plots
create_qc_plots <- function(metadata, sample_name, output_dir) {
  
  # Visualize cell counts
  p1 <- metadata %>% 
    ggplot(aes(x = sample, fill = sample)) + 
    geom_bar() +
    geom_text(aes(label = ..count..), stat = "count", vjust = -0.3, colour = "black", size = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells")
  
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
  
  # Complexity score
  p4 <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  
  # Mitochondrial gene percentage
  p5 <- metadata %>% 
    ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  
  # UMI vs gene count correlation
  p6 <- metadata %>% 
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
  
  # Save plots
  combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6)
  ggsave(file.path(output_dir, paste0(sample_name, "_qc_plots.pdf")), 
         combined_plot, width = 12, height = 16)
}

# Function to perform DoubletFinder analysis
run_doubletfinder <- function(seurat_obj, sample_name, output_dir) {
  
  cat("Running DoubletFinder for", sample_name, "\n")
  
  # Optimize pK
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # pK plot
  p_pk <- ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  ggsave(file.path(output_dir, paste0(sample_name, "_pK_plot.pdf")), p_pk)
  
  # Select optimal pK
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  # Estimate homotypic doublet proportion
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.076 * nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj, 
                              PCs = 1:20, 
                              pN = 0.25, 
                              pK = pK, 
                              nExp = nExp_poi.adj,
                              sct = TRUE)
  
  # Get doublet classification column name
  df_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  
  # Visualize doublets
  p_doublet <- DimPlot(seurat_obj, reduction = 'umap', group.by = df_col)
  ggsave(file.path(output_dir, paste0(sample_name, "_doublets_umap.pdf")), p_doublet)
  
  # Doublet statistics
  doublet_stats <- table(seurat_obj@meta.data[[df_col]])
  cat("Doublet statistics for", sample_name, ":\n")
  print(doublet_stats)
  
  # Output DoubletFinder results
  DoubletFinder_output <- data.frame(
    barcode = rownames(seurat_obj@meta.data),
    predicted_doublet = seurat_obj@meta.data[[df_col]]
  )
  write.csv(DoubletFinder_output, 
            file = file.path(output_dir, paste0(sample_name, "_DoubletFinder_output.csv")), 
            row.names = FALSE)
  
  # Separate doublets and singlets and save
  DoubletFinder_doublets <- subset(DoubletFinder_output, 
                                   DoubletFinder_output$predicted_doublet == "Doublet")
  DoubletFinder_nondoublets <- subset(DoubletFinder_output, 
                                      DoubletFinder_output$predicted_doublet == "Singlet")
  
  write.csv(DoubletFinder_doublets, 
            file = file.path(output_dir, paste0(sample_name, "_DoubletFinder_doublets.csv")), 
            row.names = FALSE)
  write.csv(DoubletFinder_nondoublets, 
            file = file.path(output_dir, paste0(sample_name, "_DoubletFinder_nondoublets.csv")), 
            row.names = FALSE)
  
  # Doublet highlight plot
  p_highlight <- DimPlot(seurat_obj, reduction = "umap", 
                         cells.highlight = as.list(DoubletFinder_doublets$barcode)) + NoLegend()
  ggsave(file.path(output_dir, paste0(sample_name, "_doublets_highlight.pdf")), p_highlight)
  
  return(seurat_obj)
}

# Main processing: Process all samples
processed_samples <- list()

for (i in 1:nrow(samples_info)) {
  sample_info <- samples_info[i, ]
  
  processed_sample <- process_sample(
    sample_id = sample_info$sample_id,
    sample_name = sample_info$sample_name,
    cell_type = sample_info$cell_type,
    clustering_resolution = sample_info$clustering_resolution
  )
  
  processed_samples[[sample_info$sample_name]] <- processed_sample
}
