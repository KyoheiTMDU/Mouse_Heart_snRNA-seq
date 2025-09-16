# GSEA Analysis Script
# Gene Set Enrichment Analysis using escape package

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(UpSetR)
  library(escape)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(rstatix)
  library(pheatmap)
  library(reshape2)
  library(corrplot)
  library(ggsci)
})

set.seed(1234)
options(future.globals.maxSize = 10000 * 1024^2)  
options(Seurat.object.assay.version = "v3")

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "annotation", "data")
output_dir <- file.path(base_dir, "GSEA")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")
results_output <- file.path(output_dir, "results")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)
dir.create(results_output, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load and prepare data
# ==============================================================================
cat("Loading annotated Heart data...\n")
Heart <- readRDS(file.path(input_dir, "annotation_seurat.rds"))

cat("Loaded Heart object summary:\n")
print(Heart)

# ==============================================================================
# 2. Gene expression analysis - UpSet plot
# ==============================================================================
cat("Creating UpSet plot for marker gene co-expression...\n")

# Use the original working approach
cells_geneA <- WhichCells(
  object = Heart,
  expression = Bmi1 > 0
)

cells_geneB <- WhichCells(
  object = Heart,
  expression = Abcg2 > 0
)

cells_geneC <- WhichCells(
  object = Heart,
  expression = Ly6a > 0
)

cells_geneD <- WhichCells(
  object = Heart,
  expression = Hopx > 0
)

cells_geneE <- WhichCells(
  object = Heart,
  expression = `Nkx2-5` > 0
)

# Create UpSet data frame
all_cells <- unique(c(cells_geneA, cells_geneB, cells_geneC, cells_geneD, cells_geneE))

# Create binary matrix
upset_df <- data.frame(
  Cell = all_cells,
  Bmi1 = as.integer(all_cells %in% cells_geneA),
  Abcg2 = as.integer(all_cells %in% cells_geneB),
  Ly6a = as.integer(all_cells %in% cells_geneC),
  Hopx = as.integer(all_cells %in% cells_geneD),
  `Nkx2-5` = as.integer(all_cells %in% cells_geneE),
  check.names = FALSE
)

# Create UpSet plot
pdf(file.path(figures_output, "upset_plot_marker_genes.pdf"), width = 10, height = 6)
upset(
  upset_df,
  sets = c("Bmi1", "Hopx", "Nkx2-5", "Abcg2", "Ly6a"),
  order.by = "freq",
  main.bar.color = "darkblue",
  sets.bar.color = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  matrix.color = "gray30",
  point.size = 3.5,
  line.size = 1.5,
  text.scale = c(1.5, 1.3, 1, 1, 1.3, 1.3),
  mb.ratio = c(0.6, 0.4)
)
dev.off()

cat("UpSet plot saved successfully\n")

# ==============================================================================
# 3. GSEA Analysis using escape
# ==============================================================================
cat("Performing GSEA analysis...\n")

# Set cell type identities
Idents(Heart) <- "CellType"

# Get Hallmark gene sets
cat("Loading Hallmark gene sets...\n")
tryCatch({
  GS.hallmark <- getGeneSets(species = "Mus musculus", library = "H")
  cat("Loaded", length(GS.hallmark), "Hallmark gene sets\n")
}, error = function(e) {
  stop("Error loading gene sets: ", e$message)
})

# Run ssGSEA
cat("Running ssGSEA analysis...\n")
tryCatch({
  Heart <- runEscape(Heart, 
                     method = "ssGSEA",
                     gene.sets = GS.hallmark,
                     min.size = 0,
                     new.assay.name = "escape.ssGSEA")
  cat("ssGSEA analysis completed\n")
}, error = function(e) {
  stop("Error in ssGSEA analysis: ", e$message)
})

# ==============================================================================
# 4. Create enrichment visualizations
# ==============================================================================
cat("Creating enrichment visualizations...\n")

# Heatmap of all pathways
pdf(file.path(figures_output, "enrichment_heatmap_all_pathways.pdf"), width = 12, height = 8)
heatmapEnrichment(Heart, 
                  group.by = "CellType",
                  assay = "escape.ssGSEA",
                  gene.set.use = "all",
                  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = FALSE)
dev.off()

# Specific pathway visualizations
pathways_of_interest <- c("HALLMARK-MYOGENESIS", "HALLMARK-NOTCH-SIGNALING")
available_pathways <- intersect(pathways_of_interest, rownames(Heart[["escape.ssGSEA"]]))

if (length(available_pathways) > 0) {
  # Geyser plots for specific pathways
  geyser_plots <- list()
  
  for (pathway in available_pathways) {
    cat("Creating geyser plot for", pathway, "\n")
    tryCatch({
      geyser_plots[[pathway]] <- geyserEnrichment(Heart, 
                                                  assay = "escape.ssGSEA",
                                                  gene.set = pathway, 
                                                  color.by = "ident",
                                                  scale = TRUE) +
        scale_color_d3(palette = "category20") +
        ggtitle(gsub("HALLMARK-", "", pathway))
    }, error = function(e) {
      cat("Error creating geyser plot for", pathway, ":", e$message, "\n")
    })
  }
  
  # Save combined geyser plots
  if (length(geyser_plots) > 1) {
    combined_geyser <- wrap_plots(geyser_plots, ncol = 1)
    ggsave(file.path(figures_output, "geyser_plots_combined.pdf"), 
           combined_geyser, width = 10, height = 12)
  } else if (length(geyser_plots) == 1) {
    ggsave(file.path(figures_output, paste0("geyser_plot_", names(geyser_plots)[1], ".pdf")), 
           geyser_plots[[1]], width = 10, height = 6)
  }
}

# ==============================================================================
# 5. Extract enrichment scores and statistical analysis
# ==============================================================================
cat("Extracting enrichment scores for statistical analysis...\n")

# Extract enrichment matrix
escape_scores <- t(as.matrix(GetAssayData(Heart, assay = "escape.ssGSEA", slot = "data")))

# Verify cell order matches
if (!all(rownames(escape_scores) == colnames(Heart))) {
  warning("Cell order mismatch between scores and metadata")
}

# Add to metadata
Heart <- AddMetaData(Heart, metadata = escape_scores)

# ==============================================================================
# 6. Statistical testing for specific pathways
# ==============================================================================
cat("Performing statistical tests...\n")

# Function to perform pathway analysis
perform_pathway_analysis <- function(pathway_name, display_name) {
  if (!pathway_name %in% colnames(Heart@meta.data)) {
    cat("Pathway", pathway_name, "not found in metadata\n")
    return(NULL)
  }
  
  cat("Analyzing", display_name, "\n")
  
  # Create metadata subset
  meta_df <- Heart@meta.data %>%
    dplyr::select(CellType, !!sym(pathway_name)) %>%
    dplyr::rename(Score = !!sym(pathway_name))
  
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(Score ~ CellType, data = meta_df)
  
  # Dunn's test for pairwise comparisons
  dunn_result <- meta_df %>%
    dunn_test(Score ~ CellType, p.adjust.method = "BH")
  
  # Save results
  safe_name <- gsub("[^A-Za-z0-9_]", "_", display_name)
  writeLines(capture.output(kruskal_result), 
             file.path(results_output, paste0("kruskal_", safe_name, ".txt")))
  write_tsv(dunn_result, 
            file.path(results_output, paste0("dunn_", safe_name, ".txt")))
  
  return(list(kruskal = kruskal_result, dunn = dunn_result))
}

# Analyze specific pathways
analysis_results <- list()
pathway_analyses <- list(
  "HALLMARK-MYOGENESIS" = "Myogenesis",
  "HALLMARK-NOTCH-SIGNALING" = "Notch_Signaling"
)

for (pathway in names(pathway_analyses)) {
  analysis_results[[pathway]] <- perform_pathway_analysis(pathway, pathway_analyses[[pathway]])
}

# ==============================================================================
# 7. Create heatmaps for pairwise comparisons
# ==============================================================================
cat("Creating pairwise comparison heatmaps...\n")

# Function to create comparison heatmap
create_comparison_heatmap <- function(dunn_data, title_suffix) {
  if (is.null(dunn_data) || nrow(dunn_data) == 0) {
    cat("No data available for heatmap:", title_suffix, "\n")
    return(NULL)
  }
  
  # Create symmetric matrix
  groups <- unique(c(dunn_data$group1, dunn_data$group2))
  n_groups <- length(groups)
  
  # Initialize matrices
  p_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
  sig_matrix <- matrix("ns", nrow = n_groups, ncol = n_groups)
  rownames(p_matrix) <- colnames(p_matrix) <- groups
  rownames(sig_matrix) <- colnames(sig_matrix) <- groups
  
  # Fill matrices
  for (i in 1:nrow(dunn_data)) {
    g1 <- dunn_data$group1[i]
    g2 <- dunn_data$group2[i]
    p_val <- dunn_data$p.adj[i]
    sig <- dunn_data$p.adj.signif[i]
    
    p_matrix[g1, g2] <- p_val
    p_matrix[g2, g1] <- p_val
    sig_matrix[g1, g2] <- sig
    sig_matrix[g2, g1] <- sig
  }
  
  # Set diagonal
  diag(p_matrix) <- 0
  diag(sig_matrix) <- "-"
  
  # Create upper triangle version
  log_p_matrix <- -log10(p_matrix + 1e-100)
  log_p_matrix[lower.tri(log_p_matrix)] <- NA
  sig_matrix[lower.tri(sig_matrix)] <- ""
  
  # Convert to long format for ggplot
  melted_data <- reshape2::melt(log_p_matrix, na.rm = TRUE)
  melted_sig <- reshape2::melt(sig_matrix, na.rm = TRUE)
  
  # Combine data
  plot_data <- merge(melted_data, melted_sig, by = c("Var1", "Var2"))
  colnames(plot_data) <- c("Group1", "Group2", "LogP", "Significance")
  
  # Filter for upper triangle only
  plot_data <- plot_data[!is.na(plot_data$LogP) & plot_data$LogP > 0, ]
  
  if (nrow(plot_data) == 0) {
    cat("No significant comparisons for", title_suffix, "\n")
    return(NULL)
  }
  
  # Create ggplot heatmap
  p <- ggplot(plot_data, aes(x = Group1, y = Group2, fill = LogP)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Significance), color = "black", size = 4, fontface = "bold") +
    scale_fill_viridis_c(name = "-log10(p)", 
                         option = "plasma",
                         trans = "sqrt") +
    labs(title = "Pairwise Comparisons (Upper Triangle)",
         subtitle = paste(title_suffix, "pathway analysis"),
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    coord_fixed()
  
  return(p)
}

# Create heatmaps for each pathway
for (pathway in names(analysis_results)) {
  if (!is.null(analysis_results[[pathway]])) {
    display_name <- pathway_analyses[[pathway]]
    heatmap_plot <- create_comparison_heatmap(analysis_results[[pathway]]$dunn, display_name)
    
    if (!is.null(heatmap_plot)) {
      safe_name <- gsub("[^A-Za-z0-9_]", "_", display_name)
      ggsave(file.path(figures_output, paste0("pairwise_heatmap_", safe_name, ".pdf")), 
             heatmap_plot, width = 10, height = 8)
    }
  }
}

# ==============================================================================
# 8. Save processed data and summary
# ==============================================================================
cat("Saving processed data and summary...\n")

# Save Heart object with GSEA results
saveRDS(Heart, file.path(data_output, "Heart_with_GSEA.rds"))

# Save enrichment scores
write.csv(escape_scores, 
          file.path(data_output, "enrichment_scores_all_pathways.csv"),
          row.names = TRUE)

# Create summary
gsea_summary <- list(
  total_cells = ncol(Heart),
  cell_types = levels(Idents(Heart)),
  n_cell_types = length(levels(Idents(Heart))),
  n_pathways = nrow(Heart[["escape.ssGSEA"]]),
  pathways_analyzed = names(pathway_analyses),
  analysis_successful = !is.null(analysis_results)
)

saveRDS(gsea_summary, file.path(data_output, "GSEA_analysis_summary.rds"))

# Print summary
cat("\n=== GSEA ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", gsea_summary$total_cells, "\n")
cat("Number of cell types:", gsea_summary$n_cell_types, "\n")
cat("Cell types:", paste(gsea_summary$cell_types, collapse = ", "), "\n")
cat("Number of pathways analyzed:", gsea_summary$n_pathways, "\n")
cat("Specific pathways tested:", paste(gsea_summary$pathways_analyzed, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")
