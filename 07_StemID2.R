# StemID2 Analysis Script
# Stemness analysis using RaceID package on cardiac cell subtypes

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
})

# Check if RaceID is available
if (!requireNamespace("RaceID", quietly = TRUE)) {
  stop("RaceID package is required but not installed. Please install from Bioconductor.")
} else {
  suppressPackageStartupMessages(library(RaceID))
}

set.seed(1234)

# Set directories
base_dir <- "/Heart"
input_dir <- file.path(base_dir, "CM_subset", "data")
output_dir <- file.path(base_dir, "StemID")
data_output <- file.path(output_dir, "data")
figures_output <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load and prepare data
# ==============================================================================
cat("Loading cardiac subset data...\n")
Heart <- readRDS(file.path(input_dir, "CM_subset.rds"))
Idents(Heart) <- "SubCellType"

cat("Loaded cardiac subset summary:\n")
print(Heart)

# Store original cluster information
original_clusters <- as.character(Idents(Heart))
names(original_clusters) <- colnames(Heart)
cluster_mapping <- sort(unique(original_clusters))

cat("Original SubCellType clusters:\n")
cluster_table <- table(original_clusters)
print(cluster_table)
cat("Number of SubCellTypes:", length(cluster_mapping), "\n")

# Convert Seurat object to SCseq object
cat("Converting Seurat to SCseq object...\n")
tryCatch({
  Heart_sc <- Seurat2SCseq(Heart, rseed = 12345)
  cat("Conversion successful\n")
}, error = function(e) {
  stop("Error in Seurat2SCseq conversion: ", e$message)
})

# ==============================================================================
# 2. Compute distance matrix (essential step)
# ==============================================================================
cat("Computing distance matrix...\n")
Heart_sc <- compdist(Heart_sc, metric = "pearson")

# ==============================================================================
# 3. Maintain original clusters (skip RaceID clustering)
# ==============================================================================
cat("Setting up original cluster information...\n")

# Check cell name correspondence
seurat_cells <- names(original_clusters)
scseq_cells <- colnames(Heart_sc@ndata)

# Find common cells
common_cells <- intersect(seurat_cells, scseq_cells)
if (length(common_cells) == 0 && length(seurat_cells) == length(scseq_cells)) {
  # If cell names differ but counts match, assume same order
  names(original_clusters) <- scseq_cells
  common_cells <- scseq_cells
  cat("Cell names adjusted for SCseq object\n")
}

# Convert clusters to numeric
cluster_to_numeric <- setNames(seq_along(cluster_mapping), cluster_mapping)
numeric_clusters <- cluster_to_numeric[original_clusters[common_cells]]

# Set cluster information in SCseq object
Heart_sc@cpart <- numeric_clusters
names(Heart_sc@cpart) <- common_cells

# Explicitly set cluster number
Heart_sc@clusterpar$cln <- length(cluster_mapping)

# Set cluster colors
Heart_sc@fcol <- rainbow(length(cluster_mapping))

# Store mapping information
cluster_name_mapping <- setNames(cluster_mapping, seq_along(cluster_mapping))
cat("Cluster mapping (numeric ID -> SubCellType):\n")
print(cluster_name_mapping)
cat("Number of clusters in Heart_sc:", length(unique(Heart_sc@cpart)), "\n")

# ==============================================================================
# 4. Dimension reduction
# ==============================================================================
cat("Computing dimension reductions...\n")

tryCatch({
  Heart_sc <- comptsne(Heart_sc, perplexity = 30)
  cat("t-SNE completed\n")
}, error = function(e) {
  cat("Warning: t-SNE failed:", e$message, "\n")
})

tryCatch({
  Heart_sc <- compfr(Heart_sc, knn = 10)
  cat("Force-directed layout completed\n")
}, error = function(e) {
  cat("Warning: Force-directed layout failed:", e$message, "\n")
})

tryCatch({
  Heart_sc <- compumap(Heart_sc, n_neighbors = 15)
  cat("UMAP completed\n")
}, error = function(e) {
  cat("Warning: UMAP failed:", e$message, "\n")
})

# ==============================================================================
# 5. StemID2 analysis using original SubCellTypes
# ==============================================================================
cat("\n=== Starting StemID2 analysis with original SubCellTypes ===\n")

# Build lineage tree
cat("Building lineage tree...\n")
ltr <- Ltree(Heart_sc)
ltr <- compentropy(ltr)

# Verify cluster count hasn't changed
cat("Clusters in lineage tree:", length(unique(ltr@sc@cpart)), "\n")

# Project cells with different threshold values
success <- FALSE
cthr_values <- c(10, 20, 30, 50)

for (cthr in cthr_values) {
  cat("Trying projcells with cthr =", cthr, "\n")
  tryCatch({
    ltr <- projcells(ltr, cthr = cthr, nmode = FALSE, fr = FALSE)
    success <- TRUE
    cat("projcells successful with cthr =", cthr, "\n")
    break
  }, error = function(e) {
    cat("Failed with cthr =", cthr, ":", e$message, "\n")
  })
}

# Continue analysis if projection was successful
stemness_scores <- NULL
if (success) {
  cat("Continuing with lineage analysis...\n")
  
  tryCatch({
    ltr <- projback(ltr, pdishuf = 50)
    ltr <- lineagegraph(ltr)
    ltr <- comppvalue(ltr, pthr = 0.05)
    cat("Lineage graph construction completed\n")
  }, error = function(e) {
    cat("Warning in lineage analysis:", e$message, "\n")
  })
  
  # Calculate stemness scores
  cat("Computing stemness scores...\n")
  tryCatch({
    stemness_result <- compscore(ltr, scthr = 0.3)
    
    # Handle different return types
    if (is.numeric(stemness_result)) {
      stemness_scores <- stemness_result
    } else if (is.list(stemness_result) && "score" %in% names(stemness_result)) {
      stemness_scores <- stemness_result$score
    } else {
      stemness_scores <- unlist(stemness_result)
    }
    
    cat("Number of stemness scores:", length(stemness_scores), "\n")
    
    # Label with SubCellType names
    if (length(stemness_scores) == length(cluster_mapping)) {
      names(stemness_scores) <- cluster_name_mapping[as.character(names(stemness_scores))]
      stemness_scores <- sort(stemness_scores, decreasing = TRUE)
      
      cat("Stemness scores by SubCellType:\n")
      print(round(stemness_scores, 4))
    } else {
      cat("WARNING: Number of stemness scores (", length(stemness_scores), 
          ") does not match number of SubCellTypes (", length(cluster_mapping), ")\n")
    }
    
  }, error = function(e) {
    cat("Error in stemness calculation:", e$message, "\n")
    stemness_scores <- NULL
  })
} else {
  cat("Lineage analysis could not be completed - projection failed\n")
}

# ==============================================================================
# 6. Visualization
# ==============================================================================
cat("Creating visualizations...\n")

# Create comprehensive results plot
pdf(file.path(figures_output, "stemid2_subcelltype_results.pdf"), width = 12, height = 10)

if (!is.null(stemness_scores) && !all(is.na(stemness_scores))) {
  par(mfrow = c(2, 1))
  
  # Stemness scores barplot
  barplot(stemness_scores, 
          las = 2, 
          main = paste("Stemness Scores by SubCellType (n =", length(stemness_scores), ")"),
          ylab = "Stemness Score",
          col = heat.colors(length(stemness_scores)),
          cex.names = 0.8)
  
  # Cell distribution
  cluster_dist <- table(Heart_sc@cpart)
  names(cluster_dist) <- cluster_name_mapping[names(cluster_dist)]
  barplot(cluster_dist,
          las = 2,
          main = "Cell Distribution by SubCellType",
          ylab = "Number of Cells",
          col = rainbow(length(cluster_dist)),
          cex.names = 0.8)
} else {
  plot(1, type = "n", main = "Stemness Analysis Results")
  text(1, 1, "Stemness scores could not be calculated", cex = 1.2)
}

dev.off()

# Lineage analysis visualization
pdf(file.path(figures_output, "lineage_analysis_results.pdf"), width = 14, height = 12)

if (success && exists("ltr")) {
  par(mfrow = c(2, 2))
  
  # Lineage graph
  tryCatch({
    plotgraph(ltr, scthr = 0.3, showCells = FALSE, showMap = TRUE)
    title("Lineage Graph - SubCellTypes")
  }, error = function(e) {
    plot(1, type = "n", main = "Lineage Graph - Error")
    text(1, 1, paste("Error:", e$message), cex = 0.8)
  })
  
  # Additional plots if available
  tryCatch({
    # Plot cluster distribution
    cluster_dist <- table(Heart_sc@cpart)
    pie(cluster_dist, 
        labels = cluster_name_mapping[names(cluster_dist)],
        main = "SubCellType Distribution")
  }, error = function(e) {
    plot(1, type = "n", main = "Distribution Plot - Error")
  })
  
} else {
  plot(1, type = "n", main = "Lineage Analysis")
  text(1, 1, "Lineage analysis could not be completed", cex = 1)
}

dev.off()

# ==============================================================================
# 7. Save results
# ==============================================================================
cat("Saving results...\n")

# Save processed objects
saveRDS(Heart_sc, file = file.path(data_output, "Heart_StemID_7subtypes.rds"))

if (exists("ltr")) {
  saveRDS(ltr, file = file.path(data_output, "Heart_lineage_tree_7subtypes.rds"))
}

# Create and save summary
summary_7subtypes <- list(
  n_subtypes = length(cluster_mapping),
  subtype_names = cluster_mapping,
  n_cells_per_subtype = as.list(table(Heart_sc@cpart)),
  stemness_scores = if (!is.null(stemness_scores)) as.list(stemness_scores) else NA,
  analysis_successful = success,
  total_cells = length(Heart_sc@cpart)
)

saveRDS(summary_7subtypes, file = file.path(data_output, "Heart_7subtypes_summary.rds"))

# Save stemness scores as CSV if available
if (!is.null(stemness_scores) && !all(is.na(stemness_scores))) {
  stemness_df <- data.frame(
    SubCellType = names(stemness_scores),
    Stemness_Score = as.numeric(stemness_scores),
    stringsAsFactors = FALSE
  )
  write.csv(stemness_df, 
            file = file.path(data_output, "stemness_scores_by_subcelltype.csv"),
            row.names = FALSE)
}

# Save cluster mapping
mapping_df <- data.frame(
  Numeric_ID = names(cluster_name_mapping),
  SubCellType = cluster_name_mapping,
  Cell_Count = as.numeric(table(Heart_sc@cpart)),
  stringsAsFactors = FALSE
)
write.csv(mapping_df, 
          file = file.path(data_output, "subcelltype_mapping.csv"),
          row.names = FALSE)

# ==============================================================================
# 8. Print summary
# ==============================================================================
cat("\n=== STEMID2 ANALYSIS SUMMARY ===\n")
cat("Analyzed SubCellTypes:", length(cluster_mapping), "\n")
cat("Total cells:", summary_7subtypes$total_cells, "\n")
cat("Analysis successful:", success, "\n")

if (!is.null(stemness_scores) && !all(is.na(stemness_scores))) {
  cat("Stemness scores calculated:", length(stemness_scores), "\n")
  cat("Top stemness scores:\n")
  top_scores <- head(stemness_scores[order(stemness_scores, decreasing = TRUE)], 3)
  for (i in seq_along(top_scores)) {
    cat(sprintf("  %d. %-8s: %.4f\n", i, names(top_scores)[i], top_scores[i]))
  }
} else {
  cat("Stemness scores: Not available\n")
}

cat("Output directory:", output_dir, "\n")
