# Figure 1 Analysis Script - Final Version
# Complete analysis pipeline for main figure with panel labels

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(ggalluvial)
  library(RColorBrewer)
})

# Check optional libraries
optional_libs <- c("scplotter", "SCP", "scCustomize")
for (lib in optional_libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    cat("Warning:", lib, "not available. Some functions may not work.\n")
  } else {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}

# Set directories
base_dir <- "/Users/fujitakyouhei/Desktop/Heart"
annotation_dir <- file.path(base_dir, "annotation", "data")
cm_dir <- file.path(base_dir, "CM_subset", "data")
output_dir <- file.path(base_dir, "Figure1")
figures_output <- file.path(output_dir, "figures")
stats_output <- file.path(output_dir, "statistics")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_output, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_output, recursive = TRUE, showWarnings = FALSE)

# Load data
cat("Loading data...\n")
Heart <- readRDS(file.path(annotation_dir, "annotation_seurat.rds"))
CM <- readRDS(file.path(cm_dir, "CM_subset.rds"))

# ==============================================================================
# Fig.1a - Heart UMAP by type and CellType
# ==============================================================================
cat("Creating Fig.1a - Heart UMAP...\n")

if (exists("CellDimPlot")) {
  CellDimPlot(
    srt = Heart, group.by = c("type", "CellType"),
    reduction = "umap", theme_use = "theme_blank"
  )
}

# ==============================================================================
# Fig.1b - Feature statistics for Heart cell types
# ==============================================================================
cat("Creating Fig.1b - Feature statistics...\n")

if (exists("FeatureStatPlot")) {
  FeatureStatPlot(Heart, 
                  stat.by = c("Ttn", "Ankrd1", "Pecam1", "Egfl7", "Pdgfra", "Gsn", 
                              "Slc39a8", "Wt1", "F13a1", "Mertk", "Ms4a1", "Rgs5", "Abcc9"),
                  group.by = "CellType", add_box = TRUE, stack = TRUE)
}

# ==============================================================================
# Fig.1c - CM UMAP by SubCellType
# ==============================================================================
cat("Creating Fig.1c - CM UMAP...\n")

if (exists("CellDimPlot")) {
  CellDimPlot(
    srt = CM, group.by = "SubCellType",
    reduction = "umap", theme_use = "theme_blank"
  )
}

# Create SubType classification
CM$SubType <- dplyr::recode(CM$SubCellType,
                            "CMC_1" = "CMC",
                            "CMC_2" = "CMC",
                            "CMB_1" = "CMB",
                            "CMB_2" = "CMB",
                            "CMB_3" = "CMB",
                            "CMB_4" = "CMB",
                            "CMB_5" = "CMB")

Idents(CM) <- CM$SubType
cat("SubType levels:", levels(CM), "\n")
cat("SubType counts:\n")
print(table(Idents(CM)))

# ==============================================================================
# Fig.1d - Three-stage alluvial plot
# ==============================================================================
cat("Creating Fig.1d - Alluvial plot...\n")

# Function to create three-stage alluvial plot
create_three_stage_alluvial <- function(seurat_obj, 
                                        stage1_col = "type", 
                                        stage2_col = "CellType",
                                        stage3_col = "SubCellType",
                                        color_palette = "custom") {
  
  # Extract metadata for three stages
  plot_data <- seurat_obj@meta.data %>%
    select(all_of(c(stage1_col, stage2_col, stage3_col))) %>%
    filter(!is.na(.data[[stage1_col]]) & 
             !is.na(.data[[stage2_col]]) & 
             !is.na(.data[[stage3_col]])) %>%
    group_by(.data[[stage1_col]], .data[[stage2_col]], .data[[stage3_col]]) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Color palette selection
  n_colors <- length(unique(plot_data[[stage1_col]]))
  
  if (color_palette == "custom") {
    colors <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#4CAF50",
                "#9C27B0", "#FF9800", "#607D8B", "#795548", "#E91E63",
                "#3F51B5", "#009688", "#8BC34A", "#FF5722", "#673AB7",
                "#FFC107", "#00BCD4", "#CDDC39", "#FF9E80", "#B39DDB")
  } else if (color_palette == "viridis") {
    colors <- viridis::viridis(n_colors, option = "D")
  } else if (color_palette == "plasma") {
    colors <- viridis::plasma(n_colors)
  } else if (color_palette == "set1") {
    colors <- RColorBrewer::brewer.pal(min(9, n_colors), "Set1")
    if (n_colors > 9) {
      colors <- c(colors, RColorBrewer::brewer.pal(min(8, n_colors-9), "Set2"))
    }
  } else {
    colors <- rainbow(n_colors)
  }
  
  # Create plot (3 stages)
  p <- ggplot(plot_data,
              aes(axis1 = .data[[stage1_col]], 
                  axis2 = .data[[stage2_col]], 
                  axis3 = .data[[stage3_col]], 
                  y = count)) +
    geom_alluvium(aes(fill = .data[[stage1_col]]), 
                  alpha = 0.8, 
                  width = 0.12,
                  color = "white",
                  size = 0.3) +
    geom_stratum(width = 0.12, 
                 fill = "grey95", 
                 color = "black",
                 size = 0.7) +
    geom_text(stat = "stratum", 
              aes(label = after_stat(stratum)), 
              size = 3.5,
              fontface = "bold",
              color = "black") +
    scale_x_discrete(limits = c(stage1_col, stage2_col, stage3_col), 
                     expand = c(0.1, 0.1)) +
    scale_fill_manual(values = colors[1:n_colors]) +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 14, face = "bold", color = "black", 
                                 margin = margin(t = 10)),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.size = unit(0.8, "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5,
                                margin = margin(b = 20)),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30",
                                   margin = margin(b = 30)),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = "Three-Stage Cell Type Classification",
      subtitle = paste0(stage1_col, " -> ", stage2_col, " -> ", stage3_col),
      fill = stage1_col
    ) +
    guides(fill = guide_legend(title.position = "top", 
                               title.hjust = 0.5,
                               nrow = 2))
  
  return(p)
}

# Create and save three-stage plot
pub_plot_3stage <- create_three_stage_alluvial(CM)
ggsave(file.path(figures_output, "Fig1d_alluvial_plot.pdf"), pub_plot_3stage, 
       width = 12, height = 12, dpi = 300, device = cairo_pdf)

# ==============================================================================
# Fig.1e - CM feature statistics by SubCellType
# ==============================================================================
cat("Creating Fig.1e - CM feature statistics...\n")

if (exists("FeatureStatPlot")) {
  FeatureStatPlot(CM, 
                  stat.by = c("Ttn", "Ankrd1", "Actn2", "Tnnt2", "Bmi1", "Hopx", "Nkx2-5", "Hand2", "Abcg2", "Ly6a"),
                  group.by = "SubCellType", add_box = TRUE, stack = TRUE)
}

# ==============================================================================
# Fig.1f - Comparative feature plots (Heart vs CM)
# ==============================================================================
cat("Creating Fig.1f - Comparative feature plots...\n")

# Define genes for comparison
genes <- c("Ttn", "Ankrd1", "Actn2", "Tnnt2", "Mybpc3", "Myh6",
           "Hopx", "Nkx2-5", "Gata4", "Tbx5", "Tbx20", "Mef2c",
           "Bmi1", "Abcg2", "Ly6a", "Pecam1", "Twist2", "Kit")

if (exists("FeaturePlot_scCustom")) {
  # Create all plots using loop
  plots <- list()
  for(i in 1:length(genes)) {
    gene <- genes[i]
    plots[[i*2-1]] <- FeaturePlot_scCustom(seurat_object = Heart, features = gene, order = TRUE)
    plots[[i*2]] <- FeaturePlot_scCustom(seurat_object = CM, features = gene, order = TRUE)
  }
  
  # Combine plots in 6x6 grid
  combined_plot <- wrap_plots(plots, ncol = 6)
  print(combined_plot)
}

# ==============================================================================
# Fig.1g - Statistical comparison between CMC and CMB (scplotter)
# ==============================================================================
cat("Creating Fig.1g - Statistical comparison...\n")
detach("package:SCP", unload = TRUE)

# Define genes for statistical analysis
genes_to_plot <- c("Tnnt2", "Ttn", "Mybpc3", "Hopx", "Nkx2-5", "Bmi1")

if (exists("FeatureStatPlot")) {
  # Create feature stat plot
  p <- FeatureStatPlot(
    CM,
    features = genes_to_plot,
    ident = "SubType",
    add_box = FALSE,             
    add_point = FALSE,
    palette = "Dark2",           
    facet_ncol = 3,
    facet_scales = "free_y",
    theme = ggplot2::theme_linedraw, 
    comparisons = TRUE
  )
  
  # Add custom boxplot
  p1 <- p + 
    ggplot2::theme(panel.grid = ggplot2::element_blank()) + 
    ggplot2::geom_boxplot(
      fill = "white",
      color = "black",
      width = 0.1,
      outlier.shape = NA,
      inherit.aes = TRUE
    )
  
  print(p1)
}

# ==============================================================================
# Statistical Analysis - Two-group comparison
# ==============================================================================
cat("Performing statistical analysis...\n")

# Initialize results dataframe
results <- data.frame()

for (gene in genes_to_plot) {
  cat("Processing gene:", gene, "\n")
  
  # Extract expression data
  df <- FetchData(CM, vars = c(gene, "SubType"))
  colnames(df) <- c("expr", "SubType")
  
  # Data preprocessing
  df <- df[complete.cases(df), ]
  df$SubType <- as.factor(df$SubType)
  
  # Check sample sizes
  group_counts <- table(df$SubType)
  cat("Gene:", gene, "\n")
  print(group_counts)
  
  # Separate CMC and CMB data
  cmc_data <- df$expr[df$SubType == "CMC"]
  cmb_data <- df$expr[df$SubType == "CMB"]
  
  # Calculate basic statistics
  cmc_mean <- mean(cmc_data)
  cmc_median <- median(cmc_data)
  cmc_sd <- sd(cmc_data)
  cmb_mean <- mean(cmb_data)
  cmb_median <- median(cmb_data)
  cmb_sd <- sd(cmb_data)
  
  # Normality testing
  if (length(cmc_data) >= 3 && length(cmb_data) >= 3) {
    shapiro_cmc <- shapiro.test(sample(cmc_data, min(5000, length(cmc_data))))
    shapiro_cmb <- shapiro.test(sample(cmb_data, min(5000, length(cmb_data))))
    normal_cmc <- shapiro_cmc$p.value > 0.05
    normal_cmb <- shapiro_cmb$p.value > 0.05
  } else {
    normal_cmc <- FALSE
    normal_cmb <- FALSE
  }
  
  # Statistical test selection and execution
  if (normal_cmc && normal_cmb) {
    # Both groups follow normal distribution: t-test
    var_test <- var.test(cmc_data, cmb_data)
    equal_var <- var_test$p.value > 0.05
    
    t_result <- t.test(cmc_data, cmb_data, var.equal = equal_var)
    test_used <- ifelse(equal_var, "Student t-test", "Welch t-test")
    statistic <- t_result$statistic
    p_value <- t_result$p.value
    
  } else {
    # Non-normal distribution: Wilcoxon rank sum test
    wilcox_result <- wilcox.test(cmc_data, cmb_data)
    test_used <- "Wilcoxon rank sum test"
    statistic <- wilcox_result$statistic
    p_value <- wilcox_result$p.value
  }
  
  # Effect size calculation
  if (normal_cmc && normal_cmb) {
    # Cohen's d
    pooled_sd <- sqrt(((length(cmc_data) - 1) * var(cmc_data) + 
                         (length(cmb_data) - 1) * var(cmb_data)) / 
                        (length(cmc_data) + length(cmb_data) - 2))
    effect_size <- abs(cmc_mean - cmb_mean) / pooled_sd
    effect_measure <- "Cohen's d"
  } else {
    # r (effect size for Mann-Whitney)
    z_score <- qnorm(p_value/2, lower.tail = FALSE)
    effect_size <- abs(z_score) / sqrt(length(cmc_data) + length(cmb_data))
    effect_measure <- "r"
  }
  
  # Add results to dataframe
  results <- rbind(results, data.frame(
    gene = gene,
    test_used = test_used,
    statistic = statistic,
    p_value = p_value,
    cmc_n = length(cmc_data),
    cmb_n = length(cmb_data),
    cmc_mean = cmc_mean,
    cmc_median = cmc_median,
    cmc_sd = cmc_sd,
    cmb_mean = cmb_mean,
    cmb_median = cmb_median,
    cmb_sd = cmb_sd,
    effect_size = effect_size,
    effect_measure = effect_measure,
    normal_cmc = normal_cmc,
    normal_cmb = normal_cmb
  ))
}

# Multiple testing correction
results$p_adjusted_bonferroni <- p.adjust(results$p_value, method = "bonferroni")
results$p_adjusted_bh <- p.adjust(results$p_value, method = "BH")

# Significance determination
results$significant_raw <- results$p_value < 0.05
results$significant_bh <- results$p_adjusted_bh < 0.05
results$significant_bonferroni <- results$p_adjusted_bonferroni < 0.05

# Results summary
cat("\n=== Statistical Analysis Results Summary ===\n")
print(results[, c("gene", "test_used", "p_value", "p_adjusted_bh", 
                  "cmc_mean", "cmb_mean", "effect_size", "effect_measure")])

cat("\n=== Significant genes (BH corrected p < 0.05) ===\n")
sig_genes <- results[results$p_adjusted_bh < 0.05, ]
if (nrow(sig_genes) > 0) {
  print(sig_genes[, c("gene", "p_value", "p_adjusted_bh", "cmc_mean", "cmb_mean")])
} else {
  cat("No significant genes found\n")
}

# Save results
write.csv(results, file.path(stats_output, "two_group_comparison_results.csv"), row.names = FALSE)

cat("\n=== FIGURE 1 ANALYSIS COMPLETE ===\n")
cat("Results saved to:", file.path(stats_output, "two_group_comparison_results.csv"), "\n")
cat("Output directory:", output_dir, "\n")