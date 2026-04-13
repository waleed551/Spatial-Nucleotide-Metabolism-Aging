setwd("/Users/freya/Desktop/cibr/paper-draft/figures/Figure5")

library(tidyverse)
library(tibble)
library(readxl)

###############################################################################
# 0. Global settings
###############################################################################
theme_set(theme_minimal(base_family = "sans"))

gene_color <- "#377EB8"
flux_color <- "#4DAF4A"
metabolite_color <- "#E6862D"

# enlarged text sizes
base_size_global <- 24
plot_title_size <- 28
axis_title_size <- 22
axis_text_size <- 19
legend_title_size <- 20
legend_text_size <- 19

###############################################################################
# 1. Normalized entropy function
###############################################################################
calc_normalized_entropy <- function(vec) {
  vec <- as.numeric(vec)
  vec <- vec[!is.na(vec)]
  
  if (length(vec) == 0) return(NA_real_)
  
  # keep Shannon input non-negative for gene layer as well
  vec <- abs(vec)
  
  if (sum(vec, na.rm = TRUE) <= 0) return(NA_real_)
  
  p <- vec / sum(vec, na.rm = TRUE)
  p <- p[p > 0]
  
  if (length(p) <= 1) return(0)
  
  h <- -sum(p * log2(p))
  h / log2(length(p))
}

get_entropy_by_age <- function(df, label) {
  df <- as.matrix(df)
  
  if (ncol(df) == 0) {
    return(tibble(
      age = character(0),
      entropy = numeric(0),
      type = character(0)
    ))
  }
  
  entropy_vec <- apply(df, 2, calc_normalized_entropy)
  
  tibble(
    age = colnames(df),
    entropy = as.numeric(entropy_vec),
    type = label
  )
}

###############################################################################
# 2. File paths
###############################################################################
figure2_root <- "/Users/freya/Desktop/cibr/paper-draft/figures/Figure2/Figure2_thalamus_integrated_k3_unique_results"
figure2_tables <- file.path(figure2_root, "01_tables")

file_gene_expr_annotated <- file.path(figure2_tables, "filtered_expression_annotated_Thalamus.csv")
file_cluster_assignment <- file.path(figure2_tables, "clusters_Thalamus_k3.csv")

file_gene_module_metabolite <- "/Users/freya/Desktop/cibr/paper-draft/figures/Figure2/Figure2_code1_whole_brain_results/01_tables/gene_module_metabolite.csv"
file_metabolite_raw <- "/Users/freya/Desktop/cibr/paper-draft/figures/Figure2/sum_normalized_data.xlsx"
file_flux_raw <- "/Users/freya/Desktop/cibr/paper-draft/figures/Figure2/region_flux_csv/all_flux_data_long.csv"

###############################################################################
# 3. Input checks
###############################################################################
required_files <- c(
  file_gene_expr_annotated,
  file_cluster_assignment,
  file_gene_module_metabolite,
  file_metabolite_raw,
  file_flux_raw
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing files:\n", paste(missing_files, collapse = "\n"))
}

###############################################################################
# 4. Read data
###############################################################################
gene_expr_annotated <- read.csv(file_gene_expr_annotated, check.names = FALSE)
cluster_assignment <- read.csv(file_cluster_assignment, check.names = FALSE)
gene_module_metabolite <- read.csv(file_gene_module_metabolite, check.names = FALSE)
df_flux_raw <- read.csv(file_flux_raw, check.names = FALSE)

met_expr_df_all <- read_excel(file_metabolite_raw)
met_expr_df_all <- as.data.frame(met_expr_df_all)

if ("rNs" %in% colnames(met_expr_df_all)) {
  colnames(met_expr_df_all)[colnames(met_expr_df_all) == "rNs"] <- "metabolite"
}

###############################################################################
# 5. Cluster settings
###############################################################################
region <- "Thalamus"
clusters <- 1:3

all_entropy_data <- data.frame()

###############################################################################
# 6. Iterate through clusters
###############################################################################
for (clus in clusters) {
  
  prefix <- paste0(region, "_cluster_", clus)
  
  cat("\n==============================\n")
  cat("Processing", prefix, "\n")
  cat("==============================\n")
  
  ###########################################################################
  # 6.1 Modules in current cluster
  ###########################################################################
  modules_in_cluster <- cluster_assignment %>%
    dplyr::filter(cluster == clus) %>%
    dplyr::pull(module)
  
  modules_in_cluster <- unique(modules_in_cluster[!is.na(modules_in_cluster) & nzchar(modules_in_cluster)])
  
  cat("Modules:", length(modules_in_cluster), "\n")
  
  if (length(modules_in_cluster) == 0) {
    message("⛔ No modules found for: ", prefix)
    next
  }
  
  ###########################################################################
  # 6.2 All genes in current cluster
  ###########################################################################
  cluster_genes <- gene_module_metabolite %>%
    dplyr::filter(module %in% modules_in_cluster) %>%
    dplyr::distinct(gene_symbol) %>%
    dplyr::pull(gene_symbol)
  
  cluster_genes <- unique(cluster_genes[!is.na(cluster_genes) & nzchar(cluster_genes)])
  
  cat("Genes:", length(cluster_genes), "\n")
  
  if (length(cluster_genes) == 0) {
    message("⛔ No genes found for: ", prefix)
    next
  }
  
  ###########################################################################
  # 6.3 Gene matrix
  ###########################################################################
  gene_sub <- gene_expr_annotated %>%
    dplyr::filter(gene_symbol %in% cluster_genes)
  
  gene_age_cols <- setdiff(colnames(gene_sub), c("gene_symbol", "metabolite", "module"))
  
  gene_long <- gene_sub %>%
    dplyr::select(gene_symbol, all_of(gene_age_cols)) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop") %>%
    tidyr::pivot_longer(
      cols = -gene_symbol,
      names_to = "sample",
      values_to = "expression"
    )
  
  sample_age_map <- sapply(gene_long$sample, function(s) tail(strsplit(s, "_")[[1]], 1))
  sample_age_map <- as.numeric(sample_age_map)
  
  gene_long <- gene_long %>%
    dplyr::mutate(
      age_num = sample_age_map,
      age = paste0(age_num, "M")
    )
  
  gene_age_mean <- gene_long %>%
    dplyr::group_by(gene_symbol, age) %>%
    dplyr::summarise(expression = mean(as.numeric(expression), na.rm = TRUE), .groups = "drop")
  
  gene_entropy_input <- gene_age_mean %>%
    dplyr::select(gene_symbol, age, expression) %>%
    tidyr::pivot_wider(names_from = age, values_from = expression)
  
  gene_entropy_input <- as.data.frame(gene_entropy_input, check.names = FALSE)
  rownames(gene_entropy_input) <- gene_entropy_input$gene_symbol
  gene_entropy_input$gene_symbol <- NULL
  
  gene_age_order <- c("3M", "12M", "15M", "18M", "21M", "26M", "28M")
  for (ag in gene_age_order) {
    if (!ag %in% colnames(gene_entropy_input)) {
      gene_entropy_input[[ag]] <- NA_real_
    }
  }
  gene_entropy_input <- gene_entropy_input[, gene_age_order, drop = FALSE]
  
  ###########################################################################
  # 6.4 Metabolite matrix
  ###########################################################################
  linked_metabolites <- gene_module_metabolite %>%
    dplyr::filter(gene_symbol %in% cluster_genes) %>%
    dplyr::distinct(metabolite) %>%
    dplyr::pull(metabolite)
  
  linked_metabolites <- unique(linked_metabolites[!is.na(linked_metabolites) & nzchar(linked_metabolites)])
  
  cat("Metabolites:", length(linked_metabolites), "\n")
  
  if (length(linked_metabolites) == 0) {
    message("⛔ No metabolites found for: ", prefix)
    next
  }
  
  met_sub <- met_expr_df_all %>%
    dplyr::filter(region == "Thalamus", metabolite %in% linked_metabolites) %>%
    dplyr::group_by(metabolite, age) %>%
    dplyr::summarise(mean_intensity = mean(mean_intensity, na.rm = TRUE), .groups = "drop")
  
  metabolite_entropy_input <- met_sub %>%
    dplyr::select(metabolite, age, mean_intensity) %>%
    tidyr::pivot_wider(names_from = age, values_from = mean_intensity) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(metabolite_entropy_input) <- metabolite_entropy_input$metabolite
  metabolite_entropy_input$metabolite <- NULL
  
  met_age_order <- c("6M", "15M", "21M")
  for (ag in met_age_order) {
    if (!ag %in% colnames(metabolite_entropy_input)) {
      metabolite_entropy_input[[ag]] <- NA_real_
    }
  }
  metabolite_entropy_input <- metabolite_entropy_input[, met_age_order, drop = FALSE]
  
  ###########################################################################
  # 6.5 Flux matrix
  ###########################################################################
  flux_sub_long <- df_flux_raw %>%
    dplyr::filter(region == "Thalamus", module %in% modules_in_cluster)
  
  if (nrow(flux_sub_long) == 0) {
    message("⛔ No flux data found for: ", prefix)
    next
  }
  
  flux_entropy_input <- flux_sub_long %>%
    dplyr::select(module, age, expression) %>%
    tidyr::pivot_wider(names_from = age, values_from = expression) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(flux_entropy_input) <- flux_entropy_input$module
  flux_entropy_input$module <- NULL
  
  if ("young" %in% colnames(flux_entropy_input)) colnames(flux_entropy_input)[colnames(flux_entropy_input) == "young"] <- "6M"
  if ("mid" %in% colnames(flux_entropy_input)) colnames(flux_entropy_input)[colnames(flux_entropy_input) == "mid"] <- "18M"
  if ("old" %in% colnames(flux_entropy_input)) colnames(flux_entropy_input)[colnames(flux_entropy_input) == "old"] <- "21M"
  
  flux_age_order <- c("6M", "18M", "21M")
  for (ag in flux_age_order) {
    if (!ag %in% colnames(flux_entropy_input)) {
      flux_entropy_input[[ag]] <- NA_real_
    }
  }
  flux_entropy_input <- flux_entropy_input[, flux_age_order, drop = FALSE]
  
  ###########################################################################
  # 6.6 Calculate normalized entropy
  ###########################################################################
  entropy_gene <- get_entropy_by_age(gene_entropy_input, "Gene")
  entropy_metabolite <- get_entropy_by_age(metabolite_entropy_input, "Metabolite")
  entropy_flux <- get_entropy_by_age(flux_entropy_input, "Flux")
  
  cat("\n--- entropy_gene ---\n")
  print(entropy_gene)
  
  cat("\n--- entropy_metabolite ---\n")
  print(entropy_metabolite)
  
  cat("\n--- entropy_flux ---\n")
  print(entropy_flux)
  
  write.csv(entropy_gene, paste0("entropy_gene_", prefix, ".csv"), row.names = FALSE)
  write.csv(entropy_metabolite, paste0("entropy_metabolite_", prefix, ".csv"), row.names = FALSE)
  write.csv(entropy_flux, paste0("entropy_flux_", prefix, ".csv"), row.names = FALSE)
  
  ###########################################################################
  # 6.7 Merge
  ###########################################################################
  full_age_levels <- c("3M", "6M", "12M", "15M", "18M", "21M", "26M", "28M")
  
  entropy_all <- bind_rows(entropy_gene, entropy_metabolite, entropy_flux) %>%
    dplyr::mutate(
      age = factor(age, levels = full_age_levels),
      age_numeric = as.numeric(gsub("M", "", as.character(age)))
    ) %>%
    dplyr::filter(!is.na(entropy)) %>%
    dplyr::arrange(type, age_numeric)
  
  if (nrow(entropy_all) == 0) {
    message("⛔ Entropy result is empty for: ", prefix)
    next
  }
  
  entropy_all$region <- region
  entropy_all$cluster <- clus
  
  all_entropy_data <- bind_rows(all_entropy_data, entropy_all)
  
  ###########################################################################
  # 6.8 Plot with larger text ratio
  ###########################################################################
  p <- ggplot(entropy_all, aes(x = age, y = entropy, color = type, group = type)) +
    geom_line(linewidth = 2.4, aes(linetype = type)) +
    geom_point(size = 5, shape = 21, fill = "white", stroke = 1.8) +
    theme_minimal(base_size = base_size_global, base_family = "sans") +
    labs(
      title = paste(region, "Cluster", clus),
      x = "Age",
      y = "Normalized Shannon Entropy",
      color = "Layer",
      linetype = "Layer"
    ) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = legend_title_size, face = "bold", family = "sans"),
      legend.text = element_text(size = legend_text_size, family = "sans"),
      plot.title = element_text(size = plot_title_size, face = "bold", color = "darkblue", family = "sans"),
      axis.title = element_text(size = axis_title_size, face = "bold", color = "darkred", family = "sans"),
      axis.text = element_text(size = axis_text_size, color = "black", family = "sans"),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.7, linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      plot.margin = margin(18, 18, 18, 18)
    ) +
    scale_color_manual(
      values = c(
        "Gene" = gene_color,
        "Metabolite" = metabolite_color,
        "Flux" = flux_color
      )
    ) +
    coord_cartesian(ylim = c(0, 1))
  
  out_file_pdf <- paste0("normalized_entropy_plot_", prefix, ".pdf")
  out_file_png <- paste0("normalized_entropy_plot_", prefix, ".png")
  
  ggsave(out_file_pdf, plot = p, device = cairo_pdf, width = 8.8, height = 6.8)
  ggsave(out_file_png, plot = p, width = 8.8, height = 6.8, dpi = 450)
  
  message("✅ Saved figure: ", out_file_pdf)
  message("✅ Saved figure: ", out_file_png)
  
  ###########################################################################
  # 6.9 Save merged entropy data
  ###########################################################################
  entropy_out_file <- paste0("normalized_entropy_all_", prefix, ".csv")
  write.csv(entropy_all, entropy_out_file, row.names = FALSE)
}

###############################################################################
# 7. Save combined entropy data
###############################################################################
final_entropy_out_file <- "all_normalized_entropy_combined.csv"
write.csv(all_entropy_data, final_entropy_out_file, row.names = FALSE)
message("✅ Saved combined normalized entropy data: ", final_entropy_out_file)
