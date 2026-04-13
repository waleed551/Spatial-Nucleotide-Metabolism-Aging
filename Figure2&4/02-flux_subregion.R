setwd("D:/Koziol_lab/Spatial-Nucleotide-Metabolism-Aging/Figure2")

library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
library(grid)
library(gridExtra)

rm(list = ls())
gc()

###############################################################################
# Global colors
###############################################################################
flux_color <- "#4DAF4A"
gene_color <- "#377EB8"
metabolite_color <- "#E6862D"
accent_text_color <- "#7A0000"

###############################################################################
# Global text sizes
###############################################################################
title_size <- 22
subtitle_size <- 15

panel_title_size <- 15
axis_title_size <- 13
axis_text_size <- 11
strip_text_size <- 11

pathway_text_size <- 5.2

###############################################################################
# 0. Output folders
###############################################################################
code1_root <- "Figure2_code1_whole_brain_results"
code1_tables <- file.path(code1_root, "01_tables")

output_root <- "Figure2_thalamus_integrated_k3_unique_results"
dir_tables <- file.path(output_root, "01_tables")
dir_figures <- file.path(output_root, "02_figures")

dirs <- c(output_root, dir_tables, dir_figures)
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

###############################################################################
# 1. Helper functions
###############################################################################
parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
}

convert_gene <- function(gene_symbols) {
  bitr(
    unique(gene_symbols),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
}

make_region_expression_annotated <- function(data_matrix, gene_module_metabolite, region) {
  gene_list <- unique(gene_module_metabolite$gene_symbol)
  
  rownames(data_matrix) <- paste(
    data_matrix$BrainRegion,
    data_matrix$mouse_id,
    data_matrix$sex,
    data_matrix$age,
    sep = "_"
  )
  
  region_data <- data_matrix[data_matrix$BrainRegion == region, , drop = FALSE]
  
  expression_only <- region_data[, !names(region_data) %in% c("BrainRegion", "mouse_id", "sex", "age"), drop = FALSE]
  expression_only <- as.matrix(expression_only)
  mode(expression_only) <- "numeric"
  
  expression_matrix <- t(expression_only)
  rownames(expression_matrix) <- colnames(expression_only)
  expression_matrix <- expression_matrix[rownames(expression_matrix) %in% gene_list, , drop = FALSE]
  
  expression_df <- as.data.frame(expression_matrix)
  expression_df$gene_symbol <- rownames(expression_df)
  
  expression_annotated <- expression_df %>%
    left_join(gene_module_metabolite, by = "gene_symbol") %>%
    relocate(gene_symbol, metabolite, module)
  
  age_info <- region_data[, c("age"), drop = FALSE]
  sample_names <- rownames(region_data)
  age_info_df <- data.frame(Sample = sample_names, age = age_info$age)
  age_info_df$age <- as.numeric(as.character(age_info_df$age))
  age_info_df <- age_info_df[order(age_info_df$age), , drop = FALSE]
  
  sorted_samples <- age_info_df$Sample
  expression_annotated <- expression_annotated[, c("gene_symbol", "metabolite", "module", sorted_samples)]
  
  list(
    expression_annotated = expression_annotated,
    age_info_df = age_info_df
  )
}

get_thalamus_flux_cluster3 <- function(df_flux, module_cluster_df_from_code1, region = "Thalamus", k = 3) {
  cluster1_module <- module_cluster_df_from_code1[module_cluster_df_from_code1$cluster == 1, "module"]
  df_region <- df_flux[df_flux$region == region & df_flux$module %in% cluster1_module, , drop = FALSE]
  
  mat <- dcast(df_region, module ~ age, value.var = "expression", fun.aggregate = mean, na.rm = TRUE)
  rownames(mat) <- mat$module
  mat$module <- NULL
  
  wanted_cols <- intersect(c("young", "mid", "old"), colnames(mat))
  mat <- mat[, wanted_cols, drop = FALSE]
  
  initial_res <- pheatmap(
    as.matrix(mat),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    silent = TRUE
  )
  
  module_clusters <- cutree(initial_res$tree_row, k = k)
  
  cluster_df <- data.frame(
    module = names(module_clusters),
    cluster = as.integer(module_clusters),
    stringsAsFactors = FALSE
  )
  
  list(
    flux_df = df_region,
    flux_matrix = mat,
    cluster_df = cluster_df
  )
}

make_flux_trend_df <- function(df_flux_region, cluster_df) {
  flux_age_labels <- c("young" = "3M", "mid" = "12M", "old" = "21M")
  
  df_flux_region %>%
    inner_join(cluster_df, by = "module") %>%
    group_by(cluster, age, module) %>%
    summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    group_by(cluster, age) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      sd_expression = sd(expression, na.rm = TRUE),
      n = dplyr::n(),
      se_expression = ifelse(n > 1, sd_expression / sqrt(n), 0),
      .groups = "drop"
    ) %>%
    mutate(
      age = factor(age, levels = c("young", "mid", "old")),
      age_label = factor(flux_age_labels[as.character(age)], levels = c("3M", "12M", "21M"))
    )
}

###############################################################################
# 1A. Unique gene helper functions
###############################################################################
make_cluster_gene_membership <- function(cluster_df, gene_module_metabolite) {
  cluster_gene_df <- cluster_df %>%
    inner_join(gene_module_metabolite[, c("module", "gene_symbol", "metabolite")], by = "module") %>%
    dplyr::select(cluster, module, gene_symbol, metabolite) %>%
    distinct()
  
  cluster_gene_presence <- cluster_gene_df %>%
    dplyr::select(cluster, gene_symbol) %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarise(
      n_clusters = n_distinct(cluster),
      cluster_membership = paste(sort(unique(cluster)), collapse = ";"),
      .groups = "drop"
    )
  
  specific_gene_list <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
    genes_cl <- unique(cluster_gene_df$gene_symbol[cluster_gene_df$cluster == cl])
    other_genes <- unique(cluster_gene_df$gene_symbol[cluster_gene_df$cluster != cl])
    setdiff(genes_cl, other_genes)
  })
  names(specific_gene_list) <- paste0("Cluster", sort(unique(cluster_df$cluster)))
  
  list(
    cluster_gene_df = cluster_gene_df,
    cluster_gene_presence = cluster_gene_presence,
    specific_gene_list = specific_gene_list
  )
}

make_gene_trend_unique_all <- function(expression_annotated,
                                       age_info_df,
                                       unique_genes) {
  if (length(unique_genes) == 0) return(NULL)
  
  sub_expr <- expression_annotated[!duplicated(expression_annotated$gene_symbol), , drop = FALSE]
  if (nrow(sub_expr) == 0) return(NULL)
  
  expr_mat <- sub_expr[, !(colnames(sub_expr) %in% c("gene_symbol", "metabolite", "module")), drop = FALSE]
  rownames(expr_mat) <- sub_expr$gene_symbol
  
  expr_long <- as.data.frame(expr_mat)
  expr_long$gene_symbol <- rownames(expr_long)
  
  expr_long <- pivot_longer(
    expr_long,
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Expression"
  )
  
  expr_long <- expr_long %>%
    left_join(age_info_df, by = c("Sample" = "Sample")) %>%
    group_by(gene_symbol, age) %>%
    summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    filter(gene_symbol %in% unique_genes) %>%
    group_by(gene_symbol) %>%
    mutate(
      z = if (sd(Expression, na.rm = TRUE) == 0 || is.na(sd(Expression, na.rm = TRUE))) {
        0
      } else {
        as.numeric(scale(Expression))
      }
    ) %>%
    ungroup()
  
  expr_long$age <- factor(expr_long$age, levels = c(3, 12, 15, 18, 21, 26, 28))
  expr_long$age_label <- factor(
    paste0(as.character(expr_long$age), "M"),
    levels = c("3M", "12M", "15M", "18M", "21M", "26M", "28M")
  )
  expr_long$gene_symbol <- factor(expr_long$gene_symbol, levels = unique_genes)
  
  expr_long
}

make_metabolite_trend_df_from_unique_genes <- function(met_expr_df_all,
                                                       gene_module_metabolite,
                                                       unique_genes,
                                                       region = "Thalamus") {
  if (length(unique_genes) == 0) return(NULL)
  
  linked_metabolites <- gene_module_metabolite %>%
    filter(gene_symbol %in% unique_genes) %>%
    distinct(gene_symbol, metabolite)
  
  if (nrow(linked_metabolites) == 0) return(NULL)
  
  met_rank <- linked_metabolites %>%
    group_by(metabolite) %>%
    summarise(n_linked_genes = n_distinct(gene_symbol), .groups = "drop") %>%
    arrange(desc(n_linked_genes), metabolite)
  
  ordered_mets <- met_rank$metabolite
  
  met_expr_df <- met_expr_df_all[met_expr_df_all$region == region, , drop = FALSE]
  
  met_expr_mean <- met_expr_df %>%
    filter(metabolite %in% ordered_mets) %>%
    group_by(metabolite, age) %>%
    summarise(mean_intensity = mean(mean_intensity, na.rm = TRUE), .groups = "drop")
  
  if (nrow(met_expr_mean) == 0) return(NULL)
  
  met_expr_mean <- met_expr_mean %>%
    left_join(met_rank, by = "metabolite") %>%
    group_by(metabolite) %>%
    mutate(
      z = if (sd(mean_intensity, na.rm = TRUE) == 0 || is.na(sd(mean_intensity, na.rm = TRUE))) {
        0
      } else {
        as.numeric(scale(mean_intensity))
      }
    ) %>%
    ungroup()
  
  met_expr_mean$age <- factor(met_expr_mean$age, levels = c("6M", "15M", "21M"))
  met_expr_mean$age_label <- factor(as.character(met_expr_mean$age), levels = c("6M", "15M", "21M"))
  met_expr_mean$metabolite <- factor(met_expr_mean$metabolite, levels = ordered_mets)
  
  met_expr_mean
}

get_top_go_terms <- function(gene_symbols, top_n = 4) {
  genes_df <- convert_gene(gene_symbols)
  if (is.null(genes_df) || nrow(genes_df) == 0) return("No significant GO terms")
  
  ego <- enrichGO(
    gene = genes_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return("No significant GO terms")
  
  ego <- simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  ego_df <- as.data.frame(ego)
  if (is.null(ego_df) || nrow(ego_df) == 0) return("No significant GO terms")
  
  ego_df <- ego_df[order(ego_df$p.adjust, -ego_df$Count), , drop = FALSE]
  terms <- head(ego_df$Description, top_n)
  paste(terms, collapse = "\n")
}

make_cluster_summary_table <- function(cluster_id,
                                       modules_in_cluster,
                                       gene_module_metabolite,
                                       gene_trend_df,
                                       met_trend_df,
                                       go_text,
                                       unique_genes) {
  genes <- unique(gene_module_metabolite$gene_symbol[gene_module_metabolite$module %in% modules_in_cluster])
  mets <- unique(gene_module_metabolite$metabolite[gene_module_metabolite$module %in% modules_in_cluster])
  
  data.frame(
    cluster = cluster_id,
    n_modules = length(unique(modules_in_cluster)),
    n_genes = length(genes),
    n_metabolites = length(mets),
    n_unique_genes = length(unique_genes),
    shown_unique_genes = paste(unique(gene_trend_df$gene_symbol), collapse = "; "),
    shown_unique_linked_metabolites = paste(unique(met_trend_df$metabolite), collapse = "; "),
    top_go_terms = gsub("\n", "; ", go_text),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 1B. Plot functions
###############################################################################
plot_gene_panel <- function(gene_df, cluster_id) {
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    return(
      ggplot() + theme_void() + ggtitle(paste0("Cluster ", cluster_id, " unique genes"))
    )
  }
  
  ggplot(gene_df, aes(x = age_label, y = z, group = gene_symbol)) +
    geom_line(color = gene_color, linewidth = 1.05) +
    geom_point(color = gene_color, size = 2.4) +
    facet_wrap(~ gene_symbol, ncol = 2, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Cluster ", cluster_id, ": Unique gene trends"),
      x = "Age",
      y = "Z-score"
    ) +
    theme(
      plot.title = element_text(size = panel_title_size, face = "bold"),
      axis.title.x = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.title.y = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = axis_text_size),
      strip.text = element_text(size = strip_text_size, face = "bold")
    )
}

plot_flux_panel <- function(flux_df, cluster_id) {
  sub_df <- flux_df[flux_df$cluster == cluster_id, , drop = FALSE]
  if (nrow(sub_df) == 0) {
    return(
      ggplot() + theme_void() + ggtitle(paste0("Cluster ", cluster_id, " flux"))
    )
  }
  
  ggplot(sub_df, aes(x = age_label, y = mean_expression, group = 1)) +
    geom_line(color = flux_color, linewidth = 1.25) +
    geom_point(color = flux_color, size = 3) +
    geom_errorbar(
      aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression),
      width = 0.15,
      color = flux_color,
      linewidth = 0.9
    ) +
    theme_bw() +
    labs(
      title = paste0("Cluster ", cluster_id, ": Flux trend"),
      x = "Age",
      y = "Mean flux"
    ) +
    theme(
      plot.title = element_text(size = panel_title_size, face = "bold"),
      axis.title.x = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.title.y = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size)
    )
}

plot_metabolite_panel <- function(met_df, cluster_id) {
  if (is.null(met_df) || nrow(met_df) == 0) {
    return(
      ggplot() + theme_void() + ggtitle(paste0("Cluster ", cluster_id, " unique-linked metabolites"))
    )
  }
  
  ggplot(met_df, aes(x = age_label, y = z, group = metabolite)) +
    geom_line(color = metabolite_color, linewidth = 1.05) +
    geom_point(color = metabolite_color, size = 2.4) +
    facet_wrap(~ metabolite, ncol = 2, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Cluster ", cluster_id, ": Metabolites linked to unique genes"),
      x = "Age",
      y = "Z-score"
    ) +
    theme(
      plot.title = element_text(size = panel_title_size, face = "bold"),
      axis.title.x = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.title.y = element_text(size = axis_title_size, color = accent_text_color, face = "bold"),
      axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = axis_text_size),
      strip.text = element_text(size = strip_text_size, face = "bold")
    )
}

plot_pathway_panel <- function(go_text, cluster_id, n_unique) {
  df_txt <- data.frame(
    x = 0,
    y = 0,
    label = paste0(
      "Cluster ", cluster_id, ": Pathway summary\n\n",
      "Unique genes: ", n_unique, "\n\n",
      go_text
    ),
    stringsAsFactors = FALSE
  )
  
  ggplot(df_txt, aes(x = x, y = y, label = label)) +
    geom_text(
      hjust = 0, vjust = 1,
      size = pathway_text_size,
      lineheight = 1.15,
      color = "black"
    ) +
    xlim(0, 1) +
    ylim(-1, 1) +
    theme_void()
}

###############################################################################
# 2. Input data
###############################################################################
gene_module_metabolite <- read.csv(
  file.path(code1_tables, "gene_module_metabolite.csv"),
  check.names = FALSE
)

module_cluster_df_from_code1 <- read.csv(
  file.path(code1_tables, "module_cluster_df.csv"),
  check.names = FALSE
)

df_flux <- read.csv("./region_flux_csv/all_flux_data_long.csv")
data_matrix <- read.csv("normalized_counts_with_metadata.csv")
met_expr_df_all <- read_excel("sum_normalized_data.xlsx")
met_expr_df_all <- as.data.frame(met_expr_df_all)

# Unify metabolite column name
if ("rNs" %in% colnames(met_expr_df_all)) {
  colnames(met_expr_df_all)[colnames(met_expr_df_all) == "rNs"] <- "metabolite"
}

###############################################################################
# 3. Build Thalamus k=3 cluster result
###############################################################################
region <- "Thalamus"

thalamus_obj <- get_thalamus_flux_cluster3(
  df_flux = df_flux,
  module_cluster_df_from_code1 = module_cluster_df_from_code1,
  region = region,
  k = 3
)

cluster_df <- thalamus_obj$cluster_df
df_flux_region <- thalamus_obj$flux_df
flux_matrix <- thalamus_obj$flux_matrix

write.csv(cluster_df, file.path(dir_tables, "clusters_Thalamus_k3.csv"), row.names = FALSE)
write.csv(flux_matrix, file.path(dir_tables, "flux_matrix_Thalamus_k3.csv"))

###############################################################################
# 4. Expression annotation
###############################################################################
expr_obj <- make_region_expression_annotated(
  data_matrix = data_matrix,
  gene_module_metabolite = gene_module_metabolite,
  region = region
)

expression_annotated <- expr_obj$expression_annotated
age_info_df <- expr_obj$age_info_df

write.csv(expression_annotated, file.path(dir_tables, "filtered_expression_annotated_Thalamus.csv"), row.names = FALSE)
write.csv(age_info_df, file.path(dir_tables, "sample_age_info_Thalamus.csv"), row.names = FALSE)

###############################################################################
# 5. Flux trend summary
###############################################################################
flux_trend_df <- make_flux_trend_df(
  df_flux_region = df_flux_region,
  cluster_df = cluster_df
)

write.csv(flux_trend_df, file.path(dir_tables, "flux_trend_Thalamus_k3.csv"), row.names = FALSE)

###############################################################################
# 6. Unique gene membership
###############################################################################
gene_membership_obj <- make_cluster_gene_membership(
  cluster_df = cluster_df,
  gene_module_metabolite = gene_module_metabolite
)

cluster_gene_df <- gene_membership_obj$cluster_gene_df
cluster_gene_presence <- gene_membership_obj$cluster_gene_presence
specific_gene_list <- gene_membership_obj$specific_gene_list

write.csv(cluster_gene_df, file.path(dir_tables, "cluster_gene_membership_Thalamus_k3.csv"), row.names = FALSE)
write.csv(cluster_gene_presence, file.path(dir_tables, "cluster_gene_presence_Thalamus_k3.csv"), row.names = FALSE)

for (nm in names(specific_gene_list)) {
  write.csv(
    data.frame(gene_symbol = specific_gene_list[[nm]]),
    file.path(dir_tables, paste0(nm, "_unique_genes.csv")),
    row.names = FALSE
  )
}

###############################################################################
# 7. Build integrated panels
###############################################################################
panel_list <- list()
summary_table_all <- list()

for (cluster_id in sort(unique(cluster_df$cluster))) {
  
  modules_in_cluster <- cluster_df$module[cluster_df$cluster == cluster_id]
  unique_genes_this <- specific_gene_list[[paste0("Cluster", cluster_id)]]
  
  gene_trend_df <- make_gene_trend_unique_all(
    expression_annotated = expression_annotated,
    age_info_df = age_info_df,
    unique_genes = unique_genes_this
  )
  
  met_trend_df <- make_metabolite_trend_df_from_unique_genes(
    met_expr_df_all = met_expr_df_all,
    gene_module_metabolite = gene_module_metabolite,
    unique_genes = unique_genes_this,
    region = region
  )
  
  genes_for_go <- unique(gene_module_metabolite$gene_symbol[gene_module_metabolite$module %in% modules_in_cluster])
  go_text <- get_top_go_terms(genes_for_go, top_n = 4)
  
  cluster_summary <- make_cluster_summary_table(
    cluster_id = cluster_id,
    modules_in_cluster = modules_in_cluster,
    gene_module_metabolite = gene_module_metabolite,
    gene_trend_df = if (is.null(gene_trend_df)) data.frame(gene_symbol = character()) else gene_trend_df,
    met_trend_df = if (is.null(met_trend_df)) data.frame(metabolite = character()) else met_trend_df,
    go_text = go_text,
    unique_genes = unique_genes_this
  )
  summary_table_all[[paste0("Cluster", cluster_id)]] <- cluster_summary
  
  p_gene <- plot_gene_panel(gene_trend_df, cluster_id)
  p_flux <- plot_flux_panel(flux_trend_df, cluster_id)
  p_met <- plot_metabolite_panel(met_trend_df, cluster_id)
  p_path <- plot_pathway_panel(go_text, cluster_id, length(unique_genes_this))
  
  row_plot <- p_gene + p_flux + p_met + p_path +
    plot_layout(widths = c(1.8, 1.0, 1.4, 1.1))
  
  panel_list[[paste0("Cluster", cluster_id)]] <- row_plot
}

summary_table_all <- bind_rows(summary_table_all)
write.csv(summary_table_all, file.path(dir_tables, "Thalamus_k3_integrated_summary_table_unique.csv"), row.names = FALSE)

###############################################################################
# 8. Combine all panels
###############################################################################
final_plot <- wrap_plots(panel_list, ncol = 1) +
  plot_annotation(
    title = "Integrated multi-layer view of Thalamus clusters (k = 3)",
    subtitle = "Gene (blue) -> Flux (green) -> Metabolite (orange) -> Pathway convergence",
    theme = theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = subtitle_size, hjust = 0.5)
    )
  )

ggsave(
  filename = file.path(dir_figures, "Integrated_Thalamus_k3_unique_multilayer.pdf"),
  plot = final_plot,
  width = 20,
  height = 18
)

ggsave(
  filename = file.path(dir_figures, "Integrated_Thalamus_k3_unique_multilayer.png"),
  plot = final_plot,
  width = 20,
  height = 18,
  dpi = 400
)

message("Done. Outputs saved in: ", output_root)
