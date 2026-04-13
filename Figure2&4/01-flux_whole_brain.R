setwd("D:/Koziol_lab/Chenyang_paper/My Updates/Code/Figure2")

library(readxl)
library(tools)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

###############################################################################
# 0. Output folders
###############################################################################
output_root <- "Figure2_code1_whole_brain_results"

dir_tables <- file.path(output_root, "01_tables")
dir_whole_brain <- file.path(output_root, "02_whole_brain")
dir_cluster_trends <- file.path(output_root, "03_cluster_trends")
dir_venn <- file.path(output_root, "04_venn")
dir_go <- file.path(output_root, "05_go")

dirs <- c(output_root, dir_tables, dir_whole_brain, dir_cluster_trends, dir_venn, dir_go)
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

###############################################################################
# 1. Helper functions
###############################################################################
standardize_gene_vector <- function(x) {
  x <- na.omit(x)
  x <- tools::toTitleCase(tolower(x))
  unique(x[nzchar(x)])
}

make_gene_metabolite_module <- function(gene_metabolite_df, module_gene_sets_list) {
  result_df <- data.frame(
    gene_symbol = character(),
    metabolite = character(),
    module = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(gene_metabolite_df))) {
    gene <- gene_metabolite_df$gene_symbol[i]
    metab <- gene_metabolite_df$metabolite[i]
    
    matched_modules <- names(module_gene_sets_list)[
      sapply(module_gene_sets_list, function(glist) gene %in% glist)
    ]
    
    if (length(matched_modules) > 0) {
      temp_df <- data.frame(
        gene_symbol = gene,
        metabolite = metab,
        module = matched_modules,
        stringsAsFactors = FALSE
      )
      result_df <- bind_rows(result_df, temp_df)
    }
  }
  
  result_df
}

make_whole_brain_flux_matrix <- function(df_flux, valid_modules) {
  df_flux <- df_flux[df_flux$module %in% valid_modules, , drop = FALSE]
  df_flux$sample <- paste(df_flux$region, df_flux$age, sep = "_")
  
  mat <- dcast(df_flux, module ~ sample, value.var = "expression")
  rownames(mat) <- mat$module
  mat <- mat[, -1, drop = FALSE]
  
  annotation_col <- data.frame(
    age = sapply(colnames(mat), function(x) strsplit(x, "_")[[1]][2]),
    region = sapply(colnames(mat), function(x) strsplit(x, "_")[[1]][1])
  )
  rownames(annotation_col) <- colnames(mat)
  
  annotation_col$age <- factor(annotation_col$age, levels = c("young", "mid", "old"))
  sorted_cols <- rownames(annotation_col)[order(annotation_col$age, annotation_col$region)]
  
  mat <- mat[, sorted_cols, drop = FALSE]
  annotation_col <- annotation_col[sorted_cols, , drop = FALSE]
  
  list(mat = mat, annotation_col = annotation_col)
}

plot_whole_brain_clustered_heatmap <- function(mat, annotation_col, outdir) {
  region_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(unique(annotation_col$region))),
    unique(annotation_col$region)
  )
  
  age_levels <- c("young", "mid", "old")
  age_colors <- setNames(brewer.pal(length(age_levels), "Dark2"), age_levels)
  
  ann_colors <- list(
    region = region_colors,
    age = age_colors
  )
  
  heatmap_colors <- rev(brewer.pal(11, "RdBu"))
  
  ph <- pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    annotation_colors = ann_colors,
    clustering_distance_rows = "euclidean",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    fontsize_col = 8,
    color = heatmap_colors,
    main = "Metabolic flux (with clustering)",
    silent = TRUE
  )
  
  row_clusters <- cutree(ph$tree_row, k = 3)
  module_cluster_df <- data.frame(
    module = names(row_clusters),
    cluster = as.factor(row_clusters)
  )
  
  row_anno <- data.frame(cluster = factor(row_clusters))
  rownames(row_anno) <- names(row_clusters)
  
  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_row = row_anno,
    cluster_cols = FALSE,
    annotation_colors = ann_colors,
    clustering_distance_rows = "euclidean",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    fontsize_col = 8,
    color = heatmap_colors,
    main = "Metabolic flux - Clustered",
    filename = file.path(outdir, "flux_clustered.pdf"),
    width = 6,
    height = 6
  )
  
  module_cluster_df
}

plot_sankey_by_cluster <- function(df, cluster_id, outdir) {
  if (nrow(df) == 0) return(NULL)
  
  df$metabolite <- factor(df$metabolite)
  df$gene_label <- factor(df$gene_label)
  df$module <- factor(df$module)
  
  p <- ggplot(df, aes(axis1 = metabolite, axis2 = gene_label, axis3 = module)) +
    geom_alluvium(aes(fill = metabolite), width = 1/12, alpha = 0.8) +
    geom_stratum(width = 1/12, fill = "gray90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, hjust = 0) +
    scale_x_discrete(limits = c("Metabolite", "Gene", "Module"), expand = c(.05, .05)) +
    theme_minimal() +
    ggtitle(paste("Sankey: Cluster", cluster_id)) +
    theme(axis.text.y = element_blank())
  
  ggsave(
    file.path(outdir, paste0("sankey_cluster_", cluster_id, ".pdf")),
    plot = p,
    width = 10,
    height = 20
  )
}

make_expression_annotated <- function(data_matrix, gene_module_metabolite) {
  gene_list <- unique(gene_module_metabolite$gene_symbol)
  
  rownames(data_matrix) <- paste(
    data_matrix$BrainRegion,
    data_matrix$mouse_id,
    data_matrix$sex,
    data_matrix$age,
    sep = "_"
  )
  
  expression_only <- data_matrix[, !names(data_matrix) %in% c("BrainRegion", "mouse_id", "sex", "age"), drop = FALSE]
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
  
  age_info <- data_matrix[, c("age"), drop = FALSE]
  sample_names <- rownames(data_matrix)
  age_info_df <- data.frame(Sample = sample_names, age = age_info$age)
  age_info_df$age <- as.numeric(as.character(age_info_df$age))
  age_info_df <- age_info_df[order(age_info_df$age), , drop = FALSE]
  
  sorted_samples <- age_info_df$Sample
  expression_annotated <- expression_annotated[, c("gene_symbol", "metabolite", "module", sorted_samples)]
  
  list(expression_annotated = expression_annotated, age_info_df = age_info_df)
}

plot_cluster_flux_heatmap <- function(cluster_df, cluster_id, outdir) {
  if (nrow(cluster_df) == 0) return(NULL)
  
  mat_cluster <- dcast(cluster_df, module ~ age, value.var = "expression", fun.aggregate = mean, na.rm = TRUE)
  rownames(mat_cluster) <- mat_cluster$module
  mat_cluster$module <- NULL
  
  wanted_cols <- intersect(c("young", "mid", "old"), colnames(mat_cluster))
  mat_cluster <- mat_cluster[, wanted_cols, drop = FALSE]
  
  if (ncol(mat_cluster) == 0 || nrow(mat_cluster) == 0) return(NULL)
  
  mat_scaled <- t(scale(t(as.matrix(mat_cluster))))
  
  annotation_col_cluster <- data.frame(age = factor(colnames(mat_scaled), levels = c("young", "mid", "old")))
  rownames(annotation_col_cluster) <- colnames(mat_scaled)
  
  write.csv(mat_scaled, file.path(outdir, paste0("flux_", cluster_id, ".csv")))
  
  pheatmap(
    as.matrix(mat_cluster),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    annotation_col = annotation_col_cluster,
    color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
    fontsize = 10,
    main = paste("Cluster", cluster_id),
    filename = file.path(outdir, paste0("flux_", cluster_id, ".pdf")),
    width = 3,
    height = 4.5
  )
}

plot_cluster_gene_heatmap <- function(expression_annotated, age_info_df, cluster_module, cluster_id, outdir) {
  rownames(age_info_df) <- age_info_df$Sample
  age_info_df$age <- as.numeric(as.character(age_info_df$age))
  
  annotation_col <- data.frame(Age = age_info_df$age)
  rownames(annotation_col) <- age_info_df$Sample
  
  sub_expr <- expression_annotated[expression_annotated$module %in% cluster_module, , drop = FALSE]
  sub_expr <- sub_expr[!duplicated(sub_expr$gene_symbol), , drop = FALSE]
  if (nrow(sub_expr) == 0) return(NULL)
  
  sub_expr_matrix <- sub_expr[, !(colnames(sub_expr) %in% c("gene_symbol", "metabolite", "module")), drop = FALSE]
  rownames(sub_expr_matrix) <- sub_expr$gene_symbol
  
  expr_with_age <- t(sub_expr_matrix)
  expr_with_age <- cbind(expr_with_age, Age = annotation_col[rownames(expr_with_age), "Age"])
  
  expr_with_age_df <- as.data.frame(expr_with_age)
  expr_with_age_df$Age <- as.numeric(as.character(expr_with_age_df$Age))
  
  expr_long <- reshape2::melt(expr_with_age_df, id.vars = "Age", variable.name = "Gene", value.name = "Expression")
  expr_avg <- aggregate(Expression ~ Gene + Age, data = expr_long, FUN = mean)
  expr_wide <- reshape2::dcast(expr_avg, Gene ~ Age, value.var = "Expression")
  
  rownames(expr_wide) <- expr_wide$Gene
  expr_wide$Gene <- NULL
  sub_expr_matrix <- as.matrix(expr_wide)
  
  sorted_sample_names <- colnames(sub_expr_matrix)
  sub_annotation_col <- data.frame(Age = as.numeric(sorted_sample_names))
  rownames(sub_annotation_col) <- sorted_sample_names
  
  age_order <- order(sub_annotation_col$Age)
  sorted_sample_names <- rownames(sub_annotation_col)[age_order]
  
  sub_expr_matrix <- sub_expr_matrix[, sorted_sample_names, drop = FALSE]
  sub_expr_matrix <- sub_expr_matrix[order(rownames(sub_expr_matrix)), , drop = FALSE]
  sub_annotation_col <- sub_annotation_col[sorted_sample_names, , drop = FALSE]
  
  sub_annotation_col$Age <- factor(
    as.character(sub_annotation_col$Age),
    levels = c("3", "12", "15", "18", "21", "26", "28"),
    labels = c("3M", "12M", "15M", "18M", "21M", "26M", "28M")
  )
  
  sub_annotation_col_colors <- list(
    Age = c(
      "3M" = "#e41a1c",
      "12M" = "#377eb8",
      "15M" = "#4daf4a",
      "18M" = "#984ea3",
      "21M" = "#ff7f00",
      "26M" = "#ffff33",
      "28M" = "#a65628"
    )
  )
  
  sub_expr_matrix <- na.omit(sub_expr_matrix)
  if (nrow(sub_expr_matrix) == 0) return(NULL)
  
  heatmap_height <- max(3, min(0.1 * nrow(sub_expr_matrix), 25))
  
  pheatmap(
    sub_expr_matrix,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    scale = "row",
    annotation_col = sub_annotation_col,
    annotation_colors = sub_annotation_col_colors,
    fontsize_row = 8,
    fontsize = 8,
    main = paste0(cluster_id),
    filename = file.path(outdir, paste0("heatmap_gene_", cluster_id, ".pdf")),
    color = rev(brewer.pal(11, "RdBu")),
    width = 3,
    height = heatmap_height
  )
  
  write.csv(sub_expr_matrix, file.path(outdir, paste0("heatmap_gene_", cluster_id, ".csv")))
}

plot_cluster_metabolite_heatmap <- function(met_expr_df, gene_module_metabolite, cluster_module, cluster_id, outdir) {
  met_expr_mean <- met_expr_df %>%
    group_by(rNs, age) %>%
    summarise(mean_intensity = mean(mean_intensity, na.rm = TRUE), .groups = "drop")
  
  gene_module_metabolite_cluster <- gene_module_metabolite[gene_module_metabolite$module %in% cluster_module, , drop = FALSE]
  sub_met_expr_mean <- met_expr_mean[met_expr_mean$rNs %in% gene_module_metabolite_cluster$metabolite, , drop = FALSE]
  
  if (nrow(sub_met_expr_mean) == 0) return(NULL)
  
  heatmap_matrix <- sub_met_expr_mean %>%
    pivot_wider(names_from = age, values_from = mean_intensity)
  
  heatmap_matrix <- as.data.frame(heatmap_matrix)
  rownames(heatmap_matrix) <- heatmap_matrix$rNs
  heatmap_matrix$rNs <- NULL
  heatmap_matrix <- as.matrix(heatmap_matrix)
  
  desired_order <- intersect(c("6M", "15M", "21M"), colnames(heatmap_matrix))
  if (length(desired_order) > 0) {
    heatmap_matrix <- heatmap_matrix[, desired_order, drop = FALSE]
  }
  
  scaled_matrix <- t(apply(heatmap_matrix, 1, function(x) {
    if (sum(!is.na(x)) >= 2) {
      as.numeric(scale(x))
    } else {
      rep(NA, length(x))
    }
  }))
  
  rownames(scaled_matrix) <- rownames(heatmap_matrix)
  colnames(scaled_matrix) <- colnames(heatmap_matrix)
  scaled_matrix <- scaled_matrix[rowSums(!is.na(scaled_matrix)) > 0, , drop = FALSE]
  
  if (nrow(scaled_matrix) == 0) return(NULL)
  
  met_annotation_col <- data.frame(age = colnames(scaled_matrix))
  rownames(met_annotation_col) <- colnames(scaled_matrix)
  
  annotation_colors <- list(
    age = setNames(
      colorRampPalette(c("skyblue", "orange", "green"))(ncol(scaled_matrix)),
      colnames(scaled_matrix)
    )
  )
  
  met_heatmap_height <- max(1, min(0.08 * nrow(scaled_matrix), 25))
  
  pheatmap(
    scaled_matrix,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    annotation_col = met_annotation_col,
    annotation_colors = annotation_colors,
    fontsize_row = 8,
    fontsize = 8,
    main = paste0(cluster_id),
    filename = file.path(outdir, paste0("heatmap_", cluster_id, "_metabolites.pdf")),
    color = rev(brewer.pal(11, "RdBu")),
    width = 3,
    height = met_heatmap_height
  )
  
  write.csv(heatmap_matrix, file.path(outdir, paste0("heatmap_", cluster_id, "_metabolites.csv")))
}

plot_venn_and_export_groups <- function(gene_csv_files, labels, out_prefix, outdir) {
  genes_list <- list()
  
  for (i in seq_along(gene_csv_files)) {
    f <- gene_csv_files[i]
    if (file.exists(f)) {
      tmp <- read.csv(f, row.names = 1, check.names = FALSE)
      genes_list[[labels[i]]] <- rownames(tmp)
    }
  }
  
  if (length(genes_list) < 2) return(NULL)
  
  all_genes <- unique(unlist(genes_list))
  presence_matrix <- sapply(genes_list, function(glist) all_genes %in% glist)
  rownames(presence_matrix) <- all_genes
  
  presence_df <- as.data.frame(presence_matrix)
  presence_df$group <- apply(presence_df, 1, function(x) {
    present_clusters <- names(which(x))
    if (length(present_clusters) == 0) return("None")
    paste(present_clusters, collapse = "_")
  })
  
  split_genes <- split(rownames(presence_df), presence_df$group)
  
  gene_group_dir <- file.path(outdir, paste0(out_prefix, "_gene_groups"))
  if (!dir.exists(gene_group_dir)) dir.create(gene_group_dir, recursive = TRUE)
  
  for (group_name in names(split_genes)) {
    gene_vec <- split_genes[[group_name]]
    if (length(gene_vec) > 0) {
      write.csv(
        data.frame(Gene = gene_vec),
        file = file.path(gene_group_dir, paste0("genes_", gsub("[^A-Za-z0-9_]", "_", group_name), ".csv")),
        row.names = FALSE
      )
    }
  }
  
  pdf(file.path(outdir, paste0(out_prefix, ".pdf")), width = 8, height = 8)
  print(
    ggVennDiagram(genes_list, label = "count", label_alpha = 0) +
      scale_fill_gradient(low = "#E0F7FA", high = "#56B4E9") +
      theme_void(base_size = 18) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      ggtitle("Gene Overlap")
  )
  dev.off()
}

parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(v) {
    as.numeric(v[1]) / as.numeric(v[2])
  })
}

convert_gene <- function(gene_symbols) {
  bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
}

plot_go_enrichment <- function(gene_symbols, pdf_file, title_prefix, top_n = 10) {
  genes_df <- convert_gene(gene_symbols)
  if (is.null(genes_df) || nrow(genes_df) == 0) {
    message("No valid genes mapped for: ", title_prefix)
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = genes_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  if (is.null(ego_df) || nrow(ego_df) == 0) {
    message("No significant GO terms for: ", title_prefix)
    return(NULL)
  }
  
  ego_df$GeneRatioNum <- parse_ratio(ego_df$GeneRatio)
  ego_df <- ego_df[order(ego_df$p.adjust, ego_df$Count, decreasing = c(FALSE, TRUE)), , drop = FALSE]
  ego_df <- head(ego_df, top_n)
  ego_df$Description <- factor(ego_df$Description, levels = rev(ego_df$Description))
  
  p1 <- ggplot(ego_df, aes(x = Count, y = Description)) +
    geom_col(fill = "#4C97C2") +
    theme_bw() +
    labs(
      title = paste0("GO BP Enrichment - ", title_prefix),
      x = "Gene Count",
      y = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  p2 <- ggplot(ego_df, aes(x = GeneRatioNum, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    theme_bw() +
    labs(
      title = paste0("GO BP Enrichment (Dotplot) - ", title_prefix),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  pdf(pdf_file, width = 8, height = 6)
  print(p1)
  print(p2)
  dev.off()
  
  invisible(ego_df)
}

plot_go_dotplot_3clusters <- function(cluster_gene_list, pdf_file, title_prefix = "Clusters 1-3", top_n = 10) {
  cluster_entrez_list <- lapply(cluster_gene_list, function(gene_symbols) {
    genes_df <- convert_gene(gene_symbols)
    if (is.null(genes_df) || nrow(genes_df) == 0) return(NULL)
    unique(genes_df$ENTREZID)
  })
  
  cluster_entrez_list <- cluster_entrez_list[!sapply(cluster_entrez_list, is.null)]
  
  if (length(cluster_entrez_list) == 0) {
    message("No valid gene sets available for combined GO enrichment.")
    return(NULL)
  }
  
  cc <- compareCluster(
    geneCluster = cluster_entrez_list,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  cc_df <- as.data.frame(cc)
  if (is.null(cc_df) || nrow(cc_df) == 0) {
    message("No significant GO terms found for combined clusters.")
    return(NULL)
  }
  
  cc_df <- cc_df %>%
    group_by(Cluster) %>%
    arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  cc_df$Description <- factor(cc_df$Description, levels = rev(unique(cc_df$Description)))
  
  p <- ggplot(cc_df, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
    scale_size(range = c(2.5, 6)) +
    theme_bw() +
    labs(
      title = paste0("GO BP Enrichment Dotplot - ", title_prefix),
      x = NULL,
      y = NULL,
      color = "Adjusted p-value",
      size = "Gene Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  ggsave(
    filename = pdf_file,
    plot = p,
    width = 6,
    height = 4.5
  )
  
  write.csv(
    cc_df,
    file.path(dirname(pdf_file), "enrichment_Clusters123_GO_dotplot_table.csv"),
    row.names = FALSE
  )
  
  invisible(cc_df)
}

###############################################################################
# 2. Input data
###############################################################################
module_gene_sets_df <- read.csv("module_gene_complete_mouse_m168.csv", check.names = FALSE)
related_gene_sets_all_df <- read_excel("Gene_Names_Results.xlsx", sheet = 1)
met_expr_df <- read_excel("sum_normalized_data.xlsx")
met_expr_df <- as.data.frame(met_expr_df)
data_matrix <- read.csv("normalized_counts_with_metadata.csv")
df_flux <- read.csv("./region_flux_csv/all_flux_data_long.csv")
gene_connectivity_labeled <- read.csv("gene_connectivity_labeled.csv", row.names = 1, check.names = FALSE)

###############################################################################
# 3. Build gene-module-metabolite relationships
###############################################################################
module_gene_sets_list <- lapply(module_gene_sets_df, standardize_gene_vector)

real_met <- unique(met_expr_df$rNs)
related_gene_sets_all_filtered_df <- related_gene_sets_all_df %>%
  dplyr::select(all_of(real_met))

related_gene_sets_all_list <- lapply(related_gene_sets_all_filtered_df, standardize_gene_vector)

gene_metabolite_df <- bind_rows(lapply(names(related_gene_sets_all_list), function(metabolite) {
  data.frame(
    gene_symbol = related_gene_sets_all_list[[metabolite]],
    metabolite = metabolite,
    stringsAsFactors = FALSE
  )
}))

result_df <- make_gene_metabolite_module(gene_metabolite_df, module_gene_sets_list)
write.csv(result_df, file.path(dir_tables, "gene_module_metabolite.csv"), row.names = FALSE)

###############################################################################
# 4. Gene connectivity
###############################################################################
gene_connectivity <- result_df %>%
  group_by(gene_symbol) %>%
  summarise(
    metabolite_degree = n_distinct(metabolite),
    module_degree = n_distinct(module),
    total_degree = metabolite_degree + module_degree,
    .groups = "drop"
  ) %>%
  arrange(desc(total_degree))

write.csv(gene_connectivity, file.path(dir_tables, "gene_connectivity.csv"), row.names = FALSE)

###############################################################################
# 5. Whole-brain clustering and clustered heatmap only
###############################################################################
whole_brain_obj <- make_whole_brain_flux_matrix(df_flux, result_df$module)
mat <- whole_brain_obj$mat
annotation_col <- whole_brain_obj$annotation_col

module_cluster_df <- plot_whole_brain_clustered_heatmap(mat, annotation_col, dir_whole_brain)
write.csv(module_cluster_df, file.path(dir_tables, "module_cluster_df.csv"), row.names = FALSE)

###############################################################################
# 6. Sankey data and plots
###############################################################################
df_alluvial <- result_df %>%
  dplyr::select(metabolite, gene_symbol, module) %>%
  na.omit()

df_alluvial_labeled <- df_alluvial %>%
  left_join(gene_connectivity_labeled %>% dplyr::select(gene_symbol, gene), by = "gene_symbol") %>%
  mutate(gene_label = ifelse(is.na(gene), as.character(gene_symbol), gene)) %>%
  left_join(module_cluster_df, by = "module") %>%
  filter(!is.na(cluster))

write.csv(df_alluvial_labeled, file.path(dir_tables, "df_alluvial_labeled.csv"), row.names = FALSE)

df_alluvial_split <- split(df_alluvial_labeled, df_alluvial_labeled$cluster)
for (cluster_id in names(df_alluvial_split)) {
  plot_sankey_by_cluster(df_alluvial_split[[cluster_id]], cluster_id, dir_whole_brain)
}

###############################################################################
# 7. Cluster trend expression matrix
###############################################################################
expr_obj <- make_expression_annotated(data_matrix, result_df)
expression_annotated <- expr_obj$expression_annotated
age_info_df <- expr_obj$age_info_df

write.csv(expression_annotated, file.path(dir_tables, "filtered_expression_annotated_whole_regions.csv"), row.names = FALSE)
write.csv(age_info_df, file.path(dir_tables, "sample_age_info_whole_regions.csv"), row.names = FALSE)

###############################################################################
# 8. Cluster-level age trend plots
###############################################################################
cluster_ids <- c("1", "2", "3")

for (cluster_id in cluster_ids) {
  cluster_module <- module_cluster_df[module_cluster_df$cluster == cluster_id, "module"]
  cluster_df <- df_flux[df_flux$module %in% cluster_module, , drop = FALSE]
  
  plot_cluster_flux_heatmap(cluster_df, cluster_id, dir_cluster_trends)
  plot_cluster_gene_heatmap(expression_annotated, age_info_df, cluster_module, cluster_id, dir_cluster_trends)
  plot_cluster_metabolite_heatmap(met_expr_df, result_df, cluster_module, cluster_id, dir_cluster_trends)
}

###############################################################################
# 9. Cluster trend venn
###############################################################################
plot_venn_and_export_groups(
  gene_csv_files = c(
    file.path(dir_cluster_trends, "heatmap_gene_1.csv"),
    file.path(dir_cluster_trends, "heatmap_gene_2.csv"),
    file.path(dir_cluster_trends, "heatmap_gene_3.csv")
  ),
  labels = c("Cluster_1", "Cluster_2", "Cluster_3"),
  out_prefix = "venn_clusters_1_to_3_cluster_trends",
  outdir = dir_venn
)

###############################################################################
# 10. Whole-brain GO enrichment only
###############################################################################
genes_cluster1 <- c("Adk", "Adpgk", "Ahcy", "Ahcyl1", "Ahcyl2", "Ak1", "Ak2", "Ak3", "Ak4", "Ak5", "Ak7", "Ak8",
                    "Asns", "Atic", "B4galt1", "Dnmt1", "Dnmt3a", "Dnmt3b", "Gne", "Nudt2", "Nudt5")

genes_cluster2 <- c("Ada", "Adsl", "Gmps", "Rnf2", "Slc27a5", "Uck1", "Uck2", "Uckl1", "Upp1", "Upp2", "Uprt",
                    "Aacs", "Ampd1", "Ampd2", "Ampd3", "Carns1", "Pygb", "Xdh")

genes_cluster3 <- c("Prps1", "Prps1l1", "Prps2", "Acsbg1", "Acsbg2", "Acsl1", "Acsl3", "Acsl4", "Acsl5", "Acsl6",
                    "Aprt", "Cant1", "Cmpk1", "Cmpk2", "Dck", "Dctpp1", "Entpd4", "Entpd5", "Entpd6")

plot_go_enrichment(genes_cluster1, file.path(dir_go, "enrichment_Cluster1_GO.pdf"), "Cluster1")
plot_go_enrichment(genes_cluster2, file.path(dir_go, "enrichment_Cluster2_GO.pdf"), "Cluster2")
plot_go_enrichment(genes_cluster3, file.path(dir_go, "enrichment_Cluster3_GO.pdf"), "Cluster3")

###############################################################################
# 11. Standalone GO dotplot across 3 clusters
###############################################################################
cluster_gene_list <- list(
  Cluster1 = genes_cluster1,
  Cluster2 = genes_cluster2,
  Cluster3 = genes_cluster3
)

plot_go_dotplot_3clusters(
  cluster_gene_list = cluster_gene_list,
  pdf_file = file.path(dir_go, "enrichment_Clusters123_GO_dotplot.pdf"),
  title_prefix = "Clusters 1-3",
  top_n = 10
)

