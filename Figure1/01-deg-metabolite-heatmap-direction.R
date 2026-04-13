setwd("D:/Koziol_lab/Chenyang_paper/My Updates/Code/Figure1")
load("dds_BulkSeq_Aging.bin")

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

start_time <- Sys.time()

cat("========== START ==========\n")
cat("[1/9] Working directory:\n")
cat(getwd(), "\n\n")

cat("[2/9] Objects loaded from bin:\n")
print(ls())
cat("\n")

##==============================
## 1. 读取 metabolite-gene 关联表
##==============================
cat("[3/9] Reading metabolite-gene table: Gene_Names_Results.xlsx\n")

gene_sets_df <- read_excel("Gene_Names_Results.xlsx", sheet = 1)

cat("    Table dimension:", nrow(gene_sets_df), "rows x", ncol(gene_sets_df), "columns\n\n")

format_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

gene_sets_list <- apply(gene_sets_df, 2, function(col) unique(format_gene(col)))
gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) > 0]

cat("    Number of metabolites:", length(gene_sets_list), "\n")
cat("    Metabolite names:\n")
print(names(gene_sets_list))
cat("\n")

cat("    Gene counts per metabolite:\n")
print(sapply(gene_sets_list, length))
cat("\n")

##==============================
## 2. comparison 统一方向
##==============================
normalize_comparison <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "_vs_")[[1]]
  
  if (length(parts) != 2) return(x)
  
  a <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", parts[1])))
  b <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", parts[2])))
  
  if (is.na(a) || is.na(b)) return(x)
  
  paste0(max(a, b), "_vs_", min(a, b))
}

##==============================
## 3. 检查 bin 结构
##==============================
cat("[4/9] Checking results_list_CA1 structure...\n")

if (!exists("results_list_CA1")) {
  stop("Object 'results_list_CA1' not found in dds_BulkSeq_Aging.bin")
}

if (!"Both" %in% names(results_list_CA1)) {
  stop("results_list_CA1 does not contain 'Both'")
}

cat("    Regions found:\n")
print(names(results_list_CA1$Both))
cat("\n")

##==============================
## 4. 提取显著基因的 log2FC
##==============================
cat("[5/9] Extracting significant log2FoldChange from DE results...\n")

plot_list <- list()
block_count <- 0

all_regions <- names(results_list_CA1$Both)
total_regions <- length(all_regions)

for (i in seq_along(all_regions)) {
  
  region <- all_regions[i]
  region_obj <- results_list_CA1$Both[[region]]
  
  cat("--------------------------------------------------\n")
  cat("Region", i, "/", total_regions, ":", region, "\n")
  
  all_comparisons <- names(region_obj)
  total_comparisons <- length(all_comparisons)
  
  cat("Comparisons in this region:", total_comparisons, "\n")
  print(all_comparisons)
  cat("\n")
  
  for (j in seq_along(all_comparisons)) {
    
    comparison_raw <- all_comparisons[j]
    comparison <- normalize_comparison(comparison_raw)
    res_obj <- region_obj[[comparison_raw]]
    
    cat("  -> Comparison", j, "/", total_comparisons, ":", comparison_raw,
        " -> normalized as ", comparison, "\n", sep = "")
    
    if (!"resOrdered" %in% names(res_obj)) {
      cat("     [Skip] no resOrdered found\n\n")
      next
    }
    
    res_df <- as.data.frame(res_obj$resOrdered)
    if (nrow(res_df) == 0) {
      cat("     [Skip] resOrdered is empty\n\n")
      next
    }
    
    if (!"gene_symbol" %in% colnames(res_df)) {
      cat("     [Skip] gene_symbol column not found\n\n")
      next
    }
    
    res_df$gene <- res_df$gene_symbol
    res_df$gene <- format_gene(res_df$gene)
    
    lfc_tbl <- res_df %>%
      dplyr::select(gene, log2FoldChange, padj, pvalue) %>%
      filter(
        !is.na(gene), gene != "",
        !is.na(log2FoldChange),
        !is.na(padj), padj < 0.05
      ) %>%
      group_by(gene) %>%
      slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
      ungroup()
    
    cat("     Significant genes (padj < 0.05):", nrow(lfc_tbl), "\n")
    
    if (nrow(lfc_tbl) == 0) {
      cat("     [Skip] no significant genes after padj filtering\n\n")
      next
    }
    
    hit_metabolite_count <- 0
    
    for (metabolite in names(gene_sets_list)) {
      
      gene_set <- gene_sets_list[[metabolite]]
      hit_genes <- intersect(gene_set, lfc_tbl$gene)
      
      if (length(hit_genes) == 0) next
      
      hit_metabolite_count <- hit_metabolite_count + 1
      
      cat("       *", metabolite, ":", length(hit_genes), "matched significant genes\n")
      
      tmp <- lfc_tbl %>%
        filter(gene %in% hit_genes) %>%
        mutate(
          metabolite = metabolite,
          comparison = comparison,
          region = region,
          row_id = paste(region, comparison, metabolite, sep = " | ")
        )
      
      plot_list[[length(plot_list) + 1]] <- tmp
      block_count <- block_count + 1
    }
    
    cat("     Matched metabolites in this comparison:", hit_metabolite_count, "\n")
    cat("     Accumulated blocks so far:", block_count, "\n\n")
  }
}

if (length(plot_list) == 0) {
  stop("No matched significant genes found between Gene_Names_Results.xlsx and dds_BulkSeq_Aging.bin")
}

plot_df <- bind_rows(plot_list)

##==============================
## 5. 去掉反向 comparison 重复
##==============================
cat("[6/9] Removing reversed comparison duplicates...\n")

plot_df <- plot_df %>%
  mutate(
    comparison = sapply(comparison, normalize_comparison),
    row_id = paste(region, comparison, metabolite, sep = " | ")
  ) %>%
  distinct(row_id, gene, .keep_all = TRUE)

cat("    Combined long table rows:", nrow(plot_df), "\n")
cat("    Unique row_id:", length(unique(plot_df$row_id)), "\n")
cat("    Unique genes:", length(unique(plot_df$gene)), "\n\n")

##==============================
## 6. 生成 heatmap 矩阵
##==============================
cat("[7/9] Building heatmap matrix...\n")

heat_df <- plot_df %>%
  dplyr::select(row_id, gene, log2FoldChange) %>%
  distinct() %>%
  pivot_wider(
    names_from = gene,
    values_from = log2FoldChange,
    values_fill = 0
  )

mat <- as.data.frame(heat_df)
rownames(mat) <- mat$row_id
mat$row_id <- NULL
mat <- as.matrix(mat)

cat("    Matrix dimension before filtering:", nrow(mat), "x", ncol(mat), "\n")

threshold_lfc <- 1
mat[abs(mat) < threshold_lfc] <- 0
cat("    Applied abs(log2FC) <", threshold_lfc, " -> set to 0\n")

mat <- mat[
  rowSums(mat != 0, na.rm = TRUE) > 0,
  colSums(mat != 0, na.rm = TRUE) > 0,
  drop = FALSE
]
cat("    Matrix dimension after removing all-zero rows/cols:", nrow(mat), "x", ncol(mat), "\n")

min_nonzero_per_row <- 2
min_nonzero_per_col <- 2

mat <- mat[
  rowSums(mat != 0, na.rm = TRUE) >= min_nonzero_per_row,
  colSums(mat != 0, na.rm = TRUE) >= min_nonzero_per_col,
  drop = FALSE
]

## 过滤完成后再保留 1 位小数，避免 0.96 -> 1.0 这种阈值误差
mat <- round(mat, 1)
plot_df$log2FoldChange <- round(plot_df$log2FoldChange, 1)

cat("    Matrix dimension after stricter filtering:", nrow(mat), "x", ncol(mat), "\n")
cat("    Rule: keep rows with >=", min_nonzero_per_row,
    "non-zero values; keep cols with >=", min_nonzero_per_col, "non-zero values\n\n")

if (nrow(mat) == 0 || ncol(mat) == 0) {
  stop("Matrix is empty after filtering. Try reducing threshold_lfc or min_nonzero_per_row/col.")
}

##==============================
## 7. 行注释 + row split
##==============================
cat("[8/9] Building row annotations...\n")

annotation_row <- plot_df %>%
  dplyr::select(row_id, metabolite, comparison, region) %>%
  distinct() %>%
  as.data.frame()

rownames(annotation_row) <- annotation_row$row_id
annotation_row$row_id <- NULL
annotation_row <- annotation_row[rownames(mat), , drop = FALSE]

wanted_metabolite_order <- c(
  "AdG",
  "AMPdGMP",
  "dAMP",
  "GMP",
  "GTP",
  "ITP",
  "UMPpseudoUMP",
  "UpseudoU"
)

present_metabolites <- wanted_metabolite_order[wanted_metabolite_order %in% annotation_row$metabolite]

annotation_row <- annotation_row %>%
  mutate(
    region = as.character(region),
    comparison = as.character(comparison),
    metabolite = factor(as.character(metabolite), levels = present_metabolites)
  ) %>%
  arrange(metabolite, region, comparison)

mat <- mat[rownames(annotation_row), , drop = FALSE]

cat("    Annotation rows:", nrow(annotation_row), "\n")
cat("    Unique metabolites kept:", length(unique(annotation_row$metabolite)), "\n")
cat("    Unique comparisons kept:", length(unique(annotation_row$comparison)), "\n")
cat("    Unique regions kept:", length(unique(annotation_row$region)), "\n\n")

make_named_colors <- function(x, palette_name = "Set3") {
  ux <- unique(as.character(x))
  n <- length(ux)
  if (n <= 12) {
    cols <- brewer.pal(max(3, n), palette_name)[1:n]
  } else {
    cols <- colorRampPalette(brewer.pal(12, palette_name))(n)
  }
  names(cols) <- ux
  cols
}

annotation_colors <- list(
  metabolite = make_named_colors(annotation_row$metabolite, "Dark2"),
  comparison = make_named_colors(annotation_row$comparison, "Pastel1"),
  region = make_named_colors(annotation_row$region, "Set3")
)

##==============================
## legend label 映射
##==============================
region_label_map <- c(
  cc  = "Corpus callosum",
  cer = "Cerebellum",
  cor = "Cortex",
  cp  = "Caudoputamen",
  ent = "Entorhinal cortex",
  hi  = "Hippocampus",
  hy  = "Hypothalamus",
  med = "Medulla",
  olf = "Olfactory area",
  plx = "Plexiform layer",
  pon = "Pons",
  svz = "Subventricular zone",
  th  = "Thalamus"
)

metabolite_label_map <- c(
  AdG          = "AdG",
  AMPdGMP      = "AMPdGMP",
  dAMP         = "dAMP",
  GMP          = "GMP",
  GTP          = "GTP",
  ITP          = "ITP",
  UMPpseudoUMP = "UMPpseudoUMP",
  UpseudoU     = "UpseudoU"
)

##==============================
## 8. 按基因方向分 panel
##==============================
cat("[9/9] Defining gene panels and drawing heatmap...\n")

classify_gene_panel <- function(x) {
  x_nz <- x[x != 0]
  if (length(x_nz) == 0) return(NA_character_)
  if (all(x_nz >= 0) && any(x_nz > 0)) return("Up")
  if (all(x_nz <= 0) && any(x_nz < 0)) return("Down")
  return("Mixed")
}

gene_panel <- apply(mat, 2, classify_gene_panel)
gene_panel <- factor(gene_panel, levels = c("Up", "Down", "Mixed"))

keep_gene <- !is.na(gene_panel)
mat <- mat[, keep_gene, drop = FALSE]
gene_panel <- gene_panel[keep_gene]

cat("    Genes in Up panel   :", sum(gene_panel == "Up"), "\n")
cat("    Genes in Down panel :", sum(gene_panel == "Down"), "\n")
cat("    Genes in Mixed panel:", sum(gene_panel == "Mixed"), "\n\n")

##==============================
## 9. 图形美化参数
##==============================
n_row <- nrow(mat)
n_col <- ncol(mat)

cat("    Final matrix size:", n_row, "rows x", n_col, "columns\n")

max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- circlize::colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#2C7BB6", "#F7F7F7", "#D73027")
)

pdf_width  <- max(12, min(20, 8 + n_col * 0.22))
pdf_height <- max(8,  min(14, 4 + n_row * 0.25))

png_width  <- pdf_width * 180
png_height <- pdf_height * 180

cat("    PDF size:", pdf_width, "x", pdf_height, "\n")
cat("    PNG size:", png_width, "x", png_height, "\n\n")

## 左侧 annotation：保留 metabolite 色条，只显示注释名，不显示具体 metabolite 名字
row_ha <- rowAnnotation(
  metabolite = annotation_row$metabolite,
  region = annotation_row$region,
  comparison = annotation_row$comparison,
  col = list(
    metabolite = annotation_colors$metabolite,
    region = annotation_colors$region,
    comparison = annotation_colors$comparison
  ),
  show_annotation_name = c(metabolite = TRUE, region = TRUE, comparison = TRUE),
  show_legend = c(metabolite = FALSE, region = FALSE, comparison = FALSE),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  annotation_name_rot = 90,
  annotation_width = unit(c(0.35, 0.4, 0.5), "cm"),
  gap = unit(0.8, "mm"),
  gp = gpar(col = NA)
)

## 手工图例
region_lgd <- Legend(
  title = "region",
  at = names(annotation_colors$region),
  labels = region_label_map[names(annotation_colors$region)],
  legend_gp = gpar(fill = annotation_colors$region),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

comparison_lgd <- Legend(
  title = "comparison",
  at = names(annotation_colors$comparison),
  labels = names(annotation_colors$comparison),
  legend_gp = gpar(fill = annotation_colors$comparison),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

metabolite_lgd <- Legend(
  title = "metabolite",
  at = names(annotation_colors$metabolite),
  labels = metabolite_label_map[names(annotation_colors$metabolite)],
  legend_gp = gpar(fill = annotation_colors$metabolite),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

log2fc_lgd <- Legend(
  title = "Log2FC",
  col_fun = col_fun,
  at = c(-max_abs, 0, max_abs),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 10),
  legend_height = unit(3.8, "cm"),
  grid_width = unit(4, "mm")
)

ht <- Heatmap(
  mat,
  name = "Log2FC",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_method_columns = "complete",
  
  column_split = gene_panel,
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  row_split = annotation_row$metabolite,
  cluster_row_slices = FALSE,
  row_title = "",
  row_title_gp = gpar(fontsize = 1),
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = colnames(mat),
  column_names_rot = 90,
  column_names_centered = TRUE,
  column_names_gp = gpar(
    fontsize = 8,
    fontface = "bold.italic"
  ),
  column_names_max_height = unit(4.5, "cm"),
  
  left_annotation = row_ha,
  
  row_gap = unit(3, "mm"),
  column_gap = unit(3, "mm"),
  
  rect_gp = gpar(col = "#F2F2F2", lwd = 0.4),
  border = FALSE,
  
  use_raster = TRUE,
  raster_quality = 3,
  
  width = unit(min(max(ncol(mat) * 5.5, 90), 220), "mm"),
  height = unit(min(max(nrow(mat) * 6.5, 85), 220), "mm"),
  
  show_heatmap_legend = FALSE
)

pdf(
  "metabolite_gene_log2fc_heatmap_gene_metabolite_panels_final.pdf",
  width = pdf_width,
  height = pdf_height
)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE,
  annotation_legend_list = list(log2fc_lgd, region_lgd, comparison_lgd, metabolite_lgd),
  padding = unit(c(5, 14, 24, 5), "mm")
)
dev.off()

png(
  "metabolite_gene_log2fc_heatmap_gene_metabolite_panels_final.png",
  width = png_width,
  height = png_height,
  res = 180
)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE,
  annotation_legend_list = list(log2fc_lgd, region_lgd, comparison_lgd, metabolite_lgd),
  padding = unit(c(5, 14, 24, 5), "mm")
)
dev.off()

##==============================
## 10. 导出表格
##==============================
write.csv(plot_df, "metabolite_gene_log2fc_long_table_significant.csv", row.names = FALSE)
write.csv(
  as.data.frame(mat) %>% rownames_to_column("row_id"),
  "metabolite_gene_log2fc_matrix_gene_metabolite_panels.csv",
  row.names = FALSE
)
write.csv(
  data.frame(gene = colnames(mat), panel = as.character(gene_panel)),
  "metabolite_gene_gene_panels.csv",
  row.names = FALSE
)
write.csv(
  annotation_row %>% rownames_to_column("row_id"),
  "metabolite_gene_heatmap_annotation_significant.csv",
  row.names = FALSE
)

cat("\n========== DONE ==========\n")
cat("Output files generated:\n")
cat(" - metabolite_gene_log2fc_heatmap_gene_metabolite_panels_final.pdf\n")
cat(" - metabolite_gene_log2fc_heatmap_gene_metabolite_panels_final.png\n")
cat(" - metabolite_gene_log2fc_long_table_significant.csv\n")
cat(" - metabolite_gene_log2fc_matrix_gene_metabolite_panels.csv\n")
cat(" - metabolite_gene_gene_panels.csv\n")
cat(" - metabolite_gene_heatmap_annotation_significant.csv\n")
cat("Total runtime:", Sys.time() - start_time, "\n")
