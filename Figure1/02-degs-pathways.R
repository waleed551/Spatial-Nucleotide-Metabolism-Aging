setwd("D:/Koziol_lab/Spatial-Nucleotide-Metabolism-Aging/Figure1")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(igraph)
library(ggraph)

start_time <- Sys.time()

cat("========== START HMDB PLOTTING ==========\n")
cat("[1/6] Working directory:\n")
cat(getwd(), "\n\n")

##==============================
## 1. 读取 HMDB 富集结果
##==============================
cat("[2/6] Reading HMDB enrichment table...\n")

input_file <- "HMDB_Metabolites_table.txt"

if (!file.exists(input_file)) {
  stop("File not found: HMDB_Metabolites_table.txt")
}

hmdb_df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

cat("    Rows:", nrow(hmdb_df), "\n")
cat("    Columns:\n")
print(colnames(hmdb_df))
cat("\n")

##==============================
## 2. 清理字段
##==============================
cat("[3/6] Cleaning columns...\n")

colnames(hmdb_df) <- c(
  "Term", "Overlap", "P_value", "Adjusted_P_value",
  "Old_P_value", "Old_Adjusted_P_value",
  "Odds_Ratio", "Combined_Score", "Genes"
)

safe_text <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x[is.na(x)] <- ""
  x
}

extract_overlap_numbers <- function(x) {
  x <- safe_text(x)
  x <- gsub("[^0-9]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x
}

get_overlap_num <- function(x, pos = 1) {
  x <- extract_overlap_numbers(x)
  sapply(strsplit(x, " "), function(v) {
    if (length(v) < pos) return(NA_real_)
    suppressWarnings(as.numeric(v[pos]))
  })
}

hmdb_df <- hmdb_df %>%
  dplyr::mutate(
    Term = safe_text(Term),
    Overlap_raw = safe_text(Overlap),
    Genes = safe_text(Genes),
    overlap_num = get_overlap_num(Overlap_raw, 1),
    overlap_den = get_overlap_num(Overlap_raw, 2),
    Overlap_clean = paste0(overlap_num, "/", overlap_den),
    GeneRatio = overlap_num / overlap_den,
    P_value = suppressWarnings(as.numeric(P_value)),
    Adjusted_P_value = suppressWarnings(as.numeric(Adjusted_P_value)),
    Odds_Ratio = suppressWarnings(as.numeric(Odds_Ratio)),
    Combined_Score = suppressWarnings(as.numeric(Combined_Score)),
    neglog10_FDR = -log10(Adjusted_P_value),
    Term_clean = stringr::str_replace(Term, " \\(HMDB[0-9]+\\)$", "")
  ) %>%
  dplyr::filter(
    !is.na(overlap_num),
    !is.na(overlap_den),
    overlap_den > 0,
    !is.na(Adjusted_P_value)
  )

cat("    Rows after cleaning:", nrow(hmdb_df), "\n\n")

##==============================
## 3. 过滤冗余项
##==============================
cat("[4/6] Filtering redundant terms...\n")

hmdb_plot_df <- hmdb_df %>%
  dplyr::filter(
    !stringr::str_detect(Term_clean, "^PI\\("),
    !stringr::str_detect(Term_clean, "^C[0-9]+H[0-9]+O[0-9]+P[0-9]+"),
    !stringr::str_detect(Term_clean, "1,2-dihexadecanoyl"),
    !stringr::str_detect(Term_clean, "sn-glycero-3-phospho"),
    !stringr::str_detect(Term_clean, "^C11H22O22P4$"),
    !stringr::str_detect(Term_clean, "^C41H"),
    !stringr::str_detect(Term_clean, "^C43H"),
    !stringr::str_detect(Term_clean, "^C45H"),
    !stringr::str_detect(Term_clean, "^C47H"),
    !stringr::str_detect(Term_clean, "^C48H"),
    !stringr::str_detect(Term_clean, "^C50H")
  ) %>%
  dplyr::arrange(Adjusted_P_value, dplyr::desc(overlap_num))

sig_df <- hmdb_plot_df %>%
  dplyr::filter(Adjusted_P_value < 0.05)

if (nrow(sig_df) >= 8) {
  top_df <- sig_df %>% dplyr::slice_head(n = 15)
} else {
  top_df <- hmdb_plot_df %>% dplyr::slice_head(n = 15)
}

top_df <- top_df %>%
  dplyr::mutate(
    Term_clean = ifelse(is.na(Term_clean) | Term_clean == "", Term, Term_clean),
    Term_plot = stringr::str_wrap(Term_clean, width = 24),
    Term_plot = forcats::fct_reorder(Term_plot, GeneRatio)
  )

write.csv(hmdb_df, "HMDB_Metabolites_table_cleaned.csv", row.names = FALSE)
write.csv(top_df, "HMDB_Metabolites_top_terms_for_plot.csv", row.names = FALSE)

cat("    Terms kept for plotting:", nrow(top_df), "\n")
cat("    Output tables written:\n")
cat("    - HMDB_Metabolites_table_cleaned.csv\n")
cat("    - HMDB_Metabolites_top_terms_for_plot.csv\n\n")

##==============================
## 4. 绘制小尺寸气泡图 + 大字体
##==============================
cat("[5/6] Drawing HMDB bubble plot...\n")

p_bubble <- ggplot(
  top_df,
  aes(
    x = GeneRatio,
    y = Term_plot,
    size = overlap_num,
    color = neglog10_FDR
  )
) +
  geom_point(alpha = 0.95) +
  scale_size(range = c(4, 10)) +
  scale_color_gradient(
    low = "#9ECAE1",
    high = "#D73027",
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    title = "HMDB metabolite enrichment",
    x = "Gene ratio",
    y = NULL,
    size = "Hits"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#EAEAEA")
  )

pdf("HMDB_metabolite_bubbleplot_small.pdf", width = 8, height = 6.5)
print(p_bubble)
dev.off()

png("HMDB_metabolite_bubbleplot_small.png", width = 1600, height = 1300, res = 220)
print(p_bubble)
dev.off()

##==============================
## 5. 绘制小尺寸 network 图 + 大字体
##==============================
cat("[6/6] Drawing metabolite-gene network...\n")

network_terms <- top_df %>%
  dplyr::mutate(
    Term_clean = ifelse(is.na(Term_clean) | Term_clean == "", Term, Term_clean)
  ) %>%
  dplyr::slice_head(n = 8)

edge_df <- network_terms %>%
  dplyr::select(Term_clean, Genes, Adjusted_P_value) %>%
  tidyr::separate_rows(Genes, sep = ";") %>%
  dplyr::rename(
    metabolite = Term_clean,
    gene = Genes
  ) %>%
  dplyr::mutate(
    gene = stringr::str_trim(gene)
  ) %>%
  dplyr::filter(gene != "")

node_metabolite <- edge_df %>%
  dplyr::distinct(metabolite, Adjusted_P_value) %>%
  dplyr::transmute(
    name = metabolite,
    type = "Metabolite",
    score = -log10(Adjusted_P_value)
  )

node_gene <- edge_df %>%
  dplyr::distinct(gene) %>%
  dplyr::transmute(
    name = gene,
    type = "Gene",
    score = 2
  )

node_df <- dplyr::bind_rows(node_metabolite, node_gene)

g <- igraph::graph_from_data_frame(
  d = edge_df %>% dplyr::select(metabolite, gene),
  vertices = node_df,
  directed = FALSE
)

p_net <- ggraph::ggraph(g, layout = "fr") +
  ggraph::geom_edge_link(alpha = 0.25, colour = "grey65") +
  ggraph::geom_node_point(
    aes(size = score, color = type)
  ) +
  ggraph::geom_node_text(
    aes(
      label = name,
      fontface = ifelse(type == "Metabolite", "bold", "plain")
    ),
    size = 4.5,
    repel = TRUE
  ) +
  scale_color_manual(
    values = c("Metabolite" = "#D73027", "Gene" = "#2C7BB6")
  ) +
  scale_size_continuous(range = c(4, 10), guide = "none") +
  theme_void(base_size = 16) +
  labs(
    title = "Metabolite-gene network"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

pdf("HMDB_metabolite_gene_network_small.pdf", width = 8.5, height = 7)
print(p_net)
dev.off()

png("HMDB_metabolite_gene_network_small.png", width = 1700, height = 1400, res = 220)
print(p_net)
dev.off()

cat("\n========== DONE ==========\n")
cat("Output files generated:\n")
cat(" - HMDB_Metabolites_table_cleaned.csv\n")
cat(" - HMDB_Metabolites_top_terms_for_plot.csv\n")
cat(" - HMDB_metabolite_bubbleplot_small.pdf\n")
cat(" - HMDB_metabolite_bubbleplot_small.png\n")
cat(" - HMDB_metabolite_gene_network_small.pdf\n")
cat(" - HMDB_metabolite_gene_network_small.png\n")
cat("Total runtime:", Sys.time() - start_time, "\n")

