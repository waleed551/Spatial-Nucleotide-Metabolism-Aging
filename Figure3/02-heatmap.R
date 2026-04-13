# 加载所需的库
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(writexl)
library(openxlsx)

setwd("D:/Koziol_lab/Spatial-Nucleotide-Metabolism-Aging/Figure3")

file_path <- "sum_normalized_data.xlsx"
data <- read_excel(file_path)

# 转换为宽格式
wide_data <- data %>%
  pivot_wider(names_from = rNs, values_from = mean_intensity)

write.csv(wide_data,"wide_data.csv")
write.xlsx(wide_data, file = "wide_data.xlsx")


# 将sheet列转换为行名
heatmap_data <- as.data.frame(wide_data)

# 数据转换为矩阵
heatmap_matrix <- as.matrix(heatmap_data[4:ncol(heatmap_data)])
rownames(heatmap_matrix)<-rownames(heatmap_data)

# 对数据取对数（自然对数或对数2，根据需要选择）
heatmap_matrix_log <- log10(heatmap_matrix + 1)  # +1 避免取对数0的问题
# 创建一个新的数据框，包含前两列的分组标签
annotations <- heatmap_data[, 1:2]
annotations$age <- factor(annotations$age, levels = c("6M", "15M", "21M"))

colnames(annotations) <- c("age", "region")  # 给注释列命名

# 2. 生成按 age 顺序排列的列索引
sorted_cols <- order(annotations$age)

# 3. 重排矩阵列顺序和注释
heatmap_matrix_log_sorted <- heatmap_matrix_log[sorted_cols, ]
annotations_sorted <- annotations[sorted_cols, ]


# 为 region 和 age 设置手动颜色
annotation_colors <- list(
  region = c(
    "olfactory bulbs" = "lightblue",
    "cortex" = "yellow",
    "corpus callosum" = "orange",
    
    "anterior olfactory" = "blue",
    "basal forebrain" = "green",
    "ventral striatum" = "darkblue",
    "caudate putamen" = "red",
  
    "fornix" = "cyan",
    "hippocampus" = "magenta",
    "Thalamus" = "darkred",
    "Hypothalamus" = "brown",
    
    "midbrain" = "pink",
    "cerebellum" = "purple",
    "pons&medulla" = "darkgreen"
    
    
  ),
  age = c(
    "6M" = "seagreen",
    "15M" = "darkorange",
    "21M" = "darkviolet"
    
  )
)

# 检查annotations数据框的列名
print(colnames(annotations))





# 确保 pheatmap 中的 `annotation_row` 参数指向正确的列
pheatmap_obj <- pheatmap(t(heatmap_matrix_log_sorted),
                         #scale = "column",
                         cluster_rows = TRUE,  # 对行进行聚类
                         cluster_cols = FALSE,  # 不对列进行聚类
                         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                         color = colorRampPalette(c("white", "firebrick3"))(50),
                         
                         main = "Intensity of nucleosides in mouse brain region",
                         fontsize_row = 8,
                         fontsize_col = 5,
                         angle_col = 45,
                         cellwidth = 1,  # 设置格子宽度
                         cellheight = 10,  # 设置格子高度
                         show_colnames =FALSE,
                         show_rownames = TRUE,
                         border = F,
                         cutree_rows = 4,
                         cutree_cols = 1,
                         annotation_col  = annotations_sorted,  # 添加行注释
                         annotation_colors = annotation_colors,
                         legend_position = "top" )  # 使用手动设置的颜色







# 定义保存pheatmap图像的函数
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(inherits(x, 'pheatmap'))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(pheatmap_obj, "heatmap.pdf")

