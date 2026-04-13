setwd("D:/Koziol_lab/Spatial-Nucleotide-Metabolism-Aging/Figure3")
# 加载必要的包
library(ggplot2)
library(reshape2)
library(dplyr)
library(ropls)
library(cowplot)
library(readxl)
library(ggrepel)
library(tidyr)
library(tidyverse)

# 读取 Excel 数据
file_path <- "sum_normalized_data.xlsx"
data <- read_excel(file_path)

# 转换为宽格式
wide_data <- data %>%
  pivot_wider(names_from = rNs, values_from = mean_intensity)
write.csv(wide_data,"wide_data.csv")
# 查看转换后的数据
print(wide_data)

# 设置年龄顺序
age_levels <- c("6M", "15M", "21M")
wide_data$age <- factor(wide_data$age, levels = age_levels)
wide_data <- wide_data[, !names(wide_data) %in% "bio"]
#######去除一些脑区

# 检查数据

# 提取独特的区域
regions <- unique(wide_data$region)
# 定义图例顺序
age_order <- c("6M", "15M", "21M")
colors <- c("red", "blue", "green")  # 颜色与顺序对
############################################################t-test和p-value
head(wide_data)


# 假设wide_data是您的数据集
# 筛选出我们关心的年龄组：6M, 15M, 21M
age_groups_of_interest <- c("6M", "15M", "21M")
region_data_filtered <- wide_data %>%
  filter(age %in% age_groups_of_interest)  # 仅保留6M, 15M, 21M这三组
# 存储最终结果的列表
result_list <- list()
# 按照region分组
regions <- unique(region_data_filtered$region)

# 循环处理每个region
for (reg in regions) {
  # 对每个region进行子集化
  region_data <- region_data_filtered %>% filter(region == reg)
  
  # 获取代谢物的列名（从第3列开始是代谢物数据）
  metabolite_columns <- colnames(region_data)[3:ncol(region_data)]
  
  # 进行t检验并计算logFC
  for (metabolite in metabolite_columns) {
    
    # 获取代谢物数据
    metabolite_data <- region_data %>%
      dplyr::select(age, !!sym(metabolite)) %>%
      drop_na()  # 删除缺失数据
    
    # 获取年龄组
    age_groups <- unique(metabolite_data$age)
    
    # 只计算指定的三组比较：6M/15M, 6M/21M, 15M/21M
    if (length(age_groups) > 1) {
      # 执行 pairwise t-test 计算
      t_test_results <- pairwise.t.test(metabolite_data[[metabolite]], metabolite_data$age, 
                                        p.adjust.method = "BH", 
                                        levels = age_groups_of_interest)  # 指定只计算这三组比较
      
      # 获取p值矩阵并转换为数据框格式
      p_values <- as.data.frame(t_test_results$p.value)
      p_values$age_group_1 <- rownames(p_values)
      p_values_long <- gather(p_values, key = "age_group_2", value = "p_value", -age_group_1)
      
      # 去除 age_group_1 == age_group_2 的行和 NA 值
      #p_values_long <- p_values_long %>%
      #  filter(age_group_1 != age_group_2 & !is.na(p_value))
      
      # 计算logFC（fold change的对数转换）
      p_values_long$log2FC <- NA
      for (i in 1:nrow(p_values_long)) {
        group1_data <- region_data %>% filter(age == p_values_long$age_group_1[i]) %>% pull(metabolite)
        group2_data <- region_data %>% filter(age == p_values_long$age_group_2[i]) %>% pull(metabolite)
        group1_mean<- mean(group1_data) 
        group2_mean<- mean(group2_data)
      
        # 计算fold change
        if(group1_mean!=0 & group2_mean!=0){
          fold_change <- group1_mean/group2_mean
        }
        else{
          fold_change<-NA
        }
        
        p_values_long$log2FC[i] <- log2(fold_change)
      }
      
      # 添加region和metabolite信息
      p_values_long$region <- reg
      p_values_long$metabolite <- metabolite
      
      # 将结果存储在结果列表中
      result_list[[paste(reg, metabolite, sep = "_")]] <- p_values_long
    }
  }
}

# 将所有结果合并为一个数据框
final_results <- bind_rows(result_list)
final_results <- na.omit(final_results)
# 筛选出符合条件的行：age_group_1 = 6M 且 (age_group_2 = 15M 或 age_group_2 = 21M)
final_results_filtered <- final_results %>%
  filter((age_group_1 == "15M" & age_group_2 == "6M") |
           (age_group_1 == "21M" & age_group_2 == "6M") |
           (age_group_1 == "21M" & age_group_2 == "15M"))


# 显示最终结果
head(final_results_filtered)
#=================================================分脑区多组火山差异图
plots <- list()

for (ron in regions) {
  
  region_results_filtered <- final_results_filtered[final_results_filtered$region == ron, ]
  
  region_results_filtered <- region_results_filtered %>%
    mutate(group = case_when(
      age_group_1 == "15M" & age_group_2 == "6M" ~ "15M/6M",
      age_group_1 == "21M" & age_group_2 == "6M" ~ "21M/6M",
      age_group_1 == "21M" & age_group_2 == "15M" ~ "21M/15M",
      TRUE ~ NA_character_
    )) %>%
    filter(is.finite(log2FC) & abs(log2FC) > 1) %>%
    mutate(sig_color = case_when(
      p_value < 0.05 & log2FC > 1 ~ "Up",
      p_value < 0.05 & log2FC < 1 ~ "Down",
      TRUE ~ "Not Significant"
    ))
  
  top10sig <- bind_rows(
    region_results_filtered %>% filter(group == "15M/6M", p_value < 0.05) %>% distinct(metabolite, .keep_all = TRUE) %>% top_n(10, abs(log2FC)),
    region_results_filtered %>% filter(group == "21M/6M", p_value < 0.05) %>% distinct(metabolite, .keep_all = TRUE) %>% top_n(10, abs(log2FC)),
    region_results_filtered %>% filter(group == "21M/15M", p_value < 0.05) %>% distinct(metabolite, .keep_all = TRUE) %>% top_n(10, abs(log2FC))
  )
  
  region_results_filtered$size <- ifelse(region_results_filtered$metabolite %in% top10sig$metabolite, 2, 1)
  
  dt <- filter(region_results_filtered, size == 1)
  
  dfcol <- data.frame(
    x = c(1:3),
    y = 6,
    label = c("15M/6M", "21M/15M", "21M/6M")
  )
  
  mycol <- c("#E64B357F", "#4DBBD57F", "#00A0877F")
  
  plot <- ggplot() +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey50", size = 0.8) +
    
    geom_jitter(data = dt, aes(x = group, y = log2FC, color = sig_color), size = 5, width = 0.4, alpha = 0.75) +
    geom_jitter(data = top10sig, aes(x = group, y = log2FC, color = sig_color), size = 6, width = 0.4, alpha = 0.95) +
    
    geom_tile(data = dfcol, aes(x = x, y = y), height = 0.8, width = 0.9, fill = mycol, color = "black", alpha = 0.6) +
    geom_text(data = dfcol, aes(x = x, y = y, label = label), size = 9, color = "black", fontface = "bold") +
    
    geom_text_repel(data = top10sig, aes(x = group, y = log2FC, label = metabolite),
                    size = 7, box.padding = 0.3, max.overlaps = 50,
                    arrow = arrow(length = unit(0.015, "npc"), type = "open", ends = "last")) +
    
    scale_color_manual(
      name = "Significance",
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "Not Significant" = "grey60"
      )
    ) +
    
    labs(x = NULL, y = "log2 Fold Change", title = ron) +
    
    theme_minimal(base_size = 18) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 22),
      axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black", size = 1),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = margin(5, 5, 5, 5),
      legend.box.margin = margin(10, 10, 10, 10),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_x_discrete()
  
  plots[[ron]] <- plot
}

# 合并子图
panel_plot <- plot_grid(plotlist = plots, ncol = 4)

# 保存高清图
ggsave("plsda_panel_with_volcano.pdf", panel_plot, width = 28, height = 22)




############################################################vip
# 存储 PLS-DA 子图
plots <- list()
# 遍历每个区域，计算 VIP 值并绘制柱状图
for (region in regions) {
  # 针对每个区域提取子集
  subset_data <- wide_data[wide_data$region == region, ]
  
  # 确保特征和标签格式正确
  feature_columns <- colnames(subset_data)[3:(ncol(wide_data))]  # 动态选择特征列
  X <- as.matrix(subset_data[, feature_columns])
  Y <- as.factor(subset_data$age)
  
  # 执行 PLS-DA 分析
  plsda_model <- opls(X, Y, predI = 2, permI = 0)
  
  # 提取 VIP 值
  vip_values <- plsda_model@vipVn
  vip_df <- data.frame(Feature = names(vip_values), VIP = vip_values)
  
  # 排序并选择前 10 个特征
  top_vip_df <- vip_df[order(-vip_df$VIP), ][1:10, ]
  
  # 创建每个区域的 ggplot 图，背景为白色
  plot <- ggplot(top_vip_df, aes(x = reorder(Feature, VIP), y = VIP)) +
    geom_bar(stat = "identity", fill = "steelblue",width = 0.6) +
    coord_flip() +
    ggtitle(paste("", region)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加黑色边框
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),  # 增加标题字体大小
      axis.title = element_text(size = 20),  # 增加坐标轴标题字体大小
      axis.text = element_text(size = 20),  # 增加坐标轴刻度字体大小
      strip.text = element_text(size = 20)) +  # 调整分面标题字体大小
    labs(x = "Features", y = "VIP Score")
  
  plots[[region]] <- plot
}

# 使用 cowplot 包合并所有子图
panel_plot <- plot_grid(plotlist = plots, ncol = 4)  # 调整列数以适配显示

# 保存结果
ggsave("vip_panel_plot_ropls.pdf", panel_plot, width = 20, height = 15)

####################################################################################plsda


plots <- list()
model_summary <- data.frame()  # ⬅️ 用于存储所有指标

for (region in regions) {
  # 提取数据
  subset_data <- wide_data[wide_data$region == region, ]
  feature_columns <- colnames(subset_data)[3:ncol(wide_data)]
  X <- as.matrix(subset_data[, feature_columns])
  Y <- as.factor(subset_data$age)
  
  # 拟合 PLS-DA
  plsda_model <- opls(X, Y, predI = 2, permI = 0)
  
  # 得分矩阵
  score_df <- as.data.frame(plsda_model@scoreMN)
  colnames(score_df) <- c("Comp1", "Comp2")
  score_df$age <- Y
  
  # 提取模型指标
  R2Y <- round(plsda_model@summaryDF$`R2Y(cum)`[1], 3)
  Q2  <- round(plsda_model@summaryDF$`Q2(cum)`[1], 3)
  
  # ⬅️ 保存到汇总表
  model_summary <- rbind(model_summary, data.frame(
    Region = region,
    R2Y = R2Y,
    Q2 = Q2
  ))
  
  # 绘图
  plot <- ggplot(score_df, aes(x = Comp1, y = Comp2, color = age)) +
    geom_point(size = 3, alpha = 0.9) +
    stat_ellipse(aes(group = age), linetype = "dashed", size = 1.2) +
    ggtitle(paste("PLS-DA -", region, "\nR² =", R2Y, ", Q² =", Q2)) +
    labs(x = "Component 1", y = "Component 2", color = "Age Group") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  plots[[region]] <- plot
}

# 合并子图
panel_plot <- plot_grid(plotlist = plots, ncol = 4)

# 保存图像
ggsave("plsda_score_panel_plot_with_ellipse.pdf", panel_plot, width = 20, height = 15)

# ⬅️ 导出模型指标表格
write.csv(model_summary, "PLSDA_model_summary.csv", row.names = FALSE)


#######################################################绝对丰度变化

averaged_data <- wide_data %>%
  group_by(age, region) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) 

# 检查结果
head(averaged_data)

# 转换为长格式
melted_data <- melt(averaged_data, id.vars = c("region", "age"), 
                    variable.name = "Metabolite", 
                    value.name = "Intensity")

# 获取所有代谢物的名字
metabolites <- unique(melted_data$Metabolite)
metabolites<-top_vip_df$Feature
# 循环绘制每个代谢物的分面图
# 设置年龄顺序
age_levels <- c("6M", "15M", "21M")

# 调整年龄为因子，并指定顺序
melted_data$age <- factor(melted_data$age, levels = age_levels)

for (metabolite in metabolites) {
  # 提取当前代谢物的数据
  subset_data <- melted_data[melted_data$Metabolite == metabolite, ]
  
  # 创建分面图，按脑区（region）分面
  plot <- ggplot(subset_data, aes(x = age, y = Intensity, color = region, group = region)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~ region, ncol = 4, scales = "free_y") +  # 按 region 分面
    labs(title = paste("Metabolite:", metabolite),
         x = "Aging Group",
         y = "Intensity") +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 20), 
      legend.title = element_text(size = 20),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1), # 黑色边框
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),  # 加粗并居中
      strip.text = element_text(size = 20, face = "bold"),  # 分面标题加粗
      axis.text = element_text(size = 20),  # 坐标轴文字大小
      axis.title = element_text(size = 20)  # 坐标轴标题大小
    )
  
  # 保存图形
  #file_name <- paste0("metabolite_", gsub(" ", "_", metabolite), "_facet_plot.pdf")
  #ggsave(file_name, plot, width = 18, height = 10, dpi = 300)
  
  # 显示当前代谢物的图形（可选）
  #print(plot)
}

# 筛选数据中属于 VIP 前 10 的代谢物
filtered_data <- melted_data[melted_data$Metabolite %in% top_vip_df$Feature, ]


# 创建分面图，按脑区（region）分面，仅显示 VIP 前 10 的代谢物
plot <- ggplot(filtered_data, aes(x = age, y = log10(Intensity+1), color = Metabolite, group = Metabolite)) +
  geom_line(size = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ region, ncol = 7, scales = "free_y") +  # 每个脑区一个分面
  labs(title = "Intensity Changes for Top 10 VIP Metabolites Across Aging Groups",
       x = "Aging Group",
       y = "Intensity(log10)") +
  theme_minimal(base_size = 20) +  # Increase base font size
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 黑色边框
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),  # 图标题加粗
    strip.text = element_text(size = 20, face = "bold"),  # 分面标题加粗
    axis.text = element_text(size = 20),  # 坐标轴文字大小
    axis.title = element_text(size = 20),  # 坐标轴标题大小
    legend.position = "bottom",
    legend.title = element_blank()
  )

# 保存图形
ggsave("top_10_vip_metabolites_by_region_facet_plot.pdf", plot, width = 20, height = 15, dpi = 300)




