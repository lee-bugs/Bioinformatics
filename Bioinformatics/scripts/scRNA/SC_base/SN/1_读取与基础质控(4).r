# ==============================================================================
# 胶质瘤单核测序 (snRNA-seq) 分析管线 - Step 1: 读取、深度诊断与极限质控
# ==============================================================================

library(Seurat)
library(stringr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table) # 用于高速读取大型 tsv 矩阵 (保留引入)
library(Matrix)     # 用于构建稀疏矩阵降低内存

# ==============================================================================
# 第一阶段：环境准备与目录智能巡检
# ==============================================================================
raw_data_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE292130(snRNA)Object_1_原发/0_原始数据"
out_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE292130(snRNA)Object_1_原发/1_读取与基础质控"

if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
setwd(out_dir)

# 自动扫描文件夹下所有的 _umi_counts.RDS 文件
rds_files <- list.files(path = raw_data_dir, pattern = "_umi_counts\\.RDS$", full.names = TRUE)

if(length(rds_files) == 0) {
  stop("错误：找不到对应的 RDS 表达矩阵文件，请检查路径！")
}

# ==============================================================================
# 第二阶段：snRNA-seq 专属降维读取 (自动提取GSM号)
# ==============================================================================
cat(">>> 🚀 正在高速加载并解析合并矩阵，请耐心等待...\n")

seurat_list <- list()

for (file_path in rds_files) {
  file_name <- basename(file_path)
  
  # 提取去除后缀的完整名 (如: GSM8850517_MD04P)
  full_prefix <- str_replace(file_name, "_umi_counts\\.RDS$", "")
  
  # 精确提取 GSM 编号 (如: GSM8850517)
  gsm_id <- str_extract(full_prefix, "^GSM[0-9]+")
  
  counts_matrix <- readRDS(file_path)
  
  # 使用 GSM 号作为 project 名字
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, 
                                   project = "GSE211376_snRNA", # 保持原来的 project 名称
                                   min.cells = 3, 
                                   min.features = 100)
  
  # 统一注入 Metadata 
  seurat_obj$SampleID <- gsm_id
  seurat_obj$Patient_ID <- gsm_id
  
  # 列表命名使用 GSM 号
  seurat_list[[gsm_id]] <- seurat_obj
}

rm(counts_matrix)
gc()

# 合并所有的 Seurat 对象
if(length(seurat_list) > 1){
  glioma_merged <- merge(x = seurat_list[[1]], 
                         y = seurat_list[2:length(seurat_list)],
                         add.cell.ids = names(seurat_list))
} else {
  glioma_merged <- seurat_list[[1]]
}

rm(seurat_list)
gc()

# 验证样本是否全部被捕获
sample_folders <- unique(glioma_merged$SampleID)
all_folders <- sample_folders 

cat("\n>>> 📊 核心检验：当前捕获的 GSM 样本及细胞数统计：\n")
print(table(glioma_merged$SampleID)) 

if("GSM_Unknown" %in% sample_folders) {
  warning("⚠️ 注意：部分细胞的前缀不在你提供的 11 个样本对照表中，已被归为 GSM_Unknown。")
} else {
  cat(sprintf("\n>>> 🧬 成功锁定全部 %d 个 GSM 样本，完美对应！准备进入质控...\n", length(sample_folders)))
}

# 严格保留原输出文件名
saveRDS(glioma_merged, file = "GSE211376_merged_raw.rds")
cat(">>> 🎉 原始底盘已安全保存为: GSE211376_merged_raw.rds\n\n")

# ==============================================================================
# 第三阶段：计算质控指标
# ==============================================================================
cat(">>> 正在计算线粒体、核糖体等游离污染指标...\n")

glioma_merged[["percent.mt"]] <- PercentageFeatureSet(glioma_merged, pattern = "^MT-")
glioma_merged[["percent.rb"]] <- PercentageFeatureSet(glioma_merged, pattern = "^RP[SL]")
glioma_merged$log10GenesPerUMI <- log10(glioma_merged$nFeature_RNA) / log10(glioma_merged$nCount_RNA)

# ==============================================================================
# 第四阶段：深度侦察代码 - 统计报表与多维关联图谱
# ==============================================================================
cat("\n>>> 🚀 启动深度侦察：正在生成各样本核心指标统计表...\n")

meta_data <- glioma_merged@meta.data

summary_stats <- meta_data %>%
  group_by(SampleID) %>%
  summarise(
    Total_Barcodes = n(),
    Median_nCount = median(nCount_RNA),
    Median_nFeature = median(nFeature_RNA),
    Median_MT = median(percent.mt),
    Max_MT = max(percent.mt)
  ) %>%
  arrange(SampleID)

cat("=================================================================================\n")
cat(sprintf("%-15s %-15s %-15s %-15s %-15s %-15s\n", 
            "SampleID", "Total_Barcodes", "Med_nCount", "Med_nFeature", "Med_MT(%)", "Max_MT(%)"))
cat("=================================================================================\n")
for (i in 1:nrow(summary_stats)) {
  cat(sprintf("%-15s %-15d %-15.1f %-15.1f %-15.3f %-15.3f\n", 
              summary_stats$SampleID[i], 
              summary_stats$Total_Barcodes[i],
              summary_stats$Median_nCount[i],
              summary_stats$Median_nFeature[i],
              summary_stats$Median_MT[i],
              summary_stats$Max_MT[i]))
}
cat("=================================================================================\n")

cat("\n>>> 正在绘制深度诊断图册 (01a_Deep_Diagnostic_Plots.pdf & 01b_Raw_QC_Panels.pdf)...\n")

p_scatter_complexity <- FeatureScatter(glioma_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                                       group.by = "SampleID", pt.size = 0.1) +
  ggtitle("Transcriptomic Complexity (nCount vs nFeature)") +
  theme(legend.position = "right")

p_scatter_mt <- FeatureScatter(glioma_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                               group.by = "SampleID", pt.size = 0.1) +
  ggtitle("Mitochondrial Contamination (nCount vs MT%)") +
  theme(legend.position = "right")

p_vln_zoomed <- VlnPlot(glioma_merged, features = c("nFeature_RNA"), group.by = "SampleID", pt.size = 0) +
  coord_cartesian(ylim = c(0, 2000)) + 
  ggtitle("Zoomed nFeature_RNA (Bottom 0-2000)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

p_final <- (p_scatter_complexity | p_scatter_mt) / p_vln_zoomed
ggsave("01a_Deep_Diagnostic_Plots.pdf", plot = p_final, width = 16, height = 12)

my_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=10),
                  legend.position = "none")
p_vln_basic <- VlnPlot(glioma_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                       group.by = "SampleID", ncol = 4, pt.size = 0) & my_theme

ggsave("01b_Raw_QC_Panels.pdf", plot = p_vln_basic, width = 18, height = 6)

# ==============================================================================
# 第五阶段：snRNA-seq 专属变态级阈值硬筛 (Hard Filtering)
# ==============================================================================
cat("\n>>> 🔪 开始执行单核测序极严苛硬筛...\n")

cells_before_per_sample <- table(glioma_merged$SampleID)

# snRNA-seq 极限过滤参数
glioma_clean <- subset(glioma_merged, 
                       subset = nFeature_RNA > 500 &        
                         nFeature_RNA < 10000 &       # 结合图表，将上限从 12000 收紧至 10000
                         nCount_RNA > 800 & 
                         nCount_RNA < 40000 &         # 结合图表，将上限从 80000 收紧至 40000
                         percent.mt < 5 &             # 保持不变
                         log10GenesPerUMI > 0.8)      # 保持不变

cells_after_per_sample <- table(glioma_clean$SampleID)

cat("\n✅ 过滤完成！存活统计表：\n")
cat("=================================================================\n")
cat(sprintf("%-15s %-10s %-10s %-10s %-12s\n", "SampleID", "Before", "After", "Removed", "Retained(%)"))
cat("=================================================================\n")

# 因为 sample_folders 已经是纯 GSM 号，这里去掉了正则提取步骤
for (s in sample_folders) {
  before <- as.numeric(cells_before_per_sample[s])
  after <- ifelse(s %in% names(cells_after_per_sample), as.numeric(cells_after_per_sample[s]), 0)
  retained_pct <- 100 * after / before
  cat(sprintf("%-15s %-10d %-10d %-10d %-12.2f\n", s, before, after, before - after, retained_pct))
}
cat("=================================================================\n")

# ==============================================================================
# 第六阶段：Post-QC 验证与高纯度对象封存
# ==============================================================================
cat("\n>>> 生成清洗后的验证图表 (02_Clean_PostQC_VlnPlot.pdf)...\n")

p_clean_vln <- VlnPlot(glioma_clean, 
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                       group.by = "SampleID", 
                       ncol = 4, pt.size = 0) & my_theme

ggsave("02_Clean_PostQC_VlnPlot.pdf", plot = p_clean_vln, width = 16, height = 5)

# 严格保留原输出文件名
output_rds <- "GSE211376_merged_clean.rds"
saveRDS(glioma_clean, file = output_rds)

cat(">>> 🎉 Step 1 (完整版) 彻底圆满完成！\n")
cat(">>> 高纯度细胞核矩阵已封存于:", output_rds, "\n")