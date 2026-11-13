# ==============================================================================
# 胶质瘤单细胞分析管线 - Step 1: 读取原始数据与多维基础质控 (适于 10X MTX 三件套)
# ==============================================================================

# 加载必要的包
library(Seurat)
library(stringr)
library(ggplot2)
library(patchwork)
library(Matrix)      
library(data.table)  

# ==============================================================================
# 第一阶段：环境准备与 10X 三件套智能读取
# ==============================================================================
# 设置绝对路径（请将 GSE223063 替换为您当前的实际编号）
raw_data_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/0_原始数据"
out_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/1_读取与基础质控"

if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
setwd(out_dir)

# 1. 获取 0_原始数据 下的所有样本子文件夹路径
sample_dirs <- list.dirs(raw_data_dir, full.names = TRUE, recursive = FALSE)

if(length(sample_dirs) == 0) {
  stop("错误：在 0_原始数据 目录下找不到任何样本子文件夹，请检查路径。")
}

cat(">>> 开始读取 10X 三件套并创建 Seurat 对象列表...\n")

# 2. 批量读取 (动态匹配后缀，忽略前缀)
seurat_list <- lapply(sample_dirs, function(dir_path) {
  
  # 提取文件夹名称作为样本名，例如从路径中提取 "GSM5518629"
  raw_id <- basename(dir_path)
  
  print(paste("========================================"))
  print(paste("正在处理原始样本:", raw_id))
  
  # 在当前文件夹中，利用正则表达式匹配特定后缀的文件
  mtx_path <- list.files(dir_path, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE, ignore.case = TRUE)
  bc_path <- list.files(dir_path, pattern = "barcodes\\.tsv\\.gz$", full.names = TRUE, ignore.case = TRUE)
  feat_path <- list.files(dir_path, pattern = "features\\.tsv\\.gz$", full.names = TRUE, ignore.case = TRUE)
  
  # 检查是否三件套齐全
  if (length(mtx_path) == 0 || length(bc_path) == 0 || length(feat_path) == 0) {
    stop(paste("\n错误: 样本", raw_id, "缺少核心的 mtx, barcodes 或 features 文件！\n"))
  }
  
  # 使用 ReadMtx 读取，取匹配到的第一个文件（避免同一文件夹下有多余文件干扰）
  raw_counts <- ReadMtx(mtx = mtx_path[1], cells = bc_path[1], features = feat_path[1])
  
  # 创建 Seurat 对象，直接使用文件夹名 raw_id 作为 project 名称
  srt <- CreateSeuratObject(counts = raw_counts, 
                            project = raw_id, 
                            min.cells = 3, 
                            min.features = 200)
  
  # 注入元数据
  srt$SampleID <- raw_id
  srt$Disease_Status <- "Glioma"
  
  rm(raw_counts)
  gc()
  
  print(paste("✅ 样本", raw_id, "读取成功！"))
  return(srt)
})

# 提取所有样本名用于合并时的防重复条形码前缀
unique_raw_names <- basename(sample_dirs)
names(seurat_list) <- unique_raw_names

cat("\n>>> 正在合并所有样本...\n")
# 合并所有的 Seurat 对象
glioma_merged <- merge(x = seurat_list[[1]], 
                       y = seurat_list[-1], 
                       add.cell.ids = unique_raw_names, 
                       project = "Glioma_Merged")

# Seurat v5 层连接
glioma_merged <- JoinLayers(glioma_merged)

rm(seurat_list)
gc()

saveRDS(glioma_merged, file = "Glioma_merged_raw.rds")
cat(">>> 🎉 所有样本已合并！原始底盘已保存为: Glioma_merged_raw.rds\n\n")


# ==============================================================================
# 第二阶段：计算核心质控指标 (MT, RB, HB, Novelty)
# ==============================================================================
cat(">>> 正在计算线粒体、核糖体、红细胞基因比例及转录本复杂度...\n")

glioma_merged[["percent.mt"]] <- PercentageFeatureSet(glioma_merged, pattern = "^MT-")
glioma_merged[["percent.rb"]] <- PercentageFeatureSet(glioma_merged, pattern = "^RP[SL]")
glioma_merged[["percent.hb"]] <- PercentageFeatureSet(glioma_merged, pattern = "^HB[^(P)]")
glioma_merged$log10GenesPerUMI <- log10(glioma_merged$nFeature_RNA) / log10(glioma_merged$nCount_RNA)


# ==============================================================================
# 第三阶段：预质控 (Pre-QC) 可视化输出
# ==============================================================================
cat(">>> 正在生成 Pre-QC 图表并保存为 PDF...\n")

my_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                  legend.position = "none")

# 图表将自动以原始样本名为单位进行汇总展示
p1 <- VlnPlot(glioma_merged, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "SampleID", 
              ncol = 3, pt.size = 0) & my_theme

p2 <- VlnPlot(glioma_merged, 
              features = c("percent.rb", "percent.hb", "log10GenesPerUMI"), 
              group.by = "SampleID", 
              ncol = 3, pt.size = 0) & my_theme

p3 <- FeatureScatter(glioma_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "SampleID")
p4 <- FeatureScatter(glioma_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "SampleID")

pdf("01_PreQC_Metrics_Basic.pdf", width = 14, height = 6)
print(p1)
dev.off()

pdf("02_PreQC_Metrics_Extra.pdf", width = 14, height = 6)
print(p2)
dev.off()

pdf("03_PreQC_Scatter.pdf", width = 14, height = 6)
print(p3 + p4)
dev.off()


# ==============================================================================
# 第四阶段：执行硬性阈值过滤并逐样本统计细胞量
# ==============================================================================
cat(">>> 开始执行细胞级硬性阈值过滤...\n")

cells_before_per_sample <- table(glioma_merged$SampleID)

# 胶质瘤特征严格过滤 
glioma_clean <- subset(glioma_merged, 
                       subset = nFeature_RNA > 300 & 
                         nFeature_RNA < 8000 &     # 建议调高
                         nCount_RNA >= 500 &
                         nCount_RNA <= 80000 &     # 建议调高
                         percent.mt < 15 & 
                         percent.hb < 5 & 
                         log10GenesPerUMI > 0.8)

cells_after_per_sample <- table(glioma_clean$SampleID)

integrated_patients <- sort(unique(glioma_merged$SampleID))

cat("\n✅ 过滤完成！各样本细胞量变化明细如下：\n")
cat("================================================================\n")
cat(sprintf("%-25s %-12s %-12s %-10s\n", "SampleID", "Before_QC", "After_QC", "Removed"))
cat("================================================================\n")

for (s in integrated_patients) {
  before <- as.numeric(cells_before_per_sample[s])
  after <- ifelse(s %in% names(cells_after_per_sample), as.numeric(cells_after_per_sample[s]), 0)
  removed <- before - after
  cat(sprintf("%-25s %-12d %-12d %-10d\n", s, before, after, removed))
}

cat("================================================================\n")
total_before <- sum(cells_before_per_sample)
total_after <- sum(cells_after_per_sample)
cat(sprintf("%-25s %-12d %-12d %-10d\n", "TOTAL", total_before, total_after, total_before - total_after))
cat("================================================================\n\n")

cat(sprintf(
  "%-20s %-10s %-10s %-10s %-12s\n",
  "SampleID", "Before", "After", "Removed", "Retained(%)"
))

for (s in integrated_patients) {
  before <- as.numeric(cells_before_per_sample[s])
  after <- if (s %in% names(cells_after_per_sample)) {
    as.numeric(cells_after_per_sample[s])
  } else {
    0
  }
  
  retained_pct <- 100 * after / before
  
  cat(sprintf(
    "%-20s %-10d %-10d %-10d %-12.2f\n",
    s, before, after, before - after, retained_pct
  ))
}


# ==============================================================================
# 第五阶段：过滤后 (Post-QC) 验证与对象保存
# ==============================================================================
cat(">>> 生成 Post-QC 验证图表...\n")

if (!dir.exists("Figures_QC")) { dir.create("Figures_QC") }

p_vln_post <- VlnPlot(glioma_clean, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "log10GenesPerUMI"), 
                      group.by = "SampleID",
                      ncol = 5, pt.size = 0) & my_theme

ggsave("Figures_QC/04_QC_Violin_Post.pdf", plot = p_vln_post, width = 18, height = 5)

output_clean_rds <- "Glioma_merged_clean.rds"
saveRDS(glioma_clean, file = output_clean_rds)

cat(">>> 高质量数据已成功保存为:", output_clean_rds, "\n")
cat(">>> 第一步 (读取与基础质控) 彻底完成！\n")