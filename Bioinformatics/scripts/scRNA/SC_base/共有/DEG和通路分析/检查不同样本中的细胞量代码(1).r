# ==============================================================================
# 样本细胞数统计模块 (对齐样本顺序 & 自动补零)
# ==============================================================================
rm(list = ls())
gc()

# 加载依赖包
library(Seurat)
library(dplyr)

# 1. 路径配置
rds_file <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/Glioma_3B_Fully_Annotated.rds"
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/通路分析/PPAT/6_DEG差异性分析_富集"

# 如果工作目录不存在则自动创建
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# 2. 读取完整注释的 RDS 文件
cat(">>> 正在加载单细胞数据对象...\n")
sc_obj <- readRDS(rds_file)

# ==============================================================================
# 统计方式 A：统计【所有细胞】并确立统一的样本排序模板
# ==============================================================================
cat(">>> 正在统计全集数据的样本细胞数，并建立基准排序模板...\n")

# 提取所有样本并转换为字符型（防止 Factor 导致合并时的类型冲突）
sample_counts_all <- as.data.frame(table(as.character(sc_obj$SampleID)), stringsAsFactors = FALSE)
colnames(sample_counts_all) <- c("SampleID", "Total_Cell_Count")

# 按总细胞量降序排列，这一顺序将作为后续所有 CSV 表格的基准行顺序
sample_counts_all <- sample_counts_all %>% arrange(desc(Total_Cell_Count))

print("=== 全数据集各样本细胞数 ===")
print(sample_counts_all)

# 导出全数据集 CSV
write.csv(sample_counts_all, file.path(work_dir, "All_Cells_Sample_Counts.csv"), row.names = FALSE)

# 提取排序好的 SampleID 列作为后续子集比对的“模板”
master_samples <- sample_counts_all %>% select(SampleID)

# ==============================================================================
# 统计方式 B：循环统计【多个指定亚群】（对齐模板排序，缺失补 0）
# ==============================================================================
# 确保当前的 Idents 设置为包含亚群信息的列（请根据实际 metadata 列名确认）
Idents(sc_obj) <- "CellType" 

# 在此处指定你需要单独统计的细胞类型列表
target_cell_types <- c("AC-like", "NPC-like")

cat("\n>>> 正在批量提取指定亚群并统计样本分布...\n")

for (ctype in target_cell_types) {
  cat(sprintf("\n-> 正在处理 [%s] 亚群...\n", ctype))
  
  if (ctype %in% levels(Idents(sc_obj))) {
    
    # 获取子集
    sub_obj <- subset(sc_obj, idents = ctype)
    
    # 统计该亚群中的样本细胞数
    sample_counts_sub <- as.data.frame(table(as.character(sub_obj$SampleID)), stringsAsFactors = FALSE)
    
    # 处理列名，替换连字符避免命名报错
    ctype_safe <- gsub("-", "_", ctype)
    count_col_name <- paste0(ctype_safe, "_Count")
    colnames(sample_counts_sub) <- c("SampleID", count_col_name)
    
    # 【核心逻辑】：使用 left_join 将该亚群结果合并到 master_samples 模板中
    # 这样可严格保持原有的总细胞数降序排列，并带出所有样本名称
    aligned_counts <- left_join(master_samples, sample_counts_sub, by = "SampleID")
    
    # 将匹配不到的缺失样本（NA）用 0 填充
    aligned_counts[[count_col_name]][is.na(aligned_counts[[count_col_name]])] <- 0
    
    print(paste0("=== ", ctype, " 亚群各样本细胞数 ==="))
    print(aligned_counts)
    
    # 导出 CSV，生成的表格行顺序会与 All_Cells 报表完全一致
    out_filename <- paste0(ctype_safe, "_Sample_Counts.csv")
    write.csv(aligned_counts, file.path(work_dir, out_filename), row.names = FALSE)
    
  } else {
    cat(sprintf("⚠️ 警告：当前数据的 Idents 中未找到名为 '%s' 的细胞类型，已跳过。\n", ctype))
  }
}

cat("\n>>> 🎉 所有统计完毕，包含所有样本且顺序一致的报表已生成！\n")