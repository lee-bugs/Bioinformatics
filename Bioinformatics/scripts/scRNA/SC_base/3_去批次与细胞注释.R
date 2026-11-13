# ==============================================================================
# 胶质瘤单细胞 SOP: Step 3 - 纯净流形定型、极简风降维图谱与细胞注释
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(patchwork)
library(clustree)
library(scales)      
library(ggrepel)     # 用于绘制 Style 3 (虚线牵引)

work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/3.注释新图形"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# ------------------------------------------------------------------------------
# 预定义 1：SCI 级配色面板与颜色扩展函数
# ------------------------------------------------------------------------------
my_color_palette <- c("#e31a1c", "#fd8d3c", "#fecc5c", "#bd0026", "#3182bd", "#9e9ac8", "#31a354", "#74c476", "#bdbdbd")
expand_pal <- function(pal, n) { colorRampPalette(pal)(n) }

# ------------------------------------------------------------------------------
# 预定义 2：图形 2 专属极简坐标轴主题 (theme_sci_arrow)
# ------------------------------------------------------------------------------
theme_sci_arrow <- theme_classic() +
  theme(
    axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ------------------------------------------------------------------------------
# 预定义 3：图形 3 虚线牵引绘图函数 (强制推开距离并绘制虚线)
# ------------------------------------------------------------------------------
plot_dashed_repel <- function(seurat_obj, reduction_name, group_var, pal, title) {
  coords <- as.data.frame(Embeddings(seurat_obj, reduction_name))
  colnames(coords) <- c("DIM_1", "DIM_2")
  coords$Group <- seurat_obj@meta.data[[group_var]]
  
  centroids <- aggregate(cbind(DIM_1, DIM_2) ~ Group, data = coords, FUN = median)
  global_x <- median(coords$DIM_1)
  global_y <- median(coords$DIM_2)
  centroids$dir_x <- centroids$DIM_1 - global_x
  centroids$dir_y <- centroids$DIM_2 - global_y
  dist <- sqrt(centroids$dir_x^2 + centroids$dir_y^2)
  dist[dist == 0] <- 1 
  
  range_x <- max(coords$DIM_1) - min(coords$DIM_1)
  range_y <- max(coords$DIM_2) - min(coords$DIM_2)
  centroids$nudge_x <- (centroids$dir_x / dist) * range_x * 0.15
  centroids$nudge_y <- (centroids$dir_y / dist) * range_y * 0.15
  
  ggplot(coords, aes(x = DIM_1, y = DIM_2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.8) + 
    scale_color_manual(values = pal) +
    geom_text_repel(data = centroids, aes(x = DIM_1, y = DIM_2, label = Group), 
                    nudge_x = centroids$nudge_x, nudge_y = centroids$nudge_y,
                    segment.linetype = "dashed", segment.color = "black", 
                    color = "black", size = 5, fontface = "bold",
                    min.segment.length = 0, box.padding = 0.5) +
    theme_classic() + ggtitle(title) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    guides(color = guide_legend(override.aes = list(size = 4)))
}


cat(">>> [步骤 1] 加载切除病灶后的绝对纯净底盘...\n")
path_3a_input <- "Glioma_3A_ReadyForIntegration.rds"
if (!file.exists(path_3a_input)) stop("找不到输入文件，请确认路径: ", path_3a_input)

glioma_clean <- readRDS(path_3a_input)
# SeuratObject::DefaultAssay(glioma_clean) <- "decontX"
#=========================== 修改默认 assays 为 RNA ============================###
SeuratObject::DefaultAssay(glioma_clean) <- "RNA"
glioma_clean = NormalizeData(glioma_clean)

# ==============================================================================
# [步骤 2] 在纯净底盘上重新寻找高变基因 (HVG) 并执行终极黑名单过滤
# ==============================================================================
cat(">>> [步骤 2] 正在寻找 6000 个高变基因备选池，为 3000 纯净池提供充足冗余...\n")
# 🚨 更改点 1：将 nfeatures 从 4000 扩大到 6000
glioma_clean <- FindVariableFeatures(glioma_clean, selection.method = "vst", nfeatures = 6000, verbose = FALSE)

all_genes <- rownames(glioma_clean)

# (保持原有的黑名单列表不变...)
cycle_genes <- unique(c(Seurat::cc.genes.updated.2019$s.genes, Seurat::cc.genes.updated.2019$g2m.genes, "MKI67", "TOP2A", "UBE2C", "BIRC5", "CENPF", "TK1", "PCNA"))
mito_ribo_hb <- grep("^MT-|^RP[SL]|^MRP[SL]|^HBA|^HBB|^HBD|^HBE1|^HBG|^HBM|^HBQ|^HBZ", all_genes, value = TRUE)
stress_genes <- c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL8", "DUSP1", "DUSP2", "EGR1", "EGR2", "EGR3", "EGR4", "FOS", "FOSB", "FOSL2", "GADD45A", "GADD45B", "GADD45G", "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "IER2", "IER3", "IER5", "IL6", "JUN", "JUNB", "JUND", "NFKBIA", "NR4A1", "NR4A2", "NR4A3", "MYC", "PPP1R15A", "SOCS3", "UBC", "ZFP36")
hsp_dnaj <- grep("^HSP|^DNAJ", all_genes, value = TRUE)
sex_genes <- c("XIST", "TSIX", "RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "USP9Y", "UTY", "ZFY")
platelet_genes <- c("PPBP", "PF4", "TUBB1")
sticky_lncrna <- c("MALAT1", "NEAT1", "KCNQ1OT1")
myelin_genes <- c("MBP", "PLP1", "MOBP", "MOG", "MAG", "CNP")

# 🚨 更改点 2：增强非编码与假基因的正则匹配，增加环境污染与看家基因黑名单
nc_pseudo_genes <- grep("^[A-Z]{2}[0-9]+\\.[0-9]+|^LINC[0-9]+|^MIR[0-9]+|^SNORD[0-9]+|^SNOU[0-9]+|^RN7[A-Z]+|^RNU[0-9]+", all_genes, value = TRUE)
ambient_housekeeping <- c("ACTB", "GAPDH", "B2M", "TMSB4X", "TMSB10", "FTL", "FTH1", "EEF1A1", "EEF2")

# 将新增的 ambient_housekeeping 加入终极黑名单列表
ultimate_blacklist <- intersect(
  unique(c(cycle_genes, mito_ribo_hb, stress_genes, hsp_dnaj, sex_genes, platelet_genes, sticky_lncrna, myelin_genes, nc_pseudo_genes, ambient_housekeeping)), 
  all_genes
)

raw_hvgs <- VariableFeatures(glioma_clean)
clean_hvgs <- setdiff(raw_hvgs, ultimate_blacklist)

# 🚨 更改点 3：截取前 3000 个纯净高变基因
top3000_pure_hvgs <- head(clean_hvgs, 3000)
VariableFeatures(glioma_clean) <- top3000_pure_hvgs

cat(sprintf(">>> 已锁定用于纯净空间重塑的 %d 个高质量 HVG (已过滤看家基因与更多非编码RNA)。\n", length(VariableFeatures(glioma_clean))))
hvg_info <- HVFInfo(glioma_clean)
hvg_info_clean <- hvg_info[top3000_pure_hvgs, ]
hvg_df <- data.frame(Rank = 1:length(top3000_pure_hvgs), Gene = top3000_pure_hvgs, Mean_Expression = hvg_info_clean$mean, Standardized_Variance = hvg_info_clean$variance.standardized)
write.table(hvg_df, file = "00_Purified_Top3000_HVGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# ==============================================================================
# [步骤 3] 缩放、降维与生成拐点图 (ElbowPlot)
# ==============================================================================
cat(">>> [步骤 3] 正在执行 ScaleData 与 PCA...\n")
glioma_clean <- ScaleData(glioma_clean, features = top3000_pure_hvgs, verbose = FALSE)
glioma_clean <- RunPCA(glioma_clean, features = top3000_pure_hvgs, npcs = 50, seed.use = 20260808, verbose = FALSE)

p_elbow <- ElbowPlot(glioma_clean, ndims = 50) + 
  theme_bw() + 
  ggtitle("PCA Elbow Plot (Purified Dataset)") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), panel.grid.minor = element_blank())
ggsave("01_Purified_PCA_ElbowPlot.pdf", plot = p_elbow, width = 8, height = 6)


# ==============================================================================
# 胶质瘤单细胞 SOP: Step 3 (中) - Harmony 去批次与聚类探库
# ==============================================================================
cat("\n=================================================================\n")
cat(">>> [步骤 4] 启动 Harmony2 整合引擎 (N_PCS = 30)...\n")
glioma_clean <- RunHarmony(glioma_clean, group.by.vars = "SampleID", reduction.use = "pca", reduction.save = "harmony", dims.use = 1:30, verbose = FALSE)

cat(">>> [步骤 5] 正在生成双套底层降维坐标 (基础版与参数重塑版)...\n")
# 基础降维系 (供 Style 1 和 Style 3 使用)
glioma_clean <- RunUMAP(glioma_clean, reduction = "harmony", dims = 1:30, seed.use = 20260808, verbose = FALSE)
glioma_clean <- RunTSNE(glioma_clean, reduction = "harmony", dims = 1:30, seed.use = 20260808, verbose = FALSE)

# 🚨 重塑降维系 (专门供图形 2 使用，彻底改变随机种子与超参数，颠覆空间拓扑)
glioma_clean <- RunUMAP(glioma_clean, reduction = "harmony", dims = 1:30, 
                        seed.use = 999, min.dist = 0.8, n.neighbors = 40, 
                        reduction.name = "umap_alt", reduction.key = "UMAPALT_", verbose = FALSE)
glioma_clean <- RunTSNE(glioma_clean, reduction = "harmony", dims = 1:30, 
                        seed.use = 999, perplexity = 50, 
                        reduction.name = "tsne_alt", reduction.key = "tSNEALT_", verbose = FALSE)

# ------------------------------------------------------------------------------
# 批次效应核查图 - 输出 3 种风格 (涵盖 UMAP 与 t-SNE)
# ------------------------------------------------------------------------------
cat(">>> 绘制去批次后的 SampleID 分布图...\n")
sample_pal <- expand_pal(my_color_palette, length(unique(glioma_clean$SampleID)))

# 【批次 UMAP - 3 种风格】
ggsave("02_Integrated_Batch_UMAP_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "umap", group.by = "SampleID", cols = sample_pal, pt.size = 2.5, raster = TRUE, raster.dpi = c(1200, 1200)) + ggtitle("Style 1: Standard UMAP Batch"), width = 9, height = 7)
# 🚨 Style 2 应用极简绘图参数 (去除标签、附带箭头主题、透明度设置)
ggsave("02_Integrated_Batch_UMAP_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "umap_alt", group.by = "SampleID", cols = sample_pal, pt.size = 0.1, alpha = 0.8) + ggtitle("Style 2: Reshaped Minimal UMAP Batch") + theme_sci_arrow + labs(x = "UMAP 1", y = "UMAP 2"), width = 8.5, height = 7)
ggsave("02_Integrated_Batch_UMAP_Style3_Dashed.pdf", plot = plot_dashed_repel(glioma_clean, "umap", "SampleID", sample_pal, "Style 3: Dashed Repel Batch"), width = 9, height = 7)

# 【批次 tSNE - 3 种风格】
ggsave("02_Integrated_Batch_tSNE_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "tsne", group.by = "SampleID", cols = sample_pal, pt.size = 2.5, raster = TRUE, raster.dpi = c(1200, 1200)) + ggtitle("Style 1: Standard t-SNE Batch"), width = 9, height = 7)
# 🚨 Style 2 应用极简绘图参数 
ggsave("02_Integrated_Batch_tSNE_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "tsne_alt", group.by = "SampleID", cols = sample_pal, pt.size = 0.1, alpha = 0.8) + ggtitle("Style 2: Reshaped Minimal t-SNE Batch") + theme_sci_arrow + labs(x = "t-SNE 1", y = "t-SNE 2"), width = 8.5, height = 7)
ggsave("02_Integrated_Batch_tSNE_Style3_Dashed.pdf", plot = plot_dashed_repel(glioma_clean, "tsne", "SampleID", sample_pal, "Style 3: Dashed Repel Batch"), width = 9, height = 7)


# ==============================================================================
# [步骤 6] 构建相邻网络并执行多分辨率扫描
# ==============================================================================
cat("\n>>> [步骤 6] 正在构建全新相邻网络...\n")
glioma_clean <- FindNeighbors(glioma_clean, reduction = "harmony", dims = 1:30, verbose = FALSE)
res_range <- seq(0.1, 1.2, by = 0.1)
glioma_clean <- FindClusters(glioma_clean, resolution = res_range, random.seed = 20260808, verbose = FALSE)

p_clustree <- clustree(glioma_clean, prefix = "RNA_snn_res.") + theme(legend.position = "right") + ggtitle("Clustree on Integrated Data (Res 0.1 - 1.2)")
ggsave("03_Integrated_Clustree.pdf", plot = p_clustree, width = 10, height = 12)


# ==============================================================================
# [步骤 7] 锁定黄金分辨率并生成 [未注释] 聚类概览 (UMAP & tSNE 均 3 种形式)
# ==============================================================================
target_res <- 0.5                                   # 人工可改 。
res_col <- paste0("RNA_snn_res.", target_res)  
Idents(glioma_clean) <- res_col

current_clusters <- as.numeric(as.character(glioma_clean@meta.data[[res_col]]))
correct_levels <- as.character(sort(unique(current_clusters)))
glioma_clean[[res_col]] <- factor(as.character(current_clusters), levels = correct_levels)
Idents(glioma_clean) <- res_col

pal_unanno <- expand_pal(my_color_palette, length(correct_levels))
cat(sprintf(">>> 已锁定目标分辨率 %s，正在绘制未注释聚类图...\n", target_res))

# 【未注释 UMAP 3 种风格】
ggsave("04a_Unannotated_UMAP_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "umap", cols = pal_unanno, label = TRUE, repel = TRUE) + ggtitle("Style 1: Standard UMAP"), width = 9, height = 7)
# 🚨 Style 2：无标签、改用新坐标系、特定的点大小、加箭头主题
ggsave("04a_Unannotated_UMAP_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "umap_alt", cols = pal_unanno, pt.size = 0.3) + ggtitle("Style 2: Reshaped Minimal UMAP Clusters") + theme_sci_arrow + labs(x = "UMAP 1", y = "UMAP 2"), width = 8.5, height = 7)
ggsave("04a_Unannotated_UMAP_Style3_Dashed.pdf", plot = plot_dashed_repel(glioma_clean, "umap", res_col, pal_unanno, "Style 3: Dashed Repel UMAP"), width = 9, height = 7)

# 【未注释 tSNE 3 种风格】
ggsave("04b_Unannotated_tSNE_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "tsne", cols = pal_unanno, label = TRUE, repel = TRUE) + ggtitle("Style 1: Standard t-SNE"), width = 9, height = 7)
ggsave("04b_Unannotated_tSNE_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "tsne_alt", cols = pal_unanno, pt.size = 0.3) + ggtitle("Style 2: Reshaped Minimal t-SNE Clusters") + theme_sci_arrow + labs(x = "t-SNE 1", y = "t-SNE 2"), width = 8.5, height = 7)
ggsave("04b_Unannotated_tSNE_Style3_Dashed.pdf", plot = plot_dashed_repel(glioma_clean, "tsne", res_col, pal_unanno, "Style 3: Dashed Repel t-SNE"), width = 9, height = 7)


# ==============================================================================
# [步骤 8] 胶质瘤专属探索性气泡图 (4色出图) & 空间连续投影图
# ==============================================================================
cat("\n>>> [步骤 8] 绘制探索性 Marker 气泡图 (输出4种经典配色)...\n")
Glioma_SCI_Markers <- list(
  "NPC-like"    = c("SOX4", "SOX2", "SOX11", "DCX", "DLL3", "CD24"),
  "OPC-like"    = c("OLIG2", "TNR", "ALCAM", "OMG", "PLLP", "PDGFRA"),
  "AC-like"     = c("GFAP", "S100B", "HOPX", "SLC1A3", "MLC1"),
  "MES-like"    = c("CD44", "CHI3L1", "ANXA1", "ANXA2", "VIM"),
  "TAM"         = c("C1QA", "C1QB", "C1QC", "TYROBP", "FCER1G", "LST1", "AIF1"),
  "Lymphocytes" = c("CD3D", "CD3E", "TRAC", "NKG7", "KLRD1", "CCL5"),
  "EC"          = c("CLDN5", "VWF", "PECAM1", "CDH5", "FLT1", "PDGFRB")
)
Glioma_SCI_Markers_clean <- lapply(Glioma_SCI_Markers, function(x) x[x %in% rownames(glioma_clean)])
Glioma_SCI_Markers_clean <- Glioma_SCI_Markers_clean[lengths(Glioma_SCI_Markers_clean) > 0]

p_dot_base <- DotPlot(glioma_clean, features = Glioma_SCI_Markers_clean) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold.italic", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"), strip.text.x = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 6, t = 6)), strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + labs(x = "Marker Categories", y = "Cluster")

ggsave("05_Exploratory_DotPlot_Color1_Classic.pdf", plot = p_dot_base + scale_colour_gradient2(low = "#4575b4", mid = "white", high = "#d73027"), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color2_Magma.pdf", plot = p_dot_base + scale_colour_viridis_c(option = "magma"), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color3_Spectral.pdf", plot = p_dot_base + scale_colour_gradientn(colors = c("#00FF00", "black", "#FF0000")), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color4_Plasma.pdf", plot = p_dot_base + scale_colour_viridis_c(option = "plasma"), width = 20, height = 9)

cat("\n>>> 正在生成胶质瘤核心 Marker 的 空间投影图...\n")
core_features <- c("SOX2", "OLIG2", "GFAP", "CD44", "CD45", "CLDN5")
core_features <- core_features[core_features %in% rownames(glioma_clean)]
if(length(core_features) > 0){
  ggsave("06_Spatial_FeaturePlot_UMAP_Standard.pdf", plot = FeaturePlot(glioma_clean, features = core_features, ncol = 3, pt.size = 0.1, cols = c("lightgrey", "red"), reduction = "umap", order = TRUE), width = 16, height = 10)
  ggsave("06_Spatial_FeaturePlot_UMAP_Reshaped.pdf", plot = FeaturePlot(glioma_clean, features = core_features, ncol = 3, pt.size = 0.1, cols = c("lightgrey", "red"), reduction = "umap_alt", order = TRUE), width = 16, height = 10)
  ggsave("06_Spatial_FeaturePlot_tSNE_Standard.pdf", plot = FeaturePlot(glioma_clean, features = core_features, ncol = 3, pt.size = 0.1, cols = c("lightgrey", "red"), reduction = "tsne", order = TRUE), width = 16, height = 10)
  ggsave("06_Spatial_FeaturePlot_tSNE_Reshaped.pdf", plot = FeaturePlot(glioma_clean, features = core_features, ncol = 3, pt.size = 0.1, cols = c("lightgrey", "red"), reduction = "tsne_alt", order = TRUE), width = 16, height = 10)
}
saveRDS(glioma_clean, "Glioma_3A_ReadyForAnnotation.rds")


# ==============================================================================
# 胶质瘤单细胞 SOP: Step 3 (下) - 终极生物学注释与 SCI 级图谱生成
# ==============================================================================
cat(">>> [步骤 9] 加载待注释的纯净底盘...\n")
if (!exists("glioma_clean")) {
  glioma_clean <- readRDS("Glioma_3A_ReadyForAnnotation.rds")
}
target_res <- 0.7                                  
res_col <- paste0("RNA_snn_res.", target_res)  
Idents(glioma_clean) <- res_col

# ==============================================================================
# [步骤 10] 执行核心亚群的生物学身份加冕
# ==============================================================================
cat(">>> [步骤 10] 正在执行精准重映射，清空历史数字簇号...\n")
expert_assignments <- list(
  "NPC-like"    = c(17, 9, 16),
  "OPC-like"    = c(5, 13, 8),
  "AC-like"     = c(6, 12),
  "MES-like"    = c(1, 2, 4, 7),
  "TAM"         = c(10, 3, 0),
  "Lymphocytes" = c(11),
  "EC"          = c(14, 15)
)
raw_clusters <- as.character(glioma_clean@meta.data[[res_col]])
new_idents <- raw_clusters
for (cell_type in names(expert_assignments)) {
  clusters_to_match <- as.character(expert_assignments[[cell_type]])
  if (length(clusters_to_match) > 0) {
    new_idents[raw_clusters %in% clusters_to_match] <- cell_type
  }
}
cell_type_order <- names(expert_assignments) 
unassigned <- setdiff(unique(new_idents), cell_type_order)
final_levels <- c(cell_type_order, unassigned)
glioma_clean$CellType <- factor(new_idents, levels = final_levels)
Idents(glioma_clean) <- glioma_clean$CellType

# ==============================================================================
# [步骤 11] 绘制 [已注释] 空间降维图谱 (UMAP & tSNE 均 3 种风格)
# ==============================================================================
cat(">>> [步骤 11] 正在绘制并输出最终注释版的 3 种形式 UMAP 与 t-SNE 图谱...\n")

# 获取总细胞数量，并格式化加上千位分隔符 (例如: 189,647)
total_cells <- format(ncol(glioma_clean), big.mark = ",")
cell_count_label <- paste0("Number of cells = ", total_cells)

# 【已注释 UMAP 3 种风格】
ggsave("07a_Annotated_UMAP_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, cols = my_color_palette) + ggtitle("Style 1: Standard UMAP"), width = 9, height = 7)

ggsave("07a_Annotated_UMAP_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "umap_alt", cols = my_color_palette, pt.size = 0.3) + ggtitle("Style 2: Reshaped Minimal UMAP Annotated") + theme_sci_arrow + labs(x = "UMAP 1", y = "UMAP 2"), width = 8.5, height = 7)

# 🚨 图形 3 修改点：通过 labs 添加 subtitle，并利用 theme 将其固定在右上角
p_umap_style3 <- plot_dashed_repel(glioma_clean, "umap", "CellType", my_color_palette, "Style 3: Dashed Repel UMAP") + 
  labs(subtitle = cell_count_label) + 
  theme(plot.subtitle = element_text(hjust = 1, face = "bold", size = 12, color = "black"))
ggsave("07a_Annotated_UMAP_Style3_Dashed.pdf", plot = p_umap_style3, width = 9, height = 7)


# 【已注释 tSNE 3 种风格】
ggsave("07b_Annotated_tSNE_Style1_Standard.pdf", plot = DimPlot(glioma_clean, reduction = "tsne", label = TRUE, label.size = 5, repel = TRUE, cols = my_color_palette) + ggtitle("Style 1: Standard t-SNE"), width = 9, height = 7)

ggsave("07b_Annotated_tSNE_Style2_Minimal.pdf", plot = DimPlot(glioma_clean, reduction = "tsne_alt", cols = my_color_palette, pt.size = 0.3) + ggtitle("Style 2: Reshaped Minimal t-SNE Annotated") + theme_sci_arrow + labs(x = "t-SNE 1", y = "t-SNE 2"), width = 8.5, height = 7)

# 🚨 图形 3 修改点：通过 labs 添加 subtitle，并利用 theme 将其固定在右上角
p_tsne_style3 <- plot_dashed_repel(glioma_clean, "tsne", "CellType", my_color_palette, "Style 3: Dashed Repel t-SNE") + 
  labs(subtitle = cell_count_label) + 
  theme(plot.subtitle = element_text(hjust = 1, face = "bold", size = 12, color = "black"))
ggsave("07b_Annotated_tSNE_Style3_Dashed.pdf", plot = p_tsne_style3, width = 9, height = 7)

# ==============================================================================
# [步骤 12] 绘制细胞构成双向比例柱状图 (维持源代码格式)
# ==============================================================================
cat(">>> [步骤 12] 正在生成样本内 / 细胞类型内双向比例图...\n")
meta_df <- glioma_clean@meta.data

plot_data_type <- meta_df %>% group_by(CellType, SampleID) %>% summarise(Count = n(), .groups = 'drop') %>% group_by(CellType) %>% mutate(Proportion = Count / sum(Count))
p_bar_type <- ggplot(plot_data_type, aes(x = CellType, y = Proportion, fill = SampleID)) + geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.5) + scale_y_continuous(labels = percent_format()) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 11), axis.text.y = element_text(face = "bold", color = "black", size = 11), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), legend.position = "right", legend.title = element_text(face = "bold")) + labs(x = "Cell Types", y = "Fraction of Cells", fill = "Sample Source")
ggsave("08a_CellType_Sample_Distribution.pdf", plot = p_bar_type, width = 10, height = 6)

plot_data_sample <- meta_df %>% group_by(SampleID, CellType) %>% summarise(Count = n(), .groups = 'drop') %>% group_by(SampleID) %>% mutate(Proportion = Count / sum(Count))
p_bar_sample <- ggplot(plot_data_sample, aes(x = SampleID, y = Proportion, fill = CellType)) + geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.5) + scale_y_continuous(labels = percent_format()) + scale_fill_manual(values = my_color_palette) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 11), axis.text.y = element_text(face = "bold", color = "black", size = 11), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), legend.position = "right", legend.title = element_text(face = "bold")) + labs(x = "Sample ID", y = "Fraction of Cells", fill = "Cell Type")
ggsave("08b_Sample_CellType_Distribution.pdf", plot = p_bar_sample, width = 10, height = 6)


# ==============================================================================
# [步骤 13] 绘制顶级 SCI 阶梯状验证气泡图 (4色版本)
# ==============================================================================
cat(">>> [步骤 13] 正在构建阶梯状特征气泡图 (4种配色，并翻转Y轴让NPC排在最上方)...\n")
Glioma_SCI_Markers_Display <- list("NPC-like"=c("SOX4","DCX"), "OPC-like"=c("OLIG2","PDGFRA"), "AC-like"=c("GFAP","MLC1"), "MES-like"=c("CD44","VIM"), "TAM"=c("C1QA","AIF1"), "Lymphocytes"=c("CD3D","NKG7"), "EC"=c("PECAM1","VWF"))
Glioma_SCI_Markers_Display_clean <- lapply(Glioma_SCI_Markers_Display, function(x) x[x %in% rownames(glioma_clean)])
Glioma_SCI_Markers_Display_clean <- Glioma_SCI_Markers_Display_clean[lengths(Glioma_SCI_Markers_Display_clean) > 0]

# 🚨 增加 scale_y_discrete(limits = rev(levels(Idents(glioma_clean)))) 以翻转Y轴顺序
p_dotplot_anno_base <- DotPlot(glioma_clean, features = Glioma_SCI_Markers_Display_clean) + 
  scale_y_discrete(limits = rev(levels(Idents(glioma_clean)))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold.italic", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"), strip.text.x = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 6, t = 6)), strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.spacing = unit(0.1, "lines")) + 
  labs(x = "Core Marker Genes", y = "Annotated Cell Types") 

ggsave("09_Annotated_Marker_DotPlot_Color1_Classic.pdf", plot = p_dotplot_anno_base + scale_colour_gradient2(low = "#4575b4", mid = "white", high = "#d73027"), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color2_Magma.pdf", plot = p_dotplot_anno_base + scale_colour_viridis_c(option = "magma"), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color3_Spectral.pdf", plot = p_dotplot_anno_base + scale_colour_gradientn(colors = c("#00FF00", "black", "#FF0000")), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color4_Plasma.pdf", plot = p_dotplot_anno_base + scale_colour_viridis_c(option = "plasma"), width = 14, height = 7)





# ==============================================================================
# 附加图 1 & 2: 堆叠小提琴图与高变基因热图 (维持源代码格式)
# ==============================================================================
cat("\n>>> 附加图 1: 正在绘制原版堆叠小提琴图...\n")
top_1_markers <- sapply(Glioma_SCI_Markers_Display_clean, function(x) x[1])
top_1_markers <- unname(top_1_markers)

p_vln_stacked <- VlnPlot(glioma_clean, features = top_1_markers, stack = TRUE, flip = TRUE, cols = my_color_palette, pt.size = 0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"), axis.title.x = element_blank(), strip.text.y.right = element_text(angle = 0, face = "bold.italic", size = 11)) + ggtitle("Expression Density")
ggsave("10_Annotated_Stacked_VlnPlot_Original.pdf", plot = p_vln_stacked, width = 8, height = 12)

cat("\n>>> 附加图 2: 正在绘制原版数据驱动热图...\n")
cluster_markers <- FindAllMarkers(glioma_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
top5_hvgs_per_cluster <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
hvg_genes <- unique(top5_hvgs_per_cluster$gene)

glioma_clean <- ScaleData(glioma_clean, features = hvg_genes, verbose = FALSE)
set.seed(20260808)
glioma_downsampled <- subset(glioma_clean, downsample = 150)

p_heatmap <- DoHeatmap(glioma_downsampled, features = hvg_genes, group.colors = my_color_palette, size = 4, angle = 90) + scale_fill_gradientn(colors = c("#4575b4", "white", "#d73027")) + theme(axis.text.y = element_text(face = "bold.italic", size = 10))
ggsave("11_Annotated_HVG_Heatmap_Original.pdf", plot = p_heatmap, width = 16, height = 12)



# ==============================================================================
# 附加图 3: 绘制每个细胞亚型的绝对数量统计柱状图
# ==============================================================================
cat("\n>>> 附加图 3: 正在统计并绘制各细胞亚型数量柱状图...\n")

# 1. 统计各细胞类型数量
plot_data_count <- glioma_clean@meta.data %>% 
  group_by(CellType) %>% 
  summarise(Count = n(), .groups = 'drop')

# 2. 绘制带有数值标签的柱状图
p_cell_count <- ggplot(plot_data_count, aes(x = CellType, y = Count, fill = CellType)) + 
  geom_bar(stat = "identity", color = "black", width = 0.7, linewidth = 0.5) + 
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 4, color = "black") + 
  scale_fill_manual(values = my_color_palette) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # 顶部留出 10% 空间给标签
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 11), 
    axis.text.y = element_text(face = "bold", color = "black", size = 11), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none" # 横轴已有类别名字，隐藏图例使画面更清爽
  ) + 
  labs(x = "Annotated Cell Types", y = "Number of Cells", title = "Cell Count per Cell Type")

# 3. 保存图片
ggsave("12_CellType_Counts_BarPlot.pdf", plot = p_cell_count, width = 8, height = 6)
cat(">>> ✅ 细胞亚型数量柱状图已生成: 12_CellType_Counts_BarPlot.pdf\n")


# ==============================================================================
# [步骤 14] 固化终极 RDS 对象
# ==============================================================================
cat(">>> [步骤 14] 正在保存全副武装的终局对象...\n")
saveRDS(glioma_clean, file = "Glioma_3B_Fully_Annotated.rds")

cat("\n=================================================================\n")
cat(">>> 🎉 全线清理完毕！\n")
cat(">>> ✅ Style 4 (边框图) 已彻底移除；\n")
cat(">>> ✅ Style 2 已 100% 严格启用新随机种子/参数计算出的全新拓扑坐标；\n")
cat(">>> ✅ Style 2 的视觉完美继承无标签、半透明散点与箭头坐标轴风格；\n")
cat(">>> ✅ 其余所有质控分析图回归最严谨的原始分析 SOP 格式！\n")
cat("=================================================================\n")
















# ==============================================================================
# 附加脚本：使用 OLIG2 重新生成探索性气泡图与阶梯状注释气泡图 (8图连发)
# ==============================================================================

library(Seurat)
library(ggplot2)

# 1. 设置工作目录 (与您源代码保持一致)
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/3.注释新图形"
setwd(work_dir)

# 2. 读取已经注释好的终局对象
cat(">>> 正在加载已注释对象 Glioma_3B_Fully_Annotated.rds ...\n")
if (!file.exists("Glioma_3B_Fully_Annotated.rds")) stop("未找到注释对象，请确认路径！")
glioma_clean <- readRDS("Glioma_3B_Fully_Annotated.rds")

# 确认默认的 Idents 为注释后的 CellType
Idents(glioma_clean) <- "CellType"

# ==============================================================================
# 模块 A：重新生成 [未注释] 探索性气泡图 (4 种配色)
# 对应原始代码的 [步骤 8]
# ==============================================================================
cat(">>> 正在重新生成 探索性 Marker 气泡图 (已将 OLIG1 替换为 OLIG2)...\n")

# 🚨 将 OPC-like 的第一个 marker 从 OLIG1 修改为 OLIG2
Glioma_SCI_Markers <- list(
  "NPC-like"    = c("SOX4", "SOX2", "SOX11", "DCX", "DLL3", "CD24"),
  "OPC-like"    = c("OLIG2", "TNR", "ALCAM", "OMG", "PLLP", "PDGFRA"), 
  "AC-like"     = c("GFAP", "S100B", "HOPX", "SLC1A3", "MLC1"),
  "MES-like"    = c("CD44", "CHI3L1", "ANXA1", "ANXA2", "VIM"),
  "TAM"         = c("C1QA", "C1QB", "C1QC", "TYROBP", "FCER1G", "LST1", "AIF1"),
  "Lymphocytes" = c("CD3D", "CD3E", "TRAC", "NKG7", "KLRD1", "CCL5"),
  "EC"          = c("CLDN5", "VWF", "PECAM1", "CDH5", "FLT1", "PDGFRB")
)

Glioma_SCI_Markers_clean <- lapply(Glioma_SCI_Markers, function(x) x[x %in% rownames(glioma_clean)])
Glioma_SCI_Markers_clean <- Glioma_SCI_Markers_clean[lengths(Glioma_SCI_Markers_clean) > 0]

# 临时将 Idents 切回数字簇，以还原“未注释状态”的气泡图外观
# 注意：这里使用您源代码中的分辨率 0.6
res_col <- "RNA_snn_res.0.6"
Idents(glioma_clean) <- res_col

p_dot_base <- DotPlot(glioma_clean, features = Glioma_SCI_Markers_clean) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold.italic", color = "black"), 
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
    strip.text.x = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 6, t = 6)), 
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) + labs(x = "Marker Categories", y = "Cluster")

ggsave("05_Exploratory_DotPlot_Color1_Classic.pdf", plot = p_dot_base + scale_colour_gradient2(low = "#4575b4", mid = "white", high = "#d73027"), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color2_Magma.pdf", plot = p_dot_base + scale_colour_viridis_c(option = "magma"), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color3_Spectral.pdf", plot = p_dot_base + scale_colour_gradientn(colors = c("#00FF00", "black", "#FF0000")), width = 20, height = 9)
ggsave("05_Exploratory_DotPlot_Color4_Plasma.pdf", plot = p_dot_base + scale_colour_viridis_c(option = "plasma"), width = 20, height = 9)


# ==============================================================================
# 模块 B：重新生成 [已注释] SCI 阶梯状气泡图 (4 种配色)
# 对应原始代码的 [步骤 13]
# ==============================================================================
cat(">>> 正在重新生成 阶梯状验证气泡图 (已将 OLIG1 替换为 OLIG2)...\n")

# 切回正式的细胞类型注释
Idents(glioma_clean) <- "CellType"

# 🚨 将 OPC-like 的第一个核心展示 marker 从 OLIG1 修改为 OLIG2
Glioma_SCI_Markers_Display <- list(
  "NPC-like"    = c("SOX4","DCX"), 
  "OPC-like"    = c("OLIG2","PDGFRA"), 
  "AC-like"     = c("GFAP","MLC1"), 
  "MES-like"    = c("CD44","VIM"), 
  "TAM"         = c("C1QA","AIF1"), 
  "Lymphocytes" = c("CD3D","NKG7"), 
  "EC"          = c("PECAM1","VWF")
)

Glioma_SCI_Markers_Display_clean <- lapply(Glioma_SCI_Markers_Display, function(x) x[x %in% rownames(glioma_clean)])
Glioma_SCI_Markers_Display_clean <- Glioma_SCI_Markers_Display_clean[lengths(Glioma_SCI_Markers_Display_clean) > 0]

p_dotplot_anno_base <- DotPlot(glioma_clean, features = Glioma_SCI_Markers_Display_clean) + 
  scale_y_discrete(limits = rev(levels(Idents(glioma_clean)))) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold.italic", color = "black"), 
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
    strip.text.x = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 6, t = 6)), 
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
    panel.spacing = unit(0.1, "lines")
  ) + labs(x = "Core Marker Genes", y = "Annotated Cell Types") 

ggsave("09_Annotated_Marker_DotPlot_Color1_Classic.pdf", plot = p_dotplot_anno_base + scale_colour_gradient2(low = "#4575b4", mid = "white", high = "#d73027"), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color2_Magma.pdf", plot = p_dotplot_anno_base + scale_colour_viridis_c(option = "magma"), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color3_Spectral.pdf", plot = p_dotplot_anno_base + scale_colour_gradientn(colors = c("#00FF00", "black", "#FF0000")), width = 14, height = 7)
ggsave("09_Annotated_Marker_DotPlot_Color4_Plasma.pdf", plot = p_dotplot_anno_base + scale_colour_viridis_c(option = "plasma"), width = 14, height = 7)

cat(">>> 🎉 8 张气泡图已全部成功覆盖重写！请前往文件夹查看最新生成的 OLIG2 版本。\n")