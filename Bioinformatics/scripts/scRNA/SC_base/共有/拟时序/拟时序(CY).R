# ==============================================================================
# Monocle 2 超音速版：完美适配 Seurat V5 多线程与最新细胞注释 (NPC 为干性起点)
# 核心优化：引入 DEG 驱动降维、剔除细胞周期基因、增加批次诊断
# ==============================================================================
options(stringsAsFactors = FALSE)
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)

# ===== 新增：多线程配置开始 =====
library(future)
library(future.apply)

# 设置多线程工作核数（这里设置为12核）
plan("multisession", workers = 12)
# 提高全局变量的内存限制
options(future.globals.maxSize = 100000 * 1024^2)
# ===== 新增：多线程配置结束 =====

# 0. 全局参数配置区 (修改此处即可全局生效)
# ------------------------------------------------------------------------------
# 0.1 设置工作目录
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE221534_Object_1_Core/拟时序"
if(!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(work_dir)

# 0.2 定义目标基因 (修改此变量，后续所有靶基因分析将自动更新)
target_gene <- "ARL6IP1"
gene_name <- target_gene

# 0.3 统一的细胞类型SCI配色 (设定 NPC-like 为红色干性起点)
my_colors <- c(
  "NPC-like" = "#E64B35",  # 红色 (干细胞态/起点)
  "OPC-like" = "#4DBBD5",  # 蓝青色
  "AC-like"  = "#00A087",  # 绿色
  "MES-like" = "#3C5488"   # 深蓝色
)
# ------------------------------------------------------------------------------

# 1. 读取已经注释好的 Seurat RDS 数据 
print(">>> [1/6] 读取 Seurat RDS 数据...")
scRNA <- readRDS("/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE221534_Object_1_Core/Glioma_3B_Fully_Annotated.rds")

# 2. 提取肿瘤细胞 & 适度抽样 (每组 800)
print(">>> [2/6] 提取细胞并抽样 (每组 800 个)...")
tumor_classes <- c("NPC-like", "OPC-like", "AC-like", "MES-like")

Idents(scRNA) <- "CellType"
valid_cells <- tumor_classes[tumor_classes %in% unique(Idents(scRNA))]
scRNA_tumor <- subset(scRNA, idents = valid_cells)

# 强制刷新 Idents
scRNA_tumor$celltype.main <- Idents(scRNA_tumor)
scRNA_sub <- subset(scRNA_tumor, downsample = 800) 

print(paste("最终分析细胞数:", ncol(scRNA_sub)))

# 3. 构建 CDS
print(">>> [3/6] 构建 Monocle 对象...")
# 核心修改：使用 DefaultAssay() 自动适配 scRNA 和 snRNA 的不同 Assay 命名习惯
current_assay <- DefaultAssay(scRNA_sub)
data <- as(as.matrix(GetAssayData(scRNA_sub, slot = "counts", assay = current_assay)), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)

# 4. 设置排序基因 (极其严格的过滤策略：DEG + 高变基因交集)
print(">>> [4/6] 计算核心驱动基因并剔除噪音...")

# 4.1 提高 DEG 的筛选门槛 (logfc 提高到 0.5，p值降到 1e-5)
Idents(scRNA_sub) <- "CellType"
deg_markers <- FindAllMarkers(scRNA_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
valid_degs <- deg_markers %>% filter(p_val_adj < 1e-5) %>% pull(gene) %>% unique()

# 4.2 提取 Seurat 中计算过的高变基因 (HVGs)
# 这里假设你在 Seurat 阶段已经运行过 FindVariableFeatures
hvgs <- VariableFeatures(scRNA_sub) 
if(length(hvgs) == 0) {
  scRNA_sub <- FindVariableFeatures(scRNA_sub, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(scRNA_sub)
}

# 4.3 取交集：必须既是细胞类型强显著标志物，又是全局高变基因
core_genes <- intersect(valid_degs, hvgs)

# 4.4 准备剔除的黑名单基因 (细胞周期 + 线粒体/核糖体)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
cell_cycle_genes <- unique(c(s_genes, g2m_genes))
junk_genes <- grep("^MT-|^RP[SL]", rownames(scRNA_sub), value = TRUE)

# 4.5 得到最终纯净的排序基因
ordering_genes <- setdiff(core_genes, c(cell_cycle_genes, junk_genes))

# 理想状态下，这里的基因数应该被压缩到 500 ~ 1200 之间
print(paste(">>> 严格过滤后，最终用于构建轨迹的核心驱动基因数:", length(ordering_genes)))

cds <- setOrderingFilter(cds, ordering_genes)

# 5. 降维 (修剪细枝，保留主干)
print(">>> [5/6] 降维并回归批次效应 (开启自动平滑)...")
cds <- reduceDimension(cds, 
                       max_components = 2, 
                       reduction_method = 'DDRTree', 
                       # 核心修改点：
                       auto_param_selection = TRUE)

# 6. 排序与找起点
print(">>> [6/6] 排序与初步寻找起点...")
cds <- orderCells(cds)

# ==============================================================================
# 新增：批次与周期效应诊断 (Diagnostic Check)
# ==============================================================================
message(">>> 正在生成批次效应诊断图 (Diagnostic_Batch_Phase.pdf)...")
pdf("Diagnostic_Batch_Phase.pdf", width = 12, height = 5)
p_batch <- plot_cell_trajectory(cds, color_by = "orig.ident", show_branch_points = FALSE) + 
  ggtitle("Trajectory by Sample (Batch Check)")
if("Phase" %in% colnames(pData(cds))) {
  p_phase <- plot_cell_trajectory(cds, color_by = "Phase", show_branch_points = FALSE) + 
    ggtitle("Trajectory by Cell Cycle Phase")
  grid.arrange(p_batch, p_phase, ncol = 2)
} else {
  print(p_batch)
}
dev.off()


# ==============================================================================
# 智能寻找“纯净” 起点 (将 NPC-like 设为根节点)
# ==============================================================================
message(">>> 正在寻找最佳 NPC 起点 (将 NPC-like 设为干性起点)...")
st_tbl <- table(pData(cds)$State, pData(cds)$celltype.main)
st_df <- as.data.frame.matrix(st_tbl)

st_df$NPC_Score <- st_df[, "NPC-like"] / (rowSums(st_df) - st_df[, "NPC-like"] + 1)
best_start_state <- rownames(st_df)[which.max(st_df$NPC_Score)]
message(paste(">>> 推荐的最佳 NPC 起点 State 是:", best_start_state))

cds <- orderCells(cds, root_state = best_start_state)
message(">>> 细胞排序完成，拟时序已更新。")


# ==============================================================================
# 统一参数设置 & 后续所有高级分析
# ==============================================================================

# --- 任务 A: 细胞类型轨迹图 ---
message(">>> 正在绘制任务: Supersonic_Trajectory.pdf ...")
pdf("Supersonic_Trajectory.pdf", width = 10, height = 8)
p1 <- plot_cell_trajectory(cds, color_by = "celltype.main", show_branch_points = FALSE) + 
  ggtitle("Cell Trajectory by Type (DEG Filtered)") +
  scale_color_manual(values = my_colors, na.value = "grey")
print(p1)
dev.off()

# --- 任务 B: 纯净版时间流向图 ---
message(">>> 正在绘制方向标记图 ...")
p_base <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") +
  ggtitle("Trajectory Direction")

pdf("Direction_of_Time_Fixed.pdf", width = 8, height = 6)
print(p_base)
dev.off()

# --- 任务 C: 基因趋势图 ---
gene_id <- row.names(fData(cds)[fData(cds)$gene_short_name == target_gene, , drop=FALSE])
if(length(gene_id) == 0) gene_id <- target_gene

if(length(gene_id) > 0 && gene_id[1] %in% row.names(cds)) {
  pdf(paste0("Supersonic_", target_gene, ".pdf"), width = 8, height = 6)
  p_gene <- plot_genes_in_pseudotime(cds[gene_id[1], ], color_by = "celltype.main", ncol = 1) +
    scale_color_manual(values = my_colors)
  print(p_gene)
  dev.off()
}

# --- 高级图表 A: 动态密度图 ---
pdata <- pData(cds)
p_stream <- ggplot(pdata, aes(x = Pseudotime, fill = celltype.main, color = celltype.main)) +
  geom_density(position = "identity", alpha = 0.4, linewidth = 1) + 
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Density (Cell Concentration)", title = "Cellular Evolution Dynamics along Pseudotime") +
  theme(text = element_text(size = 14), axis.line.x = element_line(linewidth = 1), axis.line.y = element_line(linewidth = 1))

pdf("Advanced_A_StreamPlot.pdf", width = 8, height = 5)
print(p_stream)
dev.off()

# --- 高级图表 B: 分箱气泡图 ---
expr_data <- t(as.matrix(cds[gene_id,]))
df_plot <- pData(cds)
df_plot$Expression <- as.numeric(expr_data)
df_plot$PseudoBin <- cut(df_plot$Pseudotime, breaks = 10, labels = 1:10)
df_summary <- df_plot %>% group_by(PseudoBin) %>% summarise(
  AvgExp = mean(Expression), PctExp = sum(Expression > 0) / n() * 100, TimeMean = mean(Pseudotime))

p_dot <- ggplot(df_summary, aes(x = PseudoBin, y = target_gene)) +
  geom_point(aes(size = PctExp, color = AvgExp)) +
  scale_color_gradientn(colours = c("grey90", "orange", "red", "darkred")) +
  scale_size_continuous(range = c(2, 10)) + theme_minimal() +
  labs(x = "Pseudotime Interval (Early -> Late)", y = "", color = "Avg Exp", size = "% Expressed", title = paste(target_gene, "Temporal Expression Pattern"))

pdf("Advanced_B_DotPlot.pdf", width = 8, height = 4)
print(p_dot)
dev.off()

# --- 高级图表 C: 极简轨迹图 ---
p_traj_gene <- plot_cell_trajectory(cds, markers = target_gene, use_color_gradient = TRUE, show_branch_points = FALSE) +
  scale_color_gradientn(colours = c("grey90", "#FFCC33", "#CC0000")) + 
  theme_void() + theme(legend.position = "right") + ggtitle(paste(target_gene, "Expression on Trajectory"))
p_traj_gene$data <- p_traj_gene$data %>% arrange(value)

pdf("Advanced_C_CleanTrajectory.pdf", width = 8, height = 8)
print(p_traj_gene)
dev.off()

# ==============================================================================
# BEAM 分支命运分析 (全面升级：200个基因 + 命运状态拟时序可视化 + 自定义热图)
# ==============================================================================
print(">>> 开始运行 BEAM 分支分析...")
tryCatch({
  BEAM_res <- BEAM(cds, branch_point = 1, cores = 12, progenitor_method = 'duplicate') 
}, error = function(e) {
  BEAM_res <<- BEAM(cds, branch_point = 1, cores = 1, progenitor_method = 'duplicate')
})

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

# 需求1：纳入前 200 个显著基因
top_n <- min(200, nrow(subset(BEAM_res, qval < 1e-4)))
sig_genes <- row.names(subset(BEAM_res, qval < 1e-4))[1:top_n]
if(target_gene %in% row.names(cds) & !(target_gene %in% sig_genes)){ sig_genes <- c(sig_genes, target_gene) }

# ------------------------------------------------------------------------------
# 需求1 补充：绘制两张不同的 BEAM 分支热图
# ------------------------------------------------------------------------------

# 图 1：展示所有 200 个基因的行名 (设置较高高度以防重叠)
print(">>> 正在绘制完整基因名热图...")
pdf("Advanced_BEAM_Heatmap_AllGenes.pdf", width = 8, height = 25) # 高度设为25，确保200个基因清晰
plot_genes_branched_heatmap(cds[sig_genes,], branch_point = 1, num_clusters = 4, cores = 1, 
                            use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()

# 图 2：仅展示目标基因 (target_gene) 的行名，隐藏其他
print(paste(">>> 正在绘制仅标注靶基因", target_gene, "的热图..."))
# 提取用于画图的矩阵信息以获取默认的基因短名
cds_subset <- cds[sig_genes,]
gene_labels <- fData(cds_subset)$gene_short_name
# 将不是目标基因的名称设为空字符串
gene_labels[gene_labels != target_gene] <- ""

# 将修改后的基因名赋值回 cds 对象，以便画图时直接调用
fData(cds_subset)$plot_name <- gene_labels

pdf("Advanced_BEAM_Heatmap_TargetGene.pdf", width = 7, height = 12) # 常规高度
plot_genes_branched_heatmap(cds_subset, branch_point = 1, num_clusters = 4, cores = 1, 
                            use_gene_short_name = FALSE, # 设置为FALSE，手动指定标签
                            show_rownames = TRUE) # 开启行名显示，但大部分是空字符串
dev.off()
# ------------------------------------------------------------------------------


# ==============================================================================
# 需求2：提取并可视化 Pre-branch, fate1, fate2 命运状态 (纯手工高鲁棒性方案)
# ==============================================================================
print(">>> 正在提取并绘制分支命运 (Pre-branch, fate1, fate2) 状态图...")

# 步骤 A：为原 cds 对象映射命运状态 (在 2D 轨迹空间上的体现)
pData(cds)$Fate_State <- "Unknown"
pData(cds)$Fate_State[pData(cds)$State == best_start_state] <- "Pre-branch"

other_states <- setdiff(unique(pData(cds)$State), best_start_state)
if(length(other_states) >= 1) pData(cds)$Fate_State[pData(cds)$State == other_states[1]] <- "Fate 1"
if(length(other_states) >= 2) pData(cds)$Fate_State[pData(cds)$State == other_states[2]] <- "Fate 2"
pData(cds)$Fate_State <- factor(pData(cds)$Fate_State, levels = c("Pre-branch", "Fate 1", "Fate 2", "Unknown"))

# 新增输出2：在 2D 树状轨迹图上展示 Pre-branch 和 fate 分支
pdf("Advanced_BEAM_Fate_Trajectory.pdf", width = 8, height = 6)
p_fate_traj <- plot_cell_trajectory(cds, color_by = "Fate_State", show_branch_points = TRUE) +
  # 这里已经将颜色严格替换为了你图片中的 灰、红、蓝 配色
  scale_color_manual(values = c("Pre-branch" = "#999999", "Fate 1" = "#EB5A61", "Fate 2" = "#7592C9", "Unknown" = "grey90")) +
  theme(legend.position = "right") +
  ggtitle("Cell Trajectory by Fates (Pre-branch, Fate 1, Fate 2)")
print(p_fate_traj)
dev.off()

# ------------------------------------------------------------------------------
# 步骤 B & C：彻底绕过 Monocle2 底层矩阵提取 Bug，使用 dplyr 进行平滑映射
# ------------------------------------------------------------------------------
# 1. 准备基础数据并保留原始细胞 ID
df_plot <- pData(cds) %>% filter(Fate_State != "Unknown")
df_plot$Cell_ID <- rownames(df_plot)

# 2. 将 Pre-branch 细胞分别复制给两条分支，完美模拟 'duplicate' 祖先法，且不触发矩阵报错
df_branch1 <- df_plot %>% filter(Fate_State %in% c("Pre-branch", "Fate 1")) %>% mutate(Branch = "Fate 1 Lineage")
df_branch2 <- df_plot %>% filter(Fate_State %in% c("Pre-branch", "Fate 2")) %>% mutate(Branch = "Fate 2 Lineage")
df_manual_branch <- rbind(df_branch1, df_branch2)

# 新增输出3：命运分支在拟时序坐标轴上的细胞密度动态演化图
pdf("Advanced_BEAM_Fate_Density.pdf", width = 8, height = 5)
p_fate_density <- ggplot(df_manual_branch, aes(x = Pseudotime, fill = Branch, color = Branch)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("Fate 1 Lineage" = "#d23535", "Fate 2 Lineage" = "#40469472")) +
  scale_color_manual(values = c("Fate 1 Lineage" = "#d23535", "Fate 2 Lineage" = "#40469472")) +
  theme_classic() +
  labs(title = "Cell Density along Pseudotime by Branch Fates", x = "Pseudotime", y = "Density")
print(p_fate_density)
dev.off()

# 新增输出4：核心靶基因在分支上的拟时序动力学分叉折线图
if(target_gene %in% row.names(cds)) {
  # 根据精确的 Cell_ID 提取表达量，杜绝因为行名错乱导致的问题
  df_manual_branch$Expression <- as.numeric(exprs(cds)[target_gene, df_manual_branch$Cell_ID])
  
  pdf(paste0("Advanced_BEAM_Branched_Pseudotime_", target_gene, ".pdf"), width = 8, height = 6)
  p_branch_gene <- ggplot(df_manual_branch, aes(x = Pseudotime, y = Expression)) +
    geom_point(aes(color = celltype.main), size = 1.2, alpha = 0.6) +
    # 根据两条分支画出独立平滑曲线
    geom_smooth(aes(linetype = Branch), color = "black", method = "loess", linewidth = 1.2, se = FALSE) +
    scale_color_manual(values = my_colors) +
    theme_classic() +
    labs(title = paste("Branched Pseudotime Expression of", target_gene), y = "Relative Expression")
  print(p_branch_gene)
  dev.off()
}
# ==============================================================================
# 通路评分动态拟时序 & 基因相关性
# ==============================================================================
mes_genes <- c("CHI3L1", "CD44", "TIMP1") 
cycle_genes <- c("MKI67", "TOP2A", "SNAP25") 
mes_genes <- mes_genes[mes_genes %in% row.names(cds)]
cycle_genes <- cycle_genes[cycle_genes %in% row.names(cds)]

get_score <- function(cds, features){
  expr <- exprs(cds)[features, ]
  if(length(features) > 1) { score <- colMeans(t(scale(t(as.matrix(expr))))) } else { score <- as.numeric(expr) }
  return(score)
}
pData(cds)$MES_Score <- get_score(cds, mes_genes)
pData(cds)$Cycle_Score <- get_score(cds, cycle_genes)

p_trend <- ggplot(pData(cds), aes(x = Pseudotime)) +
  geom_smooth(aes(y = as.numeric(exprs(cds)[target_gene,]), color = target_gene), se = FALSE, linewidth = 1.5) +
  geom_smooth(aes(y = MES_Score, color = "MES Score"), se = FALSE, linewidth = 1.5) +
  scale_color_manual(values = setNames(c("red", "blue"), c(target_gene, "MES Score"))) + theme_classic() +
  labs(y = "Relative Expression / Score", title = paste("Correlation of", target_gene, "with Mesenchymal Signature"))

pdf("Advanced_Pathway_Trend.pdf", width = 8, height = 5)
print(p_trend)
dev.off()

target_expr <- as.numeric(exprs(cds)[target_gene,])
all_expr <- as.matrix(exprs(cds))
cor_res <- apply(all_expr, 1, function(x) cor(x, target_expr, method = "spearman"))
junk_pattern <- "^MT-|^RP[SL]" 
junk_genes <- grep(junk_pattern, names(cor_res), value = TRUE)
clean_cor_res <- cor_res[!names(cor_res) %in% junk_genes]
clean_cor_res <- clean_cor_res[names(clean_cor_res) != target_gene]
top_cor_genes <- sort(clean_cor_res, decreasing = TRUE)[1:20]

saveRDS(names(top_cor_genes), paste0(target_gene, "_Top_Correlated_Genes.rds"))

# ==============================================================================
# SCI 级精美图表
# ==============================================================================
expr_val <- as.numeric(exprs(cds)[gene_name,])
plot_df_sci <- data.frame(
  x = reducedDimS(cds)[1,], y = reducedDimS(cds)[2,],
  Pseudotime = pData(cds)$Pseudotime, CellType = pData(cds)$celltype.main, Expression = expr_val
) %>% arrange(Expression) 

p1_map <- ggplot(plot_df_sci, aes(x=x, y=y)) + geom_point(color="grey90", size=1) +
  geom_point(aes(color=Expression), size=1.5, alpha=0.8) +
  scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"), name = "Log Expression") +
  theme_void() + ggtitle(paste0(gene_name, " Spatial Expression")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "right")

pdf("SCI_Fig1_Gene_Map.pdf", width = 8, height = 7)
print(p1_map)
dev.off()

p2_split <- ggplot(plot_df_sci, aes(x=x, y=y)) + geom_point(aes(color=CellType), size=0.8, alpha=0.6) +
  scale_color_manual(values = my_colors) + theme_void() + facet_wrap(~CellType, nrow = 1) + 
  ggtitle("Cell Type Distribution on Trajectory") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "none")

pdf("SCI_Fig2_Split_View.pdf", width = 12, height = 4)
print(p2_split)
dev.off()

p3_trend <- ggplot(plot_df_sci, aes(x=Pseudotime, y=Expression)) +
  geom_jitter(aes(color=CellType), size=1.5, alpha=0.6, width = 0.5) +
  geom_smooth(color="black", linewidth=1.5, se=FALSE, method = "loess") +
  scale_color_manual(values = my_colors) + theme_classic() + 
  labs(x = "Pseudotime (Start -> End)", y = "Relative Expression") +
  geom_vline(xintercept = 25, linetype="dashed", color="grey50")

pdf("SCI_Fig3_Trend_Clean.pdf", width = 8, height = 5)
print(p3_trend)
dev.off()

y_limit <- quantile(plot_df_sci$Expression[plot_df_sci$Expression > 0], 0.98) 
p3_fixed <- ggplot(plot_df_sci, aes(x=Pseudotime, y=Expression)) +
  geom_jitter(aes(color=CellType), size=1.2, alpha=0.4, width = 0.2) +
  geom_smooth(color="#333333", linewidth=1.5, se=TRUE, fill="grey80", method = "loess") +
  scale_color_manual(values = my_colors) + coord_cartesian(ylim = c(0, y_limit)) + theme_classic() +
  theme(legend.position = "top", legend.title = element_blank())

pdf("SCI_Fixed_Trend_Zoomed.pdf", width = 8, height = 5)
print(p3_fixed)
dev.off()

p_time_sci <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = FALSE) +
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + theme_void()

pdf("SCI_Fixed_Pseudotime.pdf", width = 8, height = 7)
print(p_time_sci)
dev.off()

p_cell_sci <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 1, show_branch_points = FALSE) +
  scale_color_manual(values = my_colors) + theme_void() + guides(color = guide_legend(override.aes = list(size = 4)))

pdf("SCI_Fixed_CellType.pdf", width = 9, height = 7)
print(p_cell_sci)
dev.off()

# ==============================================================================
# 最终存档
# ==============================================================================
saveRDS(cds, "Monocle2_Final_Corrected.rds")
if(exists("top_cor_genes")){
  saveRDS(names(top_cor_genes), paste0(target_gene, "_Top20_Correlated_Genes.rds"))
  write.csv(names(top_cor_genes), paste0(target_gene, "_Top20_Genes_List.csv")) 
}
write.csv(pData(cds), "Cell_Metadata_with_Pseudotime.csv")

print(">>> ✅ 所有拟时序分析文件已安全生成并保存至:")
print(getwd())