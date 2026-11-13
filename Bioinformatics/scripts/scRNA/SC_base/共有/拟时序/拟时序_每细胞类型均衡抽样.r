# ==============================================================================
# 🌟 全新融合数据驱动流程 - 阶段一：全局客观拟时序主干构建
# 骨架规范：YKT6拟时序hlb重构.r | 计算内核：拟时序.r
# ==============================================================================

 # 随机抽取基因1200改为随机从每个细胞类型抽取300，均衡细胞类型，
# ------------------------------------------------------------------------------
# 0. 加载必要 R 包与多线程配置
# ------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)

# 开启多核与解除内存限制
plan("multisession", workers = 12)
options(future.globals.maxSize = 100000 * 1024^2)

print(">>>>>>>>>> 开启全新阶段一：客观拟时序主干构建 <<<<<<<<<<")

# ------------------------------------------------------------------------------
# 1. 设置工作目录与 SCI 统一配色
# ------------------------------------------------------------------------------
work_dir <- "/home/shazhi11/Lin/Glioma/Data/20260925/GSE182109(scRNA)(Object0_r+nd)/8_拟时序_重出"
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

my_colors <- c("NPC-like" = "#E64B35", "OPC-like" = "#4DBBD5", "AC-like" = "#00A087", "MES-like" = "#3C5488")

# ------------------------------------------------------------------------------
# 2. 读取数据与提取肿瘤细胞 (采用 "拟时序.r" 的 downsample 800 机制)
# ------------------------------------------------------------------------------
print(">>> 读取 Seurat RDS 数据并提取肿瘤细胞...")
scRNA <- readRDS("/home/shazhi11/Lin/Glioma/Data/20260925/GSE182109(scRNA)(Object0_r+nd)/Glioma_3B_Fully_Annotated.rds")

tumor_classes <- c("NPC-like", "OPC-like", "AC-like", "MES-like")
Idents(scRNA) <- "CellType"
valid_cells <- tumor_classes[tumor_classes %in% unique(Idents(scRNA))]
scRNA_tumor <- subset(scRNA, idents = valid_cells)

# 强制刷新 Idents
Idents(scRNA_tumor) <- "CellType"
scRNA_tumor$celltype.main <- Idents(scRNA_tumor)

# ==============================================================================
# 每种细胞类型随机抽取 300 个细胞：4 × 300 = 1200
# ==============================================================================
set.seed(2026)

cells_per_type <- 300

cell_counts <- table(Idents(scRNA_tumor))
print(">>> 抽样前各细胞类型数量：")
print(cell_counts)

# 检查每种细胞是否至少有300个
if(any(cell_counts[valid_cells] < cells_per_type)) {
  stop(
    paste0(
      "以下细胞类型不足 ", cells_per_type, " 个：",
      paste(names(cell_counts[valid_cells][cell_counts[valid_cells] < cells_per_type]), collapse = ", ")
    )
  )
}

selected_cells <- unlist(lapply(valid_cells, function(ct) {
  sample(WhichCells(scRNA_tumor, idents = ct), cells_per_type)
}))

scRNA_sub <- subset(scRNA_tumor, cells = selected_cells)

print(">>> 抽样后各细胞类型数量：")
print(table(Idents(scRNA_sub)))
print(paste("最终分析细胞数:", ncol(scRNA_sub)))

# ------------------------------------------------------------------------------
# 3. 严谨筛选用于建树的排序基因并构建 CDS (CDS 离散度与找基因规则采用 "拟时序.r")
# ------------------------------------------------------------------------------
print(">>> 构建 CellDataSet 对象并筛选建树核心基因...")
current_assay <- DefaultAssay(scRNA_sub)
data <- as(as.matrix(GetAssayData(scRNA_sub, slot = "counts", assay = current_assay)), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_sub@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(data), row.names = row.names(data)))

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# 仅进行 SizeFactor 估计，不进行 Dispersion 估计（采用 拟时序.r 逻辑）
cds <- estimateSizeFactors(cds)

# 找 DEG 阈值 logfc = 0.2, p < 1e-5
Idents(scRNA_sub) <- "CellType"
deg_markers <- FindAllMarkers(scRNA_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
valid_degs <- deg_markers %>% filter(p_val_adj < 1e-5) %>% pull(gene) %>% unique()

# 找 2000 个高变基因 (HVG)
hvgs <- VariableFeatures(scRNA_sub) 
if(length(hvgs) == 0) {
  scRNA_sub <- FindVariableFeatures(scRNA_sub, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(scRNA_sub)
}
core_genes <- intersect(valid_degs, hvgs)

# 准备剔除的黑名单基因 (只剔除细胞周期 + 线粒体/核糖体)
cell_cycle_genes <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))
junk_genes <- grep("^MT-|^RP[SL]", rownames(scRNA_sub), value = TRUE)
ordering_genes <- setdiff(core_genes, c(cell_cycle_genes, junk_genes))

cds <- setOrderingFilter(cds, ordering_genes)

# ------------------------------------------------------------------------------
# 4. DDRTree 降维建树 (包含 SampleID 批次矫正，采用 "重构.r" 逻辑)
# ------------------------------------------------------------------------------
print(">>> DDRTree 降维建树中 (请耐心等待)...")
# 强制开启 residualModelFormulaStr 并关闭自动参数
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', 
                       residualModelFormulaStr = "~SampleID", auto_param_selection = FALSE)
cds <- orderCells(cds)

# ------------------------------------------------------------------------------
# 5. 智能定位客观起点 (将 NPC-like 富集度最高的端点设为根节点) 
# 此处及以下全部采用 "YKT6拟时序hlb重构.r" 的规范逻辑
# ------------------------------------------------------------------------------
st_tbl <- table(pData(cds)$State, pData(cds)$celltype.main)
st_df <- as.data.frame.matrix(st_tbl)
st_df$NPC_Score <- st_df[, "NPC-like"] / rowSums(st_df)
best_start_state <- rownames(st_df)[which.max(st_df$NPC_Score)]

print(paste(">>> 自动识别到 NPC-like 最高得分的 State 为:", best_start_state, "已将其设为 Root 起点。"))
cds <- orderCells(cds, root_state = as.numeric(best_start_state))

print(">>> 正在生成基础轨迹图...")

# 图1：带有黑圈 Node 节点的分支结构图 
pdf("Global_1_Objective_State_With_Nodes.pdf", width = 8, height = 6)
p_state <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.8, show_branch_points = TRUE) + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_state)
dev.off()

# 图2：纯净版 CellType 轨迹图
pdf("Global_2_Objective_CellType_Clean.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") +
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4))) 
print(p_type)
dev.off()

# 图3：拟时序方向图
pdf("Global_3_Objective_Pseudotime.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
print(p_time)
dev.off()

# ==============================================================================
# 手动强制覆写正确起点并生成正式图表
# ==============================================================================
# 【极其重要】：在这里填入你在 TEMP 图里看到的最边缘起点数字！(根据实际图表调整)
manual_root_state <- 18

print(paste(">>> 正在强制重置起点为 State", manual_root_state, "..."))
cds <- orderCells(cds, root_state = manual_root_state)

print(">>> 正在生成正式的客观骨架图...")

pdf("Global_1_Objective_State_With_Nodes1.pdf", width = 8, height = 6)
p_state <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.8, show_branch_points = TRUE) + 
  theme_void() + theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_state)
dev.off()

pdf("Global_2_Objective_CellType_Clean1.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") +
  theme_void() + theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4))) 
print(p_type)
dev.off()

pdf("Global_3_Objective_Pseudotime1.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + 
  theme_void() + theme(legend.position = "right", legend.title = element_text(face = "bold"))
print(p_time)
dev.off()

saveRDS(cds, "Master_Objective_CDS.rds")
print(">>> ✅ 阶段一执行完毕！手动校准的客观 CDS 骨架已安全保存。")


# ==============================================================================
# 🌟 全新数据驱动流程 - 阶段二：BEAM 客观分支动力学与热图鉴定
# ==============================================================================
# 【核心参数配置区】
core_branch_point <- 4           
num_heatmap_clusters <- 4            
my_target_genes <- c("YKT6") 
validation_genes <- c("CD44", "PDGFRA") 

print(">>> 读取第一阶段生成的客观拟时序骨架...")
cds <- readRDS("Master_Objective_CDS.rds")

# 1. 全局 BEAM 计算
beam_res_file <- paste0("BEAM_Results_Node", core_branch_point, ".rds")
if (file.exists(beam_res_file)) {
  BEAM_res <- readRDS(beam_res_file)
} else {
  tryCatch({
    BEAM_res <- BEAM(cds, branch_point = core_branch_point, cores = 12, progenitor_method = 'duplicate') 
  }, error = function(e) {
    BEAM_res <<- BEAM(cds, branch_point = core_branch_point, cores = 1, progenitor_method = 'duplicate')
  })
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  saveRDS(BEAM_res, beam_res_file)
}

# 2. 独立课题循环
for (target_gene in my_target_genes) {
  if (!(target_gene %in% row.names(fData(cds)))) next
  
  # 图 A：全景无字热图
  top300_genes <- row.names(subset(BEAM_res, qval < 1e-4))[1:300]
  if (!(target_gene %in% top300_genes)) top300_genes <- c(top300_genes, target_gene)
  
  cds_global <- cds[top300_genes, ]
  pdf(paste0(target_gene, "_Objective_Global_Heatmap_Top300.pdf"), width = 6, height = 8) 
  beam_heatmap_obj <- plot_genes_branched_heatmap(cds_global, branch_point = core_branch_point, 
                                                  num_clusters = num_heatmap_clusters, 
                                                  cores = 1, show_rownames = FALSE, return_heatmap = TRUE) 
  dev.off()
  
  gene_clusters <- beam_heatmap_obj$annotation_row
  gene_clusters$Gene <- rownames(gene_clusters)
  gene_clusters <- gene_clusters[, c("Gene", "Cluster")] %>% arrange(Cluster)
  write.csv(gene_clusters, paste0(target_gene, "_Objective_BEAM_Top300_Clusters.csv"), row.names = FALSE)

  # 图 B：极简裁判热图
  genes_to_plot <- unique(c(validation_genes, target_gene))
  genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(fData(cds))]
  cds_detail <- cds[genes_to_plot, ]
  
  pdf(paste0(target_gene, "_Objective_Detail_Targets_Only.pdf"), width = 6, height = 3.5) 
  p_detail <- plot_genes_branched_heatmap(cds_detail, branch_point = core_branch_point, 
                                          num_clusters = 2, cores = 1, 
                                          use_gene_short_name = TRUE, show_rownames = TRUE) 
  print(p_detail)
  dev.off()
}


# ==============================================================================
# 🌟 全新数据驱动流程 - 阶段三：生物学命运锁定与终极动力学出图
# ==============================================================================
print(">>> 读取第一阶段的客观拟时序骨架...")
cds <- readRDS("Master_Objective_CDS.rds")
dim_coords <- reducedDimS(cds)

# 1. 【核心校准区】：基于 BEAM 铁证进行终极 State 绑定 (根据你的数据调整这些数字)
pre_branch_states <- c(17, 18, 19) 
fate1_states <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)  
fate2_states <- c(1)                     

# 初始化 Fate_State，并覆写严密的生物学标签
pData(cds)$Fate_State <- "Unknown"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(pre_branch_states)] <- "Pre-branch"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(fate1_states)] <- "Fate 1"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(fate2_states)] <- "Fate 2"

pData(cds)$Fate_State <- factor(pData(cds)$Fate_State, levels = c("Pre-branch", "Fate 1", "Fate 2", "Unknown"))

pdf("Global_4_Objective_Locked_Fate_Check.pdf", width = 8, height = 6)
p_check <- plot_cell_trajectory(cds, color_by = "Fate_State", cell_size = 0.8, show_branch_points = TRUE) + 
  scale_color_manual(values = c("Pre-branch" = "#999999", "Fate 1" = "#EB5A61", "Fate 2" = "#7592C9", "Unknown" = "black")) +
  theme_void() + theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_check)
dev.off()

saveRDS(cds, "Master_Final_Mapped_CDS.rds")

# ------------------------------------------------------------------------------
# 终极动力学折线图 (引入标准化)
# ------------------------------------------------------------------------------
target_gene <- "YKT6"   
custom_y_max <- 0.5    

# 核心标准化公式提取
expr_val <- as.numeric(exprs(cds)[target_gene, ] / pData(cds)$Size_Factor)
plot_df <- data.frame(
  Pseudotime = pData(cds)$Pseudotime,
  CellType = pData(cds)$celltype.main,
  Fate_State = pData(cds)$Fate_State,
  Expression = log1p(expr_val)
) %>% arrange(Expression) 
  
if (!is.null(custom_y_max)) {
  y_max_zoom <- custom_y_max
} else {
  pos_expr <- plot_df$Expression[plot_df$Expression > 0]
  y_max_zoom <- if(length(pos_expr) > 10) quantile(pos_expr, 0.95) else max(plot_df$Expression)
  if (y_max_zoom == 0) y_max_zoom <- 0.5 
}

branch_time <- min(plot_df$Pseudotime[plot_df$Fate_State %in% c("Fate 1", "Fate 2")], na.rm = TRUE)
max_time <- max(plot_df$Pseudotime, na.rm = TRUE)
  
df_lineage1 <- plot_df %>% filter(Fate_State %in% c("Pre-branch", "Fate 1")) %>% mutate(Branch = "Trajectory to Fate 1")
df_lineage2 <- plot_df %>% filter(Fate_State %in% c("Pre-branch", "Fate 2")) %>% mutate(Branch = "Trajectory to Fate 2")
plot_df_branches <- rbind(df_lineage1, df_lineage2)
  
theme_sci_clean <- function() {
  theme_classic(base_size = 14) %+replace% 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
          axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, face = "bold"),
          axis.line = element_line(linewidth = 0.8, color = "black"), axis.ticks = element_line(linewidth = 0.8, color = "black"),
          legend.text = element_text(size = 12), legend.title = element_text(size = 13, face = "bold"),
          legend.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank())
}

pdf(paste0(target_gene, "_Final_Fig8_BlackLine_Dynamics.pdf"), width = 9.5, height = 5.5)
p_dynamics <- ggplot(plot_df_branches, aes(x = Pseudotime)) +
  geom_point(aes(y = Expression, color = CellType), size = 1.0, alpha = 0.5, stroke = 0) + 
  geom_vline(xintercept = branch_time, linetype = "dotted", color = "grey50", linewidth = 0.8) +
  annotate("text", x = branch_time + (max_time * 0.02), y = y_max_zoom * 0.9, label = "Branch Point", color = "black", hjust = 0, size = 4) +
  geom_smooth(aes(y = Expression, linetype = Branch), color = "black", se = FALSE, method = "loess", span = 0.8, linewidth = 1.2) +
  scale_color_manual(values = my_colors, name = "Cell Type") +
  scale_linetype_manual(values = c("Trajectory to Fate 1" = "solid", "Trajectory to Fate 2" = "dashed"), name = "Branch") +
  scale_x_continuous(limits = c(0, max_time), expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(0, y_max_zoom)) + 
  labs(title = paste0(target_gene, " Expression Dynamics"), x = "Exact Pseudotime", y = "Log(Expression + 1)") + 
  theme_sci_clean() + theme(legend.position = "right", legend.key.width = unit(1.5, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
print(p_dynamics)
dev.off()

# ------------------------------------------------------------------------------
# 终局：全局纯净轨迹图与基因表达空间投射 (SCI 级)
# ------------------------------------------------------------------------------
fate_colors <- c("Pre-branch" = "#999999", "Fate 1" = "#EB5A61", "Fate 2" = "#7592C9", "Unknown" = "black")

pdf("Final_Fig1_Trajectory_by_CellType_Clean.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") + theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
print(p_type)
dev.off()

pdf("Final_Fig2_Trajectory_by_Pseudotime_Clean.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12))
print(p_time)
dev.off()

pdf("Final_Fig3_Trajectory_by_FateState_Clean.pdf", width = 8, height = 6)
p_fate <- plot_cell_trajectory(cds, color_by = "Fate_State", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_manual(values = fate_colors, name = "Lineage Fate") + theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
print(p_fate)
dev.off()

# 靶基因表达空间投射图
for (target_gene in my_target_genes) {
  if (!(target_gene %in% row.names(fData(cds)))) next
  expr_val <- as.numeric(exprs(cds)[target_gene, ] / pData(cds)$Size_Factor)
  expr_log <- log1p(expr_val)
  plot_df <- data.frame(x = dim_coords[1, ], y = dim_coords[2, ], Expression = expr_log) %>% arrange(Expression) 
  pos_expr <- plot_df$Expression[plot_df$Expression > 0]
  color_max <- if(length(pos_expr) > 5) quantile(pos_expr, 0.95) else max(plot_df$Expression)
  plot_df$Expression[plot_df$Expression > color_max] <- color_max

  pdf(paste0("Final_Fig4_Projection_", target_gene, "_Clean.pdf"), width = 8, height = 6)
  p_proj <- ggplot(plot_df, aes(x = x, y = y)) + 
    geom_point(color = "grey90", size = 1.0) + geom_point(aes(color = Expression), size = 1.5, alpha = 0.8) + 
    scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"), name = "Log(Expr)") + theme_void() + 
    theme(legend.position = "right", legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
    labs(title = paste0(target_gene, " Expression Projection"))
  print(p_proj)
  dev.off()
}

# ------------------------------------------------------------------------------
# 补充验证：批次效应 (将原版 orig.ident 替换为 SampleID)
# ------------------------------------------------------------------------------
if (!("SampleID" %in% colnames(pData(cds)))) {
  print("⚠️ 未找到 SampleID 列，如果您使用了其他列名，请在数据准备阶段修改")
} else {
  pdf("Final_Fig5_Trajectory_by_Batch_Clean.pdf", width = 8, height = 6)
  p_batch <- plot_cell_trajectory(cds, color_by = "SampleID", cell_size = 1.0, show_branch_points = FALSE) + 
    theme_void() + 
    theme(legend.position = "right", legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
    labs(title = "Trajectory Colored by SampleID") +
    guides(color = guide_legend(override.aes = list(size = 5), title = "Sample ID")) 
  print(p_batch)
  dev.off()
}

print(">>> 🎉 恭喜！全套定制数据流程已顺利执行结束！")