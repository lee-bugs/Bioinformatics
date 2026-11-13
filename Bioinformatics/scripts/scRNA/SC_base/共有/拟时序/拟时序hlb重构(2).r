# ==============================================================================
# 🌟 全新数据驱动流程 - 阶段一：全局客观拟时序主干构建
# 核心逻辑：仅依靠拓扑算法建树并定位起点，不预设任何 Fate 分支标签
# ==============================================================================

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
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/拟时序_全局重构"
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

my_colors <- c("NPC-like" = "#E64B35", "OPC-like" = "#4DBBD5", "AC-like" = "#00A087", "MES-like" = "#3C5488")

# ------------------------------------------------------------------------------
# 2. 读取数据与提取肿瘤细胞 (规范抽样 1200 个细胞保证树形稳健)
# ------------------------------------------------------------------------------
print(">>> 读取 Seurat RDS 数据并提取肿瘤细胞...")
scRNA <- readRDS("/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/Glioma_3B_Fully_Annotated.rds")

tumor_classes <- c("NPC-like", "OPC-like", "AC-like", "MES-like")
Idents(scRNA) <- "CellType"
scRNA_tumor <- subset(scRNA, idents = tumor_classes)
scRNA_tumor$celltype.main <- Idents(scRNA_tumor)

set.seed(2026) 
sampled_cells <- sample(colnames(scRNA_tumor), min(1200, ncol(scRNA_tumor)))
scRNA_sub <- subset(scRNA_tumor, cells = sampled_cells)

# ------------------------------------------------------------------------------
# 3. 严谨筛选用于建树的排序基因 (Ordering Genes) 并构建 CDS
# ------------------------------------------------------------------------------
print(">>> 构建 CellDataSet 对象并筛选建树核心基因...")
current_assay <- DefaultAssay(scRNA_sub)
data <- as(as.matrix(GetAssayData(scRNA_sub, slot = "counts", assay = current_assay)), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_sub@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(data), row.names = row.names(data)))

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 提高阈值获取高阶差异基因，并剔除细胞周期与应激基因
Idents(scRNA_sub) <- "celltype.main"
deg_markers <- FindAllMarkers(scRNA_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
valid_degs <- deg_markers %>% filter(p_val_adj < 1e-5) %>% pull(gene) %>% unique()

scRNA_sub <- FindVariableFeatures(scRNA_sub, selection.method = "vst", nfeatures = 1500)
hvgs <- VariableFeatures(scRNA_sub)
core_genes <- intersect(valid_degs, hvgs)

cell_cycle_genes <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))
junk_genes <- grep("^MT-|^RP[SL]|^HSP", rownames(scRNA_sub), value = TRUE)
ordering_genes <- setdiff(core_genes, c(cell_cycle_genes, junk_genes))

cds <- setOrderingFilter(cds, ordering_genes)

# ------------------------------------------------------------------------------
# 4. DDRTree 降维建树 (包含 orig.ident 批次矫正)
# ------------------------------------------------------------------------------
print(">>> DDRTree 降维建树中 (请耐心等待)...")
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', 
                       residualModelFormulaStr = "~orig.ident", auto_param_selection = FALSE)
cds <- orderCells(cds)

# ------------------------------------------------------------------------------
# 5. 智能定位客观起点 (将 NPC-like 富集度最高的端点设为根节点)
# ------------------------------------------------------------------------------
st_tbl <- table(pData(cds)$State, pData(cds)$celltype.main)
st_df <- as.data.frame.matrix(st_tbl)
st_df$NPC_Score <- st_df[, "NPC-like"] / rowSums(st_df)
best_start_state <- rownames(st_df)[which.max(st_df$NPC_Score)]

print(paste(">>> 自动识别到 NPC-like 最高得分的 State 为:", best_start_state, "已将其设为 Root 起点。"))
cds <- orderCells(cds, root_state = as.numeric(best_start_state))

# ------------------------------------------------------------------------------
# 6. 生成基础客观轨迹图 (重点用于下一步确认核心分叉点 Node)
# ------------------------------------------------------------------------------
print(">>> 正在生成基础轨迹图...")

# 图1：带有黑圈 Node 节点的分支结构图 (极其重要，用于确定 BEAM 的 branch_point)
pdf("Global_1_Objective_State_With_Nodes.pdf", width = 8, height = 6)
p_state <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.8, show_branch_points = TRUE) + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_state)
dev.off()

# 图2：纯净版 CellType 轨迹图 (对照看走向)
pdf("Global_2_Objective_CellType_Clean.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") +
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4))) 
print(p_type)
dev.off()

# 图3：拟时序方向图 (检查深蓝色是否正确位于起点)
pdf("Global_3_Objective_Pseudotime.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
print(p_time)
dev.off()






# ==============================================================================
# 4. 手动强制覆写正确起点并生成正式图表
# ==============================================================================
# 【极其重要】：在这里填入你刚才在 TEMP 图里看到的最边缘起点数字！
manual_root_state <- 1  

print(paste(">>> 正在强制重置起点为 State", manual_root_state, "..."))
cds <- orderCells(cds, root_state = manual_root_state)

print(">>> 正在生成正式的客观骨架图...")

# 图1：带有黑圈 Node 节点的 State 结构图 (极其重要，用于确认分叉 Node)
pdf("Global_1_Objective_State_With_Nodes1.pdf", width = 8, height = 6)
p_state <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.8, show_branch_points = TRUE) + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_state)
dev.off()

# 图2：纯净版 CellType 轨迹图
pdf("Global_2_Objective_CellType_Clean1.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") +
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4))) 
print(p_type)
dev.off()

# 图3：拟时序方向图 (验证深蓝色是否完全从你的 manual_root_state 出发)
pdf("Global_3_Objective_Pseudotime1.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.8, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + 
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
print(p_time)
dev.off()

# ------------------------------------------------------------------------------
# 5. 保存完美骨架
# ------------------------------------------------------------------------------
saveRDS(cds, "Master_Objective_CDS.rds")
print(">>> ✅ 阶段一执行完毕！手动校准的客观 CDS 骨架已安全保存。")










# ==============================================================================
# 🌟 全新数据驱动流程 - 阶段二：BEAM 客观分支动力学与热图鉴定 (极简裁判版)
# ==============================================================================

options(stringsAsFactors = FALSE)
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)

# 开启多核并行加速 BEAM 计算
plan("multisession", workers = 12)
options(future.globals.maxSize = 100000 * 1024^2)

work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/拟时序_全局重构"
setwd(work_dir)

# ==============================================================================
# 【核心参数配置区】
# ==============================================================================
core_branch_point <- 9            
num_heatmap_clusters <- 4            # 热图聚类模块数（推荐 4，专门用于图 A 的富集大表）
my_target_genes <- c("PPAT", "TBCE") # 真正需要出 Top 300 富集大表的独立课题主角
validation_genes <- c("CD44", "PDGFRA") # 每次画极简明细图时，都会强制拉过来当“标尺”的裁判

print(">>> 读取第一阶段生成的客观拟时序骨架...")
cds <- readRDS("Master_Objective_CDS.rds")

# ------------------------------------------------------------------------------
# 1. 全局 BEAM 计算 (智能缓存，整个数据集只需算一次)
# ------------------------------------------------------------------------------
print(paste(">>>>>>>>>> 开启基于 Node", core_branch_point, "的全局 BEAM 分析 <<<<<<<<<<"))
beam_res_file <- paste0("BEAM_Results_Node", core_branch_point, ".rds")

if (file.exists(beam_res_file)) {
  print(">>> 检测到已存在的 BEAM 结果，直接加载以节省时间...")
  BEAM_res <- readRDS(beam_res_file)
} else {
  print(">>> 正在启动多核 BEAM 计算 (约需 5-15 分钟，请耐心等待)...")
  tryCatch({
    BEAM_res <- BEAM(cds, branch_point = core_branch_point, cores = 12, progenitor_method = 'duplicate') 
  }, error = function(e) {
    print("⚠️ 多核遇到瓶颈，降级单核运行...")
    BEAM_res <<- BEAM(cds, branch_point = core_branch_point, cores = 1, progenitor_method = 'duplicate')
  })
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  saveRDS(BEAM_res, beam_res_file)
}

# ------------------------------------------------------------------------------
# 2. 独立课题循环 (仅循环 PPAT 和 TBCE)
# ------------------------------------------------------------------------------
for (target_gene in my_target_genes) {
  print("----------------------------------------------------------------------")
  print(paste(">>> 🚀 正在处理独立课题主角：", target_gene))
  
  if (!(target_gene %in% row.names(fData(cds)))) next
  
  # ==============================================================================
  # 图 A：全景无字热图 & 提取 Top 300 富集分析大表
  # ==============================================================================
  top300_genes <- row.names(subset(BEAM_res, qval < 1e-4))[1:300]
  if (!(target_gene %in% top300_genes)) top300_genes <- c(top300_genes, target_gene)
  
  cds_global <- cds[top300_genes, ]
  
  pdf(paste0(target_gene, "_Objective_Global_Heatmap_Top300.pdf"), width = 6, height = 8) 
  beam_heatmap_obj <- plot_genes_branched_heatmap(cds_global, 
                                                  branch_point = core_branch_point, 
                                                  num_clusters = num_heatmap_clusters, # 调用开头设定的 4 类
                                                  cores = 1, 
                                                  show_rownames = FALSE, 
                                                  return_heatmap = TRUE) 
  dev.off()
  
  gene_clusters <- beam_heatmap_obj$annotation_row
  gene_clusters$Gene <- rownames(gene_clusters)
  gene_clusters <- gene_clusters[, c("Gene", "Cluster")] %>% arrange(Cluster)
  write.csv(gene_clusters, paste0(target_gene, "_Objective_BEAM_Top300_Clusters.csv"), row.names = FALSE)

  # ==============================================================================
  # 图 B：极简裁判热图 (仅展示 裁判基因 + 当前课题主角)
  # ==============================================================================
  # 将标尺基因 (CD44, PDGFRA) 强制加入，并去重
  genes_to_plot <- unique(c(validation_genes, target_gene))
  genes_to_plot <- genes_to_plot[genes_to_plot %in% row.names(fData(cds))]
  
  cds_detail <- cds[genes_to_plot, ]
  
  # 高度压低至 3.5，防止行数太少导致色块难看
  pdf(paste0(target_gene, "_Objective_Detail_Targets_Only.pdf"), width = 6, height = 3.5) 
  p_detail <- plot_genes_branched_heatmap(cds_detail, 
                                          branch_point = core_branch_point, 
                                          num_clusters = 2, # 强制锁死 2 类，防止 3 个基因聚 4 类报错
                                          cores = 1, 
                                          use_gene_short_name = TRUE, 
                                          show_rownames = TRUE) 
  print(p_detail)
  dev.off()
  
  print(paste(">>> ✅ [", target_gene, "] 的极简裁判热图出图完毕！"))
}
print(">>> 🎉 阶段二：客观分支动力学鉴定已全部完成！")







# ==============================================================================
# 🌟 全新数据驱动流程 - 阶段三：生物学命运锁定、视觉核验与终极动力学出图
# ==============================================================================

options(stringsAsFactors = FALSE)
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/拟时序_全局重构"
setwd(work_dir)

print(">>> 读取第一阶段的客观拟时序骨架...")
cds <- readRDS("Master_Objective_CDS.rds")
dim_coords <- reducedDimS(cds)

# ==============================================================================
# 1. 【核心校准区】：基于 BEAM 铁证进行终极 State 绑定
# ⚠️ 提示：请核对 Global_1_Objective_State_With_Nodes.pdf 中的数字
# ==============================================================================
# 处于起点与 Node 9 之间的树干阶段
pre_branch_states <- c(1, 2, 3, 4, 21) 

# 走向 Cell fate 1 (左侧，证实为 CD44+ 的 MES-like) 的分支
fate1_states <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)  

# 走向 Cell fate 2 (右侧，证实为 PDGFRA+ 的 OPC-like) 的分支
fate2_states <- c(16, 17, 18, 19, 20)                     
# ==============================================================================

# 初始化 Fate_State，并覆写严密的生物学标签
pData(cds)$Fate_State <- "Unknown"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(pre_branch_states)] <- "Pre-branch"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(fate1_states)] <- "Fate 1"
pData(cds)$Fate_State[as.character(pData(cds)$State) %in% as.character(fate2_states)] <- "Fate 2"

# 转换为 Factor 锁定顺序
pData(cds)$Fate_State <- factor(pData(cds)$Fate_State, 
                                levels = c("Pre-branch", "Fate 1", "Fate 2", "Unknown"))

# ------------------------------------------------------------------------------
# 2. 生成检查图：视觉核验 State 映射是否无遗漏
# ------------------------------------------------------------------------------
print(">>> 正在生成视觉核验图，请检查是否有黑色的 Unknown 细胞...")

pdf("Global_4_Objective_Locked_Fate_Check.pdf", width = 8, height = 6)
p_check <- plot_cell_trajectory(cds, color_by = "Fate_State", cell_size = 0.8, show_branch_points = TRUE) + 
  scale_color_manual(values = c("Pre-branch" = "#999999", 
                                "Fate 1" = "#EB5A61", 
                                "Fate 2" = "#7592C9", 
                                "Unknown" = "black")) +
  theme_void() + 
  theme(legend.position = "right", legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
print(p_check)
dev.off()

# 保存这份终极完美的 CDS（包含正确的起点、Node 和 Fate_State）
saveRDS(cds, "Master_Final_Mapped_CDS.rds")
print(">>> ✅ 分支身份已客观锁定！终极 CDS 已保存。核验图 Global_4_Objective_Locked_Fate_Check.pdf 已生成。")











# ------------------------------------------------------------------------------
# 3. 单基因高定版：生成终极动力学折线图 (方便手动反复微调 Y 轴)
# ------------------------------------------------------------------------------
print(">>> 开始生成单基因终极动力学折线图...")

my_colors <- c("NPC-like" = "#E64B35", "OPC-like" = "#4DBBD5", "AC-like" = "#00A087", "MES-like" = "#3C5488")

# ==============================================================================
# 🛠️ 手动调节区：每次运行前，只需修改这里！
# ==============================================================================
target_gene <- "TBCE"   # <--- ① 修改你要出图的基因 (比如先做 "PPAT"，做完再改成 "TBCE")

# ② 修改 Y 轴上限：越小，趋势线的起伏波动就越夸张、越明显。
# 比如你可以尝试 0.35, 0.30, 甚至 0.25。如果填 NULL，算法会自动计算。
custom_y_max <- 0.15    
# ==============================================================================

print(paste(">>> 正在生成定制折线图：", target_gene, " | Y轴上限:", custom_y_max))

if (!(target_gene %in% row.names(fData(cds)))) {
  stop(paste("⚠️ 错误：基因", target_gene, "不在表达矩阵中！请检查拼写。"))
}
  
# 1. 提取表达数据
expr_val <- as.numeric(exprs(cds)[target_gene, ] / pData(cds)$Size_Factor)
plot_df <- data.frame(
  Pseudotime = pData(cds)$Pseudotime,
  CellType = pData(cds)$celltype.main,
  Fate_State = pData(cds)$Fate_State,
  Expression = log1p(expr_val)
) %>% arrange(Expression) 
  
# 2. 应用你指定的 Y 轴上限
if (!is.null(custom_y_max)) {
  y_max_zoom <- custom_y_max
} else {
  pos_expr <- plot_df$Expression[plot_df$Expression > 0]
  y_max_zoom <- if(length(pos_expr) > 10) quantile(pos_expr, 0.95) else max(plot_df$Expression)
  if (y_max_zoom == 0) y_max_zoom <- 0.5 
}

# 3. 提取分支时间点 (适配 "Fate 1" / "Fate 2" 命名)
branch_time <- min(plot_df$Pseudotime[plot_df$Fate_State %in% c("Fate 1", "Fate 2")], na.rm = TRUE)
max_time <- max(plot_df$Pseudotime, na.rm = TRUE)
  
# 4. 分离两条命运分支的数据
df_lineage1 <- plot_df %>% filter(Fate_State %in% c("Pre-branch", "Fate 1")) %>% mutate(Branch = "Trajectory to Fate 1")
df_lineage2 <- plot_df %>% filter(Fate_State %in% c("Pre-branch", "Fate 2")) %>% mutate(Branch = "Trajectory to Fate 2")
plot_df_branches <- rbind(df_lineage1, df_lineage2)
  
# SCI 极简主题
theme_sci_clean <- function() {
  theme_classic(base_size = 14) %+replace% 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(linewidth = 0.8, color = "black"), 
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank()
    )
}

# 5. 绘图并输出
pdf_name <- paste0(target_gene, "_Final_Fig8_BlackLine_Dynamics.pdf")
pdf(pdf_name, width = 9.5, height = 5.5)

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
  theme_sci_clean() +
  theme(legend.position = "right", legend.key.width = unit(1.5, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
print(p_dynamics)
dev.off()

print(paste(">>> 🎉", target_gene, "图表已成功保存！(文件：", pdf_name, ")"))











# ==============================================================================
# 🌟 全新数据驱动流程 - 终局：全局纯净轨迹图与基因表达空间投射 (SCI 级)
# ==============================================================================

options(stringsAsFactors = FALSE)
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/拟时序_全局重构"
setwd(work_dir)

print(">>> 正在加载终极完美的 CDS 骨架...")
cds <- readRDS("Master_Final_Mapped_CDS.rds")

# SCI 配色方案
my_colors <- c("NPC-like" = "#E64B35", "OPC-like" = "#4DBBD5", "AC-like" = "#00A087", "MES-like" = "#3C5488")
fate_colors <- c("Pre-branch" = "#999999", "Fate 1" = "#EB5A61", "Fate 2" = "#7592C9", "Unknown" = "black")

# ------------------------------------------------------------------------------
# 1. 基础属性纯净轨迹图 (去除所有 Node 节点)
# ------------------------------------------------------------------------------
print(">>> 正在生成全局属性纯净轨迹图...")

# 图 1：细胞类型纯净轨迹 (CellType)
pdf("Final_Fig1_Trajectory_by_CellType_Clean.pdf", width = 8, height = 6)
p_type <- plot_cell_trajectory(cds, color_by = "celltype.main", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_manual(values = my_colors, name = "Cell Type") +
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
print(p_type)
dev.off()

# 图 2：伪时间纯净轨迹 (Pseudotime) - 展现从深蓝到深红的演化
pdf("Final_Fig2_Trajectory_by_Pseudotime_Clean.pdf", width = 8, height = 6)
p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_gradientn(colours = c("#3366CC", "#FFCC33", "#CC0000"), name = "Pseudotime") + 
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12))
print(p_time)
dev.off()

# 图 3：命运分支纯净轨迹 (Fate State) - 展现你精准锁定的三大阶段
pdf("Final_Fig3_Trajectory_by_FateState_Clean.pdf", width = 8, height = 6)
p_fate <- plot_cell_trajectory(cds, color_by = "Fate_State", cell_size = 1.0, show_branch_points = FALSE) + 
  scale_color_manual(values = fate_colors, name = "Lineage Fate") +
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
print(p_fate)
dev.off()

# ------------------------------------------------------------------------------
# 2. 靶基因表达空间投射图 (Gene Spatial Projection)
# ------------------------------------------------------------------------------
print(">>> 正在生成靶基因的空间表达投射图...")

# 如果你需要跑其他基因，直接在这个数组里加即可
target_genes <- c("PPAT", "TBCE") 
dim_coords <- reducedDimS(cds)

for (target_gene in target_genes) {
  if (!(target_gene %in% row.names(fData(cds)))) next
  
  # 提取标准化表达量
  expr_val <- as.numeric(exprs(cds)[target_gene, ] / pData(cds)$Size_Factor)
  expr_log <- log1p(expr_val)
  
  # 构建绘图数据框，按表达量排序，让高表达的红点浮在灰点最上方
  plot_df <- data.frame(
    x = dim_coords[1, ],
    y = dim_coords[2, ],
    Expression = expr_log
  ) %>% arrange(Expression) 
  
  # 强制去除极值点对颜色映射的干扰，使红色更加鲜艳饱满
  pos_expr <- plot_df$Expression[plot_df$Expression > 0]
  color_max <- if(length(pos_expr) > 5) quantile(pos_expr, 0.95) else max(plot_df$Expression)
  plot_df$Expression[plot_df$Expression > color_max] <- color_max

  pdf(paste0("Final_Fig4_Projection_", target_gene, "_Clean.pdf"), width = 8, height = 6)
  
  # 采用 ggplot2 自定义绘制，不画底层的树枝黑线，视觉极其干净
  p_proj <- ggplot(plot_df, aes(x = x, y = y)) + 
    geom_point(color = "grey90", size = 1.0) + # 铺一层灰色底
    geom_point(aes(color = Expression), size = 1.5, alpha = 0.8) + # 叠加上表达量
    scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"), name = "Log(Expr)") + 
    theme_void() + 
    theme(legend.position = "right", 
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
    labs(title = paste0(target_gene, " Expression Projection"))
  
  print(p_proj)
  dev.off()
  print(paste(">>> ✅", target_gene, "投射图已保存！"))
}

print(">>> 🎉 恭喜！课题的全部数据分析与核心图表制作已圆满收官！")





# ==============================================================================
# 🌟 补充验证：批次效应 (Batch Effect) 轨迹分布图
# 目的：向审稿人证明轨迹分支是由细胞亚型驱动，而非样本批次分离所致
# ==============================================================================

options(stringsAsFactors = FALSE)
library(monocle)
library(ggplot2)

work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core/拟时序_全局重构"
setwd(work_dir)

print(">>> 正在加载终极 CDS 骨架...")
cds <- readRDS("Master_Final_Mapped_CDS.rds")

# 检查批次信息是否存在于 pData 中 (通常 Seurat 默认将批次存在 orig.ident 里)
if (!("orig.ident" %in% colnames(pData(cds)))) {
  print("⚠️ 未找到 orig.ident 列，如果您使用了其他列名（如 Sample, patient_id），请在代码中替换 'orig.ident'")
} else {
  print(">>> 正在生成批次效应验证图 (按 orig.ident 染色)...")
  
  # 图 5：批次效应纯净轨迹 (Batch / orig.ident)
  pdf("Final_Fig5_Trajectory_by_Batch_Clean.pdf", width = 8, height = 6)
  
  p_batch <- plot_cell_trajectory(cds, color_by = "orig.ident", cell_size = 1.0, show_branch_points = FALSE) + 
    # 使用 ggplot2 默认的离散色系，自动适配不同样本数量
    theme_void() + 
    theme(legend.position = "right", 
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
    labs(title = "Trajectory Colored by Sample/Batch") +
    guides(color = guide_legend(override.aes = list(size = 5), title = "Sample ID")) 
    
  print(p_batch)
  dev.off()
  
  print(">>> ✅ 批次效应验证图已成功保存为：Final_Fig5_Trajectory_by_Batch_Clean.pdf")
}