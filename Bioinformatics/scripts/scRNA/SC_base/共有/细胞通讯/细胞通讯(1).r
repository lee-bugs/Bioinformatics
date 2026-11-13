# ==============================================================================
# 0. 环境配置与全局参数设置 (综合合并版)
# 功能：
#   1. Seurat / CellChat 计算，导出显著通讯表格
#   2. 筛选 High > Low (Diff>0.02) 或 High 特有通路
#   3. 生成 CellChat 原始气泡图，以及筛选补齐对比后的气泡图
# ==============================================================================
rm(list = ls())
gc()

# 加载所有依赖包
if (!require("mixtools", quietly = TRUE)) install.packages("mixtools")
library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(mixtools)
library(future)

# --- 并行计算设置 ---
# 提升多线程允许导出的最大内存上限 (8GB)
options(future.globals.maxSize = 8000 * 1024^2)
# 全局开启 12 个线程进行并行运算
plan("multisession", workers = 12)

# --- 路径与全局参数设置 ---
# 🚨 定义目标基因 (这也会影响工作路径的名称)
target_gene <- "ARL6IP1"

# 设定新的工作路径
work_dir <- file.path("/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE221534_Object_1_Core/细胞通讯", target_gene)

# 如果工作路径不存在则自动创建
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
  message(paste(">>> 已创建新的工作目录:", work_dir))
} else {
  message(paste(">>> 工作目录已存在:", work_dir))
}

setwd(work_dir)

# 新的 rds 路径
rds_file <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE221534_Object_1_Core/Glioma_3B_Fully_Annotated.rds"
base_out_dir  <- file.path(work_dir, "output") 
if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)


# 🚨 定义要循环的细胞类型列表
target_cell_list <- c("NPC-like", "AC-like", "OPC-like") 


# 使用命名向量为不同细胞类型设置专属 cutoff
manual_cutoff_list <- c( 
  "AC-like"  = 0.8,
  "NPC-like" = 0.75,
  "OPC-like" = 0.8) 

# 差异阈值 (挑选时使用)
diff_threshold <- 0

# ==============================================================================
# 1. 定义全局函数
# ==============================================================================
# 绘图函数 1: 原始 CellChat 气泡图
plot_safe_bubble <- function(object, sources, targets, title_prefix, filename, save_dir) {
  df_check <- subsetCommunication(object, sources.use = sources, targets.use = targets, thresh = 0.05)
  if (nrow(df_check) > 0) {
    message(paste0("      ✅ [", title_prefix, "] 检测到 ", nrow(df_check), " 条显著互作，正在绘图..."))
    p <- netVisual_bubble(object, sources.use = sources, targets.use = targets, 
                          angle.x = 45, remove.isolate = FALSE, thresh = 0.05) + 
      ggtitle(paste0(title_prefix, " (Unadjusted P < 0.05)"))
    ggsave(file.path(save_dir, filename), p, width = 12, height = 10)
  } else {
    message(paste0("      ⚠️ [", title_prefix, "] 在 P<0.05 下未检测到显著互作，跳过绘图。"))
  }
}

# 辅助函数 2: 提取基名
get_base_name <- function(name) {
  name <- str_remove(name, "_High$")
  name <- str_remove(name, "_Low$")
  return(name)
}

# 绘图函数 3: 挑选后的自定义对比气泡图
plot_single_cell <- function(data, group_name, save_dir) {
  if(nrow(data) == 0) return(NULL)
  
  labels_df <- data %>%
    dplyr::select(pair_label, My_Direction, Partner_Cell, Condition) %>%
    distinct() %>%
    arrange(desc(My_Direction), Partner_Cell, Condition)
  
  data$pair_label <- factor(data$pair_label, levels = labels_df$pair_label)
  
  p <- ggplot(data, aes(x = pair_label, y = interaction_name)) +
    geom_point(data = subset(data, prob > 0), aes(size = pval_cat, color = prob)) +
    scale_color_gradientn(colors = c("blue", "white", "red"), name = "Comm Prob") +
    scale_size_manual(values = c("p < 0.01" = 2, "0.01 < p < 0.05" = 1), 
                      name = "P-value", breaks = c("p < 0.01", "0.01 < p < 0.05")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 6),
      axis.text.y = element_text(color = "black", size = 6),
      panel.grid.major = element_line(color = "grey92", linetype = "dashed"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.6), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10)
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "", y = "Interaction Name")
  
  safe_name <- gsub("/", "_", group_name)
  filename <- paste0("Bubble_", safe_name, "_Selected.pdf")
  
  h <- max(4, length(unique(data$interaction_name)) * 0.15)
  w <- max(5, length(unique(data$pair_label)) * 0.35)
  
  ggsave(file.path(save_dir, filename), p, width = w, height = h, limitsize = FALSE)
  print(paste("已保存:", filename))
}

# ==============================================================================
# 2. 加载数据并开始循环分析
# ==============================================================================
message("\n>>> [全局] 加载 Seurat 对象 (仅执行一次)...")
if (!file.exists(rds_file)) stop("找不到 RDS 文件！")
sc_obj_original <- readRDS(rds_file)

for (current_cell in target_cell_list) {
  
  message(paste0("\n", paste(rep("=", 60), collapse = "")))
  message(paste0(">>> 正在处理细胞类型: ", current_cell))
  message(paste0(paste(rep("=", 60), collapse = "")))
  
  # 独立输出目录
  out_dir <- file.path(base_out_dir, current_cell)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 使用副本，防止污染
  sc_obj <- sc_obj_original
  
  # ----------------------------------------------------------------------------
  # [STEP 1] 人工 Cutoff 分组与 GMM 绘图
  # ----------------------------------------------------------------------------
  message("   >>> [1/6] 运行人工 Cutoff 进行分组...")
  current_cutoff <- manual_cutoff_list[current_cell]
  
  cell_ids <- WhichCells(sc_obj, idents = current_cell)
  if(length(cell_ids) == 0) {
    message(paste("   ⚠️️ 跳过：未在数据中找到", current_cell))
    next
  }
  
  tryCatch({
    expr_data <- GetAssayData(sc_obj, layer = "data")[target_gene, cell_ids]
  }, error = function(e) {
    expr_data <<- GetAssayData(sc_obj, slot = "data")[target_gene, cell_ids]
  })
  
  high_ids <- names(expr_data)[expr_data > current_cutoff]
  low_ids <- names(expr_data)[expr_data <= current_cutoff]
  
  message(paste("      High 组细胞数:", length(high_ids)))
  message(paste("      Low  组细胞数:", length(low_ids)))
  
  # GMM 分布图
  set.seed(123)
  fit_data <- expr_data + runif(length(expr_data), min = 0, max = 1e-5)
  mixmdl <- tryCatch({ normalmixEM(fit_data, k = 2, verb = FALSE) }, error = function(e) NULL)
  
  if(!is.null(mixmdl)) {
    high_comp <- which.max(mixmdl$mu); low_comp <- which.min(mixmdl$mu)
    mu_high <- mixmdl$mu[high_comp]; sigma_high <- mixmdl$sigma[high_comp]; lambda_high <- mixmdl$lambda[high_comp]
    mu_low <- mixmdl$mu[low_comp]; sigma_low <- mixmdl$sigma[low_comp]; lambda_low <- mixmdl$lambda[low_comp]
    
    p_gmm <- ggplot(data.frame(Expression = expr_data), aes(x = Expression)) +
      geom_histogram(aes(y = ..density..), bins = 60, fill = "grey85", color = "black", alpha = 0.7) +
      stat_function(geom = "line", fun = function(x) lambda_low * dnorm(x, mean = mu_low, sd = sigma_low), color = "#3376B2", linewidth = 1.2) +
      stat_function(geom = "line", fun = function(x) lambda_high * dnorm(x, mean = mu_high, sd = sigma_high), color = "#D2352C", linewidth = 1.2) +
      geom_vline(xintercept = current_cutoff, linetype = "dashed", color = "black", linewidth = 1) +
      theme_classic() + coord_cartesian(ylim = c(0, max(lambda_high * dnorm(mu_high, mu_high, sigma_high)) * 2)) + 
      labs(title = paste0(target_gene, " GMM Reference - ", current_cell), x = "LogNormalize Expression", y = "Density")
    
    ggsave(file.path(out_dir, "GMM_Distribution_Check.pdf"), p_gmm, width = 8, height = 6)
  }
  
  # 修改分组
  sc_obj$cellchat_split <- as.character(Idents(sc_obj))
  sc_obj$cellchat_split[colnames(sc_obj) %in% high_ids] <- paste0(current_cell, "_High")
  sc_obj$cellchat_split[colnames(sc_obj) %in% low_ids]  <- paste0(current_cell, "_Low")
  Idents(sc_obj) <- sc_obj$cellchat_split
  group_high <- paste0(current_cell, "_High")
  group_low  <- paste0(current_cell, "_Low")
  
  # ----------------------------------------------------------------------------
  # [STEP 2] CellChat 计算流程
  # ----------------------------------------------------------------------------
  message("   >>> [2/6] 运行 CellChat 分析...")
  cellchat <- createCellChat(object = sc_obj, group.by = "cellchat_split")
  cellchat@DB <- CellChatDB.human 
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 5) 
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, file.path(out_dir, "CellChat_Log_Replicated.RDS"))
  
  # ----------------------------------------------------------------------------
  # [STEP 3] 网络图输出
  # ----------------------------------------------------------------------------
  message("   >>> [3/6] 正在绘制网络图...")
  mat <- cellchat@net$weight
  groupSize <- as.numeric(table(cellchat@idents))
  
  pdf(file.path(out_dir, "Network_Compare_In_vs_Out.pdf"), width = 16, height = 8)
  par(mfrow = c(1, 2), xpd=TRUE) 
  if(sum(mat[, c(group_high, group_low)]) > 0) {
    netVisual_circle(mat, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                     targets.use = c(group_high, group_low), title.name = paste("Incoming:", current_cell))
  } else { plot.new(); title("No Incoming Signal") }
  netVisual_circle(mat, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                   sources.use = c(group_high, group_low), title.name = paste("Outgoing:", current_cell))
  dev.off()
  
  # ----------------------------------------------------------------------------
  # [STEP 4] 基础气泡图与导出未筛选的全数据 CSV
  # ----------------------------------------------------------------------------
  message("   >>> [4/6] 导出基础气泡图与原始显著 CSV...")
  others_use <- levels(Idents(sc_obj))[!levels(Idents(sc_obj)) %in% c(group_high, group_low)]
  plot_safe_bubble(cellchat, c(group_high, group_low), others_use, "Outgoing", "Bubble_Outgoing.pdf", out_dir)
  plot_safe_bubble(cellchat, others_use, c(group_high, group_low), "Incoming", "Bubble_Incoming.pdf", out_dir)
  
  df_out <- subsetCommunication(cellchat, sources.use = c(group_high, group_low), thresh = 0.05)
  if (nrow(df_out) > 0) df_out$Direction <- "Outgoing"
  
  df_in <- subsetCommunication(cellchat, targets.use = c(group_high, group_low), thresh = 0.05)
  if (nrow(df_in) > 0) df_in$Direction <- "Incoming"
  
  df_all <- rbind(df_out, df_in)
  
  input_csv <- file.path(out_dir, "All_Interactions_Table_P005.csv")
  if (nrow(df_all) > 0) {
    write.csv(df_all, input_csv, row.names = FALSE)
  } else {
    message("   ⚠️ 未发现显著通讯通路，跳过当前细胞类型的筛选。")
    next
  }
  
  # ----------------------------------------------------------------------------
  # [STEP 5] 数据处理与筛选 (原: 细胞通讯挑选.R)
  # ----------------------------------------------------------------------------
  message("   >>> [5/6] 正在读取并挑选细胞通讯文件...")
  output_csv <- file.path(out_dir, "All_Interactions_Table_P005 -挑选后.csv")

  df_raw <- read.csv(input_csv, stringsAsFactors = FALSE)
  original_cols <- colnames(df_raw)

  # --- 预处理：识别 High/Low 和 配对信息 ---
  df_proc <- df_raw %>%
    dplyr::mutate(
      Group_Type = ifelse(str_detect(source, "_High") | str_detect(target, "_High"), "High", "Low"),
      Partner = ifelse(Direction == "Outgoing", target, source),
      Match_ID = paste(interaction_name, ligand, receptor, Direction, Partner, sep = "@")
    ) %>%
    dplyr::filter(!str_detect(Partner, "_High|_Low"))

  # --- 核心筛选：转宽表计算差异 ---
  df_stats <- df_proc %>%
    dplyr::select(Match_ID, Group_Type, prob) %>%
    tidyr::pivot_wider(names_from = Group_Type, values_from = prob, values_fill = 0, values_fn = max) %>%
    { if(!"High" %in% names(.)) dplyr::mutate(., High = 0) else . } %>%
    { if(!"Low" %in% names(.)) dplyr::mutate(., Low = 0) else . } %>%
    dplyr::filter((Low == 0 & High > 0) | (Low > 0 & (High - Low) > diff_threshold))

  selected_ids <- df_stats$Match_ID
  message(paste("      >>> 筛选完成，共找到符合条件的通路数量:", length(selected_ids)))

  if (length(selected_ids) > 0) {
    message("      >>> 正在构建对比表格...")
    
    df_high_selected <- df_proc %>% dplyr::filter(Match_ID %in% selected_ids, Group_Type == "High")
    df_low_selected  <- df_proc %>% dplyr::filter(Match_ID %in% selected_ids, Group_Type == "Low")
    
    ids_only_high <- setdiff(df_high_selected$Match_ID, df_low_selected$Match_ID)
    
    if (length(ids_only_high) > 0) {
      df_dummy_low <- df_high_selected %>%
        dplyr::filter(Match_ID %in% ids_only_high) %>%
        dplyr::mutate(
          source = str_replace(source, "_High", "_Low"),
          target = str_replace(target, "_Low", "_High"),
          target = str_replace(target, "_High", "_Low"),
          prob = 0,
          pval = 1
        )
      df_low_final <- dplyr::bind_rows(df_low_selected, df_dummy_low)
    } else {
      df_low_final <- df_low_selected
    }
    
    # --- 合并 High 和 Low 并还原格式 ---
    df_filtered_final <- dplyr::bind_rows(df_high_selected, df_low_final) %>%
      dplyr::arrange(interaction_name, desc(source)) %>%
      dplyr::select(dplyr::all_of(original_cols))
    
    # 导出挑选后的文件
    write.csv(df_filtered_final, output_csv, row.names = FALSE)
    message(paste("      ✅ 挑选处理完成！文件已保存至:", output_csv))
    
    # ----------------------------------------------------------------------------
    # [STEP 6] 挑选后数据解析与可视化绘图 (原: 可视化.R)
    # ----------------------------------------------------------------------------
    message("\n   >>> [6/6] 开始进行可视化绘图...")
    
    # 直接读取导出的表格
    df_vis <- read.csv(output_csv, stringsAsFactors = FALSE, check.names = FALSE)
    
    df_parsed <- df_vis %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        s_base = get_base_name(source),
        t_base = get_base_name(target),
        
        Condition = case_when(
          str_detect(source, "_High$") | str_detect(target, "_High$") ~ "High",
          str_detect(source, "_Low$") | str_detect(target, "_Low$") ~ "Low",
          TRUE ~ "Other"
        ),
        
        Plot_Group = case_when(
          str_detect(source, "_High$") | str_detect(source, "_Low$") ~ s_base,
          str_detect(target, "_High$") | str_detect(target, "_Low$") ~ t_base,
          TRUE ~ "Other"
        ),
        
        My_Direction = case_when(
          str_detect(source, "_High$") | str_detect(source, "_Low$") ~ "Outgoing",
          str_detect(target, "_High$") | str_detect(target, "_Low$") ~ "Incoming",
          TRUE ~ "Other"
        ),
        
        Partner_Cell = ifelse(My_Direction == "Outgoing", target, source)
      ) %>%
      dplyr::filter(Condition %in% c("High", "Low")) %>% 
      dplyr::ungroup()
    
    # 补全缺失的 Low 数据 (容错处理)
    df_high <- df_parsed %>% dplyr::filter(Condition == "High")
    df_low <- df_parsed %>% dplyr::filter(Condition == "Low")
    
    df_high <- df_high %>% dplyr::mutate(key = paste(s_base, t_base, interaction_name, sep = "|"))
    df_low <- df_low %>% dplyr::mutate(key = paste(s_base, t_base, interaction_name, sep = "|"))
    
    missing_keys <- setdiff(df_high$key, df_low$key)
    
    if(length(missing_keys) > 0) {
      df_missing <- df_high %>%
        dplyr::filter(key %in% missing_keys) %>%
        dplyr::mutate(
          Condition = "Low",
          prob = 0, pval = 1,
          source = str_replace(source, "_High$", "_Low"),
          target = str_replace(target, "_High$", "_Low")
        ) %>%
        dplyr::select(-key)
      
      df_plot_final <- bind_rows(df_parsed, df_missing)
    } else {
      df_plot_final <- df_parsed
    }
    
    # 排序逻辑与绘图准备
    df_plot <- df_plot_final %>%
      dplyr::mutate(
        pval_cat = case_when(
          pval < 0.01 ~ "p < 0.01",
          pval >= 0.01 & pval < 0.05 ~ "0.01 < p < 0.05",
          TRUE ~ "NS" 
        ),
        pair_label = paste(source, "->", target)
      )
    
    df_plot$pval_cat <- factor(df_plot$pval_cat, levels = c("p < 0.01", "0.01 < p < 0.05", "NS"))
    
    # 执行循环出图并写入各自分组目录
    groups <- unique(df_plot$Plot_Group)
    for (grp in groups) {
      message(paste("      正在处理绘图分组:", grp))
      sub_data <- df_plot %>% dplyr::filter(Plot_Group == grp)
      # 使用集成带入 save_dir 参数的统一 plot_single_cell 函数
      plot_single_cell(sub_data, grp, save_dir = out_dir) 
    }
  } else {
    message("   ⚠️ 没有筛选到符合条件的通路，跳过图表绘制。")
  }

  # CellChat 对象备份
  save_name <- paste0("CellChat_Backup_Final_", current_cell, ".rds")
  saveRDS(cellchat, file = save_name)
  message(paste0("✅ ", current_cell, " 的流程全部完成！CellChat 对象已保存在当前工作目录: ", save_name))
}

message("\n🎉 恭喜！全部细胞类型的计算、挑选与可视化流程已运行完成！")