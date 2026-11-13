# ==============================================================================
# 1. 基础设置 (路径 & 包)
# ==============================================================================

packages <- c("Seurat", "ggplot2", "dplyr", "ggsci", "patchwork", "scales", "tidyr", "parallel")
lapply(packages, library, character.only = TRUE)

# --- 请确认路径 ---
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE223063(scRNA)/4.目标基因的投射" 
data_path <- file.path(work_dir, "Glioma_3B_Fully_Annotated.rds")
base_output_dir <- file.path(work_dir, "投射")

if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

# --- 加载数据 ---
if (!exists("scRNA")) {
  message(">>> 正在读取 RDS 文件...")
  if (file.exists(data_path)) {
    scRNA <- readRDS(data_path)
    message(">>> 数据加载成功！")
  } else {
    stop("找不到RDS文件，请检查路径。")
  }
}

# ==============================================================================
# 🚨 数据基座加固：彻底避开 Seurat v5 Layer Bug 的降维打击法 (顺序修复版)
# ==============================================================================

message(">>> 1. 先在 v5 状态下执行 JoinLayers 缝合碎片化矩阵 ...")
tryCatch({
  scRNA <- JoinLayers(scRNA)
}, error = function(e) message("    - 无需缝合或非 v5 碎片，安全跳过。"))

message(">>> 2. 强制将 decontX 转换为经典的 v3 Assay 格式 ...")
scRNA[["decontX"]] <- as(object = scRNA[["decontX"]], Class = "Assay")

message(">>> 3. 强制设置默认 Assay 为 decontX ...")
DefaultAssay(scRNA) <- "decontX"

message(">>> 4. 正在进行 LogNormalize 标准化补充 data 槽位 ...")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
message(">>> 标准化完成，底层数据已完全就绪！")

# --- 自动适配颜色 (确保14个簇不报错) ---
num_clusters <- length(unique(Idents(scRNA)))
if (num_clusters <= 10) {
  my_colors <- pal_npg("nrc")(num_clusters)
} else if (num_clusters <= 20) {
  my_colors <- pal_d3("category20")(num_clusters) 
} else {
  my_colors <- pal_igv("default")(num_clusters)
}

# --- 绘图通用设置 ---
sci_theme <- theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "right", strip.background = element_blank())

# 🚨 新增：提取自注释代码中的 Style 2 专属极简箭头坐标轴主题
theme_sci_arrow <- theme_classic() +
  theme(
    axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ==============================================================================
# 2. 基因列表 与 多核安全配置 (提取自: 基因的数据投射.xlsx)
# ==============================================================================

target_genes <- c(
"RAB13", "AK2", "CKS2", "LEPROT", "MLX", "SCPEP1", "TGIF1", "PARD6B", "RFK", "GANAB",
  "SLC30A7", "RHOC", "SHC1", "IER5L", "HIST1H2AC", "MEAF6", "ZNF593", "ULBP2", "IMMT", "HEATR2",
  "CEP135", "SLC43A3", "LEPRE1", "RHBDF1", "PMM2", "NFE2L3", "DEPDC1", "ALG8", "ZNF367", "RDH10",
  "VKORC1L1", "STK40", "PXN", "PSRC1", "ARPC5", "TMSB4X", "UGCG", "SLC44A3", "HIST1H2BH", "RAD54L",
  "TEAD3", "TROAP", "DBF4", "RPN1", "SEC61A1", "DPAGT1", "ICAM3", "RAD54B", "TNFAIP3", "KRT75",
  "HIST1H4H", "FXYD5", "KCTD19", "LACTB2", "ADM2", "ARHGEF19", "HN1L", "IGLL5", "ZNF618", "PTGER2",
  "CD300E", "BTBD19", "SYDE1", "OTP", "HSPB7", "SLC25A45", "TSKU", "HIST1H2BO", "KLF17", "TTC38",
  "NFE2", "PARVA", "ADPRH", "VWA1", "FAM181A", "ACTN1", "FAF1", "HIST1H2BG", "NRD1", "HIST1H3H",
  "LYPLA2", "HIST1H2AG", "REXO2", "EXOSC10", "HIST1H2AM", "HIST2H2BE", "TOE1", "CENPN", "HIST1H2BJ", "SF3A3",
  "APOBR", "SERPINF1", "GPR171", "MTHFS", "PARP15", "ITGAM", "CEP112", "CPA6", "SOCS2", "SP6",
  "ADPGK", "SEL1L3", "CEP89", "PODNL1", "GPBP1L1", "CALML4", "ACTA1", "HIST1H2BI", "HIST1H3E", "MYO5B",
  "IRX5", "P2RY6", "SLC6A9", "PRAM1", "SLC2A10", "SGMS2", "SAMD9L", "HS2ST1", "TRAF5", "ENPEP",
  "IQCG", "TSPAN31", "CPA4", "CTBS", "DCTD", "S100A3", "PLAT", "CASP6", "ZMYM1", "LYZ",
  "CHPF2", "B4GALT5", "COLGALT1", "TXNDC12", "ELK3", "RNFT1", "SH3PXD2B", "LRRN4CL", "SHOX2", "CAPN5",
  "EDEM2", "DENND2A", "TRAM2", "TXLNA", "SPOCD1", "CDCA2", "SHCBP1", "DPY19L1", "COL5A2", "TK1",
  "DIRAS3", "ERI1", "SHISA5", "YBX3", "SEC61B", "TTF2", "SURF4", "ERCC6L", "OSTC", "SPAG5",
  "TAF12", "STIL", "HIST1H2BE", "RHOQ", "DDOST", "MRGBP", "EIF3I", "RASSF1", "DAP", "PLIN2",
  "CCDC77", "HIST1H2BD", "AGBL2", "NEXN", "BMP5", "CHIC2", "LARP1B", "FYCO1", "KIAA1257", "LRRC66",
  "NIPSNAP3A", "ITPKC", "NCK2", "LPHN2", "LMO4", "STIMATE", "UCN2", "CA2", "CD3E", "CERS5",
  "MTMR11", "RTCA", "VGLL2", "PTCRA", "CCNJL", "ARSI", "VSIG10L", "RLF", "FHOD1", "RAB36",
  "HIVEP3", "PDLIM4", "ARL9", "TMEM59", "FRMD8", "DUSP14", "DNAJC22", "FAM120AOS", "SYNC", "CCNYL1",
  "CEMIP2", "TBC1D1", "SLC47A2", "PHTF1", "CAPZA1", "FNDC11", "OSCP1", "GLB1", "GDAP2", "MLTK",
  "PHETA2", "FAM90A1", "ARAP3", "ERP44", "FAM187A", "SLC11A1", "CC2D1B", "LRRIQ1", "BCORL1", "EYA3",
  "PHYH", "SIGLEC1", "RHBDL2", "SHROOM3", "PDZK1", "PSME4", "SZRD1", "RGS3", "ARID5A", "SLC39A1",
  "S100A5", "GSG2", "GJB2", "TAX1BP3", "PDLIM7", "SCMH1", "RFX8", "SLC35A2", "SPINT1", "SLC35C1",
  "SLC25A24", "PHF13", "ASB6", "ADCK2", "RFFL", "CD5", "PIGV", "TMEM69", "UAP1", "BOLA3",
  "SBNO2", "SLC20A1", "PPCS", "DPP9", "HIST1H2BB", "GNPNAT1", "HIST3H2BB", "EXTL2", "TBC1D10C", "FAM117A",
  "SRSF4", "THRAP3", "EEF1AKMT3", "RDM1", "HIST1H2BM", "MRPS15", "PNRC2", "FAM136A", "KPNA6", "ARHGAP11A",
  "UBE2Z", "FBLIM1", "ESCO2", "CNIH4", "AKIRIN1", "POC1A", "PTPN1", "PTMA", "SDHB", "TMSB10",
  "CYB5B", "ITPRID1", "PPIH", "ARPC3", "RPAP3", "CAP1", "TIMM22", "HIST2H2AA4", "BCAS2", "TTC4",
  "SHQ1", "DNAJC8", "SEC22B", "ALG14", "RPF1", "HIST1H2BC", "HIST1H2BN", "TUBA1C", "ZBTB8OS", "EBNA1BP2",
  "NPLOC4", "SEH1L", "SPCS3", "CAPZB", "MRPS6", "NDE1", "HIST1H2BL", "SNAPC2", "SGO1", "TUBG1",
  "COPG1", "GPN1", "TRAPPC3", "PSMD11", "KTI12", "POLR2L", "MAK16", "GNAI3", "ORMDL2", "MYOSLID",
  "RRN3P2", "BRD4", "AURKB", "MCM2", "RBBP4", "NDC1", "POLA2", "PSMC4", "PCNA", "TRIP13",
  "RPA3", "SKA3", "CDC45", "SEC13", "SRSF1", "PTTG1", "MCM3", "DSCC1", "CNOT10", "STXBP3",
  "RFC3", "SPATA6", "DCAF13", "WDHD1", "ZNF69", "ZC3HAV1", "SLC41A2", "CCDC15", "NAGLU", "NOL11",
  "RHNO1", "WDR62", "SCYL2", "TXLNB", "BST1", "CENPH", "MORC4", "POLQ", "ZWILCH", "GGCX",
  "FAM122B", "CTNS", "SGCB", "ARHGAP17", "LIMD2", "UTP20", "ATAD5", "YKT6", "TMED5", "UBE2J2",
  "TAF13", "TBCE", "TCOF1", "UTP3", "TMEM53", "UNC45A", "UXS1", "WRAP73", "PPAT", "AMZ1",
  "ARL6IP1", "DNPEP", "RTKN2", "SERTAD2", "SUSD1"
)

# 💡 安全的多核分配策略：固定上限 8 个核
max_cores <- 8
cores_to_use <- min(max_cores, parallel::detectCores() - 1)
cores_to_use <- max(1, cores_to_use) 

# ==============================================================================
# 3. 第一部分：开始自动绘图 [多核并行版]
# ==============================================================================

message("\n==============================================")
message(paste0(">>> [出图阶段] 本次任务共 ", length(target_genes), " 个基因..."))
message(paste0(">>> 🚀 正在开启多核并行加速出图 (分配安全线程数: ", cores_to_use, ")，请稍等..."))
message("==============================================")

# 💡 收集每个基因的处理状态
draw_status_list <- mclapply(target_genes, function(target_gene) {
  
  # 1. 检查基因是否存在
  if (!target_gene %in% rownames(scRNA)) {
    return("Skipped") # 未找到，记录为 Skipped
  }
  
  # 2. 创建文件夹
  gene_dir <- file.path(base_output_dir, target_gene)
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  # ============================================================================
  # 优化版降维投射图 (加入 max.cutoff 截断 + 高辨识度热力色谱 + 适度增大点径)
  # ============================================================================
  
  # 推荐配色：零表达为干净淡灰，一旦表达立即呈现明亮黄/橙，高表达呈现鲜红/深红
  hi_contrast_colors <- c("grey88", "#fee08b", "#fdae61", "#f46d43", "#d73027", "#a50026")
  
  # (1) UMAP - Style 1 (Standard)
  tryCatch({
    p_umap_s1 <- FeaturePlot(scRNA, features = target_gene, reduction = "umap", 
                             pt.size = 0.35, order = TRUE, raster = FALSE,
                             #max.cutoff = "q95"  #去掉这里的截取，结果偏差较大
							 ) + 
      scale_color_gradientn(colors = hi_contrast_colors) +
      sci_theme + theme(axis.line = element_blank(), axis.text = element_blank(), 
                        axis.ticks = element_blank(), axis.title = element_blank()) +
      ggtitle(paste0(target_gene, " (UMAP Style 1)"))
    ggsave(file.path(gene_dir, paste0("FeaturePlot_UMAP_Style1_", target_gene, ".pdf")), 
           p_umap_s1, width = 6, height = 5)
  }, error = function(e) NULL)
  
  # (2) UMAP - Style 2 (Reshaped Minimal, 使用 "umap_alt" 及箭头主题)
  tryCatch({
    p_umap_s2 <- FeaturePlot(scRNA, features = target_gene, reduction = "umap_alt", 
                             pt.size = 0.35, order = TRUE, raster = FALSE,
                             #max.cutoff = "q95"  #去掉这里的截取，结果偏差较大
							 ) + 
      scale_color_gradientn(colors = hi_contrast_colors) +
      theme_sci_arrow + labs(x = "UMAP 1", y = "UMAP 2") +
      ggtitle(paste0(target_gene, " (UMAP Style 2)"))
    ggsave(file.path(gene_dir, paste0("FeaturePlot_UMAP_Style2_", target_gene, ".pdf")), 
           p_umap_s2, width = 6, height = 5)
  }, error = function(e) NULL)
  
  # (3) t-SNE - Style 1 (Standard)
  tryCatch({
    p_tsne_s1 <- FeaturePlot(scRNA, features = target_gene, reduction = "tsne", 
                             pt.size = 0.35, order = TRUE, raster = FALSE,
                             #max.cutoff = "q95"  #去掉这里的截取，结果偏差较大
							 ) + 
      scale_color_gradientn(colors = hi_contrast_colors) +
      sci_theme + theme(axis.line = element_blank(), axis.text = element_blank(), 
                        axis.ticks = element_blank(), axis.title = element_blank()) +
      ggtitle(paste0(target_gene, " (t-SNE Style 1)"))
    ggsave(file.path(gene_dir, paste0("FeaturePlot_tSNE_Style1_", target_gene, ".pdf")), 
           p_tsne_s1, width = 6, height = 5)
  }, error = function(e) NULL)
  
  # (4) t-SNE - Style 2 (Reshaped Minimal, 使用 "tsne_alt" 及箭头主题)
  tryCatch({
    p_tsne_s2 <- FeaturePlot(scRNA, features = target_gene, reduction = "tsne_alt", 
                             pt.size = 0.35, order = TRUE, raster = FALSE,
                             #max.cutoff = "q95"  #去掉这里的截取，结果偏差较大
							 ) + 
      scale_color_gradientn(colors = hi_contrast_colors) +
      theme_sci_arrow + labs(x = "t-SNE 1", y = "t-SNE 2") +
      ggtitle(paste0(target_gene, " (t-SNE Style 2)"))
    ggsave(file.path(gene_dir, paste0("FeaturePlot_tSNE_Style2_", target_gene, ".pdf")), 
           p_tsne_s2, width = 6, height = 5)
  }, error = function(e) NULL)
  
  # --- 其他原版图表 (DotPlot, BarPlots, Violin) 保持不变 ---
  tryCatch({
    p4 <- DotPlot(scRNA, features = target_gene) +
      scale_color_gradientn(colors = c("#4DBBD5", "white", "#E64B35")) + 
      sci_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Expression Pattern", x = NULL, y = NULL)
    ggsave(file.path(gene_dir, paste0("New_DotPlot_", target_gene, ".pdf")), p4, width = 6, height = 5)
    
    dot_data <- p4$data 
    
    raw_expr <- FetchData(scRNA, vars = c("ident", target_gene), slot = "data")
    colnames(raw_expr) <- c("Cluster", "Exp_Value")
    
    df_bar_avg <- raw_expr %>%
      group_by(Cluster) %>%
      summarise(val = mean(Exp_Value, na.rm = TRUE), .groups = "drop") %>%
      rename(id = Cluster)
    
    df_bar_avg$id <- factor(df_bar_avg$id, levels = levels(Idents(scRNA)))
    
    p5 <- ggplot(df_bar_avg, aes(x = id, y = val, fill = id)) +
      geom_bar(stat = "identity", color = "black", width = 0.7) +
      scale_fill_manual(values = my_colors) + 
      labs(title = paste0(target_gene, " Avg. Exp"), y = "Expression (LogNorm)", x = NULL) +
      sci_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(gene_dir, paste0("BarPlot_AvgExp_Abs_", target_gene, ".pdf")), p5, width = 6, height = 5)
    
    p6 <- ggplot(dot_data, aes(x = id, y = pct.exp, fill = id)) +
      geom_bar(stat = "identity", color = "black", width = 0.7) +
      scale_fill_manual(values = my_colors) + 
      labs(title = paste0(target_gene, " % Expressed"), y = "Percentage (%)", x = NULL) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(dot_data$pct.exp) * 1.15)) +
      sci_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(gene_dir, paste0("BarPlot_PctExp_", target_gene, ".pdf")), p6, width = 6, height = 5)
    
  }, error = function(e) NULL)
  
  tryCatch({
    p3 <- VlnPlot(scRNA, features = target_gene, pt.size = 0) + 
      scale_fill_manual(values = my_colors) + 
      labs(x = "Cell Type", y = "Expression Level") +
      ggtitle(paste0(target_gene, " Expression Distribution")) +
      sci_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(gene_dir, paste0("ViolinPlot_Clean_", target_gene, ".pdf")), p3, width = 8, height = 6)
  }, error = function(e) NULL)
  
  return("Success") 
}, mc.cores = cores_to_use)

# 💡 解析绘图结果用于最终战报
draw_status <- unlist(draw_status_list)
num_draw_success <- sum(draw_status == "Success", na.rm = TRUE)
num_draw_skipped <- sum(draw_status == "Skipped", na.rm = TRUE)


# ==============================================================================
# 4. 第二部分：生成检测统计表 (V / B 逻辑) [多核并行版]
# ==============================================================================

message("\n==============================================")
message(paste0(">>> [统计阶段] 🚀 开启多核并行计算统计表 (分配安全线程数: ", cores_to_use, ")..."))
message("==============================================")

stats_list_raw <- mclapply(target_genes, function(target_gene) {
  
  if (!target_gene %in% rownames(scRNA)) {
    return(NULL)
  }
  
  # --- 步骤 A: V 打标逻辑 ---
  p_dot <- DotPlot(scRNA, features = target_gene) +
    scale_color_gradientn(colors = c("#4DBBD5", "white", "#E64B35"))
  plot_data <- p_dot$data 
  
  plot_data$Mark_V <- ""
  val_min <- min(plot_data$avg.exp.scaled, na.rm = TRUE)
  val_max <- max(plot_data$avg.exp.scaled, na.rm = TRUE)
  val_mid <- (val_min + val_max) / 2
  
  red_indices <- which(plot_data$avg.exp.scaled > val_mid)
  
  if (length(red_indices) > 0) {
    red_subset <- plot_data[red_indices, ]
    sort_order <- order(red_subset$avg.exp.scaled, decreasing = TRUE)
    
    for (i in 1:length(sort_order)) {
      original_idx <- red_indices[sort_order[i]]
      if (i == 1) plot_data$Mark_V[original_idx] <- "VVV"
      else if (i == 2) plot_data$Mark_V[original_idx] <- "VV"
      else plot_data$Mark_V[original_idx] <- "V"
    }
  }
  
  # --- 步骤 B: B 打标逻辑 ---
  raw_data <- FetchData(scRNA, vars = c("ident", target_gene), slot = "data")
  colnames(raw_data) <- c("ident", "expression")
  
  cluster_checks <- raw_data %>%
    group_by(ident) %>%
    summarise(
      Detected_SD = sd(expression, na.rm = TRUE),
      Detected_Pct = sum(expression > 0) / n(),
      Detected_Max = max(expression, na.rm = TRUE),
      Count = n()
    ) %>%
    mutate(
      Detected_SD = ifelse(is.na(Detected_SD), 0, Detected_SD),
      Detected_Max = ifelse(is.na(Detected_Max), 0, Detected_Max),
      Mark_B = ifelse(
        Detected_SD  > 0.2  & 
          Detected_Pct > 0.05 & 
          Detected_Max > 1.0  &
          Count >= 3,
        "B", ""
      )
    )
  
  # --- 步骤 C: 合并与输出单行 df ---
  final_data <- left_join(plot_data, cluster_checks, by = c("id" = "ident"))
  final_data$Final_Tag <- paste0(final_data$Mark_V, final_data$Mark_B)
  
  df_gene <- data.frame(
    Gene = target_gene,
    CellType = final_data$id,
    Tag = final_data$Final_Tag
  )
  
  return(df_gene)
}, mc.cores = cores_to_use)

stats_list <- Filter(Negate(is.null), stats_list_raw)

# ==============================================================================
# 5. 最终保存结果 (CSV) 与 执行战报
# ==============================================================================

if (length(stats_list) > 0) {
  
  long_df <- do.call(rbind, stats_list)
  
  plot_levels <- levels(long_df$CellType)
  wide_df <- long_df %>%
    pivot_wider(names_from = CellType, values_from = Tag, values_fill = "")
  
  valid_cols <- intersect(plot_levels, colnames(wide_df))
  final_cols <- c("Gene", valid_cols)
  wide_df <- wide_df[, final_cols]
  
  output_file <- file.path(base_output_dir, "Target_Genes_Matrix_StrictCheck.csv")
  write.csv(wide_df, file = output_file, row.names = FALSE, quote = FALSE)
  
}

message("\n==============================================")
message(">>> 🎉 全部任务圆满结束！执行战报如下：")
message(paste0("    - 设定分析基因总数 : ", length(target_genes), " 个"))
message(paste0("    - ✅ 成功绘图并统计 : ", num_draw_success, " 个"))
message(paste0("    - ❌ 未在矩阵中找到 : ", num_draw_skipped, " 个 (已跳过)"))
message("----------------------------------------------")
if (length(stats_list) > 0) {
  message(paste0("    - 统计矩阵已保存至 : \n      ", output_file))
} else {
  message("    - 警告: 未生成统计矩阵，请检查基因列表。")
}
message("==============================================\n")