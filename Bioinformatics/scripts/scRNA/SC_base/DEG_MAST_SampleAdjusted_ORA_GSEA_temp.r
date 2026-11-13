# ==============================================================================
# MAST DEG + Sample Adjustment + ORA + GSEA (多细胞类型循环版)
# 适配 MSigDB v2026.1.Hs official subcollections
# 移除 GMM 密度图，新增样本批次效应矫正(latent.vars)
# 包含进度条防刷屏修复与 12 线程并行
# ==============================================================================

rm(list = ls())
gc()

# 加载依赖包
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(future)
library(fgsea)
library(clusterProfiler)
library(pbapply)

if (!requireNamespace("MAST", quietly = TRUE)) {
  stop(
    "Package 'MAST' is required. Install it with: ",
    "BiocManager::install('MAST')"
  )
}

# ------------------------------------------------------------------------------
# 🚨 修复进度条满屏刷新的核心代码：强制缩短进度条长度，防止终端宽度不足导致折行
# ------------------------------------------------------------------------------
options(width = 150)
pboptions(type = "txt", style = 3, char = "=", txt.width = 50)


# ==============================================================================
# 1. 全局路径、基础参数与【Cutoff及剔除样本定义】
# ==============================================================================

target_gene <- "PPAT"

base_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core"

# 工作路径更新为 DEG_MAST_SA
work_dir_base <- paste0(
  base_dir, 
  "/通路分析/", target_gene, "/6_DEG差异性分析_富集/DEG_MAST_SA"
)

rds_file <- file.path(base_dir, "Glioma_3B_Fully_Annotated.rds")
gmt_dir  <- file.path(base_dir, "MSigDB_Sub")


# ------------------------------------------------------------------------------
# 🚨 核心配置区
# ------------------------------------------------------------------------------

cell_type_cutoffs <- list(
  "AC-like"  = 0.15,
  "NPC-like" = 0.15
)

sample_column <- "SampleID"

# 定义需要剔除的极少细胞量样本（若不需要剔除，留空即可：c()）
exclude_samples <- c("GSM5518629")

# DEG 与 GSEA 阈值参数
deg_fdr <- 0.05
deg_fc  <- 0.20

gsea_min_size <- 10
gsea_max_size <- 500

gsea_fdr_explore <- 0.25
gsea_fdr_plot <- 0.05
gsea_nes_plot <- 1.50
gsea_le_plot  <- 5

# 强制使用 12 个线程
workers <- 12


# ==============================================================================
# 2. 辅助函数定义 
# ==============================================================================

save_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
}


find_gmt <- function(pattern) {
  
  x <- list.files(
    gmt_dir,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(x) != 1) {
    stop(
      "GMT 文件匹配异常：", pattern,
      "\n找到：", paste(basename(x), collapse = ", ")
    )
  }
  
  x
}


to_term2gene <- function(pathways) {
  
  data.frame(
    term = rep(names(pathways), lengths(pathways)),
    gene = unlist(pathways, use.names = FALSE)
  )
}


parse_ratio <- function(x) {
  
  vapply(
    strsplit(as.character(x), "/"),
    function(z) as.numeric(z[1]) / as.numeric(z[2]),
    numeric(1)
  )
}


flatten_le <- function(x) {
  
  x$leadingEdge <- vapply(
    x$leadingEdge,
    paste,
    collapse = "/",
    FUN.VALUE = ""
  )
  
  x
}


gsea_title <- function(x) {
  toupper(gsub("_", " ", x))
}


save_gsea_3styles <- function(
    pathway,
    pathways,
    gene_rank,
    db,
    NES,
    padj,
    LE,
    output_dir,
    index
) {
  
  title_text <- gsea_title(pathway)
  
  safe_name <- substr(
    gsub("[/:?*<>|\\\\]", "_", pathway),
    1,
    80
  )
  
  prefix <- sprintf("%03d_%s", index, safe_name)
  
  subtitle_full <- sprintf(
    "%s | NES = %.2f | FDR = %.3g | LE = %d",
    db, NES, padj, LE
  )
  
  # Style 1 - Classic Red
  p1 <- plotEnrichment(
    pathways[[pathway]],
    gene_rank
  ) +
    labs(
      title = title_text,
      subtitle = subtitle_full,
      x = "Gene Rank",
      y = "Enrichment Score"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  
  if (length(p1$layers)) {
    p1$layers[[1]]$aes_params$colour <- "#D2352C"
    p1$layers[[1]]$aes_params$linewidth <- 1.3
  }
  
  # Style 2 - Black / White
  p2 <- plotEnrichment(
    pathways[[pathway]],
    gene_rank
  ) +
    labs(
      title = title_text,
      subtitle = sprintf("NES = %.2f | FDR = %.3g | LE = %d", NES, padj, LE),
      x = "Gene Rank",
      y = "Enrichment Score"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  
  if (length(p2$layers)) {
    p2$layers[[1]]$aes_params$colour <- "black"
    p2$layers[[1]]$aes_params$linewidth <- 1.25
  }
  
  # Style 3 - Minimal Blue
  p3 <- plotEnrichment(
    pathways[[pathway]],
    gene_rank
  ) +
    labs(
      title = title_text,
      subtitle = subtitle_full,
      x = "Ranked Genes",
      y = "Running Enrichment Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  
  if (length(p3$layers)) {
    p3$layers[[1]]$aes_params$colour <- "#3376B2"
    p3$layers[[1]]$aes_params$linewidth <- 1.4
  }
  
  # Save
  ggsave(file.path(output_dir, paste0(prefix, "_Style1_Classic.pdf")), p1, width = 6, height = 4.5)
  ggsave(file.path(output_dir, paste0(prefix, "_Style2_Black.pdf")), p2, width = 6, height = 4.5)
  ggsave(file.path(output_dir, paste0(prefix, "_Style3_Minimal.pdf")), p3, width = 6, height = 4.5)
}


# ==============================================================================
# 3. 预先加载 MSigDB 基因集与单细胞数据
# ==============================================================================

cat("\n>>> [加载阶段] 正在加载 MSigDB 基因集与全局 RDS...\n")

gmt_files <- c(
  Hallmark     = find_gmt("^h\\.all\\..*symbols\\.gmt$"),
  GO_BP        = find_gmt("^c5\\.go\\.bp\\..*symbols\\.gmt$"),
  GO_CC        = find_gmt("^c5\\.go\\.cc\\..*symbols\\.gmt$"),
  GO_MF        = find_gmt("^c5\\.go\\.mf\\..*symbols\\.gmt$"),
  Reactome     = find_gmt("^c2\\.cp\\.reactome\\..*symbols\\.gmt$"),
  KEGG_LEGACY  = find_gmt("^c2\\.cp\\.kegg_legacy\\..*symbols\\.gmt$"),
  WikiPathways = find_gmt("^c2\\.cp\\.wikipathways\\..*symbols\\.gmt$"),
  PID          = find_gmt("^c2\\.cp\\.pid\\..*symbols\\.gmt$"),
  C6           = find_gmt("^c6\\.all\\..*symbols\\.gmt$")
)

pathway_db <- lapply(gmt_files, gmtPathways)

msig_summary <- data.frame(
  Database = names(pathway_db),
  GeneSets = lengths(pathway_db),
  File = basename(gmt_files)
)

sc_obj_global <- readRDS(rds_file)
DefaultAssay(sc_obj_global) <- "decontX"
Idents(sc_obj_global) <- "CellType"


# ==============================================================================
# 4. [PHASE 1] 快速生成全局 SampleAdjusted 目录骨架与配置表
# ==============================================================================

cat("\n\n============================================================")
cat("\n>>> [PHASE 1] 启动前置工作：搭建全套 SampleAdjusted 目录结构...\n")
cat("============================================================\n")

for (ctype in names(cell_type_cutoffs)) {
  
  group_cutoff <- cell_type_cutoffs[[ctype]]
  if (!ctype %in% levels(Idents(sc_obj_global))) {
    next
  }
  
  # 创建带 SampleAdjusted 后缀的所有目录
  out_dir  <- file.path(work_dir_base, ctype)
  deg_dir  <- file.path(out_dir, "01_DEG_MAST_SampleAdjusted")
  ora_dir  <- file.path(out_dir, "02_ORA_MAST_SampleAdjusted")
  gsea_dir <- file.path(out_dir, "03_GSEA_MAST_SampleAdjusted")
  plot_dir <- file.path(gsea_dir, "Plots")
  
  dir.create(deg_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(ora_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 写入带有 SampleAdjusted 标识的配置文件
  save_csv(
    msig_summary, 
    file.path(out_dir, "00_MAST_SampleAdjusted_MSigDB_Subcollection_Summary.csv")
  )
  
  parameters <- data.frame(
    Parameter = c(
      "TargetGene", "CellType", "GroupCutoff", "DEG_Method", 
      "Sample_Adjustment", "Sample_Column", "DEG_FDR", "DEG_Log2FC", 
      "GSEA_MinSize", "GSEA_MaxSize", "GSEA_Exploratory_FDR", 
      "GSEA_Plot_FDR", "GSEA_Plot_NES", "GSEA_Plot_LeadingEdge"
    ),
    Value = c(
      target_gene, ctype, group_cutoff, "MAST", 
      "Fixed-effect covariate via latent.vars", sample_column, 
      deg_fdr, deg_fc, gsea_min_size, gsea_max_size, 
      gsea_fdr_explore, gsea_fdr_plot, gsea_nes_plot, gsea_le_plot
    )
  )
  save_csv(parameters, file.path(out_dir, "00_MAST_SampleAdjusted_Analysis_Parameters.csv"))
  
  cat(sprintf("\n    -> [%s] 专属 SampleAdjusted 目录结构及配置表创建完成。", ctype))
}
cat("\n\n>>> [PHASE 1 完毕] 结构已就绪，准备进入 MAST 样本矫正运算！\n")


# ==============================================================================
# 5. [PHASE 2] 执行 MAST DEG (Sample Adjusted) 与下游富集分析
# ==============================================================================

cat("\n\n============================================================")
cat("\n>>> [PHASE 2] 正式执行 MAST DEG ( latent.vars ) 与 ORA/GSEA...\n")
cat("============================================================\n")

for (ctype in names(cell_type_cutoffs)) {
  
  group_cutoff <- cell_type_cutoffs[[ctype]]
  if (!ctype %in% levels(Idents(sc_obj_global))) {
    next
  }
  
  cat(sprintf("\n>>> 正在处理 [%s] 的 MAST (Sample Adjusted) 差异分析...\n", ctype))
  
  # A. 路径重新指派
  out_dir  <- file.path(work_dir_base, ctype)
  deg_dir  <- file.path(out_dir, "01_DEG_MAST_SampleAdjusted")
  ora_dir  <- file.path(out_dir, "02_ORA_MAST_SampleAdjusted")
  gsea_dir <- file.path(out_dir, "03_GSEA_MAST_SampleAdjusted")
  plot_dir <- file.path(gsea_dir, "Plots")
  
  # B. 提取亚群并依据 exclude_samples 剔除指定样本
  obj <- subset(sc_obj_global, idents = ctype)
  
  if (length(exclude_samples) > 0) {
    cells_keep <- rownames(obj@meta.data)[!obj@meta.data[[sample_column]] %in% exclude_samples]
    obj <- subset(obj, cells = cells_keep)
    cat(sprintf("    -> 已排除样本: %s\n", paste(exclude_samples, collapse = ", ")))
  }
  
  # 校验 Sample Column 是否规范
  if (!sample_column %in% colnames(obj@meta.data)) {
    stop("Sample column not found: ", sample_column)
  }
  
  obj@meta.data[[sample_column]] <- droplevels(factor(obj@meta.data[[sample_column]]))
  
  if (nlevels(obj@meta.data[[sample_column]]) < 2) {
    warning(sprintf("    -> [%s] 剩余样本数不足2个，无法进行 latent.vars 矫正，跳过本亚群。", ctype))
    next
  }
  
  # C. 高低表达分组
  expr <- tryCatch(
    GetAssayData(obj, layer = "data")[target_gene, ], 
    error = function(e) GetAssayData(obj, slot = "data")[target_gene, ]
  )
  
  high_name <- paste0(target_gene, "_High")
  low_name  <- paste0(target_gene, "_Low")
  group_col <- paste0(target_gene, "_Group")
  
  obj[[group_col]] <- factor(
    ifelse(expr > group_cutoff, high_name, low_name), 
    levels = c(low_name, high_name)
  )
  Idents(obj) <- group_col
  
  group_counts <- data.frame(
    Group = names(table(Idents(obj))), 
    Cells = as.integer(table(Idents(obj)))
  )
  save_csv(group_counts, file.path(deg_dir, "Cell_Group_Counts.csv"))
  
  # D. Sample × Group 交叉验证 (QC)
  sample_group_tab <- table(
    Sample = obj@meta.data[[sample_column]], 
    Group = obj@meta.data[[group_col]]
  )
  
  sample_group_df <- as.data.frame(sample_group_tab) %>% 
    group_by(Sample) %>% 
    mutate(
      Sample_Total = sum(Freq), 
      Percentage = Freq / Sample_Total
    ) %>% 
    ungroup()
    
  save_csv(sample_group_df, file.path(deg_dir, "Sample_Group_Distribution.csv"))
  
  if (any(sample_group_tab == 0)) {
    bad_samples <- rownames(sample_group_tab)[apply(sample_group_tab == 0, 1, any)]
    warning(sprintf(
      "\n⚠️ 警告: [%s] 的某些样本缺乏完整的高低分组 (%s)。这会导致 latent.vars 报错！已跳过该亚群分析，请检查 Sample_Group_Distribution.csv 并将其加入 exclude_samples。\n", 
      ctype, paste(bad_samples, collapse = ", ")
    ))
    next
  }
  
  # E. MAST DEG 计算 (带入 latent.vars)
  cat(sprintf("    -> 执行 FindMarkers ( latent.vars = '%s' )...\n", sample_column))
  
  plan("multisession", workers = workers)
  options(future.globals.maxSize = 200 * 1024^3)
  
  deg <- FindMarkers(
    obj, 
    ident.1 = high_name, 
    ident.2 = low_name, 
    test.use = "MAST", 
    latent.vars = sample_column, 
    logfc.threshold = 0, 
    min.pct = 0.05, 
    verbose = TRUE
  )
  
  plan("sequential")
  
  deg$gene <- rownames(deg)
  deg$FDR <- p.adjust(deg$p_val, method = "BH")
  
  first_cols <- c("gene", "p_val", "FDR", "avg_log2FC")
  deg <- deg[, c(first_cols, setdiff(colnames(deg), first_cols))]
  
  save_csv(deg, file.path(deg_dir, "DEG_Full.csv"))
  
  plot_df <- deg %>% 
    filter(gene != target_gene) %>% 
    mutate(
      logFDR_raw = -log10(pmax(FDR, .Machine$double.xmin)), 
      logFDR = logFDR_raw,
      change = case_when(
        FDR < deg_fdr & avg_log2FC >= deg_fc ~ "UP", 
        FDR < deg_fdr & avg_log2FC <= -deg_fc ~ "DOWN", 
        TRUE ~ "STABLE"
      )
    )
  
  deg_summary <- plot_df %>% count(change)
  save_csv(deg_summary, file.path(deg_dir, "DEG_Summary.csv"))
  
  # F. 绘制全部三种风格的 Volcano 图
  junk_pattern <- paste0(
    "^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^AF[0-9]+\\.|^LINC[0-9]+|",
    "^MT-|^RP[SL]|^ENSG[0-9]+|-AS[0-9]+$|-IT[0-9]+$|^MIR|^SNO|^SNORD"
  )
  
  label_pool <- plot_df %>% 
    filter(change != "STABLE", !grepl(junk_pattern, gene)) %>% 
    mutate(score = abs(avg_log2FC) * logFDR)
  
  labels <- bind_rows(
    label_pool %>% filter(change == "UP") %>% slice_max(score, n = 8, with_ties = FALSE),
    label_pool %>% filter(change == "DOWN") %>% slice_max(score, n = 6, with_ties = FALSE)
  )
  
  x_max <- ceiling(max(abs(plot_df$avg_log2FC), na.rm = TRUE) * 10) / 10
  
  common_volcano <- list(
    geom_vline(xintercept = c(-deg_fc, deg_fc), linetype = "dashed", linewidth = 0.4),
    geom_hline(yintercept = -log10(deg_fdr), linetype = "dashed", linewidth = 0.4),
    coord_cartesian(xlim = c(-x_max, x_max)), 
    labs(x = "Log2 Fold Change", y = "-Log10(FDR)")
  )
  
  # Volcano Style 1
  p_vol1 <- ggplot(plot_df, aes(avg_log2FC, logFDR, color = change)) + 
    geom_point(alpha = 0.55, size = 1.3) + 
    scale_color_manual(values = c(UP = "#D2352C", DOWN = "#3376B2", STABLE = "grey85")) + 
    common_volcano + 
    geom_text_repel(data = labels, aes(label = gene), color = "black", size = 3.5, max.overlaps = Inf) + 
    theme_classic() + 
    theme(legend.position = "none")
    
  ggsave(file.path(deg_dir, paste0("Volcano_", target_gene, "_High_vs_Low_Style1_Classic.pdf")), p_vol1, width = 6.5, height = 5.5)
  
  # Volcano Style 2
  p_vol2 <- ggplot() + 
    geom_point(data = filter(plot_df, change == "STABLE"), aes(avg_log2FC, logFDR), color = "grey88", size = 1, alpha = 0.22) + 
    geom_point(data = filter(plot_df, change != "STABLE"), aes(avg_log2FC, logFDR, color = change), size = 1.4, alpha = 0.70) + 
    scale_color_manual(values = c(UP = "#D2352C", DOWN = "#3376B2")) + 
    common_volcano + 
    geom_text_repel(data = labels, aes(avg_log2FC, logFDR, label = gene), color = "black", size = 3.5, max.overlaps = Inf) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "none")
    
  ggsave(file.path(deg_dir, paste0("Volcano_", target_gene, "_High_vs_Low_Style2_Clean.pdf")), p_vol2, width = 6.5, height = 5.5)
  
  # Volcano Style 3
  fc_limit <- quantile(abs(plot_df$avg_log2FC), 0.95, na.rm = TRUE)
  
  p_vol3 <- ggplot(plot_df, aes(avg_log2FC, logFDR)) + 
    geom_point(aes(color = avg_log2FC), size = 1.3, alpha = 0.65) + 
    scale_color_gradient2(low = "#3376B2", mid = "grey88", high = "#D2352C", midpoint = 0, limits = c(-fc_limit, fc_limit), oob = scales::squish) + 
    common_volcano + 
    geom_point(data = labels, shape = 21, size = 2.8, fill = "white", color = "black") + 
    geom_text_repel(data = labels, aes(label = gene), color = "black", size = 3.5, max.overlaps = Inf) + 
    theme_classic()
    
  ggsave(file.path(deg_dir, paste0("Volcano_", target_gene, "_High_vs_Low_Style3_Gradient.pdf")), p_vol3, width = 7, height = 5.5)
  
  # G. Pathway Size & ORA
  rank_genes <- deg$gene
  
  size_summary <- bind_rows(
    lapply(
      names(pathway_db), 
      function(db) {
        sizes <- vapply(pathway_db[[db]], function(x) { sum(x %in% rank_genes) }, integer(1))
        data.frame(
          Database = db, 
          Total = length(sizes), 
          Below10 = sum(sizes < gsea_min_size), 
          Retained_10_500 = sum(sizes >= gsea_min_size & sizes <= gsea_max_size), 
          Above500 = sum(sizes > gsea_max_size), 
          RetainedPercent = round(mean(sizes >= gsea_min_size & sizes <= gsea_max_size) * 100, 2)
        )
      }
    )
  )
  save_csv(size_summary, file.path(gsea_dir, "GSEA_Pathway_Size_Filter_Summary.csv"))
  
  up_df <- plot_df %>% filter(change == "UP")
  save_csv(up_df, file.path(ora_dir, "ORA_Input_UP_DEGs.csv"))
  
  ora_all <- bind_rows(
    lapply(
      names(pathway_db), 
      function(db) {
        res <- enricher(
          gene = up_df$gene, 
          TERM2GENE = to_term2gene(pathway_db[[db]]), 
          pAdjustMethod = "BH", 
          pvalueCutoff = 0.05, 
          qvalueCutoff = 0.25, 
          minGSSize = 10, 
          maxGSSize = 500
        )
        
        df <- as.data.frame(res)
        if (!nrow(df)) return(data.frame())
        
        df$FoldEnrichment <- parse_ratio(df$GeneRatio) / parse_ratio(df$BgRatio)
        
        df %>% 
          filter(FoldEnrichment > 1) %>% 
          mutate(Database = db, .before = 1) %>% 
          arrange(p.adjust, desc(FoldEnrichment))
      }
    )
  )
  save_csv(ora_all, file.path(ora_dir, "ORA_All_Databases.csv"))
  
  ora_summary <- data.frame(Database = names(pathway_db)) %>% 
    left_join(
      ora_all %>% count(Database, name = "Significant_ORA"), 
      by = "Database"
    ) %>% 
    mutate(Significant_ORA = ifelse(is.na(Significant_ORA), 0, Significant_ORA))
    
  save_csv(ora_summary, file.path(ora_dir, "ORA_Subcollection_Summary.csv"))
  
  # H. GSEA 排序与计算
  deg$gsea_metric <- sign(deg$avg_log2FC) * -log10(pmax(deg$p_val, 1e-300))
  
  rank_df <- data.frame(gene = deg$gene, score = deg$gsea_metric) %>% 
    filter(is.finite(score)) %>% 
    arrange(gene, desc(abs(score))) %>% 
    distinct(gene, .keep_all = TRUE)
    
  gene_rank <- sort(setNames(rank_df$score, rank_df$gene), decreasing = TRUE)
  
  save_csv(
    data.frame(Gene = names(gene_rank), GSEA_Metric = gene_rank), 
    file.path(gsea_dir, "GSEA_Gene_Ranking.csv")
  )
  
  gsea_all <- bind_rows(
    lapply(
      names(pathway_db), 
      function(db) {
        set.seed(123)
        res <- tryCatch(
          fgsea(
            pathways = pathway_db[[db]], 
            stats = gene_rank, 
            minSize = gsea_min_size, 
            maxSize = gsea_max_size, 
            eps = 0, 
            scoreType = "std", 
            nPermSimple = 100000, 
            nproc = workers
          ), 
          error = function(e) NULL
        )
        
        if (is.null(res)) return(data.frame())
        
        as.data.frame(res) %>% 
          mutate(Database = db, LeadingEdge_Count = lengths(leadingEdge), .before = 1)
      }
    )
  )
  
  gsea_all <- gsea_all %>% arrange(Database, padj, desc(NES))
  save_csv(flatten_le(gsea_all), file.path(gsea_dir, "GSEA_All_Databases.csv"))
  
  gsea_sig <- gsea_all %>% 
    filter(!is.na(padj), padj < gsea_fdr_explore, NES > 0) %>% 
    arrange(padj, desc(NES))
  save_csv(flatten_le(gsea_sig), file.path(gsea_dir, "GSEA_Significant_Positive.csv"))
  
  gsea_plot <- gsea_all %>% 
    filter(!is.na(padj), padj < gsea_fdr_plot, NES >= gsea_nes_plot, LeadingEdge_Count >= gsea_le_plot) %>% 
    arrange(padj, desc(NES))
  save_csv(flatten_le(gsea_plot), file.path(gsea_dir, "GSEA_SCI_PlotEligible.csv"))
  
  gsea_summary <- gsea_all %>% 
    group_by(Database) %>% 
    summarise(
      Tested = n(), 
      Positive_NES = sum(NES > 0, na.rm = TRUE), 
      FDR025_Positive = sum(!is.na(padj) & padj < gsea_fdr_explore & NES > 0), 
      SCI_PlotEligible = sum(!is.na(padj) & padj < gsea_fdr_plot & NES >= gsea_nes_plot & LeadingEdge_Count >= gsea_le_plot), 
      NA_Padj = sum(is.na(padj)), 
      .groups = "drop"
    )
  save_csv(gsea_summary, file.path(gsea_dir, "GSEA_Subcollection_Summary.csv"))
  
  # I. 绘制全部三种风格的 GSEA 图
  for (db in names(pathway_db)) {
    
    db_res <- gsea_plot %>% filter(Database == db)
    if (!nrow(db_res)) next
    
    db_plot_dir <- file.path(plot_dir, db)
    dir.create(db_plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (i in seq_len(nrow(db_res))) {
      x <- db_res[i, ]
      
      save_gsea_3styles(
        pathway = x$pathway, 
        pathways = pathway_db[[db]], 
        gene_rank = gene_rank, 
        db = db, 
        NES = x$NES, 
        padj = x$padj, 
        LE = x$LeadingEdge_Count, 
        output_dir = db_plot_dir, 
        index = i
      )
    }
  }
}

cat("\n============================================================\n")
cat(">>> 🎉 分析全部完成！MAST 样本效应矫正版本已运行结束。\n")
cat("============================================================\n")