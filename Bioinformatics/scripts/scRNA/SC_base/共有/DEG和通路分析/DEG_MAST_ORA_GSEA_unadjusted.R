# ==============================================================================
# MAST DEG + ORA + GSEA (多细胞类型循环 + 目录全量预建 + GMM 密度阈值图)
# 适配 MSigDB v2026.1.Hs official subcollections
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
library(mixtools)
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
# 1. 全局路径与基础参数配置
# ==============================================================================

target_gene <- "PPAT"

base_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE182109(scRNA)_ndGBM_Core"

work_dir_base <- paste0(
  base_dir, 
  "/通路分析/", target_gene, "/6_DEG差异性分析_富集/DEG_MAST"
)

rds_file <- file.path(base_dir, "Glioma_3B_Fully_Annotated.rds")
gmt_dir  <- file.path(base_dir, "MSigDB_Sub")


# ------------------------------------------------------------------------------
# 🚨 核心配置区：定义 Cutoff 与需要剔除的样本
# ------------------------------------------------------------------------------

cell_type_cutoffs <- list(
  "AC-like"  = 0.15,
  "NPC-like" = 0.15
)

# 填入需要剔除的极少细胞量样本（若不需要剔除，留空即可：c()）
samples_to_exclude <- c("") 


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
# 2. 辅助函数
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
# 4. [PHASE 1] 快速生成全局目录骨架与 GMM 密度图 (前置校验)
# ==============================================================================

cat("\n\n============================================================")
cat("\n>>> [PHASE 1] 启动前置校验：搭建全部目录结构并生成 GMM 密度图...\n")
cat("============================================================\n")

for (ctype in names(cell_type_cutoffs)) {
  
  group_cutoff <- cell_type_cutoffs[[ctype]]
  if (!ctype %in% levels(Idents(sc_obj_global))) {
    next
  }
  
  # A. 一次性创建该亚群下的所有目录结构
  out_dir  <- file.path(work_dir_base, ctype)
  deg_dir  <- file.path(out_dir, "01_DEG_MAST")
  ora_dir  <- file.path(out_dir, "02_ORA_MAST")
  gsea_dir <- file.path(out_dir, "03_GSEA_MAST")
  plot_dir <- file.path(gsea_dir, "Plots")
  
  dir.create(deg_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(ora_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 立刻将参数和 MSigDB Summary 存入对应的文件夹
  save_csv(
    msig_summary, 
    file.path(out_dir, "00_MAST_MSigDB_Subcollection_Summary.csv")
  )
  
  parameters <- data.frame(
    Parameter = c(
      "TargetGene", "CellType", "GroupCutoff", "DEG_Method", 
      "DEG_FDR", "DEG_Log2FC", "GSEA_MinSize", "GSEA_MaxSize", 
      "GSEA_Exploratory_FDR", "GSEA_Plot_FDR", "GSEA_Plot_NES", "GSEA_Plot_LeadingEdge"
    ),
    Value = c(
      target_gene, ctype, group_cutoff, "MAST", 
      deg_fdr, deg_fc, gsea_min_size, gsea_max_size, 
      gsea_fdr_explore, gsea_fdr_plot, gsea_nes_plot, gsea_le_plot
    )
  )
  save_csv(parameters, file.path(out_dir, "00_MAST_Analysis_Parameters.csv"))
  
  # B. 取出亚群并剔除样本
  obj <- subset(sc_obj_global, idents = ctype)
  
  if (length(samples_to_exclude) > 0) {
    obj <- obj[, !obj$SampleID %in% samples_to_exclude]
    cat(sprintf("\n    -> [%s] 目录已建好，并剔除了样本: %s", ctype, paste(samples_to_exclude, collapse = ", ")))
  }
  
  # C. 高低表达分组与 GMM 拟合
  cat(sprintf("\n    -> [%s] 正在生成 GMM 密度图 (Cutoff = %s)...\n", ctype, group_cutoff))
  
  expr <- tryCatch(
    GetAssayData(obj, layer = "data")[target_gene, ],
    error = function(e) {
      GetAssayData(obj, slot = "data")[target_gene, ]
    }
  )
  
  high_name <- paste0(target_gene, "_High")
  low_name  <- paste0(target_gene, "_Low")
  group_col <- paste0(target_gene, "_Group")
  
  obj[[group_col]] <- factor(
    ifelse(expr > group_cutoff, high_name, low_name), 
    levels = c(low_name, high_name)
  )
  
  group_counts <- as.data.frame(table(obj[[group_col]]))
  colnames(group_counts) <- c("Group", "Cells")
  
  # 生成初步的分组计数表
  save_csv(group_counts, file.path(deg_dir, "Cell_Group_Counts.csv"))
  
  # 拟合并作图
  set.seed(123)
  fit_data <- expr + runif(length(expr), 0, 1e-5)
  mixmdl <- tryCatch(normalmixEM(fit_data, k = 2, verb = FALSE), error = function(e) NULL)
  
  n_low  <- group_counts$Cells[group_counts$Group == low_name]
  n_high <- group_counts$Cells[group_counts$Group == high_name]
  
  if (!is.null(mixmdl)) {
    high_comp <- which.max(mixmdl$mu)
    low_comp <- which.min(mixmdl$mu)
    
    mu_high <- mixmdl$mu[high_comp]
    sigma_high <- mixmdl$sigma[high_comp]
    lambda_high <- mixmdl$lambda[high_comp]
    
    mu_low <- mixmdl$mu[low_comp]
    sigma_low <- mixmdl$sigma[low_comp]
    lambda_low <- mixmdl$lambda[low_comp]
    
    limit_y <- max(lambda_high * dnorm(mu_high, mu_high, sigma_high)) * 2
    text_y_pos <- limit_y * 0.9
    
    p_gmm <- ggplot(data.frame(Expression = expr), aes(Expression)) +
      geom_histogram(
        aes(y = after_stat(density)), 
        bins = 60, fill = "grey85", color = "black", alpha = 0.7
      ) +
      stat_function(
        fun = function(x) lambda_low * dnorm(x, mu_low, sigma_low), 
        color = "#3376B2", linewidth = 1.2
      ) +
      stat_function(
        fun = function(x) lambda_high * dnorm(x, mu_high, sigma_high), 
        color = "#D2352C", linewidth = 1.2
      ) +
      geom_vline(
        xintercept = group_cutoff, 
        linetype = "dashed", color = "black", linewidth = 1
      ) +
      annotate(
        "text", x = group_cutoff - 0.05, y = text_y_pos, 
        label = paste0(n_low, " Cells"), 
        color = "#3376B2", fontface = "bold", hjust = 1, size = 5
      ) +
      annotate(
        "text", x = group_cutoff + 0.05, y = text_y_pos, 
        label = paste0(n_high, " Cells"), 
        color = "#D2352C", fontface = "bold", hjust = 0, size = 5
      ) +
      coord_cartesian(ylim = c(0, limit_y)) + 
      theme_classic() +
      labs(
        title = paste0(ctype, " - ", target_gene, " GMM (Cutoff=", group_cutoff, ")"), 
        x = paste0(target_gene, " LogNormalize Expression"), 
        y = "Density"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    ggsave(
      file.path(out_dir, "01_TestCutoff_GMM_Density_Annotated.pdf"), 
      p_gmm, 
      width = 6.5, height = 5
    )
  }
}
cat("\n>>> [PHASE 1 完毕] 结构已就绪，所有 GMM 验证图已生成。MAST 开始进行！\n")


# ==============================================================================
# 5. [PHASE 2] 执行 MAST DEG 与下游富集分析 (后台耗时运算)
# ==============================================================================

cat("\n\n============================================================")
cat("\n>>> [PHASE 2] 启动耗时运算：正式执行 MAST DEG 与 ORA/GSEA...\n")
cat("============================================================\n")

for (ctype in names(cell_type_cutoffs)) {
  
  group_cutoff <- cell_type_cutoffs[[ctype]]
  if (!ctype %in% levels(Idents(sc_obj_global))) {
    next
  }
  
  cat(sprintf("\n>>> 正在处理 [%s] 的 MAST 差异分析...\n", ctype))
  
  # A. 路径重新指派 
  out_dir  <- file.path(work_dir_base, ctype)
  deg_dir  <- file.path(out_dir, "01_DEG_MAST")
  ora_dir  <- file.path(out_dir, "02_ORA_MAST")
  gsea_dir <- file.path(out_dir, "03_GSEA_MAST")
  plot_dir <- file.path(gsea_dir, "Plots")
  
  # B. 再次提取亚群并分组
  obj <- subset(sc_obj_global, idents = ctype)
  
  if (length(samples_to_exclude) > 0) {
    obj <- obj[, !obj$SampleID %in% samples_to_exclude]
  }
  
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
  
  # C. MAST DEG 计算
  plan("multisession", workers = workers)
  options(future.globals.maxSize = 200 * 1024^3)
  
  deg <- FindMarkers(
    obj, 
    ident.1 = high_name, 
    ident.2 = low_name, 
    test.use = "MAST", 
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
  
  # D. 绘制全部三种风格的 Volcano 图
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
  
  # E. Pathway Size & ORA
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
  
  # F. GSEA 排序与计算
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
  
  # G. 绘制全部三种风格的 GSEA 图
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
cat(">>> 🎉 分析全部完成！\n")
cat("============================================================\n")