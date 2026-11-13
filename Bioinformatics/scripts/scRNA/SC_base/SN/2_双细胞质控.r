#!/usr/bin/env Rscript

# ==============================================================================
# snRNA-seq Step 2 教学模板：只检测并去除双细胞（doublet）
# ==============================================================================
#
# 这份脚本故意只完成一件事：
#   1) 用 scDblFinder 判断每个细胞核是 singlet 还是 doublet；
#   2) 删除 scDblFinder 判定的 doublet；
#   3) 保存完整判定对象和最终 singlet 对象。
#
# 本脚本不会再次按 nFeature、nCount、percent.mt 等指标过滤细胞，因为这些工作
# 已经在第一步完成；也不会运行 SoupX、标准化、PCA、UMAP、分簇或细胞注释。
#
# 复用到其他 10x snRNA-seq 数据集时，只需要修改下面“配置区”。输入对象必须
# 保留原始、非负、整数型 UMI counts；TPM、FPKM 或标准化表达矩阵不能用于本模板。


# ==============================================================================
# 一、配置区：复用其他数据集时主要修改这里
# ==============================================================================

DATASET_NAME <- "GSE259378"

INPUT_RDS <- paste0(
  "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE292130(snRNA)Object_1_原发/1_读取与基础质控/GSE211376_merged_clean.rds"
)

OUTPUT_DIR <- paste0(
  "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE292130(snRNA)Object_1_原发/2_双细胞质控"
)

# 该列必须表示独立的 10x 捕获/样本，而不是细胞类型或实验分组。
SAMPLE_COLUMN <- "SampleID"

# 这是第一步记录的各样本细胞数，只用于防止误读其他 RDS，不参与过滤。
# 复用其他数据集时：有可靠记录就填写命名向量；没有就设为 NULL。
EXPECTED_SAMPLE_COUNTS <- NULL

# scDblFinder 使用 BiocParallel 并行。8 核对本数据集已经足够，不建议占满 128 核。
N_WORKERS <- 8L
SEED <- 20260810L

# 所有样本联合训练一个模型，再对每个样本分别确定阈值。
# 这样可以避免 GSM4119524（101 个细胞）单独训练时样本量过小。
MULTI_SAMPLE_MODE <- "singleModelSplitThres"

# FALSE：目标文件已存在时停止，防止误覆盖。
# TRUE ：明确需要重跑时允许覆盖本脚本列出的结果文件。
OVERWRITE <- FALSE


# ==============================================================================
# 二、输出文件：正常情况下无需修改
# ==============================================================================

CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "Checkpoints")
TABLE_DIR <- file.path(OUTPUT_DIR, "Tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "Figures")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

OUTPUT_FILES <- list(
  all_cells_rds = file.path(
    CHECKPOINT_DIR, "01_all_cells_with_doublet_calls.rds"
  ),
  singlets_rds = file.path(
    CHECKPOINT_DIR, "02_singlets_after_doublet_removal.rds"
  ),
  summary_tsv = file.path(
    TABLE_DIR, "01_doublet_summary_by_sample.tsv"
  ),
  all_calls_tsv = file.path(
    TABLE_DIR, "02_all_cell_doublet_calls.tsv"
  ),
  removed_tsv = file.path(
    TABLE_DIR, "03_removed_doublets.tsv"
  ),
  parameters_tsv = file.path(
    TABLE_DIR, "04_run_parameters.tsv"
  ),
  package_versions_tsv = file.path(
    TABLE_DIR, "05_package_versions.tsv"
  ),
  score_pdf = file.path(
    FIGURE_DIR, "01_scDblFinder_score_by_sample.pdf"
  ),
  rate_pdf = file.path(
    FIGURE_DIR, "02_doublet_rate_by_sample.pdf"
  ),
  qc_pdf = file.path(
    FIGURE_DIR, "03_QC_before_after_doublet_removal.pdf"
  ),
  session_info = file.path(
    LOG_DIR, "02_doublet_only_sessionInfo.txt"
  )
)


# ==============================================================================
# 三、通用辅助函数
# ==============================================================================

# 用基础 R 写 TSV，避免为了写表再增加 readr 等依赖。
write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}


# Seurat v5 使用 layer 保存 counts；较旧版本使用 slot。
# 本模板要求 RNA 中只有一个 counts layer，避免猜测或错误合并多个矩阵。
get_counts <- function(object) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    layers <- SeuratObject::Layers(object[["RNA"]])
    count_layers <- layers[layers == "counts" | grepl("^counts\\.", layers)]

    if (length(count_layers) != 1L) {
      stop(
        paste0(
          "RNA assay 必须且只能有一个 counts layer；当前发现：",
          if (length(count_layers) == 0L) {
            "0 个"
          } else {
            paste(count_layers, collapse = ", ")
          }
        ),
        call. = FALSE
      )
    }

    return(
      SeuratObject::LayerData(
        object = object,
        assay = "RNA",
        layer = count_layers[[1L]]
      )
    )
  }

  Seurat::GetAssayData(object = object, assay = "RNA", slot = "counts")
}


# 这些检查只确认输入正确，不会重算指标或过滤细胞。
check_input <- function(object, counts) {
  if (!inherits(object, "Seurat")) {
    stop("输入 RDS 不是 Seurat 对象。", call. = FALSE)
  }
  if (!"RNA" %in% SeuratObject::Assays(object)) {
    stop("输入对象缺少 RNA assay。", call. = FALSE)
  }

  metadata <- object[[]]
  if (!SAMPLE_COLUMN %in% colnames(metadata)) {
    stop(sprintf("metadata 中缺少样本列：%s", SAMPLE_COLUMN), call. = FALSE)
  }

  sample_ids <- as.character(metadata[[SAMPLE_COLUMN]])
  if (anyNA(sample_ids) || any(!nzchar(trimws(sample_ids)))) {
    stop("样本列中存在 NA 或空值。", call. = FALSE)
  }

  if (!inherits(counts, "sparseMatrix")) {
    stop("RNA counts 不是稀疏矩阵。", call. = FALSE)
  }
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop("RNA counts 维度为空。", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("RNA counts 的基因名或 barcode 缺失、为空或重复。", call. = FALSE)
  }
  if (!identical(colnames(counts), rownames(metadata)) ||
      !identical(colnames(counts), colnames(object))) {
    stop("RNA counts 与 Seurat metadata 的 barcode 顺序不一致。", call. = FALSE)
  }

  # dgCMatrix 等常用稀疏矩阵只在 @x 中保存非零值，不需要转成稠密矩阵。
  if (!"x" %in% methods::slotNames(counts)) {
    stop("RNA counts 不是可验证的数值型稀疏矩阵。", call. = FALSE)
  }
  stored_values <- methods::slot(counts, "x")
  if (length(stored_values) == 0L || !all(is.finite(stored_values)) ||
      any(stored_values < 0) ||
      any(abs(stored_values - round(stored_values)) > 1e-8)) {
    stop("RNA counts 必须是非空、有限、非负的整数型 UMI。", call. = FALSE)
  }

  observed_counts <- table(sample_ids)
  if (!is.null(EXPECTED_SAMPLE_COUNTS)) {
    expected <- EXPECTED_SAMPLE_COUNTS
    if (is.null(names(expected)) || anyDuplicated(names(expected))) {
      stop("EXPECTED_SAMPLE_COUNTS 必须是具有唯一名称的向量。", call. = FALSE)
    }

    if (!setequal(names(observed_counts), names(expected))) {
      stop(
        paste0(
          "实际样本集合与 EXPECTED_SAMPLE_COUNTS 不一致。实际：",
          paste(names(observed_counts), collapse = ", "),
          "；预期：", paste(names(expected), collapse = ", ")
        ),
        call. = FALSE
      )
    }

    observed_ordered <- as.integer(observed_counts[names(expected)])
    if (!identical(observed_ordered, as.integer(expected))) {
      mismatch <- names(expected)[observed_ordered != as.integer(expected)]
      stop(
        sprintf(
          "以下样本细胞数与第一步记录不一致：%s",
          paste(mismatch, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  cat("\n输入保护检查通过。逐样本细胞数：\n")
  print(observed_counts)
  cat(sprintf(
    "RNA counts：%d genes x %d cells；非零 UMI 项：%d\n",
    nrow(counts), ncol(counts), Matrix::nnzero(counts)
  ))

  invisible(observed_counts)
}


# scDblFinder 直接使用原始 UMI counts。
# 返回的 SingleCellExperiment 中会增加 scDblFinder.score/class 等列。
run_scdblfinder <- function(object, counts) {
  metadata <- object[[]]
  sample_ids <- as.character(metadata[[SAMPLE_COLUMN]])

  sce_coldata <- S4Vectors::DataFrame(
    sample_ids,
    row.names = colnames(counts)
  )
  colnames(sce_coldata) <- SAMPLE_COLUMN

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = sce_coldata
  )

  detected_cores <- parallel::detectCores(logical = TRUE)
  available_workers <- if (is.na(detected_cores)) {
    N_WORKERS
  } else {
    max(1L, detected_cores - 2L)
  }
  workers <- max(1L, min(as.integer(N_WORKERS), as.integer(available_workers)))

  cat(sprintf("\nscDblFinder 将使用 %d 个 CPU worker。\n", workers))
  cat(sprintf("multiSampleMode = %s\n", MULTI_SAMPLE_MODE))

  set.seed(SEED)
  bpparam <- BiocParallel::MulticoreParam(
    workers = workers,
    RNGseed = SEED,
    progressbar = TRUE
  )
  on.exit(try(BiocParallel::bpstop(bpparam), silent = TRUE), add = TRUE)

  sce <- scDblFinder::scDblFinder(
    sce,
    samples = SAMPLE_COLUMN,
    clusters = FALSE,
    dbr = NULL,
    multiSampleMode = MULTI_SAMPLE_MODE,
    BPPARAM = bpparam
  )

  sce
}


# 通过 barcode 而不是当前位置回填结果，防止顺序变化造成错配。
attach_doublet_metadata <- function(object, sce) {
  result_metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  result_metadata$cell_barcode <- rownames(result_metadata)

  sc_columns <- grep(
    "^scDblFinder\\.",
    colnames(result_metadata),
    value = TRUE
  )
  required_columns <- c("scDblFinder.score", "scDblFinder.class")
  if (!all(required_columns %in% sc_columns)) {
    stop("scDblFinder 结果中缺少 score 或 class。", call. = FALSE)
  }

  if (anyDuplicated(result_metadata$cell_barcode) ||
      !setequal(result_metadata$cell_barcode, colnames(object))) {
    stop("scDblFinder 返回的 barcode 与 Seurat 对象不一致。", call. = FALSE)
  }

  match_index <- match(colnames(object), result_metadata$cell_barcode)
  aligned_metadata <- result_metadata[match_index, sc_columns, drop = FALSE]
  rownames(aligned_metadata) <- colnames(object)

  if (anyNA(aligned_metadata$scDblFinder.score) ||
      any(!is.finite(aligned_metadata$scDblFinder.score))) {
    stop("scDblFinder.score 存在 NA 或非有限值。", call. = FALSE)
  }

  classes <- as.character(aligned_metadata$scDblFinder.class)
  if (anyNA(classes) || !all(classes %in% c("singlet", "doublet"))) {
    stop("scDblFinder.class 出现未知类别。", call. = FALSE)
  }

  SeuratObject::AddMetaData(object, metadata = aligned_metadata)
}


# 生成逐样本统计，并追加 TOTAL 行。
make_summary <- function(object) {
  metadata <- object[[]]
  sample_ids <- as.character(metadata[[SAMPLE_COLUMN]])
  classes <- as.character(metadata$scDblFinder.class)
  sample_levels <- sort(unique(sample_ids))

  rows <- lapply(sample_levels, function(sample_id) {
    selected <- sample_ids == sample_id
    n_before <- sum(selected)
    n_singlet <- sum(classes[selected] == "singlet")
    n_doublet <- sum(classes[selected] == "doublet")

    data.frame(
      SampleID = sample_id,
      n_before = n_before,
      n_singlet = n_singlet,
      n_doublet = n_doublet,
      doublet_rate = n_doublet / n_before,
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, rows)
  total_before <- sum(summary_table$n_before)
  total_singlet <- sum(summary_table$n_singlet)
  total_doublet <- sum(summary_table$n_doublet)

  total_row <- data.frame(
    SampleID = "TOTAL",
    n_before = total_before,
    n_singlet = total_singlet,
    n_doublet = total_doublet,
    doublet_rate = total_doublet / total_before,
    stringsAsFactors = FALSE
  )

  summary_table <- rbind(summary_table, total_row)

  if (any(summary_table$n_before !=
          summary_table$n_singlet + summary_table$n_doublet)) {
    stop("细胞守恒检查失败：before 不等于 singlet + doublet。", call. = FALSE)
  }

  rownames(summary_table) <- NULL
  summary_table
}


# 将 QC metadata 整理成长表，仅用于绘制“去双细胞前后”的分布对比。
make_qc_long_table <- function(before_object, after_object) {
  qc_features <- intersect(
    c(
      "nFeature_RNA", "nCount_RNA", "percent.mt",
      "percent.rb", "log10GenesPerUMI"
    ),
    colnames(before_object[[]])
  )

  before_metadata <- before_object[[]]
  after_metadata <- after_object[[]]

  parts <- lapply(qc_features, function(feature) {
    rbind(
      data.frame(
        stage = "Before doublet removal",
        SampleID = as.character(before_metadata[[SAMPLE_COLUMN]]),
        metric = feature,
        value = as.numeric(before_metadata[[feature]]),
        stringsAsFactors = FALSE
      ),
      data.frame(
        stage = "After doublet removal",
        SampleID = as.character(after_metadata[[SAMPLE_COLUMN]]),
        metric = feature,
        value = as.numeric(after_metadata[[feature]]),
        stringsAsFactors = FALSE
      )
    )
  })

  do.call(rbind, parts)
}


# ==============================================================================
# 四、主流程
# ==============================================================================

main <- function() {
  cat("============================================================\n")
  cat(sprintf("%s snRNA-seq Step 2：只检测并去除 doublet\n", DATASET_NAME))
  cat("============================================================\n")

  # --------------------------------------------------------------------------
  # 阶段0：检查运行环境和输出目录
  # --------------------------------------------------------------------------
  required_packages <- c(
    "Seurat", "SeuratObject", "Matrix", "SingleCellExperiment",
    "SummarizedExperiment", "S4Vectors", "BiocParallel",
    "scDblFinder", "ggplot2", "xgboost"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      sprintf("缺少必需 R 包：%s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!identical(Sys.info()[["sysname"]], "Linux")) {
    stop("本模板仅支持 WSL/Linux。", call. = FALSE)
  }
  if (!isTRUE(capabilities("cairo"))) {
    stop("当前 R 不支持 Cairo PDF。", call. = FALSE)
  }

  # scDblFinder 官方提示：xgboost >= 3.1 时，Bioc 3.22 应使用 > 1.24.7。
  if (utils::packageVersion("xgboost") >= "3.1.0" &&
      utils::packageVersion("scDblFinder") <= "1.24.7") {
    stop(
      "xgboost >= 3.1 需要 scDblFinder > 1.24.7；请先升级 scDblFinder。",
      call. = FALSE
    )
  }

  if (!file.exists(INPUT_RDS)) {
    stop(sprintf("找不到输入文件：%s", INPUT_RDS), call. = FALSE)
  }

  for (directory in c(CHECKPOINT_DIR, TABLE_DIR, FIGURE_DIR, LOG_DIR)) {
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
  }

  existing_outputs <- unlist(OUTPUT_FILES, use.names = FALSE)
  existing_outputs <- existing_outputs[file.exists(existing_outputs)]
  if (length(existing_outputs) > 0L && !isTRUE(OVERWRITE)) {
    stop(
      paste0(
        "以下结果已存在。若确定重跑，请把 OVERWRITE 改为 TRUE：\n",
        paste(existing_outputs, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # 阶段1：读取第一步clean对象；这里只做防读错检查，不进行二次QC
  # --------------------------------------------------------------------------
  cat("\n[阶段1/4] 读取第一步 clean RDS...\n")
  object <- readRDS(INPUT_RDS)
  
  # ---------- 新增：自动检测并合并 Layers ----------
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    current_layers <- SeuratObject::Layers(object[["RNA"]])
    count_layers <- current_layers[current_layers == "counts" | grepl("^counts\\.", current_layers)]
    
    if (length(count_layers) > 1L) {
      cat("-> 检测到", length(count_layers), "个独立的 counts 层，正在自动合并...\n")
      object <- SeuratObject::JoinLayers(object)
    } else {
      cat("-> 检测到单一 counts 层，无需合并。\n")
    }
  }
  # ---------- 新增代码结束 ----------
  
  if (!inherits(object, "Seurat") || !"RNA" %in% SeuratObject::Assays(object)) {
    stop("输入必须是含 RNA assay 的 Seurat 对象。", call. = FALSE)
  }
  SeuratObject::DefaultAssay(object) <- "RNA"
  counts <- get_counts(object)
  check_input(object, counts)
  n_input <- ncol(object)



  # --------------------------------------------------------------------------
  # 阶段2：运行scDblFinder并回填所有结果列
  # --------------------------------------------------------------------------
  cat("\n[阶段2/4] 运行 scDblFinder...\n")
  sce_result <- run_scdblfinder(object, counts)
  all_cells <- attach_doublet_metadata(object, sce_result)
  summary_table <- make_summary(all_cells)

  cat("\n逐样本 doublet 统计：\n")
  print(summary_table)

  # 保存检查点1：保留所有细胞，便于复核算法判定。
  saveRDS(all_cells, OUTPUT_FILES$all_cells_rds)
  write_tsv(summary_table, OUTPUT_FILES$summary_tsv)

  all_metadata <- all_cells[[]]
  sc_columns <- grep("^scDblFinder\\.", colnames(all_metadata), value = TRUE)
  all_calls <- data.frame(
    cell_barcode = rownames(all_metadata),
    SampleID = as.character(all_metadata[[SAMPLE_COLUMN]]),
    all_metadata[, sc_columns, drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_tsv(all_calls, OUTPUT_FILES$all_calls_tsv)

  # 图1：每个样本的score分布，并标出算法判定的doublet。
  score_plot_data <- data.frame(
    SampleID = as.character(all_metadata[[SAMPLE_COLUMN]]),
    score = as.numeric(all_metadata$scDblFinder.score),
    class = as.character(all_metadata$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  doublet_points <- score_plot_data[
    score_plot_data$class == "doublet", , drop = FALSE
  ]
  score_plot <- ggplot2::ggplot(
    score_plot_data,
    ggplot2::aes(x = SampleID, y = score)
  ) +
    ggplot2::geom_violin(fill = "#9ECAE1", color = "#3182BD", scale = "width") +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    ggplot2::geom_point(
      data = doublet_points,
      color = "#D7301F",
      size = 0.45,
      alpha = 0.55,
      position = ggplot2::position_jitter(width = 0.12, height = 0)
    ) +
    ggplot2::labs(
      title = "scDblFinder score by sample",
      subtitle = "Red points are predicted doublets",
      x = "SampleID",
      y = "scDblFinder score"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(
    OUTPUT_FILES$score_pdf,
    plot = score_plot,
    device = grDevices::cairo_pdf,
    width = 11,
    height = 6
  )

  # 图2：逐样本doublet比例。TOTAL行不参与绘图。
  rate_plot_data <- summary_table[
    summary_table$SampleID != "TOTAL", , drop = FALSE
  ]
  rate_plot <- ggplot2::ggplot(
    rate_plot_data,
    ggplot2::aes(x = SampleID, y = doublet_rate)
  ) +
    ggplot2::geom_col(fill = "#756BB1") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f%%", 100 * doublet_rate)),
      vjust = -0.35,
      size = 3.3
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.1f%%", 100 * x),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = "Predicted doublet rate by sample",
      x = "SampleID",
      y = "Doublet rate"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(
    OUTPUT_FILES$rate_pdf,
    plot = rate_plot,
    device = grDevices::cairo_pdf,
    width = 10,
    height = 6
  )

  cat(sprintf(
    "检查点1已保存：%s（%d个细胞，尚未删除doublet）\n",
    OUTPUT_FILES$all_cells_rds, ncol(all_cells)
  ))

  # --------------------------------------------------------------------------
  # 阶段3：唯一一次细胞过滤——只保留scDblFinder singlet
  # --------------------------------------------------------------------------
  cat("\n[阶段3/4] 删除 scDblFinder 判定的 doublet...\n")
  singlet_cells <- rownames(all_metadata)[
    all_metadata$scDblFinder.class == "singlet"
  ]
  singlets <- subset(all_cells, cells = singlet_cells)

  removed_doublets <- all_calls[
    as.character(all_calls$scDblFinder.class) == "doublet", , drop = FALSE
  ]
  write_tsv(removed_doublets, OUTPUT_FILES$removed_tsv)

  if (n_input != ncol(singlets) + nrow(removed_doublets)) {
    stop("全局细胞守恒检查失败。", call. = FALSE)
  }
  if (!all(as.character(singlets$scDblFinder.class) == "singlet")) {
    stop("最终对象中仍存在非singlet类别。", call. = FALSE)
  }

  # 图3仅展示去doublet前后原有QC指标分布，不产生新的过滤条件。
  qc_long <- make_qc_long_table(all_cells, singlets)
  qc_plot <- ggplot2::ggplot(
    qc_long,
    ggplot2::aes(x = stage, y = value, fill = stage)
  ) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    ggplot2::facet_grid(metric ~ SampleID, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(
      "Before doublet removal" = "#9ECAE1",
      "After doublet removal" = "#A1D99B"
    )) +
    ggplot2::labs(
      title = "QC distributions before and after doublet removal",
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  ggplot2::ggsave(
    OUTPUT_FILES$qc_pdf,
    plot = qc_plot,
    device = grDevices::cairo_pdf,
    width = 18,
    height = 12
  )

  # 检查点2就是本步骤最终对象：仅删除doublet，其余细胞全部保留。
  saveRDS(singlets, OUTPUT_FILES$singlets_rds)

  # --------------------------------------------------------------------------
  # 阶段4：记录参数和软件版本，打印最终守恒关系
  # --------------------------------------------------------------------------
  cat("\n[阶段4/4] 保存参数、版本和最终统计...\n")

  parameter_table <- data.frame(
    parameter = c(
      "DATASET_NAME", "INPUT_RDS", "OUTPUT_DIR", "SAMPLE_COLUMN",
      "N_WORKERS_requested", "SEED", "MULTI_SAMPLE_MODE", "OVERWRITE"
    ),
    value = c(
      DATASET_NAME, INPUT_RDS, OUTPUT_DIR, SAMPLE_COLUMN,
      as.character(N_WORKERS), as.character(SEED),
      MULTI_SAMPLE_MODE, as.character(OVERWRITE)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(parameter_table, OUTPUT_FILES$parameters_tsv)

  package_table <- data.frame(
    package = required_packages,
    version = vapply(
      required_packages,
      function(package) as.character(utils::packageVersion(package)),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(package_table, OUTPUT_FILES$package_versions_tsv)
  writeLines(
    capture.output(utils::sessionInfo()),
    con = OUTPUT_FILES$session_info,
    useBytes = TRUE
  )

  total_row <- summary_table[summary_table$SampleID == "TOTAL", , drop = FALSE]
  cat("\n============================================================\n")
  cat("双细胞去除完成\n")
  cat(sprintf(
    "细胞守恒：%d input = %d singlets + %d doublets\n",
    total_row$n_before, total_row$n_singlet, total_row$n_doublet
  ))
  cat(sprintf("最终 singlet RDS：%s\n", OUTPUT_FILES$singlets_rds))
  cat(sprintf("结果目录：%s\n", OUTPUT_DIR))
  cat("本步骤未执行其他细胞过滤、分簇或注释。\n")
  cat("============================================================\n")

  invisible(singlets)
}


main()
