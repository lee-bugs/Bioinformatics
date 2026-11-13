#!/usr/bin/env Rscript

# ==============================================================================
# GSE223063：RNA 污染、双细胞、濒死细胞与重分簇质控（服务器安全版）
#
# 固定流程：
#   原始 RNA counts
#     -> 每个独立捕获批次运行 DecontX
#     -> 使用 DecontX 校正 counts 按捕获批次运行 scDblFinder
#     -> 细胞级联合过滤
#     -> 使用校正 counts 保守重分簇
#
# 重要约束：
#   1. 不覆盖输入 RDS；每一步都有检查点。
#   2. SampleID 必须代表独立的上机捕获/芯片孔，而不是任意生物学分组。
#   3. 兼容服务器已经安装的 scDblFinder 1.18.0，因此使用 dbr = NULL。
#   4. 不自动删除任何 cluster；只输出需要复核的候选 cluster。
#   5. 不从高变基因中删除细胞周期、应激或髓鞘谱系基因。
# ==============================================================================

options(stringsAsFactors = FALSE)


# ==============================================================================
# 0. 用户配置：正常情况下只需要检查本区域
# ==============================================================================

INPUT_RDS <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/1_读取与基础质控/Glioma_merged_clean.rds"
WORK_DIR <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/2_RNA污染双细胞与濒死细胞质控"
SAMPLE_COLUMN <- "SampleID"

SEED <- 20260809L
N_WORKERS_SCDOUBLET <- 6L
FORCE_RERUN <- identical(Sys.getenv("GLIOMA_QC_FORCE_RERUN", "0"), "1")
PIPELINE_SCHEMA_VERSION <- "2026-08-09.2"
INPUT_FINGERPRINT <- NA_character_

# DecontX 后的细胞级过滤阈值。
LOW_COUNT_DECONTX <- 500
LOW_FEATURE_DECONTX <- 200
MAX_PERCENT_MT_DECONTX <- 10
EXTREME_CONTAMINATION <- 0.80
EXTREME_CONTAMINATION_MAX_FEATURE <- 500

# 重分簇参数。
HVG_CANDIDATES <- 4000L
HVG_FINAL <- 2000L
N_PCS <- 30L
CLUSTER_RESOLUTION <- 0.5

# 所有输出均位于 WORK_DIR 内。
DIR_LOG <- file.path(WORK_DIR, "00_logs")
DIR_DECONTX <- file.path(WORK_DIR, "01_DecontX")
DIR_DOUBLET <- file.path(WORK_DIR, "02_scDblFinder")
DIR_FILTER <- file.path(WORK_DIR, "03_Cell_Filter")
DIR_RECLUSTER <- file.path(WORK_DIR, "04_Recluster")
DIR_FIGURE <- file.path(WORK_DIR, "Figures_QC")
DIR_CHECKPOINT <- file.path(WORK_DIR, "checkpoints")

RDS_DECONTX <- file.path(WORK_DIR, "Glioma_2A_DecontX.rds")
RDS_DOUBLET <- file.path(WORK_DIR, "Glioma_2B_scDblFinder.rds")
RDS_FILTERED <- file.path(WORK_DIR, "Glioma_2C_QC_Filtered.rds")
RDS_FINAL <- file.path(WORK_DIR, "Glioma_2D_QC_Reclustered.rds")


# ==============================================================================
# 1. 纯函数：过滤决策、黑名单和统计（可独立测试）
# ==============================================================================

normalize_doublet_class <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value
}


make_cell_qc_decision <- function(
    metadata,
    low_count = LOW_COUNT_DECONTX,
    low_feature = LOW_FEATURE_DECONTX,
    max_percent_mt = MAX_PERCENT_MT_DECONTX,
    extreme_contamination = EXTREME_CONTAMINATION,
    extreme_contamination_max_feature = EXTREME_CONTAMINATION_MAX_FEATURE) {
  
  required <- c(
    SAMPLE_COLUMN, "nCount_decontX", "nFeature_decontX",
    "percent_mt_decontX", "decontX_contamination", "scDblFinder.class"
  )
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) {
    stop("质控决策缺少 metadata 列：", paste(missing, collapse = ", "))
  }
  if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("metadata 必须以唯一细胞条形码作为行名。")
  }
  
  clean_count <- as.numeric(metadata$nCount_decontX)
  clean_feature <- as.numeric(metadata$nFeature_decontX)
  clean_mt <- as.numeric(metadata$percent_mt_decontX)
  contamination <- as.numeric(metadata$decontX_contamination)
  doublet_class <- normalize_doublet_class(metadata$scDblFinder.class)
  
  remove_low_count <- !is.finite(clean_count) | clean_count < low_count
  remove_low_feature <- !is.finite(clean_feature) | clean_feature < low_feature
  remove_high_mito <- !is.finite(clean_mt) | clean_mt >= max_percent_mt
  remove_extreme_contamination <-
    is.finite(contamination) &
    contamination > extreme_contamination &
    is.finite(clean_feature) &
    clean_feature < extreme_contamination_max_feature
  remove_doublet <- doublet_class == "doublet"
  
  final_remove <-
    remove_low_count |
    remove_low_feature |
    remove_high_mito |
    remove_extreme_contamination |
    remove_doublet
  
  reason_matrix <- cbind(
    LowCount = remove_low_count,
    LowFeature = remove_low_feature,
    HighMito = remove_high_mito,
    ExtremeContamination = remove_extreme_contamination,
    Doublet = remove_doublet
  )
  reason <- apply(reason_matrix, 1L, function(x) {
    hit <- colnames(reason_matrix)[x]
    if (!length(hit)) "Keep" else paste0("Remove_", paste(hit, collapse = "+"))
  })
  
  data.frame(
    SampleID = as.character(metadata[[SAMPLE_COLUMN]]),
    Remove_LowCount = as.logical(remove_low_count),
    Remove_LowFeature = as.logical(remove_low_feature),
    Remove_HighMito = as.logical(remove_high_mito),
    Remove_ExtremeContamination = as.logical(remove_extreme_contamination),
    Remove_Doublet = as.logical(remove_doublet),
    Final_Remove = as.logical(final_remove),
    Final_QC_Action = reason,
    row.names = rownames(metadata),
    check.names = FALSE
  )
}


summarize_cell_filtering <- function(decision) {
  required <- c(
    "SampleID", "Remove_LowCount", "Remove_LowFeature", "Remove_HighMito",
    "Remove_ExtremeContamination", "Remove_Doublet", "Final_Remove"
  )
  missing <- setdiff(required, colnames(decision))
  if (length(missing)) stop("过滤决策表缺少列：", paste(missing, collapse = ", "))
  
  samples <- unique(as.character(decision$SampleID))
  rows <- lapply(samples, function(sample_id) {
    x <- decision[decision$SampleID == sample_id, , drop = FALSE]
    data.frame(
      SampleID = sample_id,
      n_before = nrow(x),
      n_removed_low_count = sum(x$Remove_LowCount),
      n_removed_low_feature = sum(x$Remove_LowFeature),
      n_removed_high_mito = sum(x$Remove_HighMito),
      n_removed_extreme_contamination = sum(x$Remove_ExtremeContamination),
      n_removed_doublet = sum(x$Remove_Doublet),
      n_removed_union = sum(x$Final_Remove),
      n_after = sum(!x$Final_Remove),
      retention_rate = mean(!x$Final_Remove),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}


technical_gene_blacklist <- function(genes) {
  genes <- as.character(genes)
  is_mito <- grepl("^MT-", genes, ignore.case = TRUE)
  is_rpl <- grepl("^RPL(P[0-9]+|[0-9]+[A-Z0-9]*)$", genes, ignore.case = TRUE)
  is_rps <- grepl("^RPS([0-9]+[A-Z0-9]*|A)$", genes, ignore.case = TRUE) &
    !grepl("^RPS6K", genes, ignore.case = TRUE)
  is_hb <- grepl(
    "^HB(A1|A2|B|D|E1|G1|G2|M|Q1|Z)$",
    genes,
    ignore.case = TRUE
  )
  genes[is_mito | is_rpl | is_rps | is_hb]
}


make_stage_signature <- function(stage, parameters) {
  values <- unlist(parameters, recursive = TRUE, use.names = TRUE)
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("检查点参数必须全部具有名称。")
  }
  values <- values[order(names(values))]
  paste(
    c(
      paste0("stage=", stage),
      paste0(names(values), "=", as.character(values))
    ),
    collapse = "|"
  )
}


checkpoint_matches <- function(
    checkpoint,
    expected_cells,
    expected_signature,
    cell_extractor,
    signature_extractor) {
  cells <- tryCatch(cell_extractor(checkpoint), error = function(e) NULL)
  signature <- tryCatch(signature_extractor(checkpoint), error = function(e) NULL)
  identical(as.character(cells), as.character(expected_cells)) &&
    identical(as.character(signature), as.character(expected_signature))
}


set_object_stage_signature <- function(object, stage, signature) {
  signatures <- object@misc$GSE223063_QC_stage_signatures
  if (is.null(signatures)) signatures <- list()
  signatures[[stage]] <- signature
  object@misc$GSE223063_QC_stage_signatures <- signatures
  object
}


get_object_stage_signature <- function(object, stage) {
  signatures <- tryCatch(
    object@misc$GSE223063_QC_stage_signatures,
    error = function(e) NULL
  )
  if (is.null(signatures) || is.null(signatures[[stage]])) NULL else signatures[[stage]]
}


seurat_checkpoint_matches <- function(object, expected_cells, stage, signature) {
  if (!inherits(object, "Seurat") || !SAMPLE_COLUMN %in% colnames(object[[]])) {
    return(FALSE)
  }
  checkpoint_matches(
    object,
    expected_cells = expected_cells,
    expected_signature = signature,
    cell_extractor = colnames,
    signature_extractor = function(x) get_object_stage_signature(x, stage)
  )
}


seurat_checkpoint_matches_samples <- function(
    object, expected_cells, expected_samples, stage, signature) {
  seurat_checkpoint_matches(object, expected_cells, stage, signature) &&
    identical(
      as.character(object[[]][[SAMPLE_COLUMN]]),
      as.character(expected_samples)
    )
}


# ==============================================================================
# 2. 通用输入输出与安全检查
# ==============================================================================

safe_sample_name <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) stop("样本名必须是非空单值。")
  sanitized <- gsub("[^A-Za-z0-9_-]", "_", x)
  codepoints <- utf8ToInt(enc2utf8(x))
  weights <- seq_along(codepoints) %% 1009L + 1L
  hash_value <- sum((as.double(codepoints) * weights) %% 2147483629) %% 2147483629
  paste0(sanitized, "_", sprintf("%08x", as.integer(hash_value)))
}


combine_named_vectors <- function(vectors) {
  if (!length(vectors)) stop("没有可合并的向量。")
  cell_names <- unlist(lapply(unname(vectors), names), use.names = FALSE)
  if (length(cell_names) != sum(lengths(vectors)) || anyNA(cell_names) ||
      any(!nzchar(cell_names)) || anyDuplicated(cell_names)) {
    stop("待合并向量的细胞条形码缺失或重复。")
  }
  values <- unlist(lapply(unname(vectors), unname), use.names = FALSE)
  names(values) <- cell_names
  values
}


validate_scdblfinder_result <- function(result, expected_cells) {
  required <- c("scDblFinder.score", "scDblFinder.class")
  missing <- setdiff(required, colnames(result))
  if (length(missing)) stop("scDblFinder 结果缺少：", paste(missing, collapse = ", "))
  if (is.null(rownames(result)) || anyDuplicated(rownames(result))) {
    stop("scDblFinder 结果的细胞条形码缺失或重复。")
  }
  if (!setequal(rownames(result), expected_cells)) {
    stop(
      "scDblFinder 条形码集合不一致；missing=",
      length(setdiff(expected_cells, rownames(result))),
      ", unexpected=", length(setdiff(rownames(result), expected_cells))
    )
  }
  if (any(!is.finite(as.numeric(result$scDblFinder.score)))) {
    stop("scDblFinder 返回非有限 score。")
  }
  classes <- normalize_doublet_class(result$scDblFinder.class)
  if (anyNA(classes) || any(!classes %in% c("singlet", "doublet"))) {
    stop("scDblFinder 返回 singlet/doublet 以外的类别。")
  }
  TRUE
}


ensure_output_directories <- function() {
  directories <- c(
    WORK_DIR, DIR_LOG, DIR_DECONTX, DIR_DOUBLET,
    DIR_FILTER, DIR_RECLUSTER, DIR_FIGURE, DIR_CHECKPOINT
  )
  for (path in directories) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) stop("无法创建输出目录：", path)
  }
  invisible(directories)
}


require_pipeline_packages <- function() {
  packages <- c(
    "Seurat", "SeuratObject", "Matrix", "decontX", "scDblFinder",
    "SingleCellExperiment", "SummarizedExperiment", "BiocParallel", "ggplot2"
  )
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (!all(installed)) {
    stop("缺少必要 R 包：", paste(packages[!installed], collapse = ", "))
  }
  
  if (utils::packageVersion("scDblFinder") < package_version("1.18.0")) {
    stop("scDblFinder 版本过低；当前脚本至少需要 1.18.0。")
  }
  if (utils::packageVersion("decontX") < package_version("1.2.0")) {
    stop("decontX 版本过低；当前脚本至少需要 1.2.0。")
  }
  if (utils::packageVersion("Seurat") < package_version("5.0.0") ||
      utils::packageVersion("SeuratObject") < package_version("5.0.0")) {
    stop("当前脚本使用 Seurat 5 layer API，需要 Seurat 与 SeuratObject >= 5.0.0。")
  }
  
  versions <- vapply(
    packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  )
  data.frame(package = names(versions), version = unname(versions))
}


get_assay_counts <- function(object, assay) {
  if (!assay %in% names(object@assays)) stop("对象缺少 assay：", assay)
  
  assay_object <- object[[assay]]
  layer_names <- SeuratObject::Layers(assay_object, search = "^counts")
  if (length(layer_names) != 1L || !identical(layer_names[[1L]], "counts")) {
    object <- SeuratObject::JoinLayers(object, assay = assay)
  }
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  list(object = object, counts = methods::as(counts, "dgCMatrix"))
}


validate_pipeline_input <- function(object) {
  if (!inherits(object, "Seurat")) stop("输入 RDS 不是 Seurat 对象。")
  metadata <- object[[]]
  if (!SAMPLE_COLUMN %in% colnames(metadata)) {
    stop("输入对象缺少样本列：", SAMPLE_COLUMN)
  }
  
  sample_values <- as.character(metadata[[SAMPLE_COLUMN]])
  if (anyNA(sample_values) || any(!nzchar(trimws(sample_values)))) {
    stop(SAMPLE_COLUMN, " 中存在 NA 或空值。")
  }
  if (anyDuplicated(colnames(object))) stop("Seurat 对象存在重复细胞条形码。")
  if (!"RNA" %in% names(object@assays)) stop("输入对象缺少 RNA assay。")
  
  extracted <- get_assay_counts(object, "RNA")
  counts <- extracted$counts
  if (!identical(colnames(counts), colnames(extracted$object))) {
    stop("RNA counts 的细胞顺序与 Seurat 对象不一致。")
  }
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("RNA counts 含有非有限值或负值。")
  }
  if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop("RNA counts 含有非整数值；请输入未经标准化的原始 UMI counts。")
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) ||
      any(!nzchar(rownames(counts))) || anyDuplicated(rownames(counts))) {
    stop("RNA counts 的基因名缺失、为空或重复。")
  }
  if (nrow(counts) < 500L || ncol(counts) < 20L) {
    stop("输入矩阵规模异常：", nrow(counts), " genes x ", ncol(counts), " cells。")
  }
  sample_sizes <- table(sample_values)
  if (any(sample_sizes < 20L)) {
    stop("以下捕获批次少于 20 个细胞：", paste(names(sample_sizes)[sample_sizes < 20L], collapse = ", "))
  }
  
  invisible(list(
    object = extracted$object,
    counts = counts,
    samples = unique(sample_values)
  ))
}


round_sparse_counts <- function(x) {
  x <- methods::as(x, "dgCMatrix")
  x@x <- pmax(0, round(x@x))
  Matrix::drop0(x)
}


atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  
  moved_previous <- FALSE
  previous <- paste0(path, ".previous")
  if (file.exists(path)) {
    if (file.exists(previous)) unlink(previous)
    if (!file.rename(path, previous)) stop("无法备份原输出：", path)
    moved_previous <- TRUE
  }
  if (!file.rename(temporary, path)) {
    if (moved_previous && file.exists(previous)) file.rename(previous, path)
    stop("无法安装新 RDS；已尝试恢复上一版本：", path)
  }
  invisible(path)
}


safe_read_rds <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      warning("检查点损坏或不可读，将重新计算：", path, "；", conditionMessage(e))
      NULL
    }
  )
}


write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(path)
}


save_pdf <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = width, height = height, family = "sans")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}


calculate_percent <- function(counts, pattern) {
  selected <- grepl(pattern, rownames(counts), ignore.case = TRUE)
  total <- Matrix::colSums(counts)
  selected_count <- if (any(selected)) Matrix::colSums(counts[selected, , drop = FALSE]) else rep(0, ncol(counts))
  result <- ifelse(total > 0, 100 * selected_count / total, NA_real_)
  names(result) <- colnames(counts)
  result
}


write_configuration <- function(package_versions) {
  config <- data.frame(
    parameter = c(
      "INPUT_RDS", "WORK_DIR", "SAMPLE_COLUMN", "SEED",
      "PIPELINE_SCHEMA_VERSION", "INPUT_FINGERPRINT",
      "LOW_COUNT_DECONTX", "LOW_FEATURE_DECONTX", "MAX_PERCENT_MT_DECONTX",
      "EXTREME_CONTAMINATION", "EXTREME_CONTAMINATION_MAX_FEATURE",
      "HVG_CANDIDATES", "HVG_FINAL", "N_PCS", "CLUSTER_RESOLUTION",
      "FORCE_RERUN"
    ),
    value = as.character(c(
      INPUT_RDS, WORK_DIR, SAMPLE_COLUMN, SEED,
      PIPELINE_SCHEMA_VERSION, INPUT_FINGERPRINT,
      LOW_COUNT_DECONTX, LOW_FEATURE_DECONTX, MAX_PERCENT_MT_DECONTX,
      EXTREME_CONTAMINATION, EXTREME_CONTAMINATION_MAX_FEATURE,
      HVG_CANDIDATES, HVG_FINAL, N_PCS, CLUSTER_RESOLUTION,
      FORCE_RERUN
    )),
    row.names = NULL
  )
  write_tsv(config, file.path(DIR_LOG, "00_pipeline_configuration.tsv"))
  write_tsv(package_versions, file.path(DIR_LOG, "00_package_versions.tsv"))
}


# ==============================================================================
# 3. DecontX：每个独立捕获批次单独校正环境 RNA
# ==============================================================================

make_decontx_preclusters <- function(counts, sample_id) {
  if (ncol(counts) < 50L) {
    warning("样本 ", sample_id, " 少于 50 个细胞；DecontX 将使用内部聚类。")
    return(NULL)
  }
  
  temporary <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = paste0("DecontX_", safe_sample_name(sample_id)),
    min.cells = 0,
    min.features = 0
  )
  temporary <- Seurat::NormalizeData(temporary, verbose = FALSE)
  temporary <- Seurat::FindVariableFeatures(
    temporary,
    selection.method = "vst",
    nfeatures = min(2000L, nrow(temporary)),
    verbose = FALSE
  )
  variable_features <- SeuratObject::VariableFeatures(temporary)
  max_pcs <- min(20L, length(variable_features) - 1L, ncol(temporary) - 1L)
  if (length(variable_features) < 50L || max_pcs < 5L) return(NULL)
  
  temporary <- Seurat::ScaleData(
    temporary,
    features = variable_features,
    verbose = FALSE
  )
  temporary <- Seurat::RunPCA(
    temporary,
    features = variable_features,
    npcs = max_pcs,
    seed.use = SEED,
    verbose = FALSE
  )
  temporary <- Seurat::FindNeighbors(
    temporary,
    dims = seq_len(max_pcs),
    verbose = FALSE
  )
  temporary <- Seurat::FindClusters(
    temporary,
    resolution = 0.4,
    algorithm = 1,
    random.seed = SEED,
    verbose = FALSE
  )
  
  z <- as.integer(factor(as.character(SeuratObject::Idents(temporary))))
  if (length(unique(z)) < 2L) NULL else z
}


run_one_decontx_sample <- function(counts, sample_id) {
  checkpoint <- file.path(
    DIR_CHECKPOINT,
    paste0("DecontX_", safe_sample_name(sample_id), ".rds")
  )
  signature <- make_stage_signature(
    "2A_DecontX_sample",
    c(
      seed = SEED,
      schema = PIPELINE_SCHEMA_VERSION,
      input_fingerprint = INPUT_FINGERPRINT,
      sample_column = SAMPLE_COLUMN,
      decontX_version = as.character(utils::packageVersion("decontX"))
    )
  )
  
  if (!FORCE_RERUN && file.exists(checkpoint)) {
    cached <- safe_read_rds(checkpoint)
    required <- c(
      "sample_id", "cells", "genes", "corrected_counts",
      "contamination", "signature"
    )
    if (is.list(cached) && all(required %in% names(cached)) &&
        identical(cached$sample_id, sample_id) &&
        identical(cached$cells, colnames(counts)) &&
        identical(cached$genes, rownames(counts)) &&
        identical(cached$signature, signature)) {
      cat("  -> 读取 DecontX 检查点：", sample_id, "\n", sep = "")
      return(cached)
    }
    warning("DecontX 检查点与当前输入不匹配，将重新计算：", sample_id)
  }
  
  cat("  -> DecontX：", sample_id, "（", ncol(counts), " cells）\n", sep = "")
  preclusters <- make_decontx_preclusters(counts, sample_id)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  
  set.seed(SEED)
  if (is.null(preclusters)) {
    result <- decontX::decontX(sce)
  } else {
    result <- decontX::decontX(sce, z = preclusters)
  }
  
  corrected <- round_sparse_counts(decontX::decontXcounts(result))
  corrected <- corrected[rownames(counts), colnames(counts), drop = FALSE]
  contamination <- as.numeric(
    SummarizedExperiment::colData(result)$decontX_contamination
  )
  names(contamination) <- colnames(result)
  contamination <- contamination[colnames(counts)]
  
  if (anyNA(contamination) || any(!is.finite(contamination))) {
    stop("DecontX 返回非有限污染比例：", sample_id)
  }
  
  output <- list(
    sample_id = sample_id,
    cells = colnames(counts),
    genes = rownames(counts),
    corrected_counts = corrected,
    contamination = contamination,
    signature = signature
  )
  atomic_save_rds(output, checkpoint, compress = TRUE)
  output
}


summarize_decontx <- function(metadata) {
  samples <- unique(as.character(metadata[[SAMPLE_COLUMN]]))
  rows <- lapply(samples, function(sample_id) {
    values <- metadata$decontX_contamination[metadata[[SAMPLE_COLUMN]] == sample_id]
    data.frame(
      SampleID = sample_id,
      n_cells = length(values),
      mean_contamination = mean(values),
      median_contamination = stats::median(values),
      p90_contamination = unname(stats::quantile(values, 0.90)),
      p95_contamination = unname(stats::quantile(values, 0.95)),
      pct_gt_20 = 100 * mean(values > 0.20),
      pct_gt_50 = 100 * mean(values > 0.50),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}


plot_decontx_qc <- function(metadata) {
  plot_data <- data.frame(
    SampleID = metadata[[SAMPLE_COLUMN]],
    contamination = metadata$decontX_contamination
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(x = SampleID, y = contamination)) +
    ggplot2::geom_violin(fill = "#4C78A8", colour = "#2F4B66", linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = c(0.20, 0.50, 0.80), linetype = c(3, 2, 1), colour = "#B22222") +
    ggplot2::labs(x = NULL, y = "DecontX contamination", title = "Ambient RNA contamination by capture") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


run_decontx_stage <- function(object) {
  signature <- make_stage_signature(
    "2A_DecontX",
    c(
      seed = SEED,
      schema = PIPELINE_SCHEMA_VERSION,
      input_fingerprint = INPUT_FINGERPRINT,
      sample_column = SAMPLE_COLUMN,
      decontX_version = as.character(utils::packageVersion("decontX")),
      clustering = "Seurat_resolution_0.4_or_internal"
    )
  )
  if (!FORCE_RERUN && file.exists(RDS_DECONTX)) {
    cached <- safe_read_rds(RDS_DECONTX)
    if (!is.null(cached) && seurat_checkpoint_matches_samples(
      cached,
      colnames(object),
      object[[]][[SAMPLE_COLUMN]],
      "2A",
      signature
    )) {
      cat(">>> 读取阶段 2A 检查点：", RDS_DECONTX, "\n", sep = "")
      return(cached)
    }
    warning("阶段 2A 检查点的细胞顺序或配置已变化，将重新计算。")
  }
  
  extracted <- get_assay_counts(object, "RNA")
  object <- extracted$object
  raw_counts <- extracted$counts
  metadata <- object[[]]
  sample_values <- as.character(metadata[[SAMPLE_COLUMN]])
  samples <- unique(sample_values)
  
  outputs <- lapply(samples, function(sample_id) {
    cells <- colnames(object)[sample_values == sample_id]
    run_one_decontx_sample(raw_counts[, cells, drop = FALSE], sample_id)
  })
  names(outputs) <- samples
  
  corrected <- do.call(cbind, lapply(outputs, `[[`, "corrected_counts"))
  corrected <- corrected[rownames(raw_counts), colnames(object), drop = FALSE]
  corrected <- methods::as(corrected, "dgCMatrix")
  
  contamination <- combine_named_vectors(lapply(outputs, `[[`, "contamination"))
  if (anyDuplicated(names(contamination)) ||
      !setequal(names(contamination), colnames(object))) {
    stop("DecontX 合并后的细胞条形码不一致。")
  }
  contamination <- contamination[colnames(object)]
  
  object[["decontX"]] <- SeuratObject::CreateAssay5Object(counts = corrected)
  object$decontX_contamination <- as.numeric(contamination)
  object$nCount_decontX <- as.numeric(Matrix::colSums(corrected))
  object$nFeature_decontX <- as.numeric(Matrix::colSums(corrected > 0))
  object$percent_mt_decontX <- as.numeric(calculate_percent(corrected, "^MT-"))
  object$percent_hb_decontX <- as.numeric(calculate_percent(
    corrected,
    "^HB(A1|A2|B|D|E1|G1|G2|M|Q1|Z)$"
  ))
  
  summary_table <- summarize_decontx(object[[]])
  write_tsv(summary_table, file.path(DIR_DECONTX, "01_DecontX_summary_by_sample.tsv"))
  save_pdf(
    plot_decontx_qc(object[[]]),
    file.path(DIR_FIGURE, "01_DecontX_contamination_by_sample.pdf"),
    width = 10,
    height = 5
  )
  object <- set_object_stage_signature(object, "2A", signature)
  atomic_save_rds(object, RDS_DECONTX, compress = FALSE)
  object
}


# ==============================================================================
# 4. scDblFinder：在 DecontX counts 上按独立捕获批次识别双细胞
# ==============================================================================

make_bpparam <- function(workers, seed) {
  workers <- max(1L, as.integer(workers))
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(
      workers = workers,
      type = "SOCK",
      RNGseed = seed,
      progressbar = TRUE
    )
  } else {
    BiocParallel::MulticoreParam(
      workers = workers,
      RNGseed = seed,
      progressbar = TRUE
    )
  }
}


summarize_doublet_calls <- function(metadata) {
  samples <- unique(as.character(metadata[[SAMPLE_COLUMN]]))
  rows <- lapply(samples, function(sample_id) {
    x <- metadata[metadata[[SAMPLE_COLUMN]] == sample_id, , drop = FALSE]
    tested <- x$scDblFinder.class %in% c("singlet", "doublet")
    is_doublet <- x$scDblFinder.class == "doublet"
    n_tested <- sum(tested)
    data.frame(
      SampleID = sample_id,
      n_total = nrow(x),
      n_tested = n_tested,
      n_not_tested = sum(!tested),
      n_doublet = sum(is_doublet),
      called_doublet_rate = if (n_tested > 0) sum(is_doublet) / n_tested else NA_real_,
      nominal_10x_rate_from_cell_number = min(1, 0.01 * n_tested / 1000),
      median_score = if (n_tested > 0) stats::median(x$scDblFinder.score[tested]) else NA_real_,
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}


plot_doublet_rate <- function(summary_table) {
  ggplot2::ggplot(
    summary_table,
    ggplot2::aes(x = SampleID, y = 100 * called_doublet_rate)
  ) +
    ggplot2::geom_col(fill = "#E45756", width = 0.72) +
    ggplot2::geom_point(
      ggplot2::aes(y = 100 * nominal_10x_rate_from_cell_number),
      colour = "#1F4E79",
      size = 2.2
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Doublet rate (%)",
      title = "scDblFinder calls by capture",
      subtitle = "Bars: called rate; blue points: nominal 10x expectation from tested cell number"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


run_scdblfinder_stage <- function(object) {
  signature <- make_stage_signature(
    "2B_scDblFinder",
    c(
      seed = SEED,
      schema = PIPELINE_SCHEMA_VERSION,
      input_fingerprint = INPUT_FINGERPRINT,
      upstream_2A = get_object_stage_signature(object, "2A"),
      sample_column = SAMPLE_COLUMN,
      scDblFinder_version = as.character(utils::packageVersion("scDblFinder")),
      clusters = FALSE,
      expected_rate = "automatic"
    )
  )
  if (!FORCE_RERUN && file.exists(RDS_DOUBLET)) {
    cached <- safe_read_rds(RDS_DOUBLET)
    if (!is.null(cached) && seurat_checkpoint_matches_samples(
      cached,
      colnames(object),
      object[[]][[SAMPLE_COLUMN]],
      "2B",
      signature
    )) {
      cat(">>> 读取阶段 2B 检查点：", RDS_DOUBLET, "\n", sep = "")
      return(cached)
    }
    warning("阶段 2B 检查点的细胞顺序或配置已变化，将重新计算。")
  }
  
  extracted <- get_assay_counts(object, "decontX")
  object <- extracted$object
  corrected <- extracted$counts
  metadata <- object[[]]
  
  eligible <- Matrix::colSums(corrected) > 0
  eligible_cells <- colnames(corrected)[eligible]
  if (length(eligible_cells) < 20L) stop("可用于双细胞识别的细胞少于 20 个。")
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = corrected[, eligible_cells, drop = FALSE])
  )
  SummarizedExperiment::colData(sce)[[SAMPLE_COLUMN]] <-
    as.character(metadata[eligible_cells, SAMPLE_COLUMN])
  
  n_samples <- length(unique(SummarizedExperiment::colData(sce)[[SAMPLE_COLUMN]]))
  workers <- min(N_WORKERS_SCDOUBLET, n_samples)
  bp <- make_bpparam(workers, SEED)
  
  set.seed(SEED)
  sce <- scDblFinder::scDblFinder(
    sce,
    samples = SAMPLE_COLUMN,
    clusters = FALSE,
    dbr = NULL,
    BPPARAM = bp,
    verbose = TRUE
  )
  
  result <- as.data.frame(SummarizedExperiment::colData(sce))
  validate_scdblfinder_result(result, eligible_cells)
  result <- result[eligible_cells, , drop = FALSE]
  
  scores <- stats::setNames(rep(NA_real_, ncol(object)), colnames(object))
  classes <- stats::setNames(rep("not_tested_zero_count", ncol(object)), colnames(object))
  scores[eligible_cells] <- as.numeric(result$scDblFinder.score)
  classes[eligible_cells] <- normalize_doublet_class(result$scDblFinder.class)
  object$scDblFinder.score <- unname(scores[colnames(object)])
  object$scDblFinder.class <- unname(classes[colnames(object)])
  
  doublet_columns <- grep("^scDblFinder\\.", colnames(result), value = TRUE)
  score_table <- data.frame(
    cell = rownames(result),
    SampleID = as.character(result[[SAMPLE_COLUMN]]),
    result[, doublet_columns, drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  write_tsv(score_table, file.path(DIR_DOUBLET, "02_scDblFinder_cell_scores.tsv"))
  
  summary_table <- summarize_doublet_calls(object[[]])
  write_tsv(summary_table, file.path(DIR_DOUBLET, "02_scDblFinder_summary_by_sample.tsv"))
  save_pdf(
    plot_doublet_rate(summary_table),
    file.path(DIR_FIGURE, "02_scDblFinder_rate_by_sample.pdf"),
    width = 10,
    height = 5
  )
  
  object <- set_object_stage_signature(object, "2B", signature)
  atomic_save_rds(object, RDS_DOUBLET, compress = FALSE)
  object
}


# ==============================================================================
# 5. 细胞级过滤：删除双细胞、低信息、濒死和极端污染细胞
# ==============================================================================

plot_filter_summary <- function(summary_table) {
  long <- rbind(
    data.frame(SampleID = summary_table$SampleID, Stage = "Before", Cells = summary_table$n_before),
    data.frame(SampleID = summary_table$SampleID, Stage = "After", Cells = summary_table$n_after)
  )
  long$Stage <- factor(long$Stage, levels = c("Before", "After"))
  ggplot2::ggplot(long, ggplot2::aes(x = SampleID, y = Cells, fill = Stage)) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::scale_fill_manual(values = c(Before = "#9ECAE1", After = "#3182BD")) +
    ggplot2::labs(x = NULL, y = "Cells", title = "Cell retention after stage-2 QC") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


run_filter_stage <- function(object) {
  signature <- make_stage_signature(
    "2C_CellFilter",
    c(
      low_count = LOW_COUNT_DECONTX,
      schema = PIPELINE_SCHEMA_VERSION,
      input_fingerprint = INPUT_FINGERPRINT,
      upstream_2B = get_object_stage_signature(object, "2B"),
      low_feature = LOW_FEATURE_DECONTX,
      max_percent_mt = MAX_PERCENT_MT_DECONTX,
      extreme_contamination = EXTREME_CONTAMINATION,
      extreme_contamination_max_feature = EXTREME_CONTAMINATION_MAX_FEATURE
    )
  )
  decision <- make_cell_qc_decision(object[[]])
  expected_keep_cells <- rownames(decision)[!decision$Final_Remove]
  if (!FORCE_RERUN && file.exists(RDS_FILTERED)) {
    cached <- safe_read_rds(RDS_FILTERED)
    expected_samples <- object[[]][expected_keep_cells, SAMPLE_COLUMN]
    if (!is.null(cached) && seurat_checkpoint_matches_samples(
      cached,
      expected_keep_cells,
      expected_samples,
      "2C",
      signature
    )) {
      cat(">>> 读取阶段 2C 检查点：", RDS_FILTERED, "\n", sep = "")
      return(cached)
    }
    warning("阶段 2C 检查点的保留细胞或过滤阈值已变化，将重新计算。")
  }
  
  flag_columns <- setdiff(colnames(decision), "SampleID")
  object <- SeuratObject::AddMetaData(
    object,
    metadata = decision[colnames(object), flag_columns, drop = FALSE]
  )
  
  summary_table <- summarize_cell_filtering(decision)
  total <- data.frame(
    SampleID = "TOTAL",
    n_before = sum(summary_table$n_before),
    n_removed_low_count = sum(summary_table$n_removed_low_count),
    n_removed_low_feature = sum(summary_table$n_removed_low_feature),
    n_removed_high_mito = sum(summary_table$n_removed_high_mito),
    n_removed_extreme_contamination = sum(summary_table$n_removed_extreme_contamination),
    n_removed_doublet = sum(summary_table$n_removed_doublet),
    n_removed_union = sum(summary_table$n_removed_union),
    n_after = sum(summary_table$n_after),
    retention_rate = sum(summary_table$n_after) / sum(summary_table$n_before)
  )
  summary_with_total <- rbind(summary_table, total)
  
  manifest <- cbind(
    data.frame(cell = rownames(decision), row.names = rownames(decision)),
    object[[]][rownames(decision), c(
      SAMPLE_COLUMN, "nCount_decontX", "nFeature_decontX",
      "percent_mt_decontX", "percent_hb_decontX",
      "decontX_contamination", "scDblFinder.score", "scDblFinder.class"
    ), drop = FALSE],
    decision[, flag_columns, drop = FALSE]
  )
  write_tsv(manifest, file.path(DIR_FILTER, "03_cell_level_QC_manifest.tsv"))
  write_tsv(summary_with_total, file.path(DIR_FILTER, "03_filter_summary_by_sample.tsv"))
  save_pdf(
    plot_filter_summary(summary_table),
    file.path(DIR_FIGURE, "03_cell_retention_by_sample.pdf"),
    width = 10,
    height = 5
  )
  
  keep_cells <- expected_keep_cells
  if (length(keep_cells) < 20L) stop("联合过滤后细胞少于 20 个，停止。")
  filtered <- subset(object, cells = keep_cells)
  filtered <- set_object_stage_signature(filtered, "2C", signature)
  atomic_save_rds(filtered, RDS_FILTERED, compress = FALSE)
  filtered
}


# ==============================================================================
# 6. 保守重分簇：仅排除技术性 HVG，不自动删除 cluster
# ==============================================================================

add_cell_cycle_scores <- function(object) {
  s_genes <- intersect(Seurat::cc.genes.updated.2019$s.genes, rownames(object))
  g2m_genes <- intersect(Seurat::cc.genes.updated.2019$g2m.genes, rownames(object))
  if (length(s_genes) >= 10L && length(g2m_genes) >= 10L) {
    object <- Seurat::CellCycleScoring(
      object,
      s.features = s_genes,
      g2m.features = g2m_genes,
      set.ident = FALSE
    )
  } else {
    warning("可匹配的细胞周期基因不足；跳过 CellCycleScoring。")
  }
  object
}


summarize_clusters_for_review <- function(metadata) {
  clusters <- unique(as.character(metadata$qc_cluster))
  rows <- lapply(clusters, function(cluster_id) {
    x <- metadata[as.character(metadata$qc_cluster) == cluster_id, , drop = FALSE]
    median_score <- if (all(is.na(x$scDblFinder.score))) NA_real_ else stats::median(x$scDblFinder.score, na.rm = TRUE)
    signals <- c(
      SmallCluster = nrow(x) < 50L,
      LowComplexity = stats::median(x$nFeature_decontX) < 600,
      HighMito = stats::median(x$percent_mt_decontX, na.rm = TRUE) > 7,
      HighContamination = stats::median(x$decontX_contamination, na.rm = TRUE) > 0.30
    )
    active <- names(signals)[signals]
    review <- if (!length(active)) "Pass" else paste0("Review_", paste(active, collapse = "+"))
    
    data.frame(
      qc_cluster = cluster_id,
      n_cells = nrow(x),
      n_samples = length(unique(x[[SAMPLE_COLUMN]])),
      median_nCount_decontX = stats::median(x$nCount_decontX),
      median_nFeature_decontX = stats::median(x$nFeature_decontX),
      median_percent_mt_decontX = stats::median(x$percent_mt_decontX, na.rm = TRUE),
      median_decontX_contamination = stats::median(x$decontX_contamination, na.rm = TRUE),
      median_scDblFinder_score = median_score,
      cluster_review_flag = review,
      row.names = NULL
    )
  })
  output <- do.call(rbind, rows)
  output[order(as.numeric(output$qc_cluster)), , drop = FALSE]
}


write_cluster_sample_composition <- function(metadata) {
  counts <- as.data.frame.matrix(table(
    qc_cluster = metadata$qc_cluster,
    SampleID = metadata[[SAMPLE_COLUMN]]
  ))
  counts <- cbind(qc_cluster = rownames(counts), counts, row.names = NULL)
  write_tsv(counts, file.path(DIR_RECLUSTER, "04_cluster_by_sample_cell_counts.tsv"))
}


run_recluster_stage <- function(object) {
  signature <- make_stage_signature(
    "2D_Recluster",
    c(
      seed = SEED,
      schema = PIPELINE_SCHEMA_VERSION,
      input_fingerprint = INPUT_FINGERPRINT,
      upstream_2C = get_object_stage_signature(object, "2C"),
      hvg_candidates = HVG_CANDIDATES,
      hvg_final = HVG_FINAL,
      n_pcs = N_PCS,
      cluster_resolution = CLUSTER_RESOLUTION,
      hvg_blacklist = "mitochondrial_ribosomal_hemoglobin_only_v1"
    )
  )
  if (!FORCE_RERUN && file.exists(RDS_FINAL)) {
    cached <- safe_read_rds(RDS_FINAL)
    if (!is.null(cached) && seurat_checkpoint_matches_samples(
      cached,
      colnames(object),
      object[[]][[SAMPLE_COLUMN]],
      "2D",
      signature
    )) {
      cat(">>> 读取阶段 2D 检查点：", RDS_FINAL, "\n", sep = "")
      return(cached)
    }
    warning("阶段 2D 检查点的细胞顺序或重分簇参数已变化，将重新计算。")
  }
  
  SeuratObject::DefaultAssay(object) <- "decontX"
  object <- Seurat::NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = min(HVG_CANDIDATES, nrow(object)),
    verbose = FALSE
  )
  
  candidates <- SeuratObject::VariableFeatures(object)
  blacklist <- technical_gene_blacklist(candidates)
  retained_hvgs <- setdiff(candidates, blacklist)
  retained_hvgs <- head(retained_hvgs, min(HVG_FINAL, length(retained_hvgs)))
  if (length(retained_hvgs) < 500L) {
    stop("排除技术基因后高变基因少于 500 个，停止重分簇。")
  }
  SeuratObject::VariableFeatures(object) <- retained_hvgs
  
  hvg_audit <- data.frame(
    gene = candidates,
    technical_blacklist = candidates %in% blacklist,
    used_for_reclustering = candidates %in% retained_hvgs,
    row.names = NULL
  )
  write_tsv(hvg_audit, file.path(DIR_RECLUSTER, "04_HVG_audit.tsv"))
  
  object <- add_cell_cycle_scores(object)
  object <- Seurat::ScaleData(object, features = retained_hvgs, verbose = FALSE)
  
  max_pcs <- min(N_PCS, length(retained_hvgs) - 1L, ncol(object) - 1L)
  if (max_pcs < 10L) stop("可用 PCA 维度少于 10。")
  dims_use <- seq_len(max_pcs)
  
  object <- Seurat::RunPCA(
    object,
    features = retained_hvgs,
    npcs = max_pcs,
    reduction.name = "pca_qc",
    reduction.key = "PCQC_",
    seed.use = SEED,
    verbose = FALSE
  )
  object <- Seurat::FindNeighbors(
    object,
    reduction = "pca_qc",
    dims = dims_use,
    graph.name = c("qc_nn", "qc_snn"),
    verbose = FALSE
  )
  object <- Seurat::FindClusters(
    object,
    graph.name = "qc_snn",
    resolution = CLUSTER_RESOLUTION,
    algorithm = 1,
    random.seed = SEED,
    verbose = FALSE
  )
  object$qc_cluster <- as.character(SeuratObject::Idents(object))
  object <- Seurat::RunUMAP(
    object,
    reduction = "pca_qc",
    dims = dims_use,
    reduction.name = "umap_qc",
    reduction.key = "UMAPQC_",
    seed.use = SEED,
    verbose = FALSE
  )
  
  cluster_qc <- summarize_clusters_for_review(object[[]])
  write_tsv(cluster_qc, file.path(DIR_RECLUSTER, "04_cluster_technical_QC_review.tsv"))
  write_cluster_sample_composition(object[[]])
  
  p_cluster <- Seurat::DimPlot(
    object,
    reduction = "umap_qc",
    group.by = "qc_cluster",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) + ggplot2::ggtitle("Post-QC clusters (review only; no cluster auto-deleted)")
  p_sample <- Seurat::DimPlot(
    object,
    reduction = "umap_qc",
    group.by = SAMPLE_COLUMN,
    raster = TRUE
  ) + ggplot2::ggtitle("Post-QC cells by capture")
  save_pdf(
    p_cluster,
    file.path(DIR_FIGURE, "04_PostQC_UMAP_clusters.pdf"),
    width = 9,
    height = 7
  )
  save_pdf(
    p_sample,
    file.path(DIR_FIGURE, "04_PostQC_UMAP_samples.pdf"),
    width = 10,
    height = 7
  )
  
  object <- set_object_stage_signature(object, "2D", signature)
  atomic_save_rds(object, RDS_FINAL, compress = FALSE)
  object
}


# ==============================================================================
# 7. 主程序
# ==============================================================================

main <- function() {
  ensure_output_directories()
  log_path <- file.path(
    DIR_LOG,
    paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
  sink(log_connection, type = "output", split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_connection)
  }, add = TRUE)
  
  cat("============================================================\n")
  cat("GSE223063 第二阶段安全质控流程\n")
  cat("开始时间：", format(Sys.time()), "\n")
  cat("输入：", INPUT_RDS, "\n")
  cat("输出：", WORK_DIR, "\n")
  cat("强制重跑：", FORCE_RERUN, "\n")
  cat("============================================================\n\n")
  
  package_versions <- require_pipeline_packages()
  if (!file.exists(INPUT_RDS)) stop("找不到输入 RDS：", INPUT_RDS)
  INPUT_FINGERPRINT <<- unname(tools::md5sum(INPUT_RDS))
  if (is.na(INPUT_FINGERPRINT) || !nzchar(INPUT_FINGERPRINT)) {
    stop("无法计算输入 RDS 的 MD5 指纹。")
  }
  write_configuration(package_versions)
  cat(">>> 读取输入 RDS...\n")
  input_object <- readRDS(INPUT_RDS)
  validation <- validate_pipeline_input(input_object)
  input_object <- validation$object
  cat(">>> 输入检查通过：", ncol(input_object), " cells；",
      length(validation$samples), " captures。\n", sep = "")
  
  cat("\n>>> 阶段 2A：DecontX 环境 RNA 校正\n")
  object_2a <- run_decontx_stage(input_object)
  rm(input_object)
  invisible(gc())
  
  cat("\n>>> 阶段 2B：scDblFinder 双细胞识别\n")
  object_2b <- run_scdblfinder_stage(object_2a)
  rm(object_2a)
  invisible(gc())
  
  cat("\n>>> 阶段 2C：细胞级联合过滤\n")
  object_2c <- run_filter_stage(object_2b)
  rm(object_2b)
  invisible(gc())
  
  cat("\n>>> 阶段 2D：保守重分簇与 cluster 复核表\n")
  object_2d <- run_recluster_stage(object_2c)
  rm(object_2c)
  invisible(gc())
  
  writeLines(
    capture.output(sessionInfo()),
    con = file.path(DIR_LOG, "sessionInfo.txt"),
    useBytes = TRUE
  )
  
  cat("\n============================================================\n")
  cat("流程完成。最终对象：", RDS_FINAL, "\n", sep = "")
  cat("最终细胞数：", ncol(object_2d), "\n", sep = "")
  cat("请优先检查：04_cluster_technical_QC_review.tsv 和两张 PostQC UMAP。\n")
  cat("脚本没有自动删除任何 cluster。\n")
  cat("结束时间：", format(Sys.time()), "\n")
  cat("============================================================\n")
  invisible(object_2d)
}


if (!identical(Sys.getenv("GLIOMA_QC_NO_AUTORUN", "0"), "1") &&
    sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message("\n[FATAL] ", conditionMessage(e))
      quit(save = "no", status = 1L)
    }
  )
}















# ==============================================================================
# 胶质瘤单细胞 Step 2 尾声 —— 【物理切除劣质簇并移交 Step 3】
# ==============================================================================
library(Seurat)

# 1. 路径与对象加载
WORK_DIR <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/2_RNA污染双细胞与濒死细胞质控"
path_2d <- file.path(WORK_DIR, "Glioma_2D_QC_Reclustered.rds")

# 创建下一阶段的工作目录
STEP3_DIR <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/GSE273274(scRNA)_瘤周区/3.注释新图形"
dir.create(STEP3_DIR, recursive = TRUE, showWarnings = FALSE)
path_3a <- file.path(STEP3_DIR, "Glioma_3A_ReadyForIntegration.rds")

cat(">>> 正在读取聚类数据...\n")
glioma_final <- readRDS(path_2d)

# 注意：新版安全管线中，聚类结果保存在 qc_cluster 列
Idents(glioma_final) <- "qc_cluster"

# 2. 锁定死刑名单：数据审查确认的单样本高污染技术废料 (Soup Blobs)
bad_clusters <- c("30") 

# 3. 执行切除
cells_before <- ncol(glioma_final)
cat(">>> 🚀 正在执行物理切除...\n")

glioma_clean <- subset(glioma_final, idents = bad_clusters, invert = TRUE)
cells_after <- ncol(glioma_clean)

cat("\n=================================================================\n")
cat(">>> 🔪 簇级终极切除汇报：\n")
cat(">>> 切除前细胞总数：", cells_before, "\n")
cat(">>> 剔除纯技术垃圾细胞数：", cells_before - cells_after, "\n")
cat(">>> 移交 Step 3 的绝对纯净细胞数：", cells_after, "\n")
cat("=================================================================\n")

# 4. 暴力抹除含有垃圾坐标的旧降维结果（强制 Step 3 重新降维）
glioma_clean@reductions <- list()
glioma_clean@graphs <- list()
glioma_clean@neighbors <- list()

# 5. 保存终极纯净底盘
saveRDS(glioma_clean, path_3a)
cat(">>> 🎉 恭喜！完美纯净底盘已保存至新目录：\n>>>", path_3a, "\n")

