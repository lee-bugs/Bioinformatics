#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# cell2location 一体化脚本（单文件 · 西柚云 cell2loc_env）
# -----------------------------------------------------------------------------
# 输入：
#   --ref-rds    单细胞 Seurat RDS（1 个；也支持 RDS 内含 Seurat 对象 list）
#   --st-rds     Visium V2 / Visium HD 的 Seurat RDS，可 1 个也可多个：
#                多个路径、通配符、目录、或列出路径的 .txt 都支持；
#                单个 RDS 内可以是「一个合并对象多切片」或「多个 Seurat 对象的 list」。
#   --labels-key 参考 RDS 中细胞类型列名（必填）
#
# 输出（--out-dir 下）：
#   input/    ref_* / st_* counts 三件套、st_coords.csv、st_key_map.csv、prep_summary.txt
#   output/   spot_annotations.csv             ★ 注释后的 spot 总表（每个 spot 一行）
#             annotated_spots_by_sample/*.csv  ★ 按样本拆分的 spot 表
#             sample_annotation_summary.csv   ★ 各样本主导细胞类型统计
#             c2l_result.h5ad / cell_abundance_q05.csv / slice_summary.csv
#             inf_aver.csv（参考签名，带一致性校验缓存）
#             abundance_overview.png / c2l_spatial_abundance.pdf
#             c2l_dominant_spatial.png / dominant_composition.png
#             seurat_annotated/*.rds（可选，回贴后的 Seurat 对象）
#             c2l_report.txt（含 RESULT: PASS）
#   _rds_helper.R / r_helper.log（自动生成的 RDS 读取器，可留作排查）
#
# 用法（西柚云）：
#   /home/shazhi11/miniforge3/envs/cell2loc_env/bin/python c2l_onefile.py \
#     --ref-rds /path/scRNA.rds \
#     --st-rds  /path/spatial_merged.rds \
#     --labels-key CellType --out-dir ~/c2l_out
#
#   多样本（多个 RDS / 目录 / 通配符，任选一种写法）：
#     --st-rds /path/S1.rds /path/S2.rds
#     --st-rds '/path/visium/*.rds'
#     --st-rds /path/visium_dir
#     --st-rds /path/rds_list.txt
#
# 先测通再正式跑（强烈建议，换数据集必做）：
#   ... --epochs-ref 2 --epochs-st 2 --max-per-type 500 --prep-only   # 只做转换
#   ... --epochs-ref 2 --epochs-st 2 --max-per-type 500                # 全链路测通
#
# 设计要点 / 已内置的坑（勿删）：
#   1. RDS 是 R 的序列化格式，Python 无法直接读 Seurat S4 对象。本脚本把 R
#      代码作为字符串内嵌、运行时落盘为 _rds_helper.R，用 miloR 的 Rscript 调用；
#      除「读 RDS / 回贴 RDS」外全部逻辑都在本 Python 文件里。
#   2. 多样本 barcode 会撞车：统一改成 "<sample>__<原barcode>"，st_key_map.csv
#      保留回贴映射；batch_key 用样本/切片，跨样本检测效率差异才能被建模。
#   3. Seurat V5 多切片 counts.N layers 先 JoinLayers；raw counts 探测优先
#      RNA / Spatial，避免误选 SCT/decontX 的 counts 层。
#   4. Visium HD：自动识别 8um/16um/square bin 的 assay（默认取最粗的 16um
#      以控制 spot 数），兼容 pxl_row_in_fullres 等坐标格式；HD 且未显式指定
#      --n-cells-per-location 时自动用 3（V2 默认 20）。
#   5. 共享 GPU 显存少：batch=512 + PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True。
#   6. scvi-tools 1.3 的 train_size=1 不能开 early_stopping；train 参数做了兼容回退。
#   7. inf_aver.csv 只在「基因列表 + 细胞类型 + 参考批次」完全一致时复用。
#   8. 空间坐标按 image 等比归一化到 0~1；HD spot 极多时散点图随机抽样并
#      rasterized，PDF 不至于爆体积。
#   9. 只有 1 类细胞类型、或共同基因 < 100 时直接报错，不做无意义映射。
#
# 超参：Kleshchevnikov et al., Nat Biotechnol 2022, doi:10.1038/s41587-021-01139-4
#       N_cells_per_location / detection_alpha 见官方 Note_on_selecting_hyperparameters.pdf
#       q05 丰度是每 spot 各类细胞的「估计数量」，不是比例。
# =============================================================================
from __future__ import annotations

import os
import sys
import gc
import re
import json
import glob
import time
import shutil
import hashlib
import argparse
import traceback
import subprocess

os.environ.setdefault("PYTHONNOUSERSITE", "True")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

# cell2location 0.1.5 依赖 pyro<=1.8 的 deep_getattr；环境已降级则无副作用
try:
    import pyro.infer.autoguide.utils as _u
    if not hasattr(_u, "deep_getattr"):
        def deep_getattr(obj, name):
            for part in name.split("."):
                obj = getattr(obj, part)
            return obj
        _u.deep_getattr = deep_getattr
except Exception:
    pass

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import torch

DEFAULT_RSCRIPT = "/home/shazhi11/miniforge3/envs/miloR/bin/Rscript"


def log(m):
    print(f"[C2L] {m}", flush=True)


# =============================================================================
# 内嵌 R 脚本：只负责 RDS -> MatrixMarket/CSV，以及结果回贴 RDS
# 约定：R 代码中一律使用双引号字符串，避免与 Python 的 r''' ''' 冲突
# =============================================================================
R_HELPER = r'''
#!/usr/bin/env Rscript
# 由 c2l_onefile.py 自动生成，勿手工编辑
suppressPackageStartupMessages({
  library(Matrix)
  if (requireNamespace("Seurat", quietly = TRUE)) library(Seurat) else library(SeuratObject)
})
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

get_all <- function(flag) {
  i <- which(args == flag)
  if (!length(i)) return(character(0))
  unlist(lapply(i, function(k) if (k < length(args)) args[k + 1L] else NULL))
}
get_one <- function(flag, default = NULL) {
  v <- get_all(flag); if (length(v)) v[1L] else default
}

mode        <- get_one("--mode", "prep")
outdir     <- path.expand(get_one("--out-dir", "./c2l_out/input"))
ref_rds    <- path.expand(get_one("--ref-rds", ""))
st_paths   <- path.expand(get_all("--st-rds"))
labels_key <- get_one("--labels-key", "")
ref_assay  <- get_one("--ref-assay", "")
st_assay   <- get_one("--st-assay", "")
ref_batch  <- get_one("--ref-batch", "")
st_batch   <- get_one("--st-batch", "")
st_bin     <- get_one("--st-bin", "")
max_per_type <- as.integer(get_one("--max-per-type", "5000"))
max_spots  <- as.integer(get_one("--max-spots-per-sample", "0"))
seed       <- as.integer(get_one("--seed", "20261106"))
ab_csv     <- path.expand(get_one("--abun-csv", ""))
key_csv    <- path.expand(get_one("--key-csv", ""))
outrds_dir <- path.expand(get_one("--out-rds-dir", ""))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
stem <- function(p) sub("\\.(rds|RDS)$", "", basename(p))
# Visium HD 分辨率字样（8um/16um/square_002um/bin/...），别用裸 "um" 以免误判 tumor 之类
BIN_PAT <- "([0-9]+um|[0-9]+µm|square|bin|visiumhd|[_.-]hd([_.-]|$))"

# ---------------------------------------------------------------- 读 RDS
# 返回 list(objs=list(Seurat...), names=character) —— 兼容单对象 / 对象 list
read_objects <- function(p) {
  x <- readRDS(p)
  if (inherits(x, "Seurat"))
    return(list(objs = list(x), names = stem(p)))
  if (is.list(x) && length(x)) {
    ok <- vapply(x, function(z) inherits(z, "Seurat"), logical(1))
    if (all(ok)) {
      nm <- names(x)
      if (is.null(nm) || any(nm == "")) nm <- paste0(stem(p), "_obj", seq_along(x))
      return(list(objs = x, names = as.character(nm)))
    }
  }
  stop(sprintf("[%s] RDS 既不是 Seurat 对象也不是 Seurat 对象列表；class=%s",
               basename(p), paste(class(x), collapse = "/")))
}

# ---------------------------------------------------------------- counts 探测
probe_counts <- function(o, tag, preferred = NULL, bin = NULL) {
  an <- tryCatch(names(o), error = function(e) NULL)
  if (is.null(an) || !length(an)) an <- tryCatch(names(o@assays), error = function(e) character(0))
  if (!length(an)) stop(sprintf("[%s] 对象没有 assay", tag))
  if (!is.null(preferred) && preferred != "") {
    if (!preferred %in% an)
      stop(sprintf("[%s] 指定 assay '%s' 不存在；可用: %s", tag, preferred,
                   paste(an, collapse = ", ")))
    an <- preferred
  } else if (!is.null(bin) && bin != "") {
    hit <- an[grepl(bin, an, ignore.case = TRUE)]
    if (!length(hit))
      stop(sprintf("[%s] --st-bin '%s' 未匹配到任何 assay；可用: %s", tag, bin,
                   paste(an, collapse = ", ")))
    an <- hit
  } else if (length(an) > 1L) {
    # Visium HD 同一个对象里可能同时有 2um/8um/16um 多个 assay：默认取最粗
    cand <- an[grepl(BIN_PAT, an, ignore.case = TRUE)]
    if (length(cand) > 1L) {
      ord <- c("16um", "016um", "16", "8um", "008um", "8", "2um", "002um", "2")
      chosen <- NULL
      for (q in ord) {
        h <- cand[grepl(paste0("^.*", q), cand, ignore.case = TRUE)]
        if (length(h)) { chosen <- h[1L]; break }
      }
      if (!is.null(chosen)) {
        cat(sprintf("  [%s][HD] 检测到多分辨率 assay，选用 '%s'（其余: %s）\n",
                    tag, chosen, paste(setdiff(an, chosen), collapse = ", ")))
        an <- chosen
      }
    }
  }
  res <- list()
  for (nm in an) {
    asy <- tryCatch(o[[nm]], error = function(e) NULL)
    if (is.null(asy)) next
    lys <- tryCatch(names(asy@layers), error = function(e) NULL)
    ncl <- if (is.null(lys)) 0L else sum(grepl("^counts", lys))
    m <- NULL
    if (ncl > 1L) {
      a2 <- tryCatch(SeuratObject::JoinLayers(asy), error = function(e) asy)
      m <- tryCatch(SeuratObject::LayerData(a2, layer = "counts"), error = function(e) NULL)
    }
    if (is.null(m)) m <- tryCatch(SeuratObject::LayerData(o, assay = nm, layer = "counts"),
                                  error = function(e) NULL)
    if (is.null(m)) m <- tryCatch(SeuratObject::LayerData(asy, layer = "counts"),
                                  error = function(e) NULL)
    if (is.null(m)) m <- tryCatch(asy@counts, error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) next
    m <- as(m, "dgCMatrix")
    xx <- m@x
    fi <- if (length(xx)) mean(abs(xx - round(xx)) < 1e-6) else 0
    res[[nm]] <- list(m = m, nnz = length(xx), int = fi)
    cat(sprintf("  [%s] assay %-18s %6d genes x %7d cells, nnz=%11d, 整数比例=%.3f\n",
                tag, nm, nrow(m), ncol(m), length(xx), fi))
  }
  if (!length(res)) stop(sprintf("[%s] 没找到任何 counts 层；可用 assay: %s",
                                 tag, paste(an, collapse = ", ")))
  ok <- res[vapply(res, function(z) z$int > 0.99, logical(1))]
  pool <- if (length(ok)) ok else res
  nnz <- vapply(pool, function(z) z$nnz, numeric(1))
  pref <- intersect(c("RNA", "Spatial", "SpatialRNA", "raw", "Spatial.016um", "Spatial.008um"),
                    names(pool))
  pref <- pref[nnz[pref] >= 0.5 * max(nnz)]
  nm <- if (length(pref)) pref[which.max(nnz[pref])] else names(pool)[which.max(nnz)]
  why <- if (length(pref)) "Seurat raw-counts 常见命名优先" else "nnz 最大"
  if (!length(ok)) cat(sprintf("  [%s][WARN] 没有整数 counts，退化选 '%s'；建议 --ref-assay/--st-assay 明确指定\n",
                               tag, nm))
  cat(sprintf("  [%s] => 选用 assay '%s' (%s)\n", tag, nm, why))
  list(m = pool[[nm]]$m, assay = nm)
}

pick_col <- function(md, key, tag) {
  if (!is.null(key) && key != "") {
    if (!key %in% colnames(md))
      stop(sprintf("[%s] 列 '%s' 不存在；可用: %s", tag, key,
                   paste(head(colnames(md), 60), collapse = ", ")))
    return(key)
  }
  for (k in c("sample", "SampleID", "sample_id", "orig.ident", "batch", "Batch", "Sample", "slice"))
    if (k %in% colnames(md)) return(k)
  cat(sprintf("  [%s] 未找到批次列，全部视为同一批次\n", tag))
  NULL
}

meta_align <- function(o, m, tag) {
  idx <- match(colnames(m), rownames(o@meta.data))
  if (anyNA(idx)) stop(sprintf("[%s] counts 列名与 meta.data 行名对不上（%d 个失配）",
                               tag, sum(is.na(idx))))
  o@meta.data[idx, , drop = FALSE]
}

# ---------------------------------------------------------------- 坐标
xy_cols <- function(cn) {
  x <- intersect(c("imagerow", "row", "x", "pxl_row_in_fullres",
                   "center_x", "centroid_x", "array_row", "X"), cn)
  y <- intersect(c("imagecol", "col", "y", "pxl_col_in_fullres",
                   "center_y", "centroid_y", "array_col", "Y"), cn)
  list(x = if (length(x)) x[1L] else NULL, y = if (length(y)) y[1L] else NULL)
}

# 把 GetTissueCoordinates 的原始 barcode 映射到 counts 列名：
#   a) 完全相等  b) 去掉 Seurat merge 的 "_N" 后缀后相等
#   c) 按最后一个 "_" 分段的短码匹配（兼容 "GBM1_AAACA...-1" 和 HD "s_008um_x_y-1"）
#      多样本撞码时优先用样本前缀消歧
map_ids <- function(raw_ids, cell_ids, image, meta_sample = NULL, taken = NULL) {
  ii <- match(raw_ids, cell_ids)
  core <- sub("_[0-9]+$", "", cell_ids)
  j <- which(is.na(ii))
  if (length(j)) ii[j] <- match(raw_ids[j], core)
  j <- which(is.na(ii))
  if (length(j)) {
    tail_cell <- sub("^.*_", "", core)
    grp <- split(seq_along(core), tail_cell)
    tail_raw <- sub("^.*_", "", raw_ids[j])
    for (q in seq_along(j)) {
      if (is.na(tail_raw[q]) || tail_raw[q] == "") next
      cand <- grp[[tail_raw[q]]]
      if (is.null(cand) || !length(cand)) next
      if (!is.null(taken)) cand <- cand[!taken[cand]]
      if (!length(cand)) next
      if (length(cand) == 1L) { ii[j[q]] <- cand[1L]; next }
      pick <- cand
      if (!is.null(meta_sample)) {
        z <- cand[!is.na(meta_sample[cand]) & meta_sample[cand] == image]
        if (length(z)) pick <- z
      }
      z <- cand[startsWith(cell_ids[cand], paste0(image, "_"))]
      if (length(z)) pick <- z
      ii[j[q]] <- pick[1L]
    }
  }
  ii
}

get_coords <- function(o, cells, tag, sample_name) {
  md <- tryCatch(o@meta.data[match(cells, rownames(o@meta.data)), , drop = FALSE],
                 error = function(e) NULL)
  ms <- NULL
  if (!is.null(md)) {
    for (k in c("orig.ident", "sample", "SampleID", "sample_id", "fov", "image"))
      if (k %in% colnames(md)) { ms <- as.character(md[[k]]); break }
  }
  ims <- tryCatch(names(o@images), error = function(e) NULL)
  co <- NULL
  taken <- rep(FALSE, length(cells))
  if (length(ims)) {
    for (im in ims) {
      cc <- tryCatch(SeuratObject::GetTissueCoordinates(o, image = im),
                     error = function(e) NULL)
      if (is.null(cc) || nrow(cc) == 0) next
      cc <- as.data.frame(cc)
      if ("in_tissue" %in% colnames(cc)) cc <- cc[cc$in_tissue == 1, , drop = FALSE]
      if (nrow(cc) == 0) next
      xy <- xy_cols(colnames(cc))
      if (is.null(xy$x) || is.null(xy$y)) next
      # HD 的坐标表常把 barcode 放在列里而不是 rownames
      bc_col <- intersect(c("barcode", "cell", "spot", "Barcode"), colnames(cc))
      raw <- if (length(bc_col)) as.character(cc[[bc_col[1L]]]) else rownames(cc)
      x <- suppressWarnings(as.numeric(cc[[xy$x]]))
      y <- suppressWarnings(as.numeric(cc[[xy$y]]))
      good <- is.finite(x) & is.finite(y) & !is.na(raw)
      cc <- cc[good, , drop = FALSE]; x <- x[good]; y <- y[good]; raw <- raw[good]
      if (!length(x)) next
      sx <- min(x); sy <- min(y)
      scale <- max(max(x) - sx, max(y) - sy, 1)
      ii <- map_ids(raw, cells, im, ms, taken)
      use <- rep(FALSE, length(ii))
      for (k in seq_along(ii)) {
        if (!is.na(ii[k]) && !taken[ii[k]]) { use[k] <- TRUE; taken[ii[k]] <- TRUE }
      }
      if (any(use)) {
        co <- rbind(co, data.frame(cell = cells[ii[use]], x = (x[use] - sx) / scale,
                                   y = (y[use] - sy) / scale, image = im,
                                   stringsAsFactors = FALSE))
      }
      cat(sprintf("  [%s] image '%s': %d 坐标, 匹配 %d spots\n",
                  tag, im, nrow(cc), sum(use)))
    }
  }
  if (is.null(co) || nrow(co) == 0) {
    # 兜底：坐标直接在 meta.data 里（HD 常见）
    if (!is.null(md)) {
      xy <- xy_cols(colnames(md))
      if (!is.null(xy$x) && !is.null(xy$y)) {
        x <- suppressWarnings(as.numeric(md[[xy$x]]))
        y <- suppressWarnings(as.numeric(md[[xy$y]]))
        good <- is.finite(x) & is.finite(y)
        if (any(good)) {
          sx <- min(x[good]); sy <- min(y[good])
          scale <- max(max(x[good]) - sx, max(y[good]) - sy, 1)
          grp <- if (is.null(ms)) rep(sample_name, length(cells)) else ms
          grp[is.na(grp)] <- sample_name
          co <- data.frame(cell = cells[good], x = (x[good] - sx) / scale,
                           y = (y[good] - sy) / scale, image = grp[good],
                           stringsAsFactors = FALSE)
          cat(sprintf("  [%s] 使用 meta.data 坐标 %d spots\n", tag, nrow(co)))
        }
      }
    }
  }
  if (is.null(co) || nrow(co) == 0) {
    cat(sprintf("  [%s][WARN] 未取得空间坐标：仍会输出注释表，但不出空间图\n", tag))
    return(NULL)
  }
  co <- co[!duplicated(co$cell), , drop = FALSE]
  cat(sprintf("  [%s] 坐标匹配合计 %d / %d spots\n", tag, nrow(co), length(cells)))
  co
}

write_triplet <- function(m, prefix, cells_df) {
  Matrix::writeMM(m, file.path(outdir, paste0(prefix, "_counts.mtx")))
  write.csv(data.frame(gene = rownames(m)),
            file.path(outdir, paste0(prefix, "_genes.csv")), row.names = FALSE)
  write.csv(cells_df, file.path(outdir, paste0(prefix, "_cells.csv")), row.names = FALSE)
  cat(sprintf("  [%s] 写出 %d genes x %d cells, nnz=%d\n",
              prefix, nrow(m), ncol(m), length(m@x)))
}

# ================================================================= PREP =====
if (mode == "prep") {
  if (ref_rds == "" || !file.exists(ref_rds)) stop("必须给存在的 --ref-rds")
  if (!length(st_paths)) stop("必须给 --st-rds")
  if (labels_key == "") stop("必须给 --labels-key（参考 RDS 里的细胞类型列名）")
  for (p in st_paths) if (!file.exists(p)) stop(sprintf("空转 RDS 不存在: %s", p))
  set.seed(seed)

  # ---- 参考 scRNA（支持对象 list，基因取交集后合并）
  cat("=== REF:", ref_rds, "===\n")
  r <- read_objects(ref_rds)
  ref_ms <- list(); ref_lab <- list(); ref_bat <- list(); ref_assays <- c()
  for (i in seq_along(r$objs)) {
    o <- r$objs[[i]]
    p <- probe_counts(o, paste0("REF", i), if (ref_assay != "") ref_assay else NULL)
    ref_assays <- c(ref_assays, p$assay)
    md <- meta_align(o, p$m, paste0("REF", i))
    if (!labels_key %in% colnames(md))
      stop(sprintf("[REF] 细胞类型列 '%s' 不存在；可用: %s", labels_key,
                   paste(head(colnames(md), 80), collapse = ", ")))
    lab <- as.character(md[[labels_key]])
    if (anyNA(lab) || any(lab == "")) stop("细胞类型列存在 NA/空值")
    bk <- pick_col(md, if (ref_batch != "") ref_batch else NULL, "REF")
    bat <- if (is.null(bk)) rep(r$names[i], length(lab)) else as.character(md[[bk]])
    bat[is.na(bat)] <- "unknown"
    ref_ms[[i]] <- p$m; ref_lab[[i]] <- lab; ref_bat[[i]] <- bat
    rm(o); invisible(gc())
  }
  genes <- Reduce(intersect, lapply(ref_ms, rownames))
  if (length(genes) < 100) stop("参考 RDS 各对象共同基因 < 100，检查是否为同一物种/命名")
  m1 <- do.call(cbind, lapply(ref_ms, function(z) z[genes, , drop = FALSE]))
  lab <- unlist(ref_lab, use.names = FALSE)
  bat <- unlist(ref_bat, use.names = FALSE)
  if (length(lab) != ncol(m1)) stop("参考 metadata 与矩阵列数不一致")
  cat(sprintf("  细胞类型 %d 类 | 参考批次 %d 个\n", length(unique(lab)), length(unique(bat))))
  keep <- unlist(lapply(split(seq_along(lab), lab), function(ii) {
    if (max_per_type > 0 && length(ii) > max_per_type) sample(ii, max_per_type) else ii
  }), use.names = FALSE)
  keep <- sort(keep)
  cat(sprintf("  分层下采样: %d / %d (每类上限 %d, seed=%d)\n",
              length(keep), ncol(m1), max_per_type, seed))
  print(table(lab[keep]))
  m1 <- m1[, keep, drop = FALSE]
  write_triplet(m1, "ref", data.frame(cell = paste0("ref_", seq_len(ncol(m1))),
                                      cell_type = lab[keep], batch = bat[keep],
                                      orig_barcode = colnames(m1)))
  n_ref <- c(nrow(m1), ncol(m1))
  rm(r, ref_ms, m1); invisible(gc())

  # ---- 空转：多个 RDS / 多个对象 / 单对象多切片
  cat("=== ST:", length(st_paths), "个 RDS ===\n")
  st_ms <- list(); st_cells <- list(); st_meta <- list(); st_coords <- list()
  st_assays <- c(); platform_hits <- c()
  sidx <- 0L
  for (p in st_paths) {
    r2 <- read_objects(p)
    for (j in seq_along(r2$objs)) {
      sidx <- sidx + 1L
      o <- r2$objs[[j]]
      # 文件名本身能区分多个 RDS；只有「一个 RDS 内含多个对象」才需要拼对象名
      sample_name <- if (length(r2$objs) > 1L)
        safe(paste0(stem(p), "__", r2$names[j])) else safe(stem(p))
      tag <- sprintf("ST%d", sidx)
      cat(sprintf("--- %s [%s] %s ---\n", tag, sample_name, basename(p)))
      pr <- probe_counts(o, tag, if (st_assay != "") st_assay else NULL,
                         if (st_bin != "") st_bin else NULL)
      st_assays <- c(st_assays, pr$assay)
      m2 <- pr$m
      cells0 <- colnames(m2)
      md <- meta_align(o, m2, tag)
      bk <- pick_col(md, if (st_batch != "") st_batch else NULL, tag)
      bat <- if (is.null(bk)) rep(sample_name, ncol(m2)) else as.character(md[[bk]])
      bat[is.na(bat)] <- sample_name

      # 多样本时 batch 必须带样本前缀，否则不同切片同名 batch 会被合并
      if (length(st_paths) > 1L || length(r2$objs) > 1L)
        bat <- paste0(sample_name, "::", bat)

      # HD 平台判定（assay / image 名带 um、square、bin、HD，或 spot 数极大）
      ims <- tryCatch(names(o@images), error = function(e) NULL)
      hd <- (ncol(m2) > 50000L ||
             any(grepl(BIN_PAT, c(pr$assay, ims), ignore.case = TRUE)))
      platform_hits <- c(platform_hits, if (hd) "VisiumHD" else "VisiumV2")

      if (max_spots > 0 && ncol(m2) > max_spots) {
        cat(sprintf("  [%s][WARN] spot 数 %d > 上限 %d，随机抽样（注释表将只含抽样 spot）\n",
                    tag, ncol(m2), max_spots))
        sel <- sort(sample(seq_len(ncol(m2)), max_spots))
        m2 <- m2[, sel, drop = FALSE]; cells0 <- cells0[sel]; bat <- bat[sel]
      }

      uid <- paste0(sample_name, "__", cells0)
      if (anyDuplicated(uid)) uid <- make.unique(uid)
      co <- get_coords(o, cells0, tag, sample_name)
      if (!is.null(co)) {
        ii <- match(co$cell, cells0)
        co$cell <- uid[ii]
        co$sample <- sample_name
        co$orig_barcode <- cells0[ii]
        co$rds <- p; co$obj <- j
        st_coords[[length(st_coords) + 1L]] <- co
      }
      st_ms[[length(st_ms) + 1L]] <- m2
      st_cells[[length(st_cells) + 1L]] <- uid
      st_meta[[length(st_meta) + 1L]] <- data.frame(
        cell = uid, batch = bat, sample = sample_name,
        orig_barcode = cells0, rds = p, obj = j,
        platform = if (hd) "VisiumHD" else "VisiumV2",
        stringsAsFactors = FALSE)
      rm(o, m2, md); invisible(gc())
    }
    rm(r2); invisible(gc())
  }

  genes2 <- Reduce(intersect, lapply(st_ms, rownames))
  if (length(genes2) < 100) stop("空转 RDS 各样本共同基因 < 100")
  cat(sprintf("  空转共同基因 %d；样本/对象 %d 个\n", length(genes2), length(st_ms)))
  M2 <- do.call(cbind, lapply(st_ms, function(z) z[genes2, , drop = FALSE]))
  cells_df <- do.call(rbind, st_meta)
  if (ncol(M2) != nrow(cells_df)) stop("空转矩阵列数与 metadata 行数不一致")
  colnames(M2) <- cells_df$cell
  write_triplet(M2, "st", cells_df)
  co_all <- if (length(st_coords)) do.call(rbind, st_coords) else NULL
  if (!is.null(co_all)) {
    write.csv(co_all[, c("cell", "x", "y", "image", "sample", "orig_barcode", "rds", "obj")],
              file.path(outdir, "st_coords.csv"), row.names = FALSE)
  }
  key <- cells_df[, c("cell", "batch", "sample", "orig_barcode", "rds", "obj", "platform")]
  if (!is.null(co_all)) {
    ii <- match(key$cell, co_all$cell)
    key$x <- co_all$x[ii]; key$y <- co_all$y[ii]; key$image <- co_all$image[ii]
  } else {
    key$x <- NA_real_; key$y <- NA_real_; key$image <- NA_character_
  }
  write.csv(key, file.path(outdir, "st_key_map.csv"), row.names = FALSE)

  platform <- if (any(platform_hits == "VisiumHD")) "VisiumHD" else "VisiumV2"
  writeLines(c(
    sprintf("ref_rds   : %s", ref_rds),
    sprintf("st_rds    : %s", paste(st_paths, collapse = " ; ")),
    sprintf("ref_assay : %s", paste(unique(ref_assays), collapse = ",")),
    sprintf("st_assay  : %s", paste(unique(st_assays), collapse = ",")),
    sprintf("labels_key: %s (%d 类)", labels_key, length(unique(lab))),
    sprintf("ref_batch : %s", if (ref_batch != "") ref_batch else "auto"),
    sprintf("st_batch  : %s", if (st_batch != "") st_batch else "auto"),
    sprintf("ref_shape : %d genes x %d cells (下采样后)", n_ref[1], n_ref[2]),
    sprintf("st_shape  : %d genes x %d spots", nrow(M2), ncol(M2)),
    sprintf("n_samples : %d", length(st_ms)),
    sprintf("n_batches : %d", length(unique(cells_df$batch))),
    sprintf("platform  : %s", platform),
    sprintf("coords    : %s", if (is.null(co_all)) "no" else "yes"),
    sprintf("seed      : %d", seed)
  ), file.path(outdir, "prep_summary.txt"))
  cat("=== 转换完成 ->", outdir, "===\n")
}

# ================================================================= ATTACH ===
if (mode == "attach") {
  if (ab_csv == "" || key_csv == "" || outrds_dir == "")
    stop("attach 模式需要 --abun-csv / --key-csv / --out-rds-dir")
  dir.create(outrds_dir, recursive = TRUE, showWarnings = FALSE)
  key <- read.csv(key_csv, check.names = FALSE)
  ab <- read.csv(ab_csv, row.names = 1, check.names = FALSE)
  colnames(ab) <- sub("^q05cell_abundance_w_sf_", "", colnames(ab))
  cat(sprintf("丰度表: %d spots x %d 细胞类型\n", nrow(ab), ncol(ab)))
  grp <- split(seq_len(nrow(key)), list(key$rds, key$obj), drop = TRUE)
  for (idx in grp) {
    rows <- key[idx, , drop = FALSE]
    p <- as.character(rows$rds[1L]); j <- as.integer(rows$obj[1L])
    cat(sprintf("--- 回贴 %s [obj %d]: %d spots ---\n", basename(p), j, nrow(rows)))
    x <- readRDS(p)
    o <- if (inherits(x, "Seurat")) x else x[[j]]
    ii <- match(rows$cell, rownames(ab))
    keep <- !is.na(ii)
    mat <- as.data.frame(ab[ii[keep], , drop = FALSE], check.names = FALSE)
    rownames(mat) <- as.character(rows$orig_barcode[keep])
    common <- intersect(rownames(mat), colnames(o))
    rate <- length(common) / max(1L, ncol(o))
    cat(sprintf("  barcode 匹配 %d / %d (%.1f%%)\n", length(common), ncol(o), 100 * rate))
    if (rate < 0.95) cat("  [WARN] 匹配率偏低，检查 barcode 是否被改动过\n")
    mat <- mat[common, , drop = FALSE]
    o <- tryCatch(SeuratObject::AddMetaData(o, metadata = mat),
                  error = function(e) Seurat::AddMetaData(o, metadata = mat))
    o$c2l_dominant <- factor(colnames(mat)[apply(as.matrix(mat), 1L, which.max)],
                             levels = colnames(mat))
    nm <- if (nrow(rows)) safe(as.character(rows$sample[1])) else safe(stem(p))
    out <- file.path(outrds_dir, paste0(nm, ".cell2location.rds"))
    saveRDS(o, out)
    cat("  已保存 ->", out, "\n")
    print(table(o$c2l_dominant))
    rm(x, o); invisible(gc())
  }
  cat("=== 回贴完成 ===\n")
}
'''


# =============================================================================
# Python 侧工具
# =============================================================================
def run_r(rscript, helper, args, log_path):
    cmd = [rscript, helper] + [str(x) for x in args]
    log("Rscript: " + " ".join(cmd[:6]) + (" ..." if len(cmd) > 6 else ""))
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in p.stdout:
            fh.write(line)
            sys.stdout.write("[R] " + line)
            sys.stdout.flush()
        rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"Rscript 失败 (exit={rc})，详见 {log_path}")


def expand_st_paths(items):
    out = []
    for it in items:
        p = os.path.expanduser(str(it))
        if os.path.isdir(p):
            cand = sorted(glob.glob(os.path.join(p, "*.rds")) +
                          glob.glob(os.path.join(p, "*.RDS")))
            if not cand:
                raise RuntimeError(f"目录里没有 RDS: {p}")
            out.extend(cand)
        elif p.lower().endswith(".txt"):
            with open(p, encoding="utf-8") as f:
                out.extend([ln.strip() for ln in f if ln.strip() and not ln.startswith("#")])
        elif any(ch in p for ch in "*?["):
            cand = sorted(glob.glob(p))
            if not cand:
                raise RuntimeError(f"通配符没匹配到文件: {p}")
            out.extend(cand)
        else:
            out.append(p)
    seen, uniq = set(), []
    for p in out:
        p = os.path.abspath(os.path.expanduser(p))
        if p not in seen:
            if not os.path.isfile(p):
                raise RuntimeError(f"空转 RDS 不存在: {p}")
            seen.add(p); uniq.append(p)
    if not uniq:
        raise RuntimeError("没有解析到任何空转 RDS")
    return uniq


def parse_prep(path):
    d = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if ":" in line:
                k, v = line.split(":", 1)
                d[k.strip()] = v.strip()
    return d


def load_triplet(input_dir, prefix):
    M = sio.mmread(os.path.join(input_dir, f"{prefix}_counts.mtx")).tocsr()
    genes = pd.read_csv(os.path.join(input_dir, f"{prefix}_genes.csv"))["gene"].astype(str).tolist()
    cells = pd.read_csv(os.path.join(input_dir, f"{prefix}_cells.csv"))
    assert M.shape == (len(genes), len(cells)), \
        f"{prefix}: 矩阵 {M.shape} 与行列文件 ({len(genes)}, {len(cells)}) 不匹配"
    X = sp.csr_matrix(M.T).astype(np.float32)
    del M; gc.collect()
    a = ad.AnnData(X=X)
    a.var_names, a.obs_names = genes, cells["cell"].astype(str).tolist()
    for c in cells.columns:
        if c != "cell":
            a.obs[c] = cells[c].astype(str).values if cells[c].dtype == object else cells[c].values
    log(f"  {prefix}: {a.n_obs} cells/spots x {a.n_vars} genes, nnz={a.X.nnz}")
    return a


def train_compat(model, **kw):
    for _ in range(4):
        try:
            return model.train(**kw)
        except TypeError as e:
            msg = str(e)
            if "accelerator" in msg or "devices" in msg:
                kw.pop("accelerator", None); kw.pop("devices", None)
                kw["use_gpu"] = torch.cuda.is_available()
            elif "early_stopping" in msg:
                kw.pop("early_stopping", None); kw.pop("early_stopping_patience", None)
            elif "use_gpu" in msg:
                kw.pop("use_gpu", None)
            else:
                raise
            log(f"train 参数回退: {msg}")
    return model.train(**kw)


def extract_inf_aver(ref):
    key = "means_per_cluster_mu_fg"
    if key not in ref.varm:
        raise RuntimeError(f"ref.varm 缺少 {key}: {list(ref.varm.keys())}")
    m = ref.varm[key]
    df = m.copy() if hasattr(m, "columns") else pd.DataFrame(np.asarray(m), index=ref.var_names)
    cats = ref.uns.get("label_to_cat")
    if cats is not None and len(list(cats)) == df.shape[1]:
        df.columns = [str(c) for c in cats]
    df.columns = [re.sub(r"^means_per_cluster_mu_fg_", "", str(c)) for c in df.columns]
    df.index = list(ref.var_names)
    return df


def clean_cols(cols):
    return [re.sub(r"^q05cell_abundance_w_sf_", "", str(c)) for c in cols]


# =============================================================================
# 可视化
# =============================================================================
def plot_overview(df, summ, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    df.mean().sort_values(ascending=False).plot(kind="bar", ax=axes[0], color="#4C72B0")
    axes[0].set_title("Mean cell-type abundance per spot")
    axes[0].tick_params(axis="x", labelrotation=45)
    im = axes[1].imshow(summ.T.values, aspect="auto", cmap="magma")
    axes[1].set_yticks(range(len(summ.columns)))
    axes[1].set_yticklabels(summ.columns, fontsize=8)
    axes[1].set_xticks(range(len(summ.index)))
    axes[1].set_xticklabels(summ.index, rotation=45, ha="right", fontsize=8)
    axes[1].set_title("Mean abundance by sample/slice")
    fig.colorbar(im, ax=axes[1], shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "abundance_overview.png"), dpi=130)
    plt.close(fig)


def plot_composition(annot, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grp = annot.groupby(["sample", "dominant_type"]).size().unstack(fill_value=0)
    frac = grp.div(grp.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(frac)), 5))
    bottom = np.zeros(len(frac))
    cmap = plt.get_cmap("tab20", max(len(frac.columns), 1))
    for j, ct in enumerate(frac.columns):
        ax.bar(range(len(frac)), frac[ct].values, bottom=bottom,
               label=ct, color=cmap(j), width=0.75)
        bottom += frac[ct].values
    ax.set_xticks(range(len(frac)))
    ax.set_xticklabels(frac.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Fraction of spots (dominant type)")
    ax.set_title("Dominant cell type composition by sample/slice")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "dominant_composition.png"), dpi=140)
    plt.close(fig)


def plot_spatial(annot, out, point_size=None, max_types=12, max_points=80000):
    """每类一页（按切片分面）的连续丰度 PDF + 主导类型总览 PNG。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    d = annot.dropna(subset=["x", "y"]).copy()
    if len(d) < 0.5 * len(annot) or d.empty:
        log("  可用空间坐标不足一半，跳过空间图")
        return False
    cts = [c for c in annot.columns if c.startswith("q05_")]
    if not cts:
        return False
    # 类型很多时只画平均丰度最高的若干类，避免生成上百页 PDF
    if len(cts) > max_types:
        keep = annot[cts].mean().nlargest(max_types).index.tolist()
        cts = [c for c in cts if c in keep]
        log(f"  细胞类型较多，空间 PDF 只画平均丰度前 {max_types} 类")
    slices = sorted(d["image"].fillna(d["sample"]).astype(str).unique())
    if not slices:
        return False
    if len(d) > max_points:
        log(f"  spot 数 {len(d)} > {max_points}，空间图随机抽样显示")
        d = d.sample(max_points, random_state=0)
    if point_size is None:
        point_size = 0.4 if len(d) > 60000 else (1.2 if len(d) > 20000 else 4.0)
    raster = len(d) > 30000
    ncols = min(3, len(slices)); nrows = int(np.ceil(len(slices) / ncols))

    pdf_path = os.path.join(out, "c2l_spatial_abundance.pdf")
    with PdfPages(pdf_path) as pdf:
        for ct in cts:
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.6 * nrows),
                                     squeeze=False)
            vmax = float(d[ct].quantile(0.99)) or 1.0
            sc_art = None
            for k, sl in enumerate(slices):
                ax = axes[k // ncols][k % ncols]
                sub = d[(d["image"].fillna(d["sample"]).astype(str) == sl)]
                if len(sub):
                    sc_art = ax.scatter(sub["x"], sub["y"], c=sub[ct], s=point_size,
                                        cmap="magma", vmin=0, vmax=vmax,
                                        linewidths=0, rasterized=raster)
                ax.set_title(f"{sl} (n={len(sub)})", fontsize=9)
                ax.set_aspect("equal"); ax.invert_yaxis()
                ax.set_xticks([]); ax.set_yticks([])
            for k in range(len(slices), nrows * ncols):
                axes[k // ncols][k % ncols].axis("off")
            fig.suptitle(f"cell2location: {ct[4:]} (q05 abundance)", fontsize=12)
            if sc_art is not None:
                fig.colorbar(sc_art, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
            pdf.savefig(fig); plt.close(fig)

    dom = np.array([c[4:] for c in cts])[d[cts].values.argmax(axis=1)]
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.3 * ncols, 5 * nrows), squeeze=False)
    cmap = plt.get_cmap("tab10" if len(cts) <= 10 else "tab20", len(cts))
    grp = d["image"].fillna(d["sample"]).astype(str).values
    for k, sl in enumerate(slices):
        ax = axes[k // ncols][k % ncols]
        m = (grp == sl)
        for j, ct in enumerate(cts):
            mm = m & (dom == ct[4:])
            if mm.any():
                ax.scatter(d["x"].values[mm], d["y"].values[mm], s=point_size,
                           color=cmap(j), linewidths=0, rasterized=raster)
        ax.set_title(f"{sl} (n={int(m.sum())})", fontsize=10)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    for k in range(len(slices), nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    handles = [plt.Line2D([], [], marker="o", ls="", color=cmap(j),
                          label=ct[4:], markersize=5) for j, ct in enumerate(cts)]
    fig.legend(handles=handles, loc="lower center", ncol=min(5, len(cts)),
               fontsize=8, frameon=False)
    fig.suptitle("cell2location: dominant cell type per spot", fontsize=13)
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(os.path.join(out, "c2l_dominant_spatial.png"), dpi=150)
    plt.close(fig)
    log("  空间图 -> c2l_spatial_abundance.pdf / c2l_dominant_spatial.png")
    return True


# =============================================================================
# 主流程
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="单文件 cell2location：scRNA RDS + Visium V2/HD RDS -> 注释 spot 表 + 可视化",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--ref-rds", required=True, help="单细胞 Seurat RDS")
    ap.add_argument("--st-rds", required=True, nargs="+",
                    help="空转 RDS，可多个路径/通配符/目录/.txt 列表")
    ap.add_argument("--labels-key", required=True, default="CellType", help="参考 RDS 里细胞类型列名")
    ap.add_argument("--out-dir", default="./c2l_out")
    ap.add_argument("--rscript", default=None, help="Rscript 路径（默认西柚云 miloR）")
    ap.add_argument("--ref-assay", default="RNA", help="手工指定参考 counts assay")
    ap.add_argument("--st-assay", default="Spatial", help="手工指定空转 counts assay")
    ap.add_argument("--ref-batch", default="SampleID", help="参考批次列名")
    ap.add_argument("--st-batch", default="", help="空转批次/切片列名")
    ap.add_argument("--st-bin", default="", help="HD 分辨率筛选，如 16um / 008um；留空自动选最粗")
    ap.add_argument("--max-per-type", type=int, default=5000, help="每类细胞下采样上限，0=不采样")
    ap.add_argument("--max-spots-per-sample", type=int, default=10000,
                    help="HD 2um 极大量 spot 时的抽样上限；0=全量（默认）")
    ap.add_argument("--seed", type=int, default=20261106)
    ap.add_argument("--epochs-ref", type=int, default=250)
    ap.add_argument("--epochs-st", type=int, default=1000)
    ap.add_argument("--n-cells-per-location", type=int, default=None,
                    help="每 spot 预期细胞数；V2 默认 20，检测到 HD 默认 3")
    ap.add_argument("--detection-alpha", type=float, default=20.0)
    ap.add_argument("--batch-size-ref", type=int, default=512)
    ap.add_argument("--batch-size-st", type=int, default=512)
    ap.add_argument("--cell-count-cutoff-ref", type=int, default=15)
    ap.add_argument("--cell-count-cutoff-st", type=int, default=5)
    ap.add_argument("--min-common-genes", type=int, default=100)
    ap.add_argument("--max-spatial-types", type=int, default=12)
    ap.add_argument("--no-gpu", action="store_true")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--no-attach", action="store_true", help="跳过回贴 Seurat RDS")
    ap.add_argument("--force-refit", action="store_true", help="忽略 inf_aver 缓存")
    ap.add_argument("--prep-only", action="store_true", help="只做 RDS->counts 转换")
    a = ap.parse_args()

    out_root = os.path.abspath(os.path.expanduser(a.out_dir))
    inp = os.path.join(out_root, "input")
    out = os.path.join(out_root, "output")
    os.makedirs(inp, exist_ok=True); os.makedirs(out, exist_ok=True)
    helper = os.path.join(out_root, "_rds_helper.R")
    rlog = os.path.join(out_root, "r_helper.log")
    report, t0 = [], time.time()

    try:
        # ---------------- 0. 环境 ----------------
        rscript = (a.rscript or os.environ.get("RSCRIPT_BIN") or DEFAULT_RSCRIPT)
        rscript = os.path.expanduser(rscript)
        if not os.path.isfile(rscript):
            alt = shutil.which("Rscript")
            if not alt:
                raise RuntimeError(f"找不到 Rscript: {rscript}；用 --rscript 指定 miloR 的 Rscript")
            rscript = alt
        ref_rds = os.path.abspath(os.path.expanduser(a.ref_rds))
        if not os.path.isfile(ref_rds):
            raise RuntimeError(f"参考 RDS 不存在: {ref_rds}")
        st_paths = expand_st_paths(a.st_rds)
        with open(helper, "w", encoding="utf-8") as f:
            f.write(R_HELPER)
        open(rlog, "w", encoding="utf-8").close()
        accel = "cpu" if a.no_gpu else ("gpu" if torch.cuda.is_available() else "cpu")
        log(f"torch {torch.__version__} | cuda {torch.version.cuda} | accelerator={accel}")
        log(f"Rscript: {rscript}")
        log(f"空转 RDS ({len(st_paths)}): " + ", ".join(os.path.basename(p) for p in st_paths))
        report.append(f"环境      : torch {torch.__version__} / {accel}")
        report.append(f"Rscript   : {rscript}")
        report.append(f"参考 RDS  : {ref_rds}")
        report.append("空转 RDS  : " + " ; ".join(st_paths))

        # ---------------- 1. RDS -> counts ----------------
        log("[1/6] RDS -> counts 三件套")
        prep_args = ["--mode", "prep", "--ref-rds", ref_rds,
                     "--out-dir", inp, "--labels-key", a.labels_key,
                     "--max-per-type", a.max_per_type,
                     "--max-spots-per-sample", a.max_spots_per_sample,
                     "--seed", a.seed]
        for flag, val in [("--ref-assay", a.ref_assay), ("--st-assay", a.st_assay),
                          ("--ref-batch", a.ref_batch), ("--st-batch", a.st_batch),
                          ("--st-bin", a.st_bin)]:
            if val:
                prep_args += [flag, val]
        for p in st_paths:
            prep_args += ["--st-rds", p]
        run_r(rscript, helper, prep_args, rlog)

        prep = parse_prep(os.path.join(inp, "prep_summary.txt"))
        for k in ("ref_assay", "st_assay", "labels_key", "n_samples", "n_batches",
                  "platform", "coords", "st_shape"):
            report.append(f"prep/{k:<12}: {prep.get(k, 'NA')}")
        platform = prep.get("platform", "VisiumV2")
        log(f"  prep: assay={prep.get('ref_assay')}/{prep.get('st_assay')}, "
            f"平台={platform}, 样本={prep.get('n_samples')}, 批次={prep.get('n_batches')}")
        log("  ★ 请核对 input/prep_summary.txt 中 assay 是否为全转录组 raw counts")

        if a.prep_only:
            report += ["", "RESULT: PREP_ONLY_OK", f"耗时: {(time.time()-t0)/60:.1f} min"]
            rc = 0
            txt = "\n".join(report)
            with open(os.path.join(out, "c2l_report.txt"), "w", encoding="utf-8") as f:
                f.write(txt + "\n")
            print("=" * 72); print(txt); print("=" * 72)
            return rc

        # ---------------- 2. 载入 + 过滤 + 交集 ----------------
        from cell2location.models import RegressionModel, Cell2location
        from cell2location.utils.filtering import filter_genes

        log("[2/6] 载入 counts 并做基因过滤/交集")
        ref, st = load_triplet(inp, "ref"), load_triplet(inp, "st")
        labels = ref.obs["cell_type"].astype(str)
        cats = sorted(labels.unique())
        if len(cats) < 2:
            raise RuntimeError(f"细胞类型只有 {len(cats)} 类，无法做细胞类型映射")
        n_batch_ref = ref.obs["batch"].nunique()
        n_batch_st = st.obs["batch"].nunique()
        report.append(f"参考      : {ref.n_obs} cells, {ref.n_vars} genes, "
                      f"{len(cats)} 种细胞类型, {n_batch_ref} 批次")
        report.append(f"空转      : {st.n_obs} spots, {st.n_vars} genes, {n_batch_st} 批次/切片")
        report.append("细胞类型  : " + ", ".join(cats))

        ref = ref[:, filter_genes(ref, cell_count_cutoff=a.cell_count_cutoff_ref,
                                  cell_percentage_cutoff2=0.05,
                                  nonz_mean_cutoff=1.12)].copy()
        st = st[:, filter_genes(st, cell_count_cutoff=a.cell_count_cutoff_st,
                                cell_percentage_cutoff2=0.03,
                                nonz_mean_cutoff=1.12)].copy()
        st_genes = set(st.var_names)
        common = [g for g in ref.var_names if g in st_genes]
        if len(common) < a.min_common_genes:
            raise RuntimeError(f"共同基因仅 {len(common)} 个，检查两边是否都是 gene symbol")
        ref, st = ref[:, common].copy(), st[:, common].copy()
        gc.collect()
        log(f"  共同基因 {len(common)}")
        report.append(f"建模基因  : {len(common)}")

        # ---------------- 3. 参考签名（带一致性校验的缓存） ----------------
        cache = os.path.join(out, "inf_aver.csv")
        meta_cache = os.path.join(out, "inf_aver.meta.json")
        fp = {
            "genes": hashlib.sha256("\n".join(common).encode()).hexdigest()[:20],
            "cats": cats,
            "ref_batches": sorted(ref.obs["batch"].astype(str).unique().tolist()),
            "ref_nobs": int(ref.n_obs),
        }
        use_cache = False
        inf_aver = None
        if os.path.exists(cache) and not a.force_refit:
            ok = False
            if os.path.exists(meta_cache):
                try:
                    with open(meta_cache, encoding="utf-8") as f:
                        ok = json.load(f) == fp
                except Exception:
                    ok = False
            else:  # 兼容旧版缓存：只校验基因列表和细胞类型
                cached = pd.read_csv(cache, index_col=0)
                ok = (list(cached.index) == list(ref.var_names)
                      and sorted(map(str, cached.columns)) == sorted(map(str, cats)))
            if ok:
                inf_aver, use_cache = pd.read_csv(cache, index_col=0), True
            else:
                log("  缓存与当前数据不一致（基因/细胞类型/批次），重新训练参考模型")
        if use_cache:
            log("[3/6] 复用缓存参考签名 inf_aver.csv")
            report.append("参考签名  : 复用缓存 inf_aver.csv")
        else:
            log(f"[3/6] RegressionModel 训练 (max_epochs={a.epochs_ref})")
            ref.obs["cell_type"] = ref.obs["cell_type"].astype("category")
            bk_ref = "batch" if n_batch_ref > 1 else None
            if bk_ref:
                ref.obs["batch"] = ref.obs["batch"].astype("category")
            RegressionModel.setup_anndata(ref, labels_key="cell_type", batch_key=bk_ref)
            mod_ref = RegressionModel(ref)
            train_compat(mod_ref, max_epochs=a.epochs_ref,
                         batch_size=a.batch_size_ref, train_size=1, accelerator=accel)
            ref = mod_ref.export_posterior(ref, sample_kwargs={"batch_size": a.batch_size_ref})
            inf_aver = extract_inf_aver(ref)
            inf_aver.to_csv(cache)
            with open(meta_cache, "w", encoding="utf-8") as f:
                json.dump(fp, f, ensure_ascii=False, indent=1)
            report.append("参考签名  : RegressionModel OK（已缓存）")
            del mod_ref; gc.collect()

        # ---------------- 4. 空间映射 ----------------
        n_cells = a.n_cells_per_location
        if n_cells is None:
            n_cells = 3 if platform == "VisiumHD" else 20
            log(f"  --n-cells-per-location 未指定，按平台 {platform} 自动设为 {n_cells}")
        log(f"[4/6] Cell2location 训练 (epochs={a.epochs_st}, "
            f"N_cells={n_cells}, alpha={a.detection_alpha})")
        bk_st = "batch" if n_batch_st > 1 else None
        if bk_st:
            st.obs["batch"] = st.obs["batch"].astype("category")
        Cell2location.setup_anndata(st, batch_key=bk_st)
        mod_st = Cell2location(st, cell_state_df=inf_aver,
                               N_cells_per_location=n_cells,
                               detection_alpha=a.detection_alpha)
        train_compat(mod_st, max_epochs=a.epochs_st, batch_size=a.batch_size_st,
                     train_size=1, accelerator=accel)
        st = mod_st.export_posterior(st, sample_kwargs={"batch_size": a.batch_size_st})
        report.append("空间映射  : Cell2location OK")
        report.append(f"训练参数  : ref={a.epochs_ref}, st={a.epochs_st}, "
                      f"batch={a.batch_size_ref}/{a.batch_size_st}, "
                      f"N_cells={n_cells}, alpha={a.detection_alpha}")
        del mod_st; gc.collect()

        # ---------------- 5. 注释 spot 表 ----------------
        log("[5/6] 导出注释 spot 表与汇总")
        abund_keys = [k for k in st.obsm.keys() if "cell_abundance" in k]
        if not abund_keys:
            raise RuntimeError(f"未生成丰度矩阵，obsm={list(st.obsm.keys())}")
        key = next((k for k in abund_keys if k.startswith("q05")), abund_keys[0])
        obj = st.obsm[key]
        cols = list(obj.columns) if hasattr(obj, "columns") else list(inf_aver.columns)
        df = pd.DataFrame(np.asarray(obj), index=st.obs_names,
                          columns=[str(c) for c in cols])
        df.index.name = "spot"
        df.to_csv(os.path.join(out, "cell_abundance_q05.csv"))
        report.append(f"丰度矩阵  : {key} {df.shape}")

        keymap = pd.read_csv(os.path.join(inp, "st_key_map.csv"))
        keymap = keymap.drop_duplicates("cell").set_index("cell")
        ab = df.copy()
        ab.columns = ["q05_" + c for c in clean_cols(ab.columns)]
        ann = keymap.reindex(df.index)
        for c in ("sample", "image", "x", "y", "orig_barcode", "batch", "platform"):
            if c not in ann.columns:
                ann[c] = np.nan
        ann["sample"] = ann["sample"].fillna(ann["batch"]).fillna("all")
        ann["image"] = ann["image"].fillna(ann["sample"])
        ann = pd.concat([ann[["sample", "image", "orig_barcode", "x", "y"]], ab], axis=1)
        ann.insert(0, "spot", ann.index)
        dom_col = ab.columns[ab.values.argmax(axis=1)]
        ann["dominant_type"] = [c[4:] for c in dom_col]
        ann["dominant_q05"] = ab.values.max(axis=1)
        ann["total_q05"] = ab.values.sum(axis=1)
        ann["dominant_fraction"] = ann["dominant_q05"] / ann["total_q05"].replace(0, np.nan)
        ann.to_csv(os.path.join(out, "spot_annotations.csv"), index=False)
        log(f"  spot_annotations.csv: {ann.shape[0]} spots x {ann.shape[1]} 列")

        by_dir = os.path.join(out, "annotated_spots_by_sample")
        os.makedirs(by_dir, exist_ok=True)
        for s, sub in ann.groupby("sample"):
            sub.to_csv(os.path.join(by_dir, f"{re.sub(r'[^A-Za-z0-9._-]', '_', str(s))}.csv"),
                       index=False)
        comp = (ann.groupby(["sample", "dominant_type"]).size()
                .unstack(fill_value=0))
        comp.to_csv(os.path.join(out, "sample_annotation_summary.csv"))
        log("  各样本主导类型统计:\n" + comp.to_string())

        summ = df.groupby(st.obs["batch"].astype(str).values).mean()
        summ.index.name = "slice"
        summ.columns = clean_cols(summ.columns)
        summ.to_csv(os.path.join(out, "slice_summary.csv"))

        # h5ad 里也带上注释列，方便 scanpy 后续分析
        for c in ("dominant_type", "dominant_q05", "total_q05", "dominant_fraction"):
            st.obs[c] = ann[c].values
        for c in ("x", "y", "image", "sample"):
            st.obs[c] = ann[c].values
        st.write(os.path.join(out, "c2l_result.h5ad"))

        # ---------------- 6. 可视化 + 回贴 ----------------
        log("[6/6] 绘图与回贴")
        if not a.no_plot:
            try:
                plot_overview(df.set_axis(clean_cols(df.columns), axis=1), summ, out)
                report.append("汇总图    : abundance_overview.png")
            except Exception as e:
                log(f"  汇总图失败（不影响结果）: {e}")
            try:
                plot_composition(ann, out)
                report.append("组成图    : dominant_composition.png")
            except Exception as e:
                log(f"  组成图失败（不影响结果）: {e}")
            try:
                if plot_spatial(ann, out, max_types=a.max_spatial_types):
                    report.append("空间图    : c2l_spatial_abundance.pdf / c2l_dominant_spatial.png")
            except Exception as e:
                log(f"  空间图失败（不影响结果）: {e}")
        else:
            log("  --no-plot：跳过全部可视化")

        if not a.no_attach:
            try:
                run_r(rscript, helper,
                      ["--mode", "attach", "--abun-csv",
                       os.path.join(out, "cell_abundance_q05.csv"),
                       "--key-csv", os.path.join(inp, "st_key_map.csv"),
                       "--out-rds-dir", os.path.join(out, "seurat_annotated")], rlog)
                report.append("回贴 RDS  : output/seurat_annotated/*.rds")
            except Exception as e:
                log(f"  回贴失败（不影响 CSV/图）: {e}")

        report.append(f"耗时      : {(time.time()-t0)/60:.1f} min")
        report.append("")
        report.append("主要输出  : output/spot_annotations.csv（注释 spot 总表）")
        report.append("           : output/annotated_spots_by_sample/*.csv（按样本拆分）")
        report.append("           : output/c2l_spatial_abundance.pdf、c2l_dominant_spatial.png")
        report.append("           : output/abundance_overview.png、dominant_composition.png")
        report.append("")
        report.append("RESULT: PASS")
        rc = 0
    except Exception as e:
        report += ["", f"RESULT: FAIL ({type(e).__name__}: {e})",
                   "---- traceback ----", traceback.format_exc()]
        rc = 1

    txt = "\n".join(report)
    with open(os.path.join(out, "c2l_report.txt"), "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    print("=" * 72); print(txt); print("=" * 72)
    return rc


if __name__ == "__main__":
    sys.exit(main())
