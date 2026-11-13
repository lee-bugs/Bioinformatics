# =====================================================================
# gene_viz_3d.R
# ---------------------------------------------------------------------
# Gene expression visualisation from an RDS file, rendered as PSEUDO-3D
# gradient-filled violins or box-and-whisker cubes (same rendering style as
# ykt6_expression_plot.R), with jittered dots and compact significance
# caps. Pure base R + grDevices for drawing.
#
# The Y axis is trimmed to the 1st-99th percentile by default: with
# single-cell data a few dropout zeros otherwise flatten every shape.
#
# ASCII-only source (safe under the Windows "C" locale).
#
# COMMAND LINE
#   Rscript gene_viz_3d.R --rds pbmc.rds --gene YKT6 --group celltype \
#       --type violin --out YKT6_violin.png \
#       --dot_size 0.8 --dot_alpha 0.6 --fill_alpha 0.30 --max_cells 1500
#
#   Rscript gene_viz_3d.R --rds obj.rds --gene YKT6,CD4 \
#       --group group,celltype --type box --show_p 1
#
#   # Text input is auto-detected; --gene and --group become optional.
#   Rscript gene_viz_3d.R --file expr.csv                  # all genes, by Type
#   Rscript gene_viz_3d.R --file expr.tsv --gene YKT6,CD4 --type box
#
#   # Nested big/small groups: name the COLUMNS explicitly, same syntax for
#   # text tables and RDS metadata. --nest is implied; do not pass --group.
#   Rscript gene_viz_3d.R --file expr.csv --gene YKT6 \
#       --outer CellLine --inner Treatment --type box
#   Rscript gene_viz_3d.R --rds obj.rds --gene YKT6 \
#       --outer group --inner celltype
#   # Optional 3rd level: --page Cohort (one image per level)
#
# INSIDE R
#   source("gene_viz_3d.R")
#   plot_gene_3d("pbmc.rds", gene = "YKT6", group_by = "celltype",
#                type = "violin", dot_size = 0.8, dot_alpha = 0.6)
#   plot_gene_3d("expr.csv", "YKT6", outer_group = "CellLine",
#                inner_group = "Treatment", type = "box")
#
# RDS FORMATS ACCEPTED
#   * Seurat object           : expr from assay/layer, meta from @meta.data
#   * SingleCellExperiment    : logcounts + colData
#   * data.frame (wide)       : gene name as a column, group columns present
#   * matrix (genes x cells)  : needs a companion meta (see list below)
#   * list                    : list(expr = <genes x cells matrix>,
#                                   meta = <cells x groups data.frame>)
#   A list is the safest portable format when Seurat is not installed.
#
# TEXT FILES (CSV / TSV / TXT, whitespace or , ; | separated, .gz allowed):
#   ID,<gene1>,<gene2>,...,Type
#   S1,7.1,3.2,...,Normal
#   ID and Type are matched case-insensitively. All other numeric columns
#   are genes; non-numeric columns can be used via --group. When --gene is
#   omitted every numeric gene column gets its own figure (capped at 20).
#
# group_by of length 1 -> one axis with all levels.
# group_by of length 2 -> first column = panels, second = groups on the axis
#                         (reproduces the two-panel layout of the original).
# =====================================================================


# ---------------------------------------------------------------------
# 0. Defaults
# ---------------------------------------------------------------------
# Muted, layered palette - low saturation reads far better than neon on a
# light canvas. Derived from the ZS palette set used in the earlier project.

DEFAULT_PAL <- c(
  "#556136", "#FEF9DE", "#B2622C",
  "#85C2E0", "#092A3D", "#F77C26",
  "#9ACBFB", "#FCA7C4", "#C5ABD3",
  "#F27F0C", "#053F5C", "#9FE7F5",
  "#1E0A26", "#F2DB66", "#9CD9D3",
  "#61589E", "#EDD7DF", "#BF5C78",
  "#161616", "#74B5B3", "#FFFFFF",
  "#273F75", "#653136", "#A98D42",
  "#2F443A", "#507550", "#D5D1C7")




C_THEME <- list(
  ink     = "#2C3A45", axis = "#84939E", median = "#2C3A45",
  grid    = "#E9EEF1", bg   = "#FCFDFE",
  shadow  = "#60788A", whisker = "#9AA8B2", muted = "#8794A0"
)

.darken <- function(hex, f = 0.68) {
  z <- col2rgb(hex) / 255
  rgb(z[1, 1] * f, z[2, 1] * f, z[3, 1] * f)
}

.lighten <- function(hex, f = 0.58) {
  z <- col2rgb(hex) / 255
  rgb(z[1, 1] + (1 - z[1, 1]) * f,
      z[2, 1] + (1 - z[2, 1]) * f,
      z[3, 1] + (1 - z[3, 1]) * f)
}

# Vertical gradient built from horizontal bands (base R has no gradient fill)
.grad_bands <- function(y0, y1, x0_fn, x1_fn, col_bot, col_top, alpha,
                        n = 28) {
  cuts <- seq(y0, y1, length.out = n + 1)
  cols <- colorRampPalette(c(col_bot, col_top))(n)
  for (i in seq_len(n)) {
    a <- cuts[i]; b <- cuts[i + 1]
    polygon(c(x0_fn(a), x1_fn(a), x1_fn(b), x0_fn(b)),
            c(a, a, b, b),
            col = adjustcolor(cols[i], alpha.f = alpha), border = NA)
  }
}


# ---------------------------------------------------------------------
# 0b. Input auto-detection: RDS vs delimited text
# ---------------------------------------------------------------------
# RDS -> unchanged path. Text -> read as a wide table and hand it to the
# existing data.frame branch of .extract(), which already treats every
# non-group numeric column as a gene.
.guess_kind <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return("rds")
  if (ext %in% c("csv", "tsv", "txt", "tab", "dat", "text")) return("table")
  # A .gz text file also starts with the gzip magic, so resolve it by name
  # before the binary sniff (otherwise it would look like compressed RDS).
  if (grepl("\\.gz$", tolower(path))) {
    inner <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
    return(if (inner == "rds") "rds" else "table")
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  h <- readBin(con, "raw", 8L)
  # gzip magic (compressed RDS) or the RDX/RDA binary header
  if (length(h) >= 2L && h[1] == as.raw(0x1f) && h[2] == as.raw(0x8b)) {
    return("rds")
  }
  vis <- h[h >= as.raw(0x20) & h <= as.raw(0x7e)]
  asc <- if (length(vis)) rawToChar(vis) else ""
  if (startsWith(asc, "RDX") || startsWith(asc, "RDA")) return("rds")
  "table"
}

.guess_sep <- function(lines) {
  cand <- c(comma = ",", tab = "\t", semi = ";", pipe = "|")
  n <- vapply(cand, function(s) {
    sum(vapply(strsplit(lines, s, fixed = TRUE), length, 1L) - 1L)
  }, 0)
  if (max(n) <= 0) return("")            # any whitespace
  cand[[names(which.max(n))]]
}

.table_con <- function(path) {
  enc <- "UTF-8-BOM"                     # harmless when no BOM is present
  if (grepl("\\.gz$", tolower(path))) {
    gzfile(path, open = "rt", encoding = enc)
  } else {
    file(path, open = "rt", encoding = enc)
  }
}

.read_table <- function(path) {
  # A single-byte (C) locale cannot even read UTF-8 bytes; read.table would
  # then fail later with a misleading "no numeric gene columns" message.
  if (!grepl("\\.gz$", tolower(path))) {
    head_bytes <- readBin(path, what = "raw", n = 4096L)
    if (any(head_bytes > as.raw(0x7F)) &&
        !isTRUE(l10n_info()[["MBCS"]]) &&
        !grepl("UTF", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
      stop("Non-ASCII bytes (e.g. Chinese column names) in a single-byte ",
           "locale. Use ASCII column names, or restart R with a UTF-8 locale.",
           call. = FALSE)
    }
  }
  con <- .table_con(path)
  ln <- readLines(con, n = 10L, warn = FALSE)
  close(con)
  ln <- trimws(ln)
  ln <- ln[nchar(ln) > 0 & !startsWith(ln, "#")]
  if (!length(ln)) stop("Text file is empty: ", path, call. = FALSE)
  sep <- .guess_sep(ln)

  rd <- function(s) {
    cc <- .table_con(path)
    on.exit(close(cc), add = TRUE)
    utils::read.table(cc, header = TRUE, sep = s, quote = "\"",
                      check.names = FALSE, stringsAsFactors = FALSE,
                      comment.char = "", fill = TRUE)
  }
  d <- tryCatch(rd(sep), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 2L) d <- tryCatch(rd(""), error = function(e) NULL)
  if (is.null(d)) stop("Could not parse text file: ", path, call. = FALSE)
  .validate_table(d)
}

.validate_table <- function(d) {
  nm <- trimws(names(d))
  if (length(nm)) nm[1] <- sub("^\ufeff", "", nm[1])
  names(d) <- nm
  up <- toupper(nm)
  i_id <- which(up == "ID")[1]
  i_ty <- which(up == "TYPE")[1]
  if (is.na(i_id)) {
    stop("Text input needs an ID column. Found: ",
         paste(nm, collapse = ", "), call. = FALSE)
  }
  d[[i_id]] <- as.character(d[[i_id]])
  has_type <- !is.na(i_ty)
  if (has_type) d[[i_ty]] <- as.character(d[[i_ty]])
  # Type is conventional but optional: a nested design may only have
  # Sex / Age / CellLine / Treatment style grouping columns.
  others <- setdiff(seq_along(d), c(i_id, if (has_type) i_ty))
  # numeric columns = genes; non-numeric ones stay as extra grouping columns
  is_gene <- vapply(others, function(j) {
    x <- d[[j]]
    if (is.numeric(x)) return(TRUE)
    z <- suppressWarnings(as.numeric(as.character(x)))
    ok <- sum(!is.na(x))
    ok > 0 && sum(!is.na(z)) >= 0.8 * ok
  }, TRUE)
  for (j in others[is_gene]) {
    if (!is.numeric(d[[j]])) {
      d[[j]] <- suppressWarnings(as.numeric(as.character(d[[j]])))
    }
  }
  genes <- nm[others[is_gene]]
  genes <- genes[nzchar(genes)]
  if (!length(genes)) {
    stop("No numeric gene columns found between ID and Type.", call. = FALSE)
  }
  extra <- nm[others[!is_gene]]
  attr(d, "gene_cols") <- genes
  attr(d, "extra_cols") <- extra
  attr(d, "default_group") <- if (has_type) "Type" else {
    if (length(extra)) extra[1] else NA_character_
  }
  d
}

.prepare_input <- function(src, gene, group_by, max_auto_genes,
                           outer_group = NULL, inner_group = NULL,
                           page_group = NULL) {
  # Explicit big/small group column names. Works identically for a text
  # table and for Seurat metadata: the names are resolved by the same
  # column checks below, so no RDS-specific code is needed.
  explicit <- !is.null(outer_group) || !is.null(inner_group) ||
              !is.null(page_group)
  if (explicit) {
    if (is.null(outer_group) || is.null(inner_group)) {
      stop("Give BOTH --outer (big group) and --inner (small group) ",
           "column names.", call. = FALSE)
    }
    if (!is.null(group_by)) {
      stop("Use --outer/--inner OR --group, not both.", call. = FALSE)
    }
    group_by <- as.character(c(page_group, outer_group, inner_group))
  }

  kind <- NULL
  if (is.character(src) && length(src) == 1L && file.exists(src)) {
    kind <- .guess_kind(src)
  } else if (is.data.frame(src)) {
    kind <- "table"
  }

  if (identical(kind, "table")) {
    d <- if (is.data.frame(src)) .validate_table(src) else .read_table(src)
    genes <- attr(d, "gene_cols")
    if (is.null(gene)) {
      if (length(genes) > max_auto_genes) {
        stop("Detected ", length(genes), " gene columns. Use --gene to pick ",
             "at most ", max_auto_genes, " at a time, or raise ",
             "--max_auto_genes.", call. = FALSE)
      }
      gene <- genes
    } else {
      gene <- as.character(gene)
      hit <- toupper(gene) %in% toupper(genes)
      if (any(!hit)) {
        stop("Gene column(s) not found: ",
             paste(gene[!hit], collapse = ", "), ". Available: ",
             paste(utils::head(genes, 12L), collapse = ", "), call. = FALSE)
      }
      gene <- vapply(gene, function(g) {
        genes[which(toupper(genes) == toupper(g))[1]]
      }, "", USE.NAMES = FALSE)
    }
    if (is.null(group_by)) {
      group_by <- attr(d, "default_group")
      if (is.na(group_by)) {
        stop("No grouping column found. Add a Type column or pass --group.",
             call. = FALSE)
      }
    }
    group_by <- as.character(group_by)
    miss <- setdiff(group_by, names(d))
    if (length(miss)) {
      stop("Group column(s) not found: ", paste(miss, collapse = ", "),
           ". Available: ", paste(names(d), collapse = ", "), call. = FALSE)
    }
    return(list(src = d, gene = gene, group_by = group_by,
                kind = "table", genes = genes, forced_nest = explicit))
  }

  # RDS path or an already-loaded object: exactly the old behaviour
  if (is.null(gene)) stop("--gene is required for RDS input", call. = FALSE)
  if (is.null(group_by)) stop("--group is required for RDS input", call. = FALSE)
  list(src = src, gene = as.character(gene),
       group_by = as.character(group_by), kind = "rds", genes = NULL,
       forced_nest = explicit)
}


# ---------------------------------------------------------------------
# 1. RDS -> expression vector + metadata
# ---------------------------------------------------------------------
.rds_read <- function(rds) {
  if (is.character(rds)) {
    if (length(rds) != 1L || !file.exists(rds)) {
      stop("RDS file not found: ", rds, call. = FALSE)
    }
    readRDS(rds)
  } else rds
}

.pick_matrix <- function(candidates) {
  for (nm in c("expr", "expression", "exprs", "logcounts", "data",
               "counts", "mat", "matrix", "counts.mat")) {
    hit <- names(candidates)[tolower(names(candidates)) == nm]
    if (length(hit)) return(candidates[[hit[1]]])
  }
  for (x in candidates) if (is.matrix(x) && is.numeric(x)) return(x)
  NULL
}

.seurat_matrix <- function(obj, assay, layer) {
  # Seurat v5 / SeuratObject
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    a <- if (is.null(assay)) "RNA" else assay
    for (ly in unique(c(layer, "data", "counts"))) {
      z <- tryCatch(
        SeuratObject::GetAssayData(obj, assay = a, layer = ly),
        error = function(e) NULL)
      if (!is.null(z) && nrow(z) > 0) return(as.matrix(z))
    }
  }
  if (requireNamespace("Seurat", quietly = TRUE)) {
    z <- tryCatch(Seurat::GetAssayData(obj, slot = "data"),
                  error = function(e) NULL)   # Seurat v4 and earlier
    if (!is.null(z)) return(as.matrix(z))
  }
  # Last resort: direct S4 slot access (works if class defs are loadable)
  z <- tryCatch({
    a <- if (is.null(assay)) names(obj@assays)[1] else assay
    as <- obj@assays[[a]]
    if (!is.null(as@data)) as.matrix(as@data) else NULL
  }, error = function(e) NULL)
  z
}

.extract <- function(obj, assay, layer) {
  cls <- paste(class(obj), collapse = ",")

  if (grepl("Seurat", cls)) {
    mat <- .seurat_matrix(obj, assay, layer)
    meta <- tryCatch(as.data.frame(obj@meta.data), error = function(e) NULL)
    if (is.null(mat) || is.null(meta)) {
      stop("Detected a Seurat object but could not read expression/meta.data. ",
           "Install SeuratObject, or re-save the RDS as ",
           "list(expr = <matrix>, meta = <data.frame>).", call. = FALSE)
    }
    return(list(mat = mat, meta = meta))
  }

  if (grepl("SingleCellExperiment", cls)) {
    mat <- tryCatch(as.matrix(SummarizedExperiment::assay(obj, "logcounts")),
                    error = function(e) NULL)
    if (is.null(mat)) mat <- tryCatch(
      as.matrix(SummarizedExperiment::assay(obj, 1)), error = function(e) NULL)
    meta <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)),
                     error = function(e) NULL)
    if (!is.null(mat) && !is.null(meta)) return(list(mat = mat, meta = meta))
  }

  if (is.data.frame(obj)) {
    # wide format: gene name is a column
    return(list(mat = NULL, meta = obj, wide = TRUE))
  }

  if (is.matrix(obj)) return(list(mat = obj, meta = NULL))

  if (is.list(obj)) {
    mat <- .pick_matrix(obj)
    meta <- NULL
    for (nm in c("meta", "metadata", "colData", "pheno", "phenoData", "pd")) {
      hit <- names(obj)[tolower(names(obj)) == tolower(nm)]
      if (length(hit)) {
        meta <- tryCatch(as.data.frame(obj[[hit[1]]]), error = function(e) NULL)
        if (!is.null(meta)) break
      }
    }
    if (!is.null(mat)) return(list(mat = mat, meta = meta))
  }
  stop("Unsupported RDS object of class '", cls, "'. Supported: Seurat, ",
       "SingleCellExperiment, data.frame (wide), matrix, or a list with ",
       "expr + meta elements.", call. = FALSE)
}

.find_gene <- function(mat, gene) {
  rn <- rownames(mat)
  if (is.null(rn)) stop("Expression matrix has no rownames (gene names).",
                        call. = FALSE)
  i <- which(rn == gene)
  if (length(i) == 0L) i <- which(toupper(rn) == toupper(gene))
  if (length(i) == 0L) {
    hint <- utils::head(agrep(gene, rn, value = TRUE, max.distance = 0.25), 5)
    extra <- if (length(hint)) paste0(" Did you mean: ", paste(hint, collapse = ", ")) else ""
    stop("Gene '", gene, "' not found in the expression matrix.", extra,
         call. = FALSE)
  }
  as.numeric(mat[i[1], ])
}

.check_cols <- function(meta, cols) {
  miss <- setdiff(cols, names(meta))
  if (length(miss)) {
    stop("Group column(s) not in metadata: ", paste(miss, collapse = ", "),
         ". Available: ", paste(names(meta), collapse = ", "), call. = FALSE)
  }
}

.levels_of <- function(x, explicit = NULL) {
  if (!is.null(explicit)) return(as.character(explicit))
  if (is.factor(x)) as.character(levels(x)) else sort(unique(as.character(x)))
}

.build_frame <- function(rds, gene, group_by, assay, layer, transform,
                         group_levels) {
  obj  <- .rds_read(rds)
  info <- .extract(obj, assay, layer)

  if (isTRUE(info$wide)) {
    meta <- info$meta
    .check_cols(meta, group_by)
    cn <- names(meta)
    i <- which(cn == gene)
    if (length(i) == 0L) i <- which(toupper(cn) == toupper(gene))
    if (length(i) == 0L) {
      stop("Gene '", gene, "' is not a column of the data.frame. ",
           "Available columns: ", paste(cn, collapse = ", "), call. = FALSE)
    }
    expr <- as.numeric(meta[[i[1]]])
  } else {
    mat <- info$mat
    meta <- info$meta
    if (is.null(meta)) {
      stop("A bare expression matrix has no grouping metadata. Save the RDS ",
           "as list(expr = mat, meta = meta) instead.", call. = FALSE)
    }
    meta <- as.data.frame(meta)
    expr <- .find_gene(mat, gene)
    if (ncol(mat) != nrow(meta)) {
      stop("Dimension mismatch: expression has ", ncol(mat), " cells, metadata ",
           "has ", nrow(meta), " rows.", call. = FALSE)
    }
  }
  .check_cols(meta, group_by)

  if (transform == "log1p") expr <- log1p(expr)
  if (transform == "log2p1") expr <- log2(expr + 1)

  d <- data.frame(expr = expr, stringsAsFactors = FALSE)
  for (cc in group_by) d[[cc]] <- as.character(meta[[cc]])
  d <- d[is.finite(d$expr), , drop = FALSE]
  for (cc in group_by) d <- d[!is.na(d[[cc]]), , drop = FALSE]
  if (nrow(d) < 3L) stop("Fewer than 3 usable cells after filtering.", call. = FALSE)

  lev <- list()
  for (j in seq_along(group_by)) {
    lev[[group_by[j]]] <- .levels_of(meta[[group_by[j]]], group_levels[[j]])
  }
  list(d = d, levels = lev, n_total = nrow(d))
}


# ---------------------------------------------------------------------
# 2. Statistics
# ---------------------------------------------------------------------
.box_stat <- function(v, whisker = c("tukey", "range")) {
  whisker <- match.arg(whisker)
  q <- as.numeric(quantile(v, c(0.25, 0.5, 0.75), type = 7))
  if (whisker == "tukey") {
    iqr <- q[3] - q[1]
    wlo <- max(min(v), q[1] - 1.5 * iqr)
    whi <- min(max(v), q[3] + 1.5 * iqr)
  } else {
    wlo <- min(v); whi <- max(v)
  }
  c(q1 = q[1], med = q[2], q3 = q[3], wlo = wlo, whi = whi)
}

.p_stars <- function(p) {
  if (is.na(p)) "ns" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else
    if (p < 0.05) "*" else "ns"
}

.test_two <- function(a, b, test) {
  p <- if (test == "t") {
    suppressWarnings(t.test(a, b)$p.value)
  } else {
    suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  }
  list(p = p, stars = .p_stars(p))
}

.test_many <- function(groups, test) {
  if (test == "t") {
    suppressWarnings(oneway.test(expr ~ g, data = groups)$p.value)
  } else {
    suppressWarnings(kruskal.test(expr ~ g, data = groups)$p.value)
  }
}


# ---------------------------------------------------------------------
# 3. Drawing helpers (shared style with ykt6_expression_plot.R)
# ---------------------------------------------------------------------
.soft_ellipse <- function(cx, cy, rx, ry, alpha = 0.06) {
  t <- seq(0, 2 * pi, length.out = 90)
  polygon(cx + rx * cos(t), cy + ry * sin(t),
          col = adjustcolor(C_THEME$shadow, alpha.f = alpha), border = NA)
}

.glossy_dots <- function(x, y, fill, edge, cex, alpha) {
  if (length(x) == 0L) return(invisible(NULL))
  # Small flat dots with a hairline white rim. Glossy spheres turn into
  # visual noise as soon as a group has more than ~100 cells.
  points(x, y, pch = 21, cex = cex, lwd = 0.65,
         col = adjustcolor("#FFFFFF", alpha.f = 0.75 * alpha),
         bg  = adjustcolor(edge, alpha.f = alpha))
}

.y_axis <- function(ticks, y_bottom, y_top, ylab, text_cex, x = 4.2) {
  segments(x, y_bottom, x, y_top,
           col = adjustcolor(C_THEME$axis, alpha.f = 0.70), lwd = 1)
  lab <- vapply(ticks, function(v) formatC(v, format = "g", digits = 4), "")
  tc <- 0.86 * text_cex
  # Place the vertical axis title to the LEFT of the widest tick label.
  # A hard-coded offset collides as soon as ticks are 3+ digits (800, 600)
  # or scientific notation.
  sw <- suppressWarnings(max(strwidth(lab, units = "user", cex = tc)))
  if (!is.finite(sw) || sw <= 0) sw <- 2.5
  for (i in seq_along(ticks)) {
    segments(x - 0.35, ticks[i], x, ticks[i], col = C_THEME$axis, lwd = 1)
    text(x - 0.75, ticks[i], lab[i], adj = 1, cex = tc,
         col = C_THEME$muted)
  }
  lab_x <- max(-4.2, x - 0.75 - sw - 1.8 * text_cex)
  text(lab_x, (y_bottom + y_top) / 2, ylab, srt = 90,
       cex = 1.0 * text_cex, font = 2, col = C_THEME$ink)
  invisible(lab_x)
}

.shorten <- function(s, n = 16) {
  s <- as.character(s)
  ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "~"), s)
}


# ---------------------------------------------------------------------
# 4. Pseudo-3D cube (boxplot)
# ---------------------------------------------------------------------
.draw_cube <- function(cx, st, v_plot, fill, edge, w, dx, dy,
                       opt, seed) {
  x0 <- cx - w / 2; x1 <- cx + w / 2
  q1 <- st["q1"]; q3 <- st["q3"]; med <- st["med"]
  xm <- cx + dx / 2
  fa <- opt$fill_alpha
  deep <- .darken(fill, 0.76)
  lite <- .lighten(fill, 0.66)

  # depth faces, deliberately quiet
  polygon(c(x1, x1 + dx, x1 + dx, x1), c(q1, q1 + dy, q3 + dy, q3),
          col = adjustcolor(deep, alpha.f = fa * 0.50), border = NA)
  polygon(c(x0, x0 + dx, x1 + dx, x1), c(q3, q3 + dy, q3 + dy, q3),
          col = adjustcolor(lite, alpha.f = fa * 0.80), border = NA)

  wc <- C_THEME$whisker
  segments(xm, q3, xm, st["whi"], col = wc, lwd = 1.1)
  segments(xm, q1, xm, st["wlo"], col = wc, lwd = 1.1)
  cap <- w * 0.16
  segments(xm - cap, st["whi"], xm + cap, st["whi"], col = wc, lwd = 1.4)
  segments(xm - cap, st["wlo"], xm + cap, st["wlo"], col = wc, lwd = 1.4)

  # front face: vertical gradient, dark at the bottom
  .grad_bands(q1, q3, function(z) x0, function(z) x1,
              fill, lite, fa)

  # median: a crisp line, not a black slab
  polygon(c(x0, x1, x1 + dx, x0 + dx), c(med, med, med + dy, med + dy),
          col = adjustcolor(deep, alpha.f = 0.14), border = NA)
  segments(x0, med, x1, med, col = .darken(fill, 0.42), lwd = 1.9)

  ec <- adjustcolor(edge, alpha.f = 0.92)
  rect(x0, q1, x1, q3, border = ec, lwd = 1.1)
  lines(c(x0, x0 + dx, x1 + dx, x1), c(q3, q3 + dy, q3 + dy, q3),
        col = ec, lwd = 1)
  lines(c(x1, x1 + dx, x1 + dx), c(q1, q1 + dy, q3 + dy), col = ec, lwd = 1)
  lines(c(x0 + dx, x1 + dx), c(q1 + dy, q1 + dy),
        col = adjustcolor(ec, alpha.f = 0.45), lwd = 0.9)

  .draw_dots(cx, v_plot, fill, edge, w, dx, dy, opt, seed)
}

# 3D violin: back silhouette + side wall + front silhouette + median slab
.draw_violin <- function(cx, v_all, v_plot, fill, edge, w, dx, dy,
                         opt, seed) {
  st <- .box_stat(v_all)
  fa <- opt$fill_alpha
  den <- tryCatch(
    density(v_all, n = 160, from = min(v_all), to = max(v_all),
            adjust = opt$bw_adjust),
    error = function(e) NULL)
  if (is.null(den) || all(!is.finite(den$y)) || max(den$y) <= 0) {
    den <- list(x = c(min(v_all), max(v_all)), y = c(0.5, 0.5))
  }
  yy <- den$x
  hw <- den$y / max(den$y) * w / 2
  # Clip the silhouette to the visible Y window, otherwise long dropout
  # tails get cut by the plot region and look like broken lines.
  if (!is.null(opt$ylo)) {
    kk <- yy >= opt$ylo & yy <= opt$yhi
    if (sum(kk) >= 3L) { yy <- yy[kk]; hw <- hw[kk] }
  }
  if (all(hw < 1e-6)) hw <- rep(w * 0.15, length(hw))
  hw[hw < max(hw) * 0.04] <- 0           # suppress hair-thin tails
  # A curved silhouette must NOT be shifted vertically: that reads as two
  # stacked violins. Keep only a tiny vertical offset, use dx for depth.
  # Depth scales with local width: a thin tail must not become a solid slab.
  dyv <- dy * 0.12
  dfac <- pmin(1, hw / max(0.15 * max(hw), 1e-6))
  dxl <- dx * dfac
  dyl <- dyv * dfac

  deep <- .darken(fill, 0.78)
  lite <- .lighten(fill, 0.68)

  # back silhouette and side wall: a whisper of depth, not a double image
  polygon(c(cx + hw + dxl, rev(cx - hw + dxl)), c(yy + dyl, rev(yy + dyl)),
          col = adjustcolor(deep, alpha.f = fa * 0.16), border = NA)
  polygon(c(cx + hw, rev(cx + hw + dxl)), c(yy, rev(yy + dyl)),
          col = adjustcolor(deep, alpha.f = fa * 0.24), border = NA)

  # front silhouette: gradient clipped to the density outline
  hw_at <- function(v) approx(yy, hw, xout = v, rule = 2)$y
  .grad_bands(min(yy), max(yy),
              function(z) cx - hw_at(z), function(z) cx + hw_at(z),
              fill, lite, fa, n = 44)
  px <- c(cx + hw, rev(cx - hw), cx + hw[1])
  py <- c(yy, rev(yy), yy[1])
  lines(px, py, col = adjustcolor(edge, alpha.f = 0.90), lwd = 1.2)

  # median + quartiles: thin marks only
  hm <- hw_at(st["med"])
  segments(cx - hm, st["med"], cx + hm, st["med"],
           col = .darken(fill, 0.42), lwd = 1.9)
  for (q in c(st["q1"], st["q3"])) {
    hq <- hw_at(q) * 0.72
    segments(cx - hq, q, cx + hq, q,
             col = adjustcolor(C_THEME$median, alpha.f = 0.42), lty = 3, lwd = 1)
  }
  # dots are confined to the local density width of the violin
  .draw_dots(cx, v_plot, fill, edge, w, dx, dyv, opt, seed, hw_at = hw_at)
}

.draw_dots <- function(cx, v_plot, fill, edge, w, dx, dy, opt, seed,
                       hw_at = NULL) {
  if (!opt$show_points || length(v_plot) == 0L) return(invisible(NULL))
  set.seed(seed)
  half <- if (is.null(hw_at)) {
    rep(w * 0.5 * opt$jitter * 0.85, length(v_plot))
  } else {
    pmax(0.02, hw_at(v_plot) * opt$jitter * 0.82)
  }
  jx <- cx + runif(length(v_plot), -1, 1) * half
  depth <- runif(length(v_plot), 0, dx * 0.8)
  jy <- v_plot + depth / max(dx, 1e-6) * dy * 0.5
  .glossy_dots(jx + depth * 0.05, jy, fill, edge,
               cex = opt$dot_size, alpha = opt$dot_alpha)
}


# ---------------------------------------------------------------------
# 5. Main plotting engine
# ---------------------------------------------------------------------
.plot_one <- function(d, gene, group_by, lev, opt) {
  # ---- panel / group structure ---------------------------------------
  if (length(group_by) == 2L) {
    panels <- lev[[group_by[1]]]
    groups <- lev[[group_by[2]]]
    d <- d[d[[group_by[1]]] %in% panels & d[[group_by[2]]] %in% groups, ]
  } else {
    panels <- ""
    groups <- lev[[group_by[1]]]
    d <- d[d[[group_by[1]]] %in% groups, ]
  }
  nest <- isTRUE(opt$nest)
  P <- length(panels); G <- length(groups)
  if (G < 1L) stop("No groups left after applying group_levels.", call. = FALSE)
  max_boxes <- if (nest) 30L else 16L
  if (P * G > max_boxes) {
    stop("Too many boxes (", P * G, "). Use group_levels to subset, or split ",
         "by gene.", call. = FALSE)
  }
  if (nest && length(group_by) != 2L) {
    stop("Nested mode needs exactly 2 grouping columns (outer, inner).",
         call. = FALSE)
  }

  getv <- function(p, g) {
    if (!nest && P == 1L && panels[1] == "") {
      sel <- d[[group_by[1]]] == g
    } else {
      sel <- d[[group_by[1]]] == p & d[[group_by[2]]] == g
    }
    as.numeric(d$expr[sel])
  }
  v_all <- list(); v_plot <- list(); n_by <- matrix(0L, P, G)
  for (p in seq_len(P)) for (g in seq_len(G)) {
    k <- (p - 1) * G + g
    v <- getv(panels[p], groups[g])
    v_all[[k]] <- v
    n_by[p, g] <- length(v)
    v_plot[[k]] <- if (length(v) > opt$max_cells) {
      set.seed(opt$seed + k); sample(v, opt$max_cells)
    } else v
  }
  keep <- vapply(v_all, length, 0L) >= 2L
  if (!any(keep)) stop("Every group has fewer than 2 cells.", call. = FALSE)

  # ---- canvas geometry ------------------------------------------------
  allv <- unlist(v_all)
  # Focus the axis on the bulk of the distribution. In single-cell data a
  # handful of dropout zeros otherwise squash 99% of cells into a thin disc.
  nz <- allv[allv > 1e-8]
  if (isTRUE(opt$y_zoom)) {
    # If a sizeable block of exact zeros exists (scRNA dropout), the 1st
    # percentile is still 0 and a plain percentile trim does nothing.
    # Anchor the view on the non-zero distribution instead.
    base <- if (length(nz) >= 0.5 * length(allv) && length(nz) > 10L) nz else allv
    q <- as.numeric(quantile(base, c(0.01, 0.99), type = 7))
    ylo <- q[1]; yhi <- q[2]
  } else {
    ylo <- min(allv); yhi <- max(allv)
  }
  # Floor must be RELATIVE to the data. A fixed 0.4 floor would blow up the
  # axis for 1e-5-scale values and flatten every shape into a line.
  sp <- max(yhi - ylo, 1e-6 * max(abs(c(ylo, yhi)), 1e-8))
  ylo <- ylo - 0.10 * sp
  yhi <- yhi + 0.10 * sp
  span2 <- yhi - ylo
  ticks <- pretty(c(ylo, yhi), n = 5)
  ticks <- ticks[ticks >= ylo & ticks <= yhi]
  if (length(ticks) < 2L) ticks <- seq(ylo, yhi, length.out = 4)
  n_clip <- sum(allv < ylo | allv > yhi)
  if (isTRUE(opt$y_zoom)) {
    for (k in seq_along(v_plot)) {
      v_plot[[k]] <- v_plot[[k]][v_plot[[k]] >= ylo & v_plot[[k]] <= yhi]
    }
  }

  y_base     <- ylo - 0.015 * span2
  y_group    <- ylo - 0.080 * span2
  y_n        <- ylo - 0.140 * span2
  y_panel    <- ylo - 0.235 * span2
  y_plot_low <- ylo - (if (nest) 0.350 else 0.300) * span2
  y_bracket  <- yhi + 0.055 * span2
  y_plot_hi  <- yhi + 0.135 * span2
  opt$ylo <- ylo; opt$yhi <- yhi

  # Keep boxes clustered and centred instead of stretching them over 0-100
  # when there are only 2-3 groups. In --nest mode an extra wide gap marks
  # the boundary between two OUTER groups (the "big group" brackets).
  if (nest) {
    gap <- min(16, 82 / (P * G + 0.9 * (P - 1L)))
    panel_gap <- if (P > 1L) gap * 1.15 else 0
  } else {
    gap <- min(19, 74 / (P * G))
    panel_gap <- if (P > 1L) gap * 0.75 else 0
  }
  total_span <- (P * G - 1) * gap + (P - 1) * panel_gap
  x_start <- max(16, (100 - total_span) / 2)
  cx <- numeric(P * G)
  for (p in seq_len(P)) for (g in seq_len(G)) {
    k <- (p - 1) * G + g
    cx[k] <- x_start + (k - 1) * gap + (p - 1) * panel_gap
  }
  x_axis <- cx[1] - gap * 0.95

  # Derive box width from the PIXEL aspect ratio, not from arbitrary x units.
  # This is what stops IQR=0.5 data from looking like a row of flying saucers.
  Wpx <- if (is.null(opt$wpx)) 1600 else opt$wpx
  Hpx <- if (is.null(opt$hpx)) 900 else opt$hpx
  px_per_x <- Wpx / 106              # xlim is c(-6, 100): room for axis title
  px_per_y <- (Hpx * 0.80) / (y_plot_hi - y_plot_low)
  # Thin tick labels that would vertically overlap on short/high-res canvases
  if (length(ticks) > 2L) {
    dpx <- min(diff(ticks)) * px_per_y
    need <- 20 * opt$text_cex
    if (dpx < need) {
      step <- max(2L, as.integer(ceiling(need / max(dpx, 1e-6))))
      ticks <- ticks[seq(1L, length(ticks), by = step)]
    }
  }
  usable <- v_all[vapply(v_all, length, 0L) >= 2L]
  iqr_med <- stats::median(vapply(usable, IQR, 0), na.rm = TRUE)
  w <- if (opt$type == "box") {
    1.20 * iqr_med * px_per_y / px_per_x
  } else {
    0.30 * span2 * px_per_y / px_per_x
  }
  w <- min(gap * 0.70, max(gap * 0.26, w))
  dx <- w * 0.16 * opt$perspective
  dy <- span2 * 0.022 * opt$perspective
  pal <- rep_len(opt$palette, G)
  fill <- rep(pal, times = P)
  edge <- vapply(fill, .darken, "", f = 0.72, USE.NAMES = FALSE)

  # ---- statistics -----------------------------------------------------
  stats <- data.frame(panel = character(), comparison = character(),
                      p_value = numeric(), stars = character(),
                      stringsAsFactors = FALSE)
  brackets <- list()
  global_p <- rep(NA_real_, P)   # one omnibus p per outer group in --nest mode
  if (opt$show_p) {
    for (p in seq_len(P)) {
      idx <- (p - 1) * G + seq_len(G)
      if (G == 2L) {
        a <- v_all[[idx[1]]]; b <- v_all[[idx[2]]]
        if (length(a) >= 2L && length(b) >= 2L) {
          tt <- .test_two(a, b, opt$test)
          stats <- rbind(stats, data.frame(
            panel = panels[p], comparison = paste(groups, collapse = " vs "),
            p_value = tt$p, stars = tt$stars, stringsAsFactors = FALSE))
          brackets[[length(brackets) + 1L]] <- list(
            x1 = cx[idx[1]], x2 = cx[idx[2]],
            y1 = max(v_all[[idx[1]]]) + 0.05 * span2,
            y2 = max(v_all[[idx[2]]]) + 0.05 * span2,
            top = y_bracket, stars = tt$stars)
        }
      } else if (G > 2L) {
        gl <- unlist(lapply(seq_len(G), function(g) {
          rep(paste0("g", g), length(v_all[[idx[g]]])) }))
        gg <- data.frame(expr = unlist(v_all[idx]), g = factor(gl))
        global_p[p] <- tryCatch(.test_many(gg, opt$test),
                                error = function(e) NA_real_)
        stats <- rbind(stats, data.frame(
          panel = panels[p], comparison = "omnibus",
          p_value = global_p[p], stars = .p_stars(global_p[p]),
          stringsAsFactors = FALSE))
      }
    }
  }

  # ---- draw ------------------------------------------------------------
  do_draw <- function() {
    # Restore ONLY what we change. A full par(no.readonly=TRUE) restore would
    # also put back fig/plt/mfg/mfrow and lock every --facet_genes panel into
    # the first cell of the grid.
    old_mar <- par("mar"); old_bg <- par("bg")
    on.exit(par(mar = old_mar, bg = old_bg), add = TRUE)
    par(mar = c(0.4, 0.4, 3.4, 0.4), bg = "white")
    plot.new()
    plot.window(c(-6, 100), c(y_plot_low, y_plot_hi))
    u <- par("usr")
    rect(u[1] - 5, u[3], u[2] + 5, u[4], col = C_THEME$bg, border = NA)

    # quiet horizontal grid + a single horizontal ground line
    for (tv in ticks) {
      segments(x_axis, tv, 98, tv, col = C_THEME$grid, lwd = 1)
    }
    segments(x_axis, y_base, 98, y_base,
             col = adjustcolor(C_THEME$axis, alpha.f = 0.55), lwd = 1)
    for (k in seq_along(cx)) {
      .soft_ellipse(cx[k] + dx / 2, y_base - 0.012 * span2,
                    w * 0.80, span2 * 0.016, alpha = 0.05)
    }

    ttl <- if (is.null(opt$title)) paste0(gene, " expression") else opt$title
    mtext(ttl, side = 3, line = 1.05, cex = opt$title_cex,
          font = 2, col = C_THEME$ink)
    if (!isTRUE(opt$compact)) {
      sub <- paste0("by ", paste(group_by, collapse = if (nest) " > " else " x "),
                    "   |   n = ", sum(n_by), " ")
      if (n_clip > 0L) {
        sub <- paste0(sub, "   |   ", n_clip, " low values clipped")
      }
      if (opt$show_p && G > 2L && !nest && !is.na(global_p[1])) {
        sub <- paste0(sub, "   |   omnibus p = ",
                      formatC(global_p[1], format = "g", digits = 3),
                      " (", .p_stars(global_p[1]), ")")
      }
      mtext(sub, side = 3, line = -0.05, cex = 0.76 * opt$text_cex,
            col = C_THEME$muted)
    }
    .y_axis(ticks, ylo, yhi, opt$ylab, opt$text_cex, x = x_axis)

    # horizontal labels: far easier to read than text rotated along a 3D axis
    glab <- .shorten(groups, opt$label_max)
    for (k in seq_along(cx)) {
      p <- (k - 1L) %/% G + 1L; g <- (k - 1L) %% G + 1L
      text(cx[k] + dx / 2, y_group, glab[g],
           cex = 0.92 * opt$text_cex, col = C_THEME$ink, adj = c(0.5, 1))
      text(cx[k] + dx / 2, y_n, paste0("n = ", n_by[p, g]),
           cex = 0.70 * opt$text_cex, col = C_THEME$muted, adj = c(0.5, 1))
    }
    if (nest) {
      # Outer ("big") group: a full-width bracket + label below the inner
      # labels, exactly like T98G / LN229 under NC / BSI1 in the reference.
      for (p in seq_len(P)) {
        idx <- (p - 1L) * G + seq_len(G)
        xm <- mean(cx[idx])
        segments(min(cx[idx]) - w * 0.65, ylo - 0.190 * span2,
                 max(cx[idx]) + w * 0.65, ylo - 0.190 * span2,
                 col = adjustcolor(C_THEME$axis, alpha.f = 0.50), lwd = 1)
        text(xm, y_panel, .shorten(panels[p], 32),
             cex = 0.94 * opt$text_cex, font = 2,
             col = C_THEME$ink, adj = c(0.5, 1))
        if (opt$show_p && G > 2L && !is.na(global_p[p])) {
          text(xm, y_panel - 0.058 * span2,
               sprintf("p = %s (%s)",
                       formatC(global_p[p], format = "g", digits = 2),
                       .p_stars(global_p[p])),
               cex = 0.68 * opt$text_cex, col = C_THEME$muted,
               adj = c(0.5, 1))
        }
      }
    } else if (!(P == 1L && panels[1] == "")) {
      for (p in seq_len(P)) {
        idx <- (p - 1L) * G + seq_len(G)
        text(mean(cx[idx]), y_panel, .shorten(panels[p], 32),
             cex = 0.80 * opt$text_cex, font = 2,
             col = C_THEME$muted, adj = c(0.5, 1))
      }
    }

    for (k in seq_along(cx)) {
      if (length(v_all[[k]]) < 2L) next
      st <- .box_stat(v_all[[k]], opt$whisker)
      args <- list(cx = cx[k], fill = fill[k], edge = edge[k], w = w,
                   dx = dx, dy = dy, opt = opt, seed = opt$seed + 500 + k)
      if (opt$type == "box") {
        do.call(.draw_cube, c(args, list(st = st, v_plot = v_plot[[k]])))
      } else {
        do.call(.draw_violin,
                c(args, list(v_all = v_all[[k]], v_plot = v_plot[[k]])))
      }
    }

    # compact brackets: a short cap close to the data, not a tall doorway
    if (length(brackets)) for (b in brackets) {
      b1 <- b$x1 + w * 0.30
      b2 <- b$x2 - w * 0.30
      top <- y_bracket
      y0 <- yhi + 0.004 * span2
      segments(b1, y0, b1, top, lwd = 1, col = C_THEME$axis)
      segments(b1, top, b2, top, lwd = 1, col = C_THEME$axis)
      segments(b2, top, b2, y0, lwd = 1, col = C_THEME$axis)
      text((b1 + b2) / 2, top + 0.018 * span2, b$stars,
           cex = 1.05 * opt$text_cex, font = 2, col = C_THEME$ink)
    }
  }

  list(do_draw = do_draw, stats = stats, n = sum(n_by),
       y_range = range(allv))
}


# ---------------------------------------------------------------------
# 6. Public API
# ---------------------------------------------------------------------
plot_gene_3d <- function(rds, gene = NULL, group_by = NULL,
                         type = c("violin", "box"),
                         out = "gene_viz_3d.png",
                         pdf = TRUE,
                         assay = NULL,
                         layer = "data",
                         transform = c("none", "log1p", "log2p1"),
                         group_levels = NULL,
                         palette = NULL,
                         dot_size = 0.65,
                         dot_alpha = 0.55,
                         fill_alpha = 0.55,
                         show_points = TRUE,
                         max_cells = 1500,
                         jitter = 0.75,
                         perspective = 1,
                         bw_adjust = 1,
                         whisker = c("tukey", "range"),
                         y_zoom = TRUE,
                         show_p = TRUE,
                         test = c("wilcox", "t"),
                         title = NULL,
                         ylab = NULL,
                         label_max = 16,
                         text_cex = 1,
                         title_cex = 1.6,
                         width = NULL,
                         height = 900,
                         res = 160,
                         max_auto_genes = 20,
                         nest = FALSE,
                         facet_genes = FALSE,
                         facet_ncol = NULL,
                         outer_group = NULL,
                         inner_group = NULL,
                         page_group = NULL,
                         seed = 20261112,
                         verbose = TRUE) {
  prep <- .prepare_input(rds, gene, group_by, max_auto_genes,
                         outer_group, inner_group, page_group)
  rds <- prep$src; gene <- prep$gene; group_by <- prep$group_by
  # Naming a big/small group column IS the nested request; --nest is implied.
  if (isTRUE(prep$forced_nest)) nest <- TRUE
  type <- match.arg(type)
  transform <- match.arg(transform)
  test <- match.arg(test)
  whisker <- match.arg(whisker)
  if (length(group_by) < 1L || length(group_by) > 3L) {
    stop("group_by must have length 1, 2 or 3", call. = FALSE)
  }
  if (nest && length(group_by) < 2L) {
    stop("--nest needs 2 grouping columns: outer (big) and inner (small).",
         call. = FALSE)
  }
  if (length(group_by) == 3L && !nest) {
    stop("3 grouping columns require nest=TRUE (page > outer > inner).",
         call. = FALSE)
  }
  # 3 columns = one PAGE per level of column 1, each page is a nested
  # outer > inner plot. Keeps the drawing engine at exactly 2 columns.
  page_col  <- if (length(group_by) == 3L) group_by[1] else NULL
  draw_cols <- if (length(group_by) == 3L) group_by[2:3] else group_by
  if (verbose) {
    cat(sprintf("Input: %s | genes: %d (%s) | groups: %s%s\n",
                prep$kind, length(gene),
                paste(utils::head(gene, 4L), collapse = ", "),
                paste(group_by, collapse = if (nest) " > " else " x "),
                if (nest) " [nested]" else ""))
    if (nest) {
      cat(sprintf("  big (outer): %s | small (inner): %s%s\n",
                  draw_cols[1], draw_cols[2],
                  if (!is.null(page_col)) paste0(" | page: ", page_col) else ""))
    }
  }
  if (is.null(palette)) palette <- DEFAULT_PAL
  ylab_user <- ylab
  if (is.null(group_levels)) group_levels <- as.list(rep(NULL, length(group_by)))
  if (!is.list(group_levels)) group_levels <- as.list(c(group_levels,
      rep(list(NULL), max(0L, length(group_by) - 1L))))
  if (length(group_levels) < length(group_by)) {
    group_levels <- c(group_levels,
                      rep(list(NULL), length(group_by) - length(group_levels)))
  }

  opt <- list(type = type, palette = palette, dot_size = dot_size,
              dot_alpha = dot_alpha, fill_alpha = fill_alpha,
              show_points = show_points, max_cells = max_cells,
              jitter = jitter, perspective = perspective,
              bw_adjust = bw_adjust, whisker = whisker, y_zoom = y_zoom,
              show_p = show_p, test = test, nest = nest, compact = FALSE,
              title = title, ylab = "expression", label_max = label_max,
              text_cex = text_cex, title_cex = title_cex, seed = seed)

  safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  results <- list()

  # ---- --facet_genes: several genes as a grid of panels in ONE image ---
  # Reproduces the multi-panel layout of the reference figure. Each panel
  # keeps its own Y scale, groups and statistics.
  if (isTRUE(facet_genes) && length(gene) > 1L) {
    if (!is.null(page_col)) {
      stop("--facet_genes cannot be combined with 3 grouping columns.",
           call. = FALSE)
    }
    nc <- if (is.null(facet_ncol)) min(3L, ceiling(sqrt(length(gene)))) else
      max(1L, as.integer(facet_ncol))
    nr <- ceiling(length(gene) / nc)
    wpx <- if (is.null(width)) 780L * nc else width
    hpx <- height * nr
    opt$compact   <- TRUE
    opt$text_cex  <- text_cex * 0.88
    opt$title_cex <- title_cex * 0.88
    opt$wpx <- wpx / nc; opt$hpx <- hpx / nr

    all_stats <- data.frame()
    draw_grid <- function() {
      par(mfrow = c(nr, nc))
      for (g in gene) {
        opt$title <- if (is.null(title)) g else title
        opt$ylab  <- if (is.null(ylab_user)) g else ylab_user
        fr <- .build_frame(rds, g, group_by, assay, layer, transform,
                           group_levels)
        one <- .plot_one(fr$d, g, draw_cols, fr$levels, opt)
        one$do_draw()
        if (nrow(one$stats)) {
          s <- one$stats; s$gene <- g
          all_stats <<- rbind(all_stats,
                              s[, c("gene", setdiff(names(s), "gene")),
                                drop = FALSE])
        }
        results[[g]] <<- list(stats = one$stats, n_cells = one$n)
      }
    }
    f_png <- out; f_pdf <- sub("\\.[^.]+$", ".pdf", f_png)
    png(f_png, width = wpx, height = hpx, res = res, bg = "white")
    draw_grid(); dev.off()
    if (pdf) {
      pdf(f_pdf, width = wpx / res, height = hpx / res)
      draw_grid(); dev.off()
    }
    if (nrow(all_stats)) {
      utils::write.csv(all_stats, sub("\\.[^.]+$", "_stats.csv", f_png),
                       row.names = FALSE)
    }
    for (g in gene) results[[g]]$png <- f_png
    if (verbose) cat(sprintf("[facet] %d panels | %s\n", length(gene), f_png))
    return(invisible(results))
  }

  # ---- classic / nested: one image per gene (and per page) ------------
  for (g in gene) {
    fr <- .build_frame(rds, g, group_by, assay, layer, transform, group_levels)
    pages <- if (!is.null(page_col)) fr$levels[[page_col]] else ""
    for (pg in pages) {
      dpg <- if (!is.null(page_col)) {
        fr$d[fr$d[[page_col]] == pg, , drop = FALSE]
      } else fr$d
      if (nrow(dpg) < 3L) {
        warning("Skipping page '", pg, "' of ", g, ": too few rows.")
        next
      }
      n_outer <- if (length(draw_cols) == 2L) {
        length(fr$levels[[draw_cols[1]]])
      } else 1L
      n_inner <- length(fr$levels[[draw_cols[length(draw_cols)]]])
      wpx <- if (is.null(width)) max(1200, 620 + 230 * n_outer * n_inner) else width
      opt$wpx <- wpx; opt$hpx <- height
      opt$ylab <- if (is.null(ylab_user)) paste0(g, " expression") else ylab_user
      # Assign through [ so a NULL title stays an (empty) list element:
      # `opt$title <- NULL` would DELETE it and `opt$title` would then
      # partially match `title_cex` (printing "1.6" as the plot title).
      opt["title"] <- list(if (!is.null(title)) title else if (!is.null(page_col)) {
        paste0(g, " (", pg, ")")
      } else NULL)
      one <- .plot_one(dpg, g, draw_cols, fr$levels, opt)

      # one gene -> use `out` verbatim; prefix with gene / page otherwise
      fname <- basename(out)
      if (length(gene) > 1L) fname <- paste0(safe(g), "_", fname)
      if (!is.null(page_col)) fname <- paste0(safe(pg), "_", fname)
      f_png <- file.path(dirname(out), fname)
      f_pdf <- sub("\\.[^.]+$", ".pdf", f_png)

      png(f_png, width = wpx, height = height, res = res, bg = "white")
      one$do_draw(); dev.off()
      if (pdf) {
        pdf(f_pdf, width = wpx / res, height = height / res)
        one$do_draw(); dev.off()
      }
      f_csv <- sub("\\.[^.]+$", "_stats.csv", f_png)
      if (nrow(one$stats)) utils::write.csv(one$stats, f_csv, row.names = FALSE)

      key <- if (!is.null(page_col)) paste0(pg, "/", g) else g
      results[[key]] <- list(png = f_png, pdf = if (pdf) f_pdf else NULL,
                             stats = one$stats, n_cells = one$n)
      if (verbose) {
        cat(sprintf("[%s] n = %d | %s | %s\n", key, one$n, type, f_png))
        if (nrow(one$stats)) print(one$stats, row.names = FALSE)
      }
    }
  }
  invisible(results)
}


# ---------------------------------------------------------------------
# 7. Command line interface
# ---------------------------------------------------------------------
.cli_help <- function() {
  cat("
gene_viz_3d.R - pseudo-3D gene expression violin / boxplot

Required:
  --rds PATH        input file; --file is an alias. RDS and delimited text
                    are detected automatically (extension + binary sniff)

  --gene NAME       gene column(s), comma-separated.
                    RDS: REQUIRED.  Text: OPTIONAL, all numeric gene
                    columns are plotted when omitted.
  --group COLS      grouping column(s), 1-3, comma-separated.
                    RDS: REQUIRED.  Text: OPTIONAL, defaults to Type.
                    1 col : simple groups
                    2 cols: default = panel x axis (faceted)
                            --nest 1 = OUTER big group > INNER small group
                    3 cols: requires --nest 1; col 1 = page (one image per
                            level), col 2 = outer, col 3 = inner

Explicit big/small group COLUMN NAMES (recommended, works for BOTH text
and RDS; --nest is implied, do NOT also pass --group):
  --outer COL       BIG / outer group column, e.g. CellLine, group, Sex
                    alias: --big_group
  --inner COL       SMALL / inner group column, e.g. Treatment, celltype, Age
                    alias: --small_group
  --page COL        optional 3rd level: one image per level of this column
  Column names must match exactly; a wrong name lists all available columns.

Text file layout (CSV / TSV / TXT / .gz, separator auto-detected):
  ID,YKT6,CD4,MS4A1,Type
  S001,7.12,3.20,1.10,Normal
  S002,10.4,3.55,1.02,Tumor
  ID and Type are case-insensitive. Every other NUMERIC column is a gene;
  non-numeric columns remain available as extra grouping columns.
  --max_auto_genes NUM  safety cap when --gene is omitted (default: 20)

Examples:
  Rscript gene_viz_3d.R --rds obj.rds --gene YKT6 --group celltype
  Rscript gene_viz_3d.R --file expr.csv                    # all genes, by Type
  Rscript gene_viz_3d.R --file expr.tsv --gene YKT6,CD4 --type box

Nested big/small groups - name the columns explicitly (text or RDS):
  Rscript gene_viz_3d.R --file data.csv --gene YKT6 \
      --outer CellLine --inner Treatment --type box
  Rscript gene_viz_3d.R --rds obj.rds --gene YKT6 \
      --outer group --inner celltype
  # optional 3rd level (one image per cohort):
  Rscript gene_viz_3d.R --file data.csv --gene YKT6 \
      --page Cohort --outer Sex --inner Age

Equivalent positional form (order = page > outer > inner):
  Rscript gene_viz_3d.R --file data.csv --gene YKT6 \
      --group Sex,Age --nest 1 --type box
  # several genes as a 2x3 panel grid in one image (reference-figure style):
  Rscript gene_viz_3d.R --file data.csv --gene A,B,C,D,E,F \
      --outer CellLine --inner Treatment --facet_genes 1

Data options:
  --assay NAME      Seurat assay (default: RNA / first assay)
  --layer NAME      data | counts | scale.data (default: data)
  --transform NAME  none | log1p | log2p1 (default: none)
  --palette HEXES   comma-separated custom colors

Plot options:
  --type NAME       violin | box (default: violin)
  --out PATH        output PNG (default: gene_viz_3d.png)
  --pdf 0|1         also write PDF (default: 1)
  --title TEXT      --ylab TEXT
  --width PX        --height PX    --res DPI

Aesthetic options:
  --dot_size NUM    point cex (default: 0.65)
  --dot_alpha NUM   point opacity 0-1 (default: 0.55)
  --fill_alpha NUM  violin/box opacity 0-1 (default: 0.55)
  --show_points 0|1 draw jittered dots (default: 1)
  --max_cells NUM   max dots drawn per group (default: 1500)
  --jitter NUM      horizontal jitter 0-1 (default: 0.75)
  --perspective NUM 0 = flat, 1 = full 3D (default: 1)
  --bw_adjust NUM   violin bandwidth multiplier (default: 1)
  --whisker NAME    tukey (1.5*IQR, default) | range (min/max)
  --y_zoom 0|1      trim Y axis to 1st-99th pct (default: 1); prevents
                    dropout zeros from flattening every box/violin
  --label_max NUM   truncate long group labels (default: 16)
  --text_cex NUM    global text scale (default: 1)

Grouping / layout:
  --nest 0|1        nested big > small groups on ONE x axis (default: 0).
                    Implied by --outer/--inner. Colors follow the INNER
                    group, brackets compare inner groups WITHIN each outer.
  --facet_genes 0|1 draw all --gene values as a grid of panels in one image
                    (default: 0); each panel keeps its own Y scale
  --facet_ncol NUM  grid columns (default: auto, max 3)

Statistics:
  --show_p 0|1      significance brackets / omnibus p (default: 1)
  --test NAME       wilcox | t (default: wilcox)

Other:
  --seed NUM        --max_auto_genes NUM   --help
")
}

.parse_cli <- function(args) {
  o <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) { i <- i + 1L; next }
    kv <- sub("^--", "", a)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv); val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      if (i < length(args) && !startsWith(args[i + 1L], "--")) {
        val <- args[i + 1L]; i <- i + 1L
      } else val <- "1"
    }
    key <- gsub("-", "_", key)
    o[[key]] <- val
    i <- i + 1L
  }
  o
}

.num <- function(x, default) {
  if (is.null(x)) return(default)
  z <- suppressWarnings(as.numeric(x))
  if (is.na(z)) default else z
}
.flag <- function(x, default) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("1", "true", "yes", "y", "t")
}
.chrvec <- function(x) {
  if (is.null(x)) return(NULL)
  trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}

.cli_main <- function(args) {
  o <- .parse_cli(args)
  if (!is.null(o$help)) { .cli_help(); return(invisible(NULL)) }
  path <- if (!is.null(o$rds)) o$rds else o$file
  if (is.null(path)) {
    stop("Missing --rds / --file. Run with --help for usage.", call. = FALSE)
  }
  plot_gene_3d(
    rds = path, gene = .chrvec(o$gene), group_by = .chrvec(o$group),
    max_auto_genes = .num(o$max_auto_genes, 20),
    nest = .flag(o$nest, FALSE),
    outer_group = if (!is.null(o$outer)) o$outer else o$big_group,
    inner_group = if (!is.null(o$inner)) o$inner else o$small_group,
    page_group  = if (!is.null(o$page)) o$page else o$page_group,
    facet_genes = .flag(o$facet_genes, FALSE),
    facet_ncol = if (is.null(o$facet_ncol)) NULL else .num(o$facet_ncol, 3),
    type = if (is.null(o$type)) "violin" else o$type,
    out = if (is.null(o$out)) "gene_viz_3d.png" else o$out,
    pdf = .flag(o$pdf, TRUE),
    assay = o$assay,
    layer = if (is.null(o$layer)) "data" else o$layer,
    transform = if (is.null(o$transform)) "none" else o$transform,
    palette = .chrvec(o$palette),
    dot_size = .num(o$dot_size, 0.65),
    dot_alpha = .num(o$dot_alpha, 0.55),
    fill_alpha = .num(o$fill_alpha, 0.55),
    show_points = .flag(o$show_points, TRUE),
    max_cells = .num(o$max_cells, 1500),
    jitter = .num(o$jitter, 0.75),
    perspective = .num(o$perspective, 1),
    bw_adjust = .num(o$bw_adjust, 1),
    whisker = if (is.null(o$whisker)) "tukey" else o$whisker,
    y_zoom = .flag(o$y_zoom, TRUE),
    show_p = .flag(o$show_p, TRUE),
    test = if (is.null(o$test)) "wilcox" else o$test,
    title = o$title, ylab = o$ylab,
    label_max = .num(o$label_max, 16),
    text_cex = .num(o$text_cex, 1),
    width = if (is.null(o$width)) NULL else .num(o$width, 1200),
    height = .num(o$height, 900),
    res = .num(o$res, 160),
    seed = .num(o$seed, 20261112)
  )
}

if (sys.nframe() == 0L) {
  .cli_main(commandArgs(trailingOnly = TRUE))
}
