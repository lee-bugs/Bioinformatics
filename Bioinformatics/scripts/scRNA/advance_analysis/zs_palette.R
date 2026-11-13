# =====================================================================
# zs_palette.R  |  ZS palette kit for R
# ---------------------------------------------------------------------
# 45 HEX colors extracted from 9 palette cards (verified pixel-by-pixel).
# Pure ASCII source: works under ANY locale, incl. Windows "C" locale.
# Core functions use base R / grDevices only. No packages required.
#
#   source("zs_palette.R")
#
# GRADIENT (choose one of 9 palettes)
#   zs_grad(1)            # the 5 base colors of palette 1
#   zs_grad(5)            # the 5 base colors of palette 5
#   zs_grad(5, n = 100)   # 100-step continuous ramp of palette 5
#   zs_grad(5, 100, interpolate = "linear")  # strict, no spline overshoot
#   zs_ramp(5)(256)       # a function(n) ramp generator for ggplot2/image
#
# DISCRETE (highly distinguishable, up to ~40 colors)
#   zs_disc(1:30)         # 30 maximally separated colors
#   zs_disc(30)           # shorthand, same result
#   zs_disc(c(3, 1, 7))   # pick the 3rd, 1st, 7th color only
#   zs_disc(1:30, as_list = TRUE)  # return a real R list
#
# GGPLOT2 (only these wrappers need ggplot2)
#   + scale_color_zs_d() / scale_fill_zs_d()
#   + scale_color_zs_c(pal = 5) / scale_fill_zs_c(pal = 5)
#
# MISC
#   zs_show()             # base-R preview of all 9 ramps + 30 discrete
#   zs_pal_names()        # ids and names
# =====================================================================


# ---------------------------------------------------------------------
# 1. The 45 colors. Group order = image reading order (L->R, T->B);
#    within a group = top to bottom in the image.
# ---------------------------------------------------------------------
.zs_pals <- list(
  "1" = c("#556136", "#242F1A", "#FEF9DE", "#D6975B", "#B2622C"),
  "2" = c("#85C2E0", "#2793B1", "#092A3D", "#FDAF31", "#F77C26"),
  "3" = c("#9ACBFB", "#B6DBFB", "#FCA7C4", "#FDC2D7", "#C5ABD3"),
  "4" = c("#F27F0C", "#F7AD19", "#053F5C", "#428EBD", "#9FE7F5"),
  "5" = c("#1E0A26", "#F2B28D", "#F2DB66", "#F2CEE6", "#9CD9D3"),
  "6" = c("#61589E", "#EB6962", "#EDD7DF", "#FCF3E7", "#BF5C78"),
  "7" = c("#161616", "#463DE8", "#74B5B3", "#D76A7B", "#FFFFFF"),
  "8" = c("#273F75", "#797784", "#653136", "#A35C4B", "#A98D42"),
  "9" = c("#2F443A", "#344F3A", "#507550", "#99A782", "#D5D1C7")
)

# Chinese names stored as Unicode escapes so the file stays ASCII-safe.
.zs_meta <- data.frame(
  id      = 1:9,
  name_en = c("forest_clay", "sea_citrus", "cotton_candy",
              "tangerine_teal", "candy_spring", "berry_plum",
              "neon_pop", "vintage_dusk", "moss_wood"),
  name_zh = c("\u68EE\u91CE\u9676\u571F", "\u6D77\u76D0\u67D1\u6A58",
              "\u68C9\u82B1\u7CD6\u7C89\u84DD", "\u6A58\u84DD\u649E\u8272",
              "\u7CD6\u679C\u6625\u5149", "\u8393\u679C\u7D2B",
              "\u7535\u5149\u9713\u8679", "\u590D\u53E4\u66AE\u8272",
              "\u82D4\u85D3\u68EE\u6797"),
  stringsAsFactors = FALSE
)

zs_pal_names <- function() .zs_meta

zs_pal <- function(pal = 1) .zs_pals[[as.character(.zs_resolve(pal))]]

.zs_resolve <- function(pal) {
  if (length(pal) != 1L || is.na(pal)) {
    stop("`pal` must be a single id (1-9) or a palette name", call. = FALSE)
  }
  key <- tolower(trimws(as.character(pal)))
  if (key %in% as.character(1:9)) return(as.integer(key))
  hit <- which(tolower(.zs_meta$name_en) == key |
                 .zs_meta$name_zh == as.character(pal))
  if (length(hit) == 1L) return(hit)
  stop("Unknown palette '", pal, "'. See zs_pal_names()", call. = FALSE)
}


# ---------------------------------------------------------------------
# 2. Color science: sRGB -> CIELAB (D65), CIE76 distance, maximin order.
#    Implemented here to avoid colorspace / farver dependencies.
# ---------------------------------------------------------------------
.zs_lab <- function(hex) {
  rgb <- t(col2rgb(hex)) / 255
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  M <- matrix(c(0.4124564, 0.3575761, 0.1804375,
                0.2126729, 0.7151522, 0.0721750,
                0.0193339, 0.1191920, 0.9503041), 3, 3, byrow = TRUE)
  xyz <- lin %*% t(M)
  f <- function(t) ifelse(t > 216 / 24389, t^(1 / 3), (841 / 108) * t + 4 / 29)
  fx <- f(xyz[, 1] / 0.95047); fy <- f(xyz[, 2]); fz <- f(xyz[, 3] / 1.08883)
  cbind(L = 116 * fy - 16, a = 500 * (fx - fy), b = 200 * (fy - fz))
}

.zs_de <- function(lab) {
  d <- matrix(0, nrow(lab), nrow(lab))
  for (j in seq_len(ncol(lab))) d <- d + outer(lab[, j], lab[, j], "-")^2
  sqrt(d)
}

# pmin/pmax drop dims, so restore the n x 3 shape explicitly.
.zs_clamp_rgb <- function(v) {
  matrix(pmax(0, pmin(1, as.vector(v) / 255)), ncol = 3)
}

# Greedy maximin: each step picks the color whose nearest already-picked
# neighbor is farthest away. So the first k entries are the k most
# mutually distinguishable colors of the whole pool.
.zs_maximin <- function(lab) {
  D <- .zs_de(lab)
  chroma <- sqrt(lab[, 2]^2 + lab[, 3]^2)
  start <- which.max(chroma - 0.15 * abs(lab[, 1] - 55))
  n <- nrow(lab)
  ord <- integer(n); ord[1L] <- start
  remain <- setdiff(seq_len(n), start)
  mind <- D[start, ]
  if (n > 1L) for (k in 2:n) {
    pick <- remain[which.max(mind[remain])]
    ord[k] <- pick
    remain <- setdiff(remain, pick)
    mind <- pmin(mind, D[pick, ])
  }
  ord
}


# ---------------------------------------------------------------------
# 3. Discrete palette: 45 colors -> drop near-whites -> maximin in Lab
# ---------------------------------------------------------------------
# x: 1:30 (returns 30 colors), a single 30 (same), or indices c(3,1,7).
# drop_light: remove colors too close to white (invisible on white bg).
# light_cut : CIE76 distance-to-white threshold; smaller = keep more lights.
# as_list   : TRUE -> list, FALSE -> named character vector (default).
zs_disc <- function(x = 1:30, drop_light = TRUE, light_cut = 20, as_list = FALSE) {
  all_hex <- unlist(.zs_pals, use.names = FALSE)
  lab <- .zs_lab(all_hex)

  if (drop_light) {
    w <- .zs_lab("#FFFFFF")
    dw <- sqrt(rowSums((lab - matrix(w, nrow(lab), 3, byrow = TRUE))^2))
    keep <- which(dw >= light_cut)
  } else {
    keep <- seq_along(all_hex)
  }

  pool <- all_hex[keep][.zs_maximin(lab[keep, , drop = FALSE])]
  names(pool) <- sprintf("D%02d", seq_along(pool))

  # A single number means "how many"; a vector means "which indices".
  if (length(x) == 1L && is.numeric(x) && x == as.integer(x) && x > 0) {
    idx <- seq_len(as.integer(x))
  } else {
    idx <- as.integer(x)
    if (anyNA(idx) || any(idx < 1L)) stop("`x` must be positive integer indices", call. = FALSE)
  }
  if (any(idx > length(pool))) {
    stop("Only ", length(pool), " distinguishable colors available (light_cut=",
         light_cut, "). Use fewer, lower light_cut, or drop_light = FALSE (45).",
         call. = FALSE)
  }

  out <- pool[idx]
  if (as_list) as.list(out) else out
}

zs_disc_max <- function(drop_light = TRUE, light_cut = 20) {
  all_hex <- unlist(.zs_pals, use.names = FALSE)
  if (!drop_light) return(length(all_hex))
  lab <- .zs_lab(all_hex)
  w <- .zs_lab("#FFFFFF")
  dw <- sqrt(rowSums((lab - matrix(w, nrow(lab), 3, byrow = TRUE))^2))
  sum(dw >= light_cut)
}


# ---------------------------------------------------------------------
# 4. Gradient palettes
# ---------------------------------------------------------------------
.zs_order_cols <- function(pal, reverse, sort_by) {
  cols <- zs_pal(pal)
  sort_by <- match.arg(sort_by, c("none", "luminance", "hue"))
  if (sort_by != "none") {
    lab <- .zs_lab(cols)
    ord <- if (sort_by == "luminance") order(lab[, 1]) else order(atan2(lab[, 3], lab[, 2]))
    cols <- cols[ord]
  }
  if (reverse) cols <- rev(cols)
  cols
}

# pal     : 1-9 or a palette name.
# n       : NULL -> the 5 base colors; integer -> n interpolated steps.
# sort_by : "none" keeps the original card order (some cards jump from
#           cool to warm); "luminance"/"hue" sorts first for a smooth ramp.
# interpolate: "spline" = smoother, may slightly overshoot the 5 colors;
#              "linear" = strictly stays between adjacent base colors.
zs_grad <- function(pal = 1, n = NULL, reverse = FALSE,
                    sort_by = c("none", "luminance", "hue"), space = "Lab",
                    interpolate = c("spline", "linear")) {
  i <- .zs_resolve(pal)
  cols <- .zs_order_cols(i, reverse, sort_by)

  if (is.null(n)) {
    names(cols) <- paste0("G", seq_along(cols))
    return(cols)
  }
  if (length(n) != 1L || !is.numeric(n) || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer, or NULL for the 5 base colors", call. = FALSE)
  }
  interpolate <- match.arg(interpolate)
  ramp <- grDevices::colorRamp(cols, space = space, interpolate = interpolate)
  v <- .zs_clamp_rgb(ramp(seq(0, 1, length.out = as.integer(n))))
  out <- grDevices::rgb(v[, 1], v[, 2], v[, 3])
  names(out) <- sprintf("C%03d", seq_along(out))
  out
}

# Returns function(n) -> n colors, for ggplot2 / image() / any ramp consumer.
zs_ramp <- function(pal = 1, reverse = FALSE,
                    sort_by = c("none", "luminance", "hue"), space = "Lab",
                    interpolate = c("spline", "linear")) {
  i <- .zs_resolve(pal)
  cols <- .zs_order_cols(i, reverse, sort_by)
  ramp <- grDevices::colorRamp(cols, space = space,
                               interpolate = match.arg(interpolate))
  function(n) {
    if (length(n) != 1L || n < 1) stop("n must be a positive integer", call. = FALSE)
    v <- .zs_clamp_rgb(ramp(seq(0, 1, length.out = as.integer(n))))
    grDevices::rgb(v[, 1], v[, 2], v[, 3])
  }
}


# ---------------------------------------------------------------------
# 5. ggplot2 wrappers (the only part that needs ggplot2)
# ---------------------------------------------------------------------
.zs_need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 first: install.packages('ggplot2')", call. = FALSE)
  }
}

scale_color_zs_d <- function(n = 30, ...) {
  .zs_need_ggplot2()
  ggplot2::scale_color_manual(values = unname(zs_disc(n)), ...)
}
scale_fill_zs_d <- function(n = 30, ...) {
  .zs_need_ggplot2()
  ggplot2::scale_fill_manual(values = unname(zs_disc(n)), ...)
}
scale_color_zs_c <- function(pal = 1, n = 256, reverse = FALSE, sort_by = "none", ...) {
  .zs_need_ggplot2()
  ggplot2::scale_color_gradientn(colours = unname(zs_grad(pal, n, reverse, sort_by)), ...)
}
scale_fill_zs_c <- function(pal = 1, n = 256, reverse = FALSE, sort_by = "none", ...) {
  .zs_need_ggplot2()
  ggplot2::scale_fill_gradientn(colours = unname(zs_grad(pal, n, reverse, sort_by)), ...)
}
scale_colour_zs_d <- scale_color_zs_d
scale_colour_zs_c <- scale_color_zs_c


# ---------------------------------------------------------------------
# 6. Preview (base R graphics only)
# ---------------------------------------------------------------------
zs_show <- function(disc = 1:30, grad_n = 100) {
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  par(mar = c(0.3, 0.3, 1.6, 0.3))
  plot.new(); plot.window(c(0, 1), c(0, 1))

  text(0.01, 0.985, "ZS PALETTES  |  9 gradients + high-contrast discrete",
       adj = c(0, 0.5), cex = 1.1, font = 2)
  text(0.01, 0.955, "Left: 5 base colors        Right: continuous ramp",
       adj = c(0, 0.5), cex = 0.75, col = "grey40")

  top <- 0.925; h <- 0.076
  for (i in 1:9) {
    y1 <- top - i * h; y0 <- y1 + h - 0.008
    text(0.012, (y0 + y1) / 2, sprintf("%d  %s", i, .zs_meta$name_en[i]),
         adj = c(0, 0.5), cex = 0.72)
    raw <- zs_grad(i)
    for (j in 1:5) {
      x0 <- 0.18 + (j - 1) * 0.048
      rect(x0, y1, x0 + 0.046, y0, col = raw[j], border = "grey85")
    }
    g <- zs_grad(i, grad_n)
    for (j in seq_along(g)) {
      x0 <- 0.44 + (j - 1) * 0.54 / length(g)
      rect(x0, y1, x0 + 0.54 / length(g) + 0.0008, y0, col = g[j], border = NA)
    }
  }

  y <- 0.185
  d <- zs_disc(disc)
  text(0.012, y + 0.022, sprintf("%d discrete colors (CIELAB maximin order)", length(d)),
       adj = c(0, 0.5), cex = 0.8, font = 2)
  nc <- ceiling(length(d) / 2); w <- 0.97 / nc
  for (k in seq_along(d)) {
    row <- (k - 1) %/% nc; col <- (k - 1) %% nc
    x0 <- 0.015 + col * w
    y1 <- y - 0.06 - row * 0.075; y0 <- y1 + 0.062
    rect(x0, y1, x0 + w - 0.002, y0, col = d[k], border = "grey85")
    lum <- sum(col2rgb(d[k]) * c(0.299, 0.587, 0.114)) / 255
    text(x0 + (w - 0.002) / 2, y1 + 0.031, k, cex = 0.6,
         col = if (lum > 0.6) "grey20" else "white")
  }
  invisible(NULL)
}
