# -------------------------------------------------------------------------------
# 3dumap.R - 通用 3D UMAP 绘图脚本（命令行传参版）
#
#   Rscript 3dumap.R <rds> [基因1 基因2 ...] [选项]
#
# 只需要两个输入：RDS 路径（必需）+ 关注基因（0~N 个，可选）。
# 分组列、assay、表达层、降维、输出目录全部自动探测，可用 -- 选项覆盖。
# 用法细节见文件末尾「执行」段，或运行 Rscript 3dumap.R --help
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# 函数
# -------------------------------------------------------------------------------
umap3d <- local({
  `%||%` <- function(x, y) if (is.null(x)) y else x
  umap3d_threads <- function() {
    n <- suppressWarnings(tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_))
    if (length(n) != 1L || is.na(n)) n <- 1L
    max(1L, min(8L, as.integer(n)))
  }
  need <- function(packages) {
    for (package in unique(packages)) {
      # 缺包安装，已有包跳过。
      if (!requireNamespace(package, quietly = TRUE)) {
        message('缺少R包，正在自动安装：', package)
        repos <- getOption('repos')
        if (is.null(repos) || !length(repos) || anyNA(repos) || any(repos == '@CRAN@'))
          repos <- c(CRAN = 'https://cloud.r-project.org')
        tryCatch(utils::install.packages(package, repos = repos), error = function(e)
          stop('自动安装 ', package, ' 失败：', conditionMessage(e),
            '。请检查网络和R包库写入权限后重试。', call. = FALSE))
        if (!requireNamespace(package, quietly = TRUE))
          stop('自动安装后仍无法加载R包：', package,
            '。请查看上方安装日志，修复安装问题后重新运行。', call. = FALSE)
        message('R包安装完成：', package)
      }
    }
    invisible(TRUE)
  }
  mkdir <- function(path) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop('无法创建目录：', path)
    normalizePath(path, winslash = '/', mustWork = TRUE)
  }
  check_output <- function(path, overwrite = FALSE) {
    if (file.exists(path) && !isTRUE(overwrite)) stop('文件已存在；请换输出目录或显式设置 overwrite=TRUE：', path)
    mkdir(dirname(path)); invisible(path)
  }
  hash <- function(value) {
    tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
    saveRDS(value, tmp, version = 2, compress = FALSE)
    unname(tools::md5sum(tmp))
  }
  with_seed <- function(seed, expr) {
    had <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
    if (had) old <- get('.Random.seed', envir = .GlobalEnv)
    on.exit(if (had) assign('.Random.seed', old, envir = .GlobalEnv) else
      if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv))
    set.seed(seed); force(expr)
  }
  as_data <- function(obj) {
    if (!inherits(obj, 'Seurat')) stop('输入仅支持Seurat对象。', call. = FALSE)
    need('SeuratObject'); cells <- colnames(obj); meta <- obj[[]]
    if (length(cells) < 10L || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells))
      stop('细胞至少需要10个，并且细胞ID必须非空且唯一。')
    list(object = obj, metadata = meta[cells, , drop = FALSE], cells = cells)
  }
  resolve_input <- function(input) {
    if (is.null(input)) stop('umap3d_input需要填写Seurat对象或Seurat RDS路径。', call. = FALSE)
    if (is.character(input)) {
      if (length(input) != 1L || is.na(input) || !nzchar(trimws(input)))
        stop('RDS路径必须是非空字符串；内存对象不加引号。', call. = FALSE)
      if (!file.exists(input) || dir.exists(input)) stop('RDS文件不存在：', input,
        '；内存对象不加引号。', call. = FALSE)
      message('读取输入RDS：', normalizePath(input, winslash = '/', mustWork = TRUE))
      input <- tryCatch(readRDS(input), error = function(e)
        stop('无法读取输入RDS：', conditionMessage(e), call. = FALSE))
    }
    if (!inherits(input, 'Seurat')) stop('输入仅支持Seurat对象；RDS必须保存Seurat对象。', call. = FALSE)
    input
  }
  inspect <- function(input) {
    d <- as_data(resolve_input(input)); reduction_ids <- reduction_names(d)
    dimensions <- vapply(reduction_ids, function(name) ncol(reduction_matrix(d, name)), integer(1))
    info <- list(cells = length(d$cells), metadata_columns = names(d$metadata),
      reductions = data.frame(name = reduction_ids, dimensions = unname(dimensions)),
      assays = SeuratObject::Assays(d$object))
    print(info); invisible(info)
  }
  reduction_names <- function(d) {
    SeuratObject::Reductions(d$object)
  }
  reduction_matrix <- function(d, name) {
    if (!name %in% reduction_names(d)) stop('缺少输入降维：', name,
      '。可用：', paste(reduction_names(d), collapse = ', '))
    x <- SeuratObject::Embeddings(d$object, name)
    if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) || anyDuplicated(rownames(x)) ||
        !setequal(rownames(x), d$cells)) stop('降维矩阵的细胞ID与 metadata 不匹配：', name)
    x <- x[d$cells, , drop = FALSE]
    if (any(!is.finite(x))) stop('降维矩阵包含非有限值：', name)
    x
  }
  input_matrix <- function(d, source_reduction = NULL, dims = NULL) {
    source_reduction <- source_reduction %||% if ('harmony' %in% reduction_names(d)) 'harmony' else 'pca'
    x <- reduction_matrix(d, source_reduction)
    dims <- dims %||% seq_len(min(30L, ncol(x)))
    if (!is.numeric(dims) || !length(dims) || anyNA(dims) || any(!is.finite(dims)) ||
        any(dims != floor(dims)) || any(dims < 1L | dims > ncol(x))) stop('dims 超出输入降维维度。')
    list(x = x[, unique(as.integer(dims)), drop = FALSE], name = source_reduction, dims = unique(as.integer(dims)))
  }
  labels_for <- function(d, column) {
    if (!is.character(column) || length(column) != 1L || anyNA(column) || !column %in% names(d$metadata))
      stop('metadata 标签列不存在：', paste(column, collapse = ','),
        '。可用：', paste(names(d$metadata), collapse = ', '))
    labels <- as.character(d$metadata[[column]])
    labels[is.na(labels) | !nzchar(trimws(labels))] <- 'Unknown'
    labels
  }
  palette <- function(n) {
    if (n <= 0L) return(character())
    grDevices::hcl(h = (seq_len(n) - 1) * 360 / n + 15, c = 70, l = 68)
  }
  read_xyz <- function(path, cells) {
    tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
    choices <- list(c('axis_1','axis_2','axis_3'), c('UMAP3D_1','UMAP3D_2','UMAP3D_3'),
      c('UMAP_1','UMAP_2','UMAP_3'))
    columns <- Filter(function(z) all(z %in% names(tab)), choices)
    if (!'cell' %in% names(tab) || !length(columns) || anyNA(tab$cell) || anyDuplicated(tab$cell) ||
        !setequal(cells, tab$cell)) stop('3D坐标CSV必须有 cell 和三轴列，且完整匹配细胞ID。')
    xyz <- as.matrix(tab[match(cells, tab$cell), columns[[1L]], drop = FALSE])
    if (!is.numeric(xyz) || any(!is.finite(xyz))) stop('3D坐标包含无效数值。')
    rownames(xyz) <- cells; xyz
  }
  normalize_methods <- function(methods) {
    if (!is.character(methods) || !length(methods) || anyNA(methods) ||
        any(!methods %in% c('umap', 'tsne')))
      stop('methods填写 "umap"、"tsne" 或 c("umap", "tsne")。', call. = FALSE)
    unique(methods)
  }
  compute_one <- function(obj, output_dir, method, source_reduction = NULL, dims = NULL,
      n_neighbors = 30L, min_dist = 0.3, metric = 'cosine', seed = 123L,
      coordinate_file = NULL, reuse_embedding = TRUE, overwrite = FALSE,
      n_epochs = NULL, verbose = TRUE, perplexity = 30, max_iter = 1000L, theta = 0.5) {
    need('jsonlite'); d <- as_data(obj)
    title <- if (method == 'umap') 'UMAP' else 't-SNE'
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max)
      stop('seed 必须是有效的非负整数。')
    xyz <- NULL; source <- NULL; input <- NULL
    if (!is.null(coordinate_file)) {
      xyz <- read_xyz(coordinate_file, d$cells)
      source <- paste0('csv:', normalizePath(coordinate_file, winslash = '/', mustWork = TRUE))
    } else if (isTRUE(reuse_embedding)) {
      candidates <- if (method == 'umap') c('umap3d', 'umap_3d', 'umap') else c('tsne3d', 'tsne_3d', 'tsne')
      for (name in intersect(candidates, reduction_names(d))) {
        x <- reduction_matrix(d, name)
        if (ncol(x) >= 3L) { xyz <- x[, 1:3, drop = FALSE]; source <- paste0('object:', name); break }
      }
    }
    if (is.null(xyz)) input <- input_matrix(d, source_reduction, dims)
    if (method == 'umap') {
      if (length(n_neighbors) != 1L || !is.numeric(n_neighbors) || !is.finite(n_neighbors) || n_neighbors < 2)
        stop('n_neighbors 必须至少为2。')
      if (length(min_dist) != 1L || !is.numeric(min_dist) || !is.finite(min_dist) || min_dist < 0)
        stop('min_dist 必须是非负数。')
      # 保持旧UMAP签名，已有匹配缓存仍可复用。
      signature <- hash(list(cells = d$cells, xyz = xyz, input = input,
        n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, seed = seed, n_epochs = n_epochs))
    } else {
      if (length(perplexity) != 1L || !is.numeric(perplexity) || !is.finite(perplexity) || perplexity <= 0)
        stop('perplexity必须是正数。')
      if (length(max_iter) != 1L || !is.numeric(max_iter) || !is.finite(max_iter) ||
          max_iter < 1 || max_iter != floor(max_iter) || max_iter > .Machine$integer.max)
        stop('max_iter必须是正整数。')
      if (length(theta) != 1L || !is.numeric(theta) || !is.finite(theta) || theta < 0 || theta > 1)
        stop('theta必须在0到1之间。')
      actual_perplexity <- min(perplexity, (length(d$cells) - 2) / 3)
      signature <- hash(list(method = method, cells = d$cells, xyz = xyz, input = input,
        perplexity = perplexity, max_iter = max_iter, theta = theta, seed = seed))
    }
    root <- mkdir(output_dir); results <- mkdir(file.path(root, 'results'))
    prefix <- paste0(method, '3d')
    cache <- file.path(results, paste0(prefix, '_coordinates.rds'))
    csv <- file.path(results, paste0(prefix, '_coordinates.csv'))
    params <- file.path(results, paste0(prefix, '_parameters.json'))
    if (file.exists(cache) && !isTRUE(overwrite)) {
      old <- readRDS(cache)
      if (!identical(old$signature, signature)) stop(title, '缓存与输入或参数不匹配；请使用新目录或 overwrite=TRUE。')
      message('复用匹配的', title, ' 3D坐标缓存：', cache, '（生成时间：', old$parameters$created, '）')
      return(invisible(old))
    }
    for (path in c(cache, csv, params)) check_output(path, overwrite)
    if (is.null(xyz)) {
      message('计算独立3D ', title, '：', length(d$cells), '个细胞，输入 ', input$name)
      if (method == 'umap') {
        need('uwot')
        arguments <- list(X = input$x, n_components = 3L, n_neighbors = min(as.integer(n_neighbors), nrow(input$x)-1L),
          min_dist = min_dist, metric = metric, n_threads = umap3d_threads(), n_sgd_threads = umap3d_threads(), verbose = verbose)
        if (!is.null(n_epochs)) arguments$n_epochs <- n_epochs
        xyz <- with_seed(seed, do.call(uwot::umap, arguments))
      } else {
        need('Rtsne')
        if (actual_perplexity < perplexity) message('perplexity按细胞数调整为：', actual_perplexity)
        xyz <- with_seed(seed, Rtsne::Rtsne(input$x, dims = 3L, pca = FALSE,
          check_duplicates = FALSE, perplexity = actual_perplexity, max_iter = as.integer(max_iter),
          theta = theta, num_threads = 1L, verbose = verbose)$Y)
      }
      rownames(xyz) <- d$cells; source <- paste0('computed:', input$name)
    } else message('使用已有真实3D ', title, '坐标：', source)
    parameters <- list(source = source, dims = input$dims, components = 3L, cells = length(d$cells),
      seed = seed, created = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), method = method)
    if (method == 'umap') {
      parameters$n_neighbors <- min(as.integer(n_neighbors), length(d$cells)-1L)
      parameters$min_dist <- min_dist; parameters$metric <- metric; parameters$n_epochs <- n_epochs
      parameters$uwot_version <- if (is.null(input)) NULL else as.character(utils::packageVersion('uwot'))
    } else {
      parameters$perplexity <- actual_perplexity; parameters$requested_perplexity <- perplexity
      parameters$max_iter <- max_iter; parameters$theta <- theta
      parameters$Rtsne_version <- if (is.null(input)) NULL else as.character(utils::packageVersion('Rtsne'))
    }
    value <- list(format = paste0('scToolbox-standalone-', prefix), version = 1L,
      xyz = xyz, signature = signature, parameters = parameters,
      paths = list(data = cache, coordinates = csv, parameters = params))
    saveRDS(value, cache)
    utils::write.csv(data.frame(cell = d$cells, axis_1 = xyz[,1], axis_2 = xyz[,2], axis_3 = xyz[,3]),
      csv, row.names = FALSE, fileEncoding = 'UTF-8')
    jsonlite::write_json(parameters, params, auto_unbox = TRUE, pretty = TRUE, null = 'null')
    message(title, ' 3D坐标已保存：', csv); invisible(value)
  }
  compute3d <- function(obj, output_dir, source_reduction = NULL, dims = NULL,
      n_neighbors = 30L, min_dist = 0.3, metric = 'cosine', seed = 123L,
      coordinate_file = NULL, reuse_embedding = TRUE, overwrite = FALSE,
      n_epochs = NULL, verbose = TRUE, methods = 'umap', perplexity = 30, max_iter = 1000L, theta = 0.5) {
    methods <- normalize_methods(methods)
    if (!is.null(coordinate_file) && length(methods) != 1L)
      stop('外部坐标CSV只用于一个method；请分别计算。')
    values <- setNames(lapply(methods, function(method) compute_one(obj, output_dir, method,
      source_reduction, dims, n_neighbors, min_dist, metric, seed, coordinate_file,
      reuse_embedding, overwrite, n_epochs, verbose, perplexity, max_iter, theta)), methods)
    invisible(if (length(values) == 1L) values[[1L]] else values)
  }
  load3d <- function(output_dir, obj = NULL, method = 'umap') {
    cache <- file.path(output_dir, 'results', paste0(method, '3d_coordinates.rds'))
    if (!file.exists(cache)) stop('缺少', method, ' 3D坐标缓存，请先运行 compute：', cache)
    value <- readRDS(cache)
    if (!identical(value$format, paste0('scToolbox-standalone-', method, '3d')) || !is.matrix(value$xyz) ||
        !is.numeric(value$xyz) || ncol(value$xyz) != 3L || any(!is.finite(value$xyz)) ||
        is.null(rownames(value$xyz)) || anyDuplicated(rownames(value$xyz))) stop('3D缓存无效。')
    if (!is.null(obj)) {
      cells <- as_data(obj)$cells
      if (!setequal(rownames(value$xyz), cells)) stop('缓存细胞ID不匹配。')
      value$xyz <- value$xyz[cells, , drop = FALSE]
    }
    value
  }
  json_script <- function(value) {
    need('jsonlite'); txt <- as.character(jsonlite::toJSON(value, auto_unbox = TRUE, digits = 10, null = 'null'))
    gsub('<', '\\u003c', txt, fixed = TRUE)
  }
  write_page <- function(template, payload, path, overwrite = TRUE, defaults_file = NULL) {
    check_output(path, overwrite)
    parts <- strsplit(template, '__SC_DATA__', fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop('内置网页数据占位符无效。')
    html <- paste0(parts[[1L]], json_script(payload), parts[[2L]])
    if (!is.null(defaults_file)) {
      need('jsonlite'); prefs <- jsonlite::fromJSON(defaults_file, simplifyVector = FALSE)
      if (!identical(prefs$format, 'scToolbox-3d-defaults') || !identical(prefs$version, 1L)) stop('显示参数JSON格式不正确。')
      lines <- strsplit(html, '\n', fixed = TRUE)[[1L]]
      index <- which(startsWith(lines, 'const embeddedDefaults='))
      if (length(index) != 1L) stop('内置显示参数占位符无效。')
      lines[index] <- paste0('const embeddedDefaults=', json_script(prefs), ';')
      html <- paste(lines, collapse = '\n')
    }
    writeLines(enc2utf8(html), path, useBytes = TRUE); normalizePath(path, winslash = '/', mustWork = TRUE)
  }
  open_page <- function(path, launch.browser) {
    if (isTRUE(launch.browser)) tryCatch(utils::browseURL(normalizePath(path, winslash = '/', mustWork = TRUE)),
      error = function(e) message('请手动打开网页：', path, '；', conditionMessage(e)))
  }
  demo_data <- function() with_seed(27L, {
    cells <- paste0('demo_cell_', seq_len(180))
    parent <- rep(c('T_cell','B_cell','Monocyte'), each = 60)
    subtype <- rep(c('CD4_T','CD8_T','Naive_B','Memory_B','CD14_Mono','FCGR3A_Mono'), each = 30)
    x <- matrix(rnorm(180*8), 180, 8, dimnames = list(cells, paste0('PC_', 1:8)))
    x[,1] <- x[,1] + rep(c(-4,0,4), each = 60); x[,2] <- x[,2] + rep(c(-1,1), 3, each = 30)
    expr <- rbind(CD3D = pmax(0, rnorm(180, ifelse(parent=='T_cell', 3, 0), .4)),
      MS4A1 = pmax(0, rnorm(180, ifelse(parent=='B_cell', 3, 0), .4)))
    colnames(expr) <- cells
    need('SeuratObject')
    obj <- SeuratObject::CreateSeuratObject(counts = round(expr * 10),
      meta.data = data.frame(celltype = parent, subtype = subtype, row.names = cells))
    if (inherits(obj[['RNA']], 'Assay5')) {
      SeuratObject::LayerData(obj, assay = 'RNA', layer = 'data') <- expr
    } else obj <- SeuratObject::SetAssayData(obj, assay = 'RNA', slot = 'data', new.data = expr)
    obj[['pca']] <- SeuratObject::CreateDimReducObject(embeddings = x, key = 'PC_', assay = 'RNA')
    obj
  })
  template3d <- paste(c(
    "<!doctype html><html lang=\"zh-CN\"><meta charset=\"utf-8\"><title>UMAP · 3D cell atlas</title>",
    "<style>",
    ":root{color-scheme:dark}*{box-sizing:border-box}html{height:100%}body{margin:0;height:100dvh;min-height:520px;display:flex;flex-direction:column;background:#000000;color:#e6eefc;font:14px system-ui}header{flex:none;padding:8px 20px 10px;border-bottom:1px solid #29354a;background:#0a1120}h2{margin:0 0 5px;font-size:21px}p{color:#a8b8cd;margin:5px 0;line-height:1.45}button,select,input[type=\"number\"],input[type=\"text\"]{background:#17253b;color:white;border:1px solid #536780;padding:6px 8px;border-radius:5px}button{cursor:pointer}button:disabled{cursor:default;opacity:.5}input[type=\"range\"]{width:132px}input[type=\"color\"]{width:34px;height:28px;padding:1px;border:1px solid #536780;border-radius:4px;background:transparent}.quick-controls,.control-panel,.preset-bar{display:flex;align-items:center;flex-wrap:wrap;gap:10px 22px}.quick-controls{padding:3px 0 8px}.quick-controls label,.control-panel label{display:flex;align-items:center;gap:8px;margin:0;white-space:nowrap}.control-tabs{display:flex;gap:2px;overflow-x:auto;border-bottom:1px solid #34445d}.control-tab{border:0;border-radius:4px 4px 0 0;background:transparent;color:#9fb0c7;padding:8px 14px;white-space:nowrap}.control-tab[aria-selected=\"true\"]{background:#1a2a43;color:#fff;box-shadow:inset 0 -2px #55dceb}.control-panel{min-height:49px;padding:10px 4px 8px}.control-panel .command-group{display:flex;align-items:center;gap:7px;flex-wrap:wrap}.preset-bar{border-top:1px solid #202d42;padding-top:8px}.section-label{color:#7f93ad;font-size:12px;font-weight:650;text-transform:uppercase}.status-inline{color:#a8b8cd}.gene-scale{color:#9ce8ee;white-space:nowrap}.gene-lookup-status{color:#f4bd6a;font-size:12px;max-width:520px}#geneSearch{width:150px}#geneFeature{min-width:110px;max-width:210px}main{display:flex;flex:1;min-height:0}canvas{display:block;position:absolute;inset:0;width:100%;height:100%;touch-action:none}#view{flex:1;min-width:0;min-height:0;position:relative;overflow:hidden}#bloom{position:absolute;inset:0;pointer-events:none;mix-blend-mode:screen;opacity:0}aside{width:220px;padding:18px;overflow:auto}i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:8px}#legend .group-row{display:flex;align-items:center;gap:7px;margin:10px 0}#legend .group-row label{flex:1;min-width:0;margin:0;overflow-wrap:anywhere}#legend input[type=\"color\"]{width:30px;height:26px;padding:0;border:1px solid #536780;border-radius:4px;background:transparent;cursor:pointer;flex:none}#legend .palette-help{font-size:12px;line-height:1.5}[hidden]{display:none!important}",
    "@media(max-width:700px){body{height:auto;min-height:100dvh}header{padding:12px}.quick-controls,.control-panel{align-items:flex-start;gap:10px 14px}.control-panel{min-height:92px}.control-panel label{white-space:normal}.control-panel label:not(.toggle){width:100%;justify-content:space-between}.control-tabs{margin-left:-12px;margin-right:-12px;padding:0 12px}input[type=\"range\"]{width:min(48vw,190px)}main{flex-direction:column;flex:none;height:auto}#view{height:65vh;flex:none}aside{box-sizing:border-box;width:100%;display:flex;flex-wrap:wrap;gap:12px;padding:12px}aside label{margin:0}}",
    "</style>",
    "<style>",
    "#headerChordSelector{margin:10px 0 1px;padding:10px 12px;border:1px solid #2b7280;border-radius:7px;background:#102938;display:flex;align-items:center;gap:12px;flex-wrap:wrap}",
    "#headerChordSelector strong{color:#9ce8ee}#headerChordSelector .target{color:#f4bd6a}#headerChordSelector button{background:#177d90;border-color:#54c9dc;font-weight:700}#headerChordSelector button.secondary{background:#17253b;border-color:#536780;font-weight:500}#headerChordStatus{color:#9ce8ee}",
    "</style>",
    "<header><h2 id=\"pageTitle\" hidden></h2><p id=\"status\" hidden></p>",
    "<div id=\"headerChordSelector\" hidden><strong id=\"headerChordProject\"></strong><span class=\"target\" id=\"headerChordTarget\"></span><button id=\"headerChordBack\" class=\"secondary\">返回上一个项目</button><button id=\"selectHeaderView\">选定此角度</button><span id=\"headerChordStatus\" role=\"status\" aria-live=\"polite\"></span></div>",
    "<div class=\"quick-controls\"><label>嵌入 <select id=\"embedding\"><option value=\"umap\">UMAP 3D</option><option value=\"tsne\" disabled>t-SNE 3D（未载入）</option></select></label><label class=\"toggle\"><input id=\"showGroupLabels\" type=\"checkbox\">显示细胞群名称</label></div>",
    "<nav class=\"control-tabs\" role=\"tablist\" aria-label=\"三维图参数分类\"><button class=\"control-tab\" id=\"geneTab\" data-panel=\"genePanel\" role=\"tab\">基因表达光效</button><button class=\"control-tab\" data-panel=\"cellPanel\" role=\"tab\">细胞粒子</button><button class=\"control-tab\" data-panel=\"membranePanel\" role=\"tab\">轮廓膜</button><button class=\"control-tab\" data-panel=\"viewPanel\" role=\"tab\">视角与旋转</button><button class=\"control-tab\" data-panel=\"exportPanel\" role=\"tab\">图片与视频</button></nav>",
    "<section class=\"control-panel\" id=\"genePanel\" role=\"tabpanel\" hidden><label id=\"geneControl\" hidden>基因 <input id=\"geneSearch\" type=\"text\" list=\"geneChoices\" autocomplete=\"off\" placeholder=\"搜索已载入基因\"><datalist id=\"geneChoices\"></datalist><select id=\"geneFeature\" aria-label=\"已载入基因\"><option value=\"\">不显示基因光效</option></select></label><span id=\"geneLookupStatus\" class=\"gene-lookup-status\" role=\"status\"></span><span id=\"geneScale\" class=\"gene-scale\" hidden></span><label id=\"geneGainControl\" hidden>光效增亮 <input title=\"仅调整光效亮度，不改变表达数值\" id=\"geneGain\" type=\"range\" min=\"0\" max=\"20\" step=\".1\" value=\"5\"></label><label id=\"geneSizeControl\" hidden>光点大小 <input id=\"geneSize\" type=\"range\" min=\".5\" max=\"3\" step=\".1\" value=\"1.2\"></label><label id=\"geneOpacityControl\" hidden>光场透明度 <input id=\"geneOpacity\" type=\"range\" min=\"0\" max=\"1\" step=\".01\" value=\".72\"></label><label id=\"genePaletteControl\" hidden>基因光色 <select id=\"genePalette\"><option value=\"#53e6fa\">冰蓝</option><option value=\"#ffc477\">暖金</option><option value=\"#ff789e\">玫红</option><option value=\"#88f0bc\">薄荷绿</option><option value=\"#bfa4ff\">淡紫</option><option value=\"#ffffff\">纯白</option><option value=\"custom\">自定义</option></select><input id=\"geneColor\" type=\"color\" value=\"#53e6fa\" aria-label=\"基因发光颜色\" title=\"自定义基因发光颜色\"></label><label id=\"geneStyleControl\" hidden>光效 <select id=\"geneStyle\"><option value=\"solid\">颗粒点</option><option value=\"soft\" selected>柔光</option></select></label><label id=\"geneSpreadControl\" hidden>光晕扩散 <input id=\"geneSpread\" type=\"range\" min=\"1\" max=\"10\" step=\".1\" value=\"4\"></label><label id=\"geneBreathControl\" class=\"toggle\" hidden><input id=\"geneBreath\" type=\"checkbox\">呼吸灯</label><label id=\"geneBreathPeriodControl\" hidden>呼吸周期 <input id=\"geneBreathPeriod\" type=\"range\" min=\"1\" max=\"10\" step=\".1\" value=\"3\">秒</label><label id=\"geneBasePointsControl\" class=\"toggle\" hidden><input id=\"geneBasePoints\" type=\"checkbox\" checked>叠加细胞类型粒子</label></section>",
    "<section class=\"control-panel\" id=\"cellPanel\" role=\"tabpanel\" hidden><span class=\"section-label\">细胞粒子</span><input id=\"mode\" type=\"hidden\" value=\"solid\"><label>细胞光晕 <input id=\"halo\" type=\"range\" min=\"0\" max=\"5\" step=\".05\" value=\".8\"></label><label>点大小 <input id=\"size\" type=\"range\" min=\".5\" max=\"3\" step=\".1\" value=\"1.2\"></label><label>点透明度 <input id=\"pointOpacity\" type=\"range\" min=\"0\" max=\"1\" step=\".01\" value=\".72\"></label><label>点亮度 <input id=\"pointBrightness\" type=\"range\" min=\"0\" max=\"5\" step=\".05\" value=\"1.8\"></label><label>点颗粒度 <input id=\"pointGrain\" type=\"range\" min=\"0\" max=\"1\" step=\".01\" value=\".65\"></label><label class=\"toggle\"><input id=\"showPoints\" type=\"checkbox\" checked>显示细胞点</label><label class=\"toggle\"><input id=\"showInside\" type=\"checkbox\" checked>显示膜内的点</label></section>",
    "<section class=\"control-panel\" id=\"membranePanel\" role=\"tabpanel\" hidden><span class=\"section-label\">轮廓膜</span><label class=\"toggle\"><input id=\"membrane\" type=\"checkbox\">显示发光透明轮廓膜</label><label>贴合度 <input id=\"membraneFit\" type=\"range\" min=\"0\" max=\"1\" step=\".01\" value=\".9\"></label><label>透明度 <input id=\"membraneOpacity\" type=\"range\" min=\"0\" max=\"1\" step=\".01\" value=\".45\"></label><label>亮度 <input id=\"membraneBrightness\" type=\"range\" min=\"0\" max=\"5\" step=\".05\" value=\"1.8\"></label></section>",
    "<section class=\"control-panel\" id=\"viewPanel\" role=\"tabpanel\" hidden><span class=\"section-label\">视角</span><label class=\"toggle\"><input id=\"autoRotate\" type=\"checkbox\" title=\"从当前视角绕屏幕竖直轴旋转\">自动旋转</label><label>旋转速度 <input id=\"rotateSpeed\" type=\"range\" min=\".1\" max=\"2\" step=\".1\" value=\".6\"></label><label>前后层次 <input id=\"depthCue\" type=\"range\" min=\"0\" max=\"1\" step=\".05\" value=\".7\"></label><label class=\"toggle\"><input id=\"showAxes\" type=\"checkbox\" checked>方向坐标轴</label><div class=\"command-group\"><button id=\"view12\">维度1–维度2正视图</button><button id=\"bestView\" title=\"减少勾选细胞群的遮挡，兼顾各群展示面积\">自动最佳视角</button><span id=\"bestViewStatus\" class=\"status-inline\" role=\"status\"></span><button id=\"reset\">重置视角</button></div></section>",
    "<section class=\"control-panel\" id=\"exportPanel\" role=\"tabpanel\" hidden><span class=\"section-label\">图片</span><label>图片比例 <select id=\"imageRatio\"><option value=\"4:3\">4:3</option><option value=\"1:1\" selected>1:1</option><option value=\"16:9\">16:9</option></select></label><label>图片清晰度 <select id=\"imageQuality\"><option value=\"standard\">标准 · 长边1600</option><option value=\"high\" selected>高清 · 长边2400</option><option value=\"ultra\">超清 · 长边3600</option></select></label><button id=\"save\">保存 PNG</button><span id=\"imageStatus\" class=\"status-inline\" role=\"status\"></span><span class=\"section-label\">视频</span><label>镜头运动 <select id=\"recordMotion\"><option value=\"rotate\" selected>绕屏幕竖直轴旋转</option><option value=\"fixed\" title=\"镜头不旋转；开启呼吸灯时光效仍会变化\">镜头静止（保留当前视角）</option></select></label><label>格式 <select id=\"videoFormat\"><option value=\"mp4\" selected>MP4（H.264）</option><option value=\"webm\">WebM（VP8 / VP9）</option></select></label><label>帧率 <select id=\"recordFps\"><option>15</option><option>24</option><option selected>30</option><option>60</option><option>120</option></select> FPS</label><label>视频清晰度 <select id=\"videoQuality\"><option value=\"standard\">标准 · 720 × 720</option><option value=\"high\" selected>高清 · 1080 × 1080</option><option value=\"ultra\">超清 · 2160 × 2160</option></select></label><label>停止条件 <select id=\"recordStop\"><option value=\"time\">按时长</option><option value=\"turns\" selected>按圈数</option></select></label><label>旋转圈数 <input id=\"recordTurns\" type=\"number\" min=\"0.25\" max=\"10\" step=\"0.25\" value=\"1\" style=\"width:60px\"></label><label>视频时长 <select id=\"recordDuration\"><option value=\"5\">5 秒</option><option value=\"8\" selected>8 秒</option><option value=\"12\">12 秒</option></select></label><div class=\"command-group\"><button id=\"record\">生成视频</button><span id=\"recordProgress\" class=\"status-inline\" role=\"status\" aria-live=\"polite\"></span></div></section>",
    "<div class=\"preset-bar\"><span class=\"section-label\">参数方案</span><button id=\"recommendedDefaults\">推荐参数</button><button id=\"saveDefaults\">保存为默认参数</button><button id=\"factoryDefaults\">恢复初始参数</button><button id=\"importDefaults\">导入参数文件</button><input id=\"defaultsFile\" type=\"file\" accept=\".json,application/json\" hidden><span id=\"defaultsStatus\" class=\"status-inline\" role=\"status\"></span></div></header>",
    "<main><div id=\"view\"><canvas id=\"cloud\"></canvas></div><aside id=\"legend\"><p class=\"palette-help\">点击每个细胞群旁的色块调整颜色。</p><button id=\"resetColors\">恢复原始配色</button><p id=\"paletteStatus\" class=\"palette-help\" role=\"status\"></p></aside></main>",
    "<script>",
    "const data=__SC_DATA__;",
    "const headerChordSelectorData=data.headerChordSelector||null;",
    "// All cells are available; visibility controls do not modify coordinates. A single uniform scale preserves embedding geometry.",
    "const embeddedDefaults={\"format\":\"scToolbox-3d-defaults\",\"version\":1,\"updatedAt\":\"2026-09-13T05:08:10.646Z\",\"controls\":{\"embedding\":\"umap\",\"halo\":0.8,\"size\":1.2,\"pointOpacity\":0.72,\"pointBrightness\":1.35,\"pointGrain\":0.77,\"membrane\":true,\"membraneFit\":0.9,\"membraneOpacity\":0.35,\"membraneBrightness\":3.4,\"showPoints\":true,\"showInside\":false,\"autoRotate\":false,\"rotateSpeed\":0.6,\"showAxes\":true,\"recordFps\":\"60\",\"recordStop\":\"turns\",\"recordTurns\":2,\"recordDuration\":\"8\"},\"view\":{\"ax\":-0.3617995098694618,\"ay\":-2.2575927756601053,\"zoom\":1}};",
    "const recommendedDefaults={\"format\":\"scToolbox-3d-defaults\",\"version\":1,\"updatedAt\":\"2026-10-12T13:56:36.462Z\",\"controls\":{\"embedding\":\"umap\",\"geneGain\":0.9,\"geneSize\":1.1,\"geneOpacity\":0.87,\"geneColor\":\"#ffa34d\",\"geneStyle\":\"soft\",\"geneBasePoints\":true,\"geneSpread\":10,\"geneBreath\":true,\"geneBreathPeriod\":1.8,\"halo\":5,\"size\":1,\"pointOpacity\":0.72,\"pointBrightness\":0.9,\"pointGrain\":0.77,\"membrane\":true,\"membraneFit\":0,\"membraneOpacity\":0.55,\"membraneBrightness\":5,\"showPoints\":true,\"showInside\":false,\"autoRotate\":false,\"rotateSpeed\":1.2,\"showAxes\":true,\"showGroupLabels\":true,\"depthCue\":0.65,\"recordMotion\":\"rotate\",\"recordFps\":\"120\",\"videoQuality\":\"ultra\",\"videoFormat\":\"mp4\",\"recordStop\":\"turns\",\"recordTurns\":2,\"recordDuration\":\"8\",\"imageRatio\":\"1:1\",\"imageQuality\":\"high\"},\"view\":{\"ax\":5.923531132677101,\"ay\":21.35833899066692,\"zoom\":0.9462407168103326},\"palette\":{\"B\":\"#eb8b86\",\"CD14 Monocyte\":\"#caa03d\",\"CD4n T\":\"#8db23b\",\"CD8m T\":\"#00be82\",\"DC\":\"#00bdc1\",\"NK\":\"#57adeb\",\"RBC\":\"#c192ed\",\"SC & Eosinophil\":\"#ed82c6\",\"B_cell\":\"#eb8b86\",\"Dendritic_cell\":\"#caa03d\",\"Erythroid_cell\":\"#8db23b\",\"Monocyte\":\"#00be82\",\"NK_cell\":\"#00bdc1\",\"Platelet\":\"#57adeb\",\"T_cell\":\"#c192ed\",\"Unknown\":\"#ed82c6\"}};",
    "const xyz=data.xyz, tsneData=Array.isArray(data.tsne)&&data.tsne.length===xyz.length?data.tsne:null, groups=Array.isArray(data.groups)?data.groups:[data.groups], labels=Array.isArray(data.labels)?data.labels:[data.labels], colors=(Array.isArray(data.colors)?data.colors:[data.colors]).slice();",
    "const initialColors=colors.slice(),colorInputs=[];let paletteMemory=Object.create(null);",
    "const geneData=data.expression||null,geneFeature=document.querySelector('#geneFeature'),geneSearch=document.querySelector('#geneSearch'),geneChoices=document.querySelector('#geneChoices'),geneLookupStatus=document.querySelector('#geneLookupStatus'),geneGain=document.querySelector('#geneGain'),geneSize=document.querySelector('#geneSize'),geneOpacity=document.querySelector('#geneOpacity'),geneBreath=document.querySelector('#geneBreath'),geneBreathPeriod=document.querySelector('#geneBreathPeriod');",
    "const geneColor=document.querySelector('#geneColor'),genePalette=document.querySelector('#genePalette'),geneStyle=document.querySelector('#geneStyle'),geneBasePoints=document.querySelector('#geneBasePoints'),geneSpread=document.querySelector('#geneSpread');",
    "const controlTabs=[...document.querySelectorAll('.control-tab')],controlPanels=[...document.querySelectorAll('.control-panel')];",
    "function selectControlPanel(id){for(const panel of controlPanels)panel.hidden=panel.id!==id;for(const tab of controlTabs){const selected=tab.dataset.panel===id;tab.setAttribute('aria-selected',String(selected));tab.tabIndex=selected?0:-1}}",
    "for(const tab of controlTabs){tab.setAttribute('aria-controls',tab.dataset.panel);tab.onclick=()=>selectControlPanel(tab.dataset.panel)}",
    "function syncGeneAppearance(){",
    " const active=geneActive();",
    " for(const el of [geneGain,geneSize,geneOpacity,geneColor,genePalette,geneStyle,geneBasePoints,geneBreath])el.disabled=!active;",
    " geneBreathPeriod.disabled=!active||!geneBreath.checked;",
    " geneSpread.disabled=!active||geneStyle.value==='solid';",
    " genePalette.value=[...genePalette.options].some(o=>o.value===geneColor.value)?geneColor.value:'custom';",
    "}",
    "genePalette.onchange=()=>{if(genePalette.value!=='custom')geneColor.value=genePalette.value;else geneColor.click();schedule()};",
    "geneColor.oninput=()=>{syncGeneAppearance();schedule()};",
    "geneStyle.onchange=()=>{syncGeneAppearance();schedule()};",
    "geneBasePoints.onchange=syncVisibility;",
    "geneBreath.onchange=()=>{syncGeneAppearance();syncGeneBreathing(true)};",
    "geneGain.oninput=geneSize.oninput=geneOpacity.oninput=geneSpread.oninput=geneBreathPeriod.oninput=schedule;",
    "let genePulse=1,geneBreathFrame=0,geneBreathClock=0;",
    "function geneBreathAt(seconds,period){return .5+.5*Math.cos(2*Math.PI*seconds/Math.max(1,period))}",
    "function geneBreathFromDarkAt(seconds,period){return 1-geneBreathAt(seconds,period)}",
    "function stopGeneBreathing(reset=true){if(geneBreathFrame){cancelAnimationFrame(geneBreathFrame);geneBreathFrame=0}geneBreathClock=0;if(reset){genePulse=1;schedule()}}",
    "function geneBreathTick(now){if(!geneBreath.checked||!geneActive()||recording){stopGeneBreathing();return}if(!geneBreathClock)geneBreathClock=now;genePulse=geneBreathAt((now-geneBreathClock)/1000,Number(geneBreathPeriod.value));schedule();geneBreathFrame=requestAnimationFrame(geneBreathTick)}",
    "function syncGeneBreathing(reset=false){if(!geneBreath.checked||!geneActive()){stopGeneBreathing();return}if(reset||!geneBreathClock){geneBreathClock=performance.now();genePulse=1}if(!geneBreathFrame&&!recording)geneBreathFrame=requestAnimationFrame(geneBreathTick);schedule()}",
    "let geneValues=null,geneMaximum=0,geneMinimum=0;",
    "function geneActive(){return geneValues!==null}",
    "function geneLevel(i){return geneActive()&&geneMaximum>0?Math.min(1,Math.max(0,geneValues[i]/geneMaximum)):0}",
    "function geneValuesFor(k){const row=geneData.values[k];if(geneData.encoding!=='sparse')return row;const values=new Float64Array(xyz.length);row.indices.forEach((i,j)=>values[i]=row.values[j]);return values}",
    "function selectGene(){",
    " const k=Number(geneFeature.value);geneValues=geneFeature.value===''?null:geneValuesFor(k);",
    " geneMaximum=geneActive()?Number(geneData.maxima[k]):0;geneMinimum=geneActive()?geneValues.reduce((low,value)=>Math.min(low,value),Infinity):0;",
    " geneSearch.value=geneActive()?geneData.names[k]:'';geneLookupStatus.textContent=`已载入 ${geneData.names.length.toLocaleString()} 个基因`;",
    " document.querySelector('#geneScale').textContent=geneActive()?`${geneData.names[k]} | ${geneData.layer} | 表达范围：${geneMinimum.toPrecision(4)} - ${geneMaximum.toPrecision(4)}`:'';",
    " geneGain.disabled=!geneActive();",
    " if(glState)updateGeneBuffers();syncVisibility();syncGeneBreathing(true);schedule();",
    "}",
    "if(geneData){",
    " if(!Array.isArray(geneData.names)||!geneData.names.length||!Array.isArray(geneData.values)||geneData.names.length!==geneData.values.length||!Array.isArray(geneData.maxima)||geneData.maxima.length!==geneData.names.length||geneData.maxima.some(v=>!Number.isFinite(v)||v<0)||geneData.values.some(a=>{if(geneData.encoding!=='sparse')return !Array.isArray(a)||a.length!==xyz.length||a.some(v=>!Number.isFinite(v)||v<0);return !a||!Array.isArray(a.indices)||!Array.isArray(a.values)||a.indices.length!==a.values.length||a.indices.some((v,i)=>!Number.isInteger(v)||v<0||v>=xyz.length||(i>0&&v<=a.indices[i-1]))||a.values.some(v=>!Number.isFinite(v)||v<0)}))throw new Error('基因表达数据与细胞坐标不匹配');",
    " for(const id of ['geneControl','geneGainControl','geneSizeControl','geneOpacityControl','geneScale','genePaletteControl','geneStyleControl','geneBasePointsControl','geneSpreadControl','geneBreathControl','geneBreathPeriodControl'])document.getElementById(id).hidden=false;",
    " geneData.names.forEach((name,k)=>{const opt=document.createElement('option'),choice=document.createElement('option');opt.value=String(k);opt.textContent=name;choice.value=name;geneFeature.append(opt);geneChoices.append(choice)});",
    " geneFeature.value='0';geneValues=geneValuesFor(0);geneMaximum=Number(geneData.maxima[0]);",
    " geneSearch.value=geneData.names[0];geneFeature.onchange=selectGene;",
    " const chooseSearchedGene=()=>{const query=geneSearch.value.trim();if(!query){geneFeature.value='';selectGene();return true}const k=geneData.names.findIndex(name=>name.toLowerCase()===query.toLowerCase());if(k<0){geneLookupStatus.textContent=query?`“${query}”未载入（共 ${geneData.names.length.toLocaleString()} 个基因）。请检查 features 并重新生成网页。`:'';return false}geneFeature.value=String(k);selectGene();return true};",
    " geneSearch.onchange=chooseSearchedGene;geneSearch.onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();chooseSearchedGene()}};",
    " document.querySelector('#geneTab').hidden=false;",
    "}",
    "selectControlPanel(geneData?'genePanel':'cellPanel');",
    "if(!geneData)geneLookupStatus.textContent='未载入表达量：gene_load设为all或selected后重绘；selected需填写features。';",
    "const canvas=document.querySelector('#cloud');let ctx=null;const bloomCanvas=document.createElement('canvas');bloomCanvas.id='bloom';document.querySelector('#view').appendChild(bloomCanvas);const bloomCtx=bloomCanvas.getContext('2d'),off=document.createElement('canvas'),oc=off.getContext('2d'),geneOff=document.createElement('canvas'),gc=geneOff.getContext('2d');",
    "const pageTitle=document.querySelector('#pageTitle'),embedding=document.querySelector('#embedding'),mode=document.querySelector('#mode'),halo=document.querySelector('#halo'),size=document.querySelector('#size'),pointOpacity=document.querySelector('#pointOpacity'),pointBrightness=document.querySelector('#pointBrightness'),pointGrain=document.querySelector('#pointGrain'),autoRotate=document.querySelector('#autoRotate'),rotateSpeed=document.querySelector('#rotateSpeed'),recordDuration=document.querySelector('#recordDuration');",
    "const videoFormat=document.querySelector('#videoFormat'),membrane=document.querySelector('#membrane'),membraneFit=document.querySelector('#membraneFit'),membraneOpacity=document.querySelector('#membraneOpacity'),membraneBrightness=document.querySelector('#membraneBrightness'),showPoints=document.querySelector('#showPoints'),showInside=document.querySelector('#showInside');",
    "if(tsneData){embedding.options[1].disabled=false;embedding.options[1].textContent='t-SNE 3D'}",
    "if(Array.isArray(data.availableEmbeddings)){for(const option of embedding.options)option.disabled=!data.availableEmbeddings.includes(option.value);embedding.value=data.availableEmbeddings[0]}",
    "const embeddingData={umap:xyz,tsne:tsneData};",
    "document.querySelector('#status').textContent=`${xyz.length.toLocaleString()} 个细胞已载入 · ${data.stage==='expression'?'基因表达':data.stage==='prepare'?'CellTypist优先候选标签（复核可选）':'已应用的审核标签'}`;",
    "let enabled=labels.map(()=>true),ax=.2,ay=-.35,az=0,zoom=1,drag=null,pending=false;",
    "let exportViewport=null;",
    "const depthCue=document.querySelector('#depthCue'),showGroupLabels=document.querySelector('#showGroupLabels');",
    "let labelAnchors=[];",
    "showGroupLabels.onchange=schedule;",
    "depthCue.oninput=schedule;",
    "let activeEmbedding='umap',points=[],members=[],interior=[],meshes=[],viewStats=[],glState=null,rotationAnchor={ax:.2,ay:-.35,az:0},autoFrame=0,autoClock=0,autoAngle=0,recording=false;",
    "function normalizeEmbedding(raw){const source=raw.map(p=>[Number(p[0]),Number(p[1]),Number.isFinite(Number(p[2]))?Number(p[2]):0]),lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];source.forEach(p=>p.forEach((v,j)=>{lo[j]=Math.min(lo[j],v);hi[j]=Math.max(hi[j],v)}));const mid=lo.map((v,j)=>(v+hi[j])/2),span=Math.max(...hi.map((v,j)=>hi[j]-lo[j]))||1;return source.map(p=>p.map((v,j)=>(v-mid[j])/span))}",
    "// Equal-sized samples keep large groups from overwhelming small groups in view scoring.",
    "function rebuildViewStats(){",
    " const directions=[];for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++)if(x||y||z)directions.push([x,y,z]);",
    " const budget=Math.max(24,Math.min(256,Math.floor(8192/Math.max(1,members.length))));",
    " viewStats=members.map(ids=>{",
    "  const center=[0,0,0],extreme=new Array(directions.length).fill(-Infinity),supportIds=new Array(directions.length).fill(-1);",
    "  for(const i of ids){const p=points[i];for(let k=0;k<3;k++)center[k]+=p[k];for(let k=0;k<directions.length;k++){const d=directions[k],v=p[0]*d[0]+p[1]*d[1]+p[2]*d[2];if(v>extreme[k]){extreme[k]=v;supportIds[k]=i}}}",
    "  if(ids.length)for(let k=0;k<3;k++)center[k]/=ids.length;",
    "  const chosen=new Set(supportIds.filter(i=>i>=0)),count=Math.min(ids.length,budget);",
    "  for(let k=0;k<count;k++)chosen.add(ids[Math.floor((k+.5)*ids.length/count)]);",
    "  return {n:ids.length,center,sample:[...chosen].map(i=>points[i])};",
    " });",
    "}",
    "function prepareEmbedding(name){const raw=embeddingData[name]||embeddingData.umap;activeEmbedding=embeddingData[name]?name:'umap';points=normalizeEmbedding(raw);members=labels.map(()=>[]);points.forEach((p,i)=>members[groups[i]].push(i));interior=points.map(()=>false);meshes=[];members.forEach((ids,j)=>buildDensityMesh(ids,j));rebuildViewStats();rebuildLabelAnchors();embedding.value=activeEmbedding;const title=activeEmbedding==='tsne'?'t-SNE 3D':'UMAP 3D';pageTitle.textContent=activeEmbedding==='tsne'?'t-SNE · 三维细胞点云':'UMAP · 三维细胞点云';document.title=`${title} · cell atlas`;membrane.disabled=false;document.querySelector('#status').textContent=`${xyz.length.toLocaleString()} 个细胞已载入 · ${data.stage==='expression'?'基因表达':data.stage==='prepare'?'CellTypist优先候选标签（复核可选）':'已应用的审核标签'} · 当前图：${title}`;if(glState)rebuildWebGL();syncVisibility();schedule()}",
    "// Smoothed 3D occupancy density + marching tetrahedra. No sphere or convex hull.",
    "// Each disconnected density island retains its own surface.",
    "function buildDensityMesh(ids,group){",
    " if(ids.length<10)return;",
    " const n=ids.length>1000?64:52,at=(x,y,z)=>(x*n+y)*n+z;",
    " const low=[Infinity,Infinity,Infinity],high=[-Infinity,-Infinity,-Infinity];",
    " ids.forEach(i=>points[i].forEach((v,k)=>{low[k]=Math.min(low[k],v);high[k]=Math.max(high[k],v)}));",
    " const step=low.map((v,k)=>Math.max(high[k]-v,.006)/(n-9));",
    " low.forEach((v,k)=>low[k]=v-4*step[k]);",
    " let field=new Float32Array(n*n*n);",
    " const cell=p=>p.map((v,k)=>Math.max(0,Math.min(n-1,Math.round((v-low[k])/step[k]))));",
    " ids.forEach(i=>{const q=cell(points[i]);field[at(...q)]++});",
    " const sortedStep=step.slice().sort((a,b)=>a-b),sigma=sortedStep[1]*1.05;",
    " const kernels=step.map(s=>{const w=[-3,-2,-1,0,1,2,3].map(v=>Math.exp(-((v*s)**2)/(2*sigma*sigma))),sum=w.reduce((a,b)=>a+b,0);return w.map(v=>v/sum)});",
    " for(let axis=0;axis<3;axis++){",
    "  const out=new Float32Array(field.length);",
    "  for(let x=0;x<n;x++)for(let y=0;y<n;y++)for(let z=0;z<n;z++){",
    "   const index=at(x,y,z),coord=axis===0?x:axis===1?y:z,stride=axis===0?n*n:axis===1?n:1;let value=0;",
    "   for(let t=-3;t<=3;t++)if(coord+t>=0&&coord+t<n)value+=field[index+t*stride]*kernels[axis][t+3];",
    "   out[at(x,y,z)]=value;",
    "  }field=out;",
    " }",
    " const pointDensity=ids.map(i=>field[at(...cell(points[i]))]).filter(v=>v>0).sort((a,b)=>a-b);",
    " const fieldDensity=Array.from(field).filter(v=>v>0).sort((a,b)=>a-b),defaultFit=.9,coverage=.98-.16*defaultFit;",
    " const pointLevel=pointDensity[Math.min(pointDensity.length-1,Math.floor((1-coverage)*Math.max(0,pointDensity.length-1)))]||0;",
    " const floorLevel=fieldDensity[Math.floor(.02*Math.max(0,fieldDensity.length-1))]||0,level=Math.max(pointLevel,floorLevel);",
    " if(!(level>0))return;",
    " ids.forEach(i=>interior[i]=field[at(...cell(points[i]))]>level*1.35);",
    " const corners=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];",
    " const tetra=[[0,5,1,6],[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6]],triangles=[];",
    " function triangle(a,b,c){const u=b.map((v,k)=>v-a[k]),v=c.map((v,k)=>v-a[k]);const normal=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],len=Math.hypot(...normal);if(len>1e-14)triangles.push({p:[a,b,c],normal:normal.map(v=>v/len)})}",
    " for(let x=0;x<n-1;x++)for(let y=0;y<n-1;y++)for(let z=0;z<n-1;z++){",
    "  const q=corners.map(c=>[x+c[0],y+c[1],z+c[2]]),values=q.map(c=>field[at(...c)]);",
    "  if(values.every(v=>v<level)||values.every(v=>v>=level))continue;",
    "  const p=q.map(c=>c.map((v,k)=>low[k]+v*step[k]));",
    "  function edge(a,b){const t=(level-values[a])/(values[b]-values[a]);return p[a].map((v,k)=>v+t*(p[b][k]-v))}",
    "  for(const t of tetra){const inside=t.filter(i=>values[i]>=level),outside=t.filter(i=>values[i]<level);",
    "   if(inside.length===1||inside.length===3){const a=inside.length===1?inside:outside,b=inside.length===1?outside:inside;triangle(...b.map(i=>edge(a[0],i)))}",
    "   else if(inside.length===2){const [a,b]=inside,[c,d]=outside,A=edge(a,c),B=edge(a,d),C=edge(b,c),D=edge(b,d);triangle(A,B,C);triangle(B,D,C)}",
    "  }",
    " }",
    " // Shared vertices eliminate cracks; Taubin pairs smooth without cumulative shrinkage.",
    " const vertices=[],lookup=new Map(),faces=[],neighbors=[];",
    " for(const t of triangles){const face=t.p.map(p=>{const key=p.map(v=>v.toFixed(9)).join(',');if(!lookup.has(key)){lookup.set(key,vertices.length);vertices.push(p.slice());neighbors.push(new Set())}return lookup.get(key)});",
    "  if(new Set(face).size<3)continue;faces.push(face);for(let k=0;k<3;k++){neighbors[face[k]].add(face[(k+1)%3]);neighbors[face[k]].add(face[(k+2)%3])}",
    " }",
    " for(let iter=0;iter<4;iter++)for(const lambda of [.45,-.47]){",
    "  const next=vertices.map((p,i)=>{const ids=Array.from(neighbors[i]);if(!ids.length)return p.slice();return p.map((v,k)=>v+lambda*(ids.reduce((a,j)=>a+vertices[j][k],0)/ids.length-v))});",
    "  next.forEach((p,i)=>vertices[i]=p);",
    " }",
    " // Density-gradient normals vary smoothly even across tetrahedron diagonals.",
    " function gradient(p){const q=cell(p);return [0,1,2].map(k=>{const a=q.slice(),b=q.slice();a[k]=Math.max(0,q[k]-1);b[k]=Math.min(n-1,q[k]+1);return (field[at(...b)]-field[at(...a)])/(2*step[k])})}",
    " const normals=vertices.map(p=>{const g=gradient(p),len=Math.hypot(...g);return len>1e-12?g.map(v=>v/len):[0,0,1]});",
    " const smooth=faces.map(face=>{const normal=[0,1,2].map(k=>face.reduce((v,i)=>v+normals[i][k],0)),len=Math.hypot(...normal);return {p:face.map(i=>vertices[i]),normal:len>1e-12?normal.map(v=>v/len):[0,0,1]}});",
    " const center=[0,1,2].map(k=>ids.reduce((sum,i)=>sum+points[i][k],0)/ids.length);",
    " meshes.push({group,center,triangles:smooth});",
    "}",
    "// Row-major world-to-view rotation; az retains rotations that cannot fit two angles.",
    "function viewMatrix(a,b,c=0){",
    " const cx=Math.cos(a),sx=Math.sin(a),cy=Math.cos(b),sy=Math.sin(b),cz=Math.cos(c),sz=Math.sin(c);",
    " return [cy*cz,-cy*sz,sy,sx*sy*cz+cx*sz,-sx*sy*sz+cx*cz,-sx*cy,-cx*sy*cz+sx*sz,cx*sy*sz+sx*cz,cx*cy];",
    "}",
    "function applyViewMatrix(p,m){return [m[0]*p[0]+m[1]*p[1]+m[2]*p[2],m[3]*p[0]+m[4]*p[1]+m[5]*p[2],m[6]*p[0]+m[7]*p[1]+m[8]*p[2]]}",
    "function multiplyViewMatrix(a,b){const m=new Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)m[r*3+c]=a[r*3]*b[c]+a[r*3+1]*b[3+c]+a[r*3+2]*b[6+c];return m}",
    "function viewPoseFromMatrix(m){",
    " const ay=Math.asin(Math.max(-1,Math.min(1,m[2]))),regular=Math.hypot(m[0],m[1])>1e-8;",
    " return {ax:regular?Math.atan2(-m[5],m[8]):Math.atan2(m[3]*Math.sign(m[2]),m[4]),ay,az:regular?Math.atan2(-m[1],m[0]):0};",
    "}",
    "function screenRotatePose(origin,yaw,pitch=0){",
    " const m=viewMatrix(origin.ax,origin.ay,origin.az||0),yawMatrix=viewMatrix(0,yaw,0),pitchMatrix=viewMatrix(pitch,0,0);",
    " return viewPoseFromMatrix(multiplyViewMatrix(pitchMatrix,multiplyViewMatrix(yawMatrix,m)));",
    "}",
    "function setViewPose(pose){ax=pose.ax;ay=pose.ay;az=pose.az||0}",
    "function rotateScreen(p,a,b,c=0){return applyViewMatrix(p,viewMatrix(a,b,c)).slice(0,2)}",
    "function viewProjection(p,a,b,c=0){return applyViewMatrix(p,viewMatrix(a,b,c))}",
    "function viewFrame(width,height){",
    " const short=Math.max(1,Math.min(width,height));return {w:width/short,h:height/short,k:.875/.72*(width<height?.85:1)};",
    "}",
    "function fitViewDistance(projected,frame){",
    " let distance=.225;for(const p of projected)distance=Math.max(distance,p[2]+.08,p[2]+frame.k*Math.abs(p[0])/(frame.w*.42),p[2]+frame.k*Math.abs(p[1])/(frame.h*.42));",
    " return Math.max(.225,Math.min(6.75,distance));",
    "}",
    "function viewScore(a,b,stats,frame,c=0){",
    " const rotation=viewMatrix(a,b,c),clouds=stats.map(s=>s.sample.map(p=>applyViewMatrix(p,rotation))),distance=fitViewDistance(clouds.flat(),frame);",
    " const width=Math.max(24,Math.round(96*frame.w/Math.max(frame.w,frame.h))),height=Math.max(24,Math.round(96*frame.h/Math.max(frame.w,frame.h))),cells=width*height;",
    " const front=new Float32Array(cells).fill(-Infinity),overlap=new Uint16Array(cells),layers=[];",
    " for(const cloud of clouds){",
    "  const depth=new Float32Array(cells).fill(-Infinity),occupied=[];",
    "  for(const p of cloud){const d=distance-p[2];if(d<=.02)continue;const px=Math.round((.5+frame.k*p[0]/d/frame.w)*(width-1)),py=Math.round((.5-frame.k*p[1]/d/frame.h)*(height-1));",
    "   for(let dy=-2;dy<=2;dy++)for(let dx=-2;dx<=2;dx++){if(dx*dx+dy*dy>5)continue;const x=px+dx,y=py+dy;if(x<0||x>=width||y<0||y>=height)continue;const at=y*width+x;if(depth[at]===-Infinity)occupied.push(at);depth[at]=Math.max(depth[at],p[2])}",
    "  }",
    "  for(const at of occupied){front[at]=Math.max(front[at],depth[at]);overlap[at]++}layers.push({depth,occupied});",
    " }",
    " let visibility=0,area=0,worst=1,visibleGroups=0,union=0,shared=0;",
    " for(let at=0;at<cells;at++)if(overlap[at]){union++;if(overlap[at]>1)shared++}",
    " for(const {depth,occupied} of layers){let exposed=0;for(const at of occupied)exposed+=depth[at]>=front[at]-.008?1/Math.sqrt(overlap[at]):.08;",
    "  const fraction=exposed/Math.max(1,occupied.length);visibility+=fraction;worst=Math.min(worst,fraction);area+=Math.sqrt(exposed/cells);if(fraction>=.3&&exposed>=Math.min(3,occupied.length))visibleGroups++;",
    " }",
    " const n=Math.max(1,stats.length),score=6*visibleGroups/n+2*visibility/n+.8*worst+4*area/n+2*union/cells-1.5*shared/Math.max(1,union);",
    " return {score,ax:a,ay:b,az:c,zoom:1.35/distance,visibleGroups,visibility:visibility/n,area:area/n,overlap:shared/Math.max(1,union)};",
    "}",
    "function findBestView(stats,frame,start){",
    " let top=[viewScore(start.ax,start.ay,stats,frame,start.az||0)];const add=(a,b,c=0)=>{top.push(viewScore(a,b,stats,frame,c));top.sort((x,y)=>y.visibleGroups-x.visibleGroups||y.score-x.score);top=top.slice(0,4)};",
    " for(let ia=-4;ia<=4;ia++)for(let ib=-8;ib<8;ib++)add(ia*Math.PI/8,ib*Math.PI/8);",
    " for(const step of [Math.PI/16,Math.PI/32,Math.PI/64]){const seeds=top.slice();for(const seed of seeds)for(let da=-1;da<=1;da++)for(let db=-1;db<=1;db++)add(seed.ax+da*step,seed.ay+db*step,seed.az)}",
    " return top[0];",
    "}",
    "function fittedViewZoom(a,b,frame,c=0){",
    " const rotation=viewMatrix(a,b,c);let distance=.225;const include=p=>{const q=applyViewMatrix(p,rotation);distance=Math.max(distance,q[2]+.08,q[2]+frame.k*Math.abs(q[0])/(frame.w*.42),q[2]+frame.k*Math.abs(q[1])/(frame.h*.42))};",
    " points.forEach((p,i)=>{if(enabled[groups[i]])include(p)});",
    " if(membrane.checked){const fit=1.025-.05*Number(membraneFit.value);for(const mesh of meshes)if(enabled[mesh.group])for(const triangle of mesh.triangles)for(const p of triangle.p)include(p.map((v,k)=>mesh.center[k]+(v-mesh.center[k])*fit))}",
    " return Math.max(.2,Math.min(6,1.35/distance));",
    "}",
    "async function bestView(){",
    " const stats=viewStats.filter((s,j)=>enabled[j]&&s.n);const button=document.querySelector('#bestView'),status=document.querySelector('#bestViewStatus');",
    " if(!stats.length){status.textContent='先勾选细胞群';return}if(recording||button.disabled)return;",
    " stopAutoRotate();autoRotate.checked=false;button.disabled=true;status.textContent='正在比较视角…';",
    " try{await new Promise(requestAnimationFrame);if(recording)return;const box=canvas.getBoundingClientRect(),frame=viewFrame(Math.max(1,box.width),Math.max(1,box.height)),best=findBestView(stats,frame,{ax,ay,az});",
    "  setViewPose(best);zoom=fittedViewZoom(ax,ay,frame,az);rotationAnchor={ax,ay,az};schedule();status.textContent=`已优化 ${stats.length} 个细胞群的遮挡与展示面积`;",
    " }finally{button.disabled=recording}",
    "}",
    "function visiblePoint(i){return enabled[groups[i]] && showPoints.checked && (geneActive() || !membrane.checked || showInside.checked || !interior[i])}",
    "function visibleCellPoint(i){return enabled[groups[i]] && showPoints.checked && (!geneActive() || geneBasePoints.checked) && (!membrane.checked || showInside.checked || !interior[i])}",
    "function syncVisibility(){syncGeneAppearance();membrane.disabled=recording;showInside.disabled=recording || !showPoints.checked || !membrane.checked || (geneActive() && !geneBasePoints.checked);showInside.title=recording?'导出完成后可调整':!showPoints.checked?'先开启显示细胞点':!membrane.checked?'先开启轮廓膜':geneActive()&&!geneBasePoints.checked?'先开启叠加细胞类型粒子':'控制细胞类型底图的膜内点，基因光点单独显示';membraneFit.disabled=recording||!membrane.checked;membraneOpacity.disabled=recording||!membrane.checked;membraneBrightness.disabled=recording||!membrane.checked;schedule()}",
    "function setGroupColor(group,hex,notify=true){",
    " if(!Number.isInteger(group)||group<0||group>=labels.length||typeof hex!=='string'||!/^#[0-9a-f]{6}$/i.test(hex))throw new Error('颜色需要#RRGGBB格式');",
    " colors[group]=hex.toLowerCase();paletteMemory[labels[group]]=colors[group];",
    " if(colorInputs[group])colorInputs[group].value=colors[group];",
    " if(glState){const rgb=hexRgb(colors[group]);glState.pointBuffers[group].color=rgb;glState.membraneBuffers[group].color=rgb}",
    " if(notify)document.querySelector('#paletteStatus').textContent=`${labels[group]}：${colors[group]}。可点击“保存为默认参数”保留配色。`;",
    " schedule();",
    "}",
    "function applyPalette(palette){paletteMemory=Object.assign(Object.create(null),palette||{});labels.forEach((label,j)=>setGroupColor(j,Object.prototype.hasOwnProperty.call(paletteMemory,label)?paletteMemory[label]:initialColors[j],false))}",
    "labels.forEach((label,j)=>{const row=document.createElement('div'),toggle=document.createElement('label'),cb=document.createElement('input'),picker=document.createElement('input');row.className='group-row';cb.type='checkbox';cb.checked=true;toggle.append(cb,document.createTextNode(label));picker.type='color';picker.value=colors[j];picker.title=`调整 ${label} 的颜色`;picker.setAttribute('aria-label',`${label} 颜色`);colorInputs[j]=picker;row.append(toggle,picker);document.querySelector('#legend').append(row);cb.onchange=()=>{enabled[j]=cb.checked;schedule()};picker.oninput=()=>setGroupColor(j,picker.value)});",
    "document.querySelector('#resetColors').onclick=()=>{initialColors.forEach((hex,j)=>setGroupColor(j,hex,false));document.querySelector('#paletteStatus').textContent='已恢复当前数据集的原始配色；点击“保存为默认参数”可保留。'};",
    "function schedule(){if(!pending){pending=true;requestAnimationFrame(()=>draw())}}",
    "let geneSpriteCache=null;",
    "function getGeneSprites(){",
    " const key=geneColor.value;if(geneSpriteCache?.key===key)return geneSpriteCache;",
    " const rgb=hexRgb(geneColor.value),makeSprite=profile=>{const sprite=document.createElement('canvas');sprite.width=sprite.height=128;const c=sprite.getContext('2d'),im=c.createImageData(128,128);",
    "  for(let y=0;y<128;y++)for(let x=0;x<128;x++){const qx=(x+.5)/64-1,qy=(y+.5)/64-1,r2=qx*qx+qy*qy,r=Math.sqrt(r2),fade=Math.max(0,Math.min(1,(1-r)/.25)),smooth=fade*fade*(3-2*fade),k=(y*128+x)*4,a=profile(qx,qy,r2)*smooth;im.data[k]=rgb[0]*255;im.data[k+1]=rgb[1]*255;im.data[k+2]=rgb[2]*255;im.data[k+3]=Math.min(1,a)*255}c.putImageData(im,0,0);return sprite};",
    " const soft=makeSprite((qx,qy,r2)=>.65*Math.exp(-28*r2)+.18*Math.exp(-5*r2));",
    " geneSpriteCache={key,soft};return geneSpriteCache;",
    "}",
    "function draw(){pending=false;if(!ctx)ctx=canvas.getContext('2d');const box=canvas.getBoundingClientRect(),dpr=exportViewport?exportViewport.dpr:Math.min(devicePixelRatio||1,2);canvas.width=exportViewport?exportViewport.width:Math.round(box.width*dpr);canvas.height=exportViewport?exportViewport.height:Math.round(box.height*dpr);off.width=geneOff.width=canvas.width;off.height=geneOff.height=canvas.height;",
    "ctx.fillStyle='#000000';ctx.fillRect(0,0,canvas.width,canvas.height);const rotation=viewMatrix(ax,ay,az);",
    "const projected=[];points.forEach((p,i)=>{if(!visiblePoint(i))return;const q=projectScenePoint(p,canvas.width,canvas.height,rotation);if(q)projected.push([q.x,q.y,q.z,groups[i],q.depth,i])});projected.sort((a,b)=>a[2]-b[2]);",
    "const grainAmount=Number(pointGrain.value),geneSprites=geneActive()&&geneStyle.value!=='solid'?getGeneSprites():null;",
    "for(const p of projected){const baseR=Math.max(.5,Number(size.value)*dpr*1.25/p[4]),geneR=Math.max(.5,Number(geneSize.value)*dpr*1.25/p[4]);",
    "const depthFade=1-Number(depthCue.value)*(.65*(1-Math.max(0,Math.min(1,(p[2]+.65)/1.3))));",
    "if(visibleCellPoint(p[5])){oc.globalCompositeOperation='source-over';oc.fillStyle=colors[p[3]];oc.globalAlpha=Number(pointOpacity.value)*depthFade;oc.beginPath();oc.arc(p[0],p[1],baseR,0,Math.PI*2);oc.fill();if(grainAmount>0&&baseR>1.5){const seed=p[0]*12.9898+p[1]*78.233+p[3]*37.719;for(let k=0;k<3;k++){const n=Math.abs(Math.sin(seed+k*17.13)*43758.5453)%1;if(n>.72){const a=Math.abs(Math.sin(seed+k*5.31)*6.283),rr=baseR*(.18+.16*(Math.abs(Math.sin(seed+k*9.17))%1))*grainAmount;oc.globalAlpha=Number(pointOpacity.value)*.42*grainAmount*depthFade;oc.beginPath();oc.arc(p[0]+Math.cos(a)*baseR*.42,p[1]+Math.sin(a)*baseR*.42,Math.max(.35,rr),0,Math.PI*2);oc.fill()}}}}",
    "if(!geneActive())continue;const e=geneLevel(p[5]);",
    "if(geneSprites){const radius=geneR*Number(geneSpread.value),gain=Number(geneGain.value);gc.globalCompositeOperation='lighter';gc.globalAlpha=genePulse*Number(geneOpacity.value)*Math.min(1,(.008+e*gain)*.3);gc.drawImage(geneSprites.soft,p[0]-radius,p[1]-radius,radius*2,radius*2);gc.globalCompositeOperation='source-over';continue;}",
    "gc.fillStyle=geneColor.value;gc.globalAlpha=genePulse*Number(geneOpacity.value)*Math.min(1,.025+e*Number(geneGain.value));gc.beginPath();gc.arc(p[0],p[1],geneR,0,Math.PI*2);gc.fill();",
    "}",
    "const pointGlow=Math.max(0,Number(pointBrightness.value)),haloStrength=Math.max(0,Number(halo.value));",
    "if(pointGlow>0 && haloStrength>0){ctx.save();ctx.globalCompositeOperation='screen';const passes=Math.max(1,Math.ceil(pointGlow));for(let pass=0;pass<passes;pass++){ctx.globalAlpha=Math.min(1,(pointGlow/passes)*(haloStrength/2));ctx.filter=`blur(${(mode.value==='glow'?5:2)*dpr}px)`;ctx.drawImage(off,0,0)}ctx.restore()}",
    "ctx.save();ctx.globalAlpha=1;ctx.drawImage(off,0,0);ctx.globalCompositeOperation='screen';ctx.drawImage(geneOff,0,0);ctx.restore();",
    "// Render translucent density surface facets into a soft offscreen layer.",
    "// No stroke, wireframe, or artificial outer ring.",
    "if(membrane.checked){",
    " oc.clearRect(0,0,off.width,off.height);oc.globalCompositeOperation='source-over';",
    " const membraneAlpha=Number(membraneOpacity.value),membraneGlow=Math.max(0,Number(membraneBrightness.value));",
    " for(const mesh of meshes){if(!enabled[mesh.group])continue;oc.fillStyle=colors[mesh.group];",
    "  const fitScale=1.025-.05*Number(membraneFit.value);",
    "  for(const t of mesh.triangles){const nz=applyViewMatrix(t.normal,rotation)[2];",
    "   oc.globalAlpha=membraneAlpha*(.045+.45*(1-Math.abs(nz))**2)*(1-Number(depthCue.value)*(nz<0?.65:.12));",
    "   oc.beginPath();t.p.forEach((p,i)=>{const q=p.map((v,k)=>mesh.center[k]+(v-mesh.center[k])*fitScale);const screen=projectScenePoint(q,canvas.width,canvas.height,rotation);if(!screen)return;if(i===0)oc.moveTo(screen.x,screen.y);else oc.lineTo(screen.x,screen.y)});oc.closePath();oc.fill();",
    "  }",
    " }",
    " ctx.save();ctx.globalCompositeOperation='screen';const passes=Math.max(1,Math.ceil(membraneGlow));",
    " for(let pass=0;pass<passes;pass++){ctx.globalAlpha=Math.min(1,membraneGlow/passes);ctx.filter=`blur(${1.2*dpr}px)`;ctx.drawImage(off,0,0);ctx.globalAlpha=Math.min(1,membraneGlow/(passes*1.8));ctx.filter=`blur(${4*dpr}px)`;ctx.drawImage(off,0,0)}ctx.restore();",
    "}",
    "if(!exportViewport){ctx.fillStyle='#d6e3f8';ctx.font=`${12*dpr}px system-ui`;ctx.fillText(`${activeEmbedding==='tsne'?'t-SNE 3D':'UMAP 3D'} | ${projected.length} / ${points.length} cells | ${mode.value} | ${data.stage} | membrane ${membrane.checked?'ON':'OFF'} | point opacity ${Number(pointOpacity.value).toFixed(2)} | point brightness ${Number(pointBrightness.value).toFixed(1)}`,15*dpr,24*dpr);}",
    "}",
    "canvas.onpointerdown=e=>{if(autoRotate.checked){autoRotate.checked=false;stopAutoRotate()}drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId)};canvas.onpointermove=e=>{if(!drag)return;setViewPose(screenRotatePose({ax,ay,az},(e.clientX-drag[0])*.006,(e.clientY-drag[1])*.006));document.querySelector('#bestViewStatus').textContent='';drag=[e.clientX,e.clientY];schedule()};canvas.onpointerup=canvas.onpointercancel=()=>drag=null;",
    "canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(.2,Math.min(6,zoom*Math.exp(-e.deltaY*.001)));schedule()};",
    "const axisCanvas=document.createElement('canvas');",
    "axisCanvas.style.cssText='position:absolute;inset:0;pointer-events:none';",
    "document.querySelector('#view').appendChild(axisCanvas);",
    "const axisCtx=axisCanvas.getContext('2d'),showAxes=document.querySelector('#showAxes');",
    "function drawDirectionAxes(c,w,h,dpr,anchor=null,showHint=true){",
    " if(!showAxes.checked)return;",
    " const ox=anchor?anchor.x:76*dpr,oy=anchor?anchor.y:h-72*dpr,len=35*dpr,rotation=viewMatrix(ax,ay,az);",
    " c.save();c.fillStyle='rgba(0,0,0,.68)';c.beginPath();c.arc(ox,oy,62*dpr,0,2*Math.PI);c.fill();",
    " c.font=`${10*dpr}px system-ui`;c.textAlign='center';c.textBaseline='middle';c.lineCap='round';",
    " const prefix=activeEmbedding==='tsne'?'tSNE':'UMAP';",
    " const axes=[[1,0,0],[0,1,0],[0,0,1]].map((p,j)=>{const [x,y,z]=applyViewMatrix(p,rotation);return {j,x,y,z}}).sort((a,b)=>a.z-b.z);",
    " for(const a of axes){const color=['#ef8799','#8bdfbb','#8eabff'][a.j];",
    "  c.strokeStyle=color;c.globalAlpha=.25;c.lineWidth=dpr;c.setLineDash([2*dpr,3*dpr]);c.beginPath();c.moveTo(ox,oy);c.lineTo(ox-a.x*len*.65,oy+a.y*len*.65);c.stroke();",
    "  c.setLineDash(a.z<-.05?[3*dpr,3*dpr]:[]);c.globalAlpha=a.z<-.05?.4:1;c.lineWidth=1.6*dpr;c.beginPath();c.moveTo(ox,oy);c.lineTo(ox+a.x*len,oy-a.y*len);c.stroke();",
    "  c.fillStyle=color;c.beginPath();c.arc(ox+a.x*len,oy-a.y*len,2.5*dpr,0,2*Math.PI);c.fill();",
    "  c.fillText(`${prefix}${a.j+1}`,ox+a.x*(len+17*dpr),oy-a.y*(len+17*dpr));",
    " }c.globalAlpha=.8;c.fillStyle='#a8b8cd';c.font=`${9*dpr}px system-ui`;if(showHint)c.fillText('亮实线近 · 暗虚线远',ox,oy+57*dpr);c.restore();",
    "}",
    "// Anchor each label to a real cell nearest its group's centroid, computed once per embedding.",
    "function projectScenePoint(p,w,h,rotation=viewMatrix(ax,ay,az)){",
    " const [x,y,z]=applyViewMatrix(p,rotation),depth=1.35/zoom-z;if(depth<=.02)return null;",
    " const frame=exportViewport,ref=frame?frame.reference:{width:w,height:h},rw=ref.width,rh=ref.height,k=Math.min(rw,rh)*.875*(rw<rh?.85:1)/.72;",
    " const X=rw/2+x*k/depth,Y=rh/2-y*k/depth;",
    " return frame?{x:(X-frame.crop.x)*w/frame.crop.width,y:(Y-frame.crop.y)*h/frame.crop.height,z,depth}:{x:X,y:Y,z,depth};",
    "}",
    "function rebuildLabelAnchors(){labelAnchors=members.map((ids,j)=>{let best=null,dist=Infinity;for(const i of ids){const p=points[i],d=p.reduce((sum,v,k)=>sum+(v-viewStats[j].center[k])**2,0);if(d<dist){dist=d;best=p}}return best})}",
    "function projectGroupLabel(p,w,h){return projectScenePoint(p,w,h)}",
    "function drawGeneScale(c,w,h,dpr){",
    " // Gene identifier only; expression range stays in the control panel.",
    " if(!geneActive()||!geneData)return;",
    " const name=String(geneData.names[Number(geneFeature.value)]||'').trim();if(!name)return;",
    " c.save();try{",
    "  c.font=`600 ${16*dpr}px system-ui`;c.globalAlpha=1;c.textAlign='left';c.textBaseline='top';",
    "  const lines=wrapLegendText(c,`Gene: ${name}`,Math.max(1,Math.min(320*dpr,w-32*dpr))),width=Math.max(...lines.map(line=>c.measureText(line).width))+12*dpr;",
    "  c.fillStyle='rgba(0,0,0,.78)';c.fillRect(10*dpr,10*dpr,width,(lines.length*21+10)*dpr);",
    "  c.fillStyle=geneColor.value;lines.forEach((line,i)=>c.fillText(line,16*dpr,(15+i*21)*dpr));",
    " }finally{c.restore()}",
    "}",
    "function drawGroupLabels(c,w,h,dpr){",
    " if(!showGroupLabels.checked)return;",
    " const candidates=labelAnchors.map((p,j)=>p&&enabled[j]?{...projectGroupLabel(p,w,h),group:j}:null).filter(p=>p&&Number.isFinite(p.x)&&p.x>=0&&p.x<=w&&p.y>=0&&p.y<=h).sort((a,b)=>b.z-a.z),occupied=[];",
    " c.save();c.font=`${13*dpr}px system-ui`;c.textAlign='center';c.textBaseline='middle';",
    " for(const p of candidates){",
    "  const lines=wrapLegendText(c,labels[p.group],Math.min(190*dpr,w-24*dpr)),bw=Math.max(...lines.map(t=>c.measureText(t).width))+16*dpr,bh=(lines.length*17+10)*dpr;",
    "  const options=[[0,-bh/2-12*dpr],[0,bh/2+12*dpr],[-bw/2-15*dpr,0],[bw/2+15*dpr,0],[0,-bh-26*dpr],[0,bh+26*dpr]];",
    "  let box=null,best=Infinity;",
    "  for(const [dx,dy] of options){const b={x:Math.max(4*dpr,Math.min(w-bw-4*dpr,p.x+dx-bw/2)),y:Math.max(4*dpr,Math.min(h-bh-4*dpr,p.y+dy-bh/2)),w:bw,h:bh};const cost=occupied.reduce((sum,a)=>sum+Math.max(0,Math.min(a.x+a.w,b.x+b.w)-Math.max(a.x,b.x))*Math.max(0,Math.min(a.y+a.h,b.y+b.h)-Math.max(a.y,b.y)),0);if(cost<best){best=cost;box=b}if(cost===0)break}",
    "  occupied.push(box);c.strokeStyle=colors[p.group];c.globalAlpha=.75;c.lineWidth=dpr;c.beginPath();c.moveTo(p.x,p.y);c.lineTo(box.x+box.w/2,box.y+box.h/2);c.stroke();",
    "  c.globalAlpha=1;c.fillStyle='rgba(0,0,0,.78)';c.fillRect(box.x,box.y,box.w,box.h);c.fillStyle=colors[p.group];c.fillRect(box.x,box.y,3*dpr,box.h);",
    "  c.fillStyle='#f4f8ff';lines.forEach((line,i)=>c.fillText(line,box.x+box.w/2,box.y+(5+8.5+i*17)*dpr));",
    " }c.restore();",
    "}",
    "function drawAxisOverlay(){const dpr=Math.min(devicePixelRatio||1,2);if(axisCanvas.width!==canvas.width||axisCanvas.height!==canvas.height){axisCanvas.width=canvas.width;axisCanvas.height=canvas.height}axisCtx.clearRect(0,0,axisCanvas.width,axisCanvas.height);drawDirectionAxes(axisCtx,axisCanvas.width,axisCanvas.height,dpr);drawGroupLabels(axisCtx,axisCanvas.width,axisCanvas.height,dpr);drawGeneScale(axisCtx,axisCanvas.width,axisCanvas.height,dpr)}",
    "showAxes.onchange=schedule;",
    "const recordFps=document.querySelector('#recordFps'),recordMotion=document.querySelector('#recordMotion'),recordStop=document.querySelector('#recordStop'),recordTurns=document.querySelector('#recordTurns'),recordProgress=document.querySelector('#recordProgress');",
    "function recordingPlan(fps,stop,turns,seconds,speed,motion='rotate'){",
    " if(![15,24,30,60,120].includes(fps)||!Number.isFinite(speed)||speed<=0)throw new Error('帧率或旋转速度无效');",
    " if(!['rotate','fixed'].includes(motion))throw new Error('视频镜头模式无效');",
    " if(!Number.isFinite(seconds)||seconds<=0)throw new Error('录制时长无效');",
    " if(motion==='fixed')return {fps,duration:seconds,angle:0,motion};",
    " if(stop==='turns'&&(!Number.isFinite(turns)||turns<.25||turns>10))throw new Error('旋转圈数请填写0.25到10之间的数字');",
    " const duration=stop==='turns'?turns*2*Math.PI/(speed*.55):seconds;",
    " return {fps,duration,angle:stop==='turns'?turns*2*Math.PI:seconds*speed*.55,motion};",
    "}",
    "function rotationPose(origin,angle){",
    " if(Math.abs(angle%(2*Math.PI))<1e-12)return {ax:origin.ax,ay:origin.ay,az:origin.az||0};",
    " return screenRotatePose(origin,angle);",
    "}",
    "function syncRecordingOptions(){const fixed=recordMotion.value==='fixed';recordStop.disabled=fixed;recordDuration.disabled=!fixed&&recordStop.value==='turns';recordTurns.disabled=fixed||recordStop.value!=='turns'}",
    "recordMotion.onchange=recordStop.onchange=syncRecordingOptions;syncRecordingOptions();",
    "function stopAutoRotate(){if(autoFrame){cancelAnimationFrame(autoFrame);autoFrame=0}autoClock=0;autoAngle=0}",
    "function autoTick(now){if(!autoRotate.checked){stopAutoRotate();return}const speed=Math.max(.1,Number(rotateSpeed.value));autoAngle+=Math.max(0,now-autoClock)/1000*speed*.55;autoClock=now;setViewPose(rotationPose(rotationAnchor,autoAngle));schedule();autoFrame=requestAnimationFrame(autoTick)}",
    "function toggleAutoRotate(){stopAutoRotate();if(autoRotate.checked){rotationAnchor={ax,ay,az};autoClock=performance.now();autoFrame=requestAnimationFrame(autoTick)}}",
    "// A fixed square composition, independent of the browser viewport.",
    "function videoQualityPlan(quality,fps){const p={standard:{width:720,min:2500000,max:12000000},high:{width:1080,min:6000000,max:30000000},ultra:{width:2160,min:20000000,max:180000000}}[quality];if(!p)throw new Error('视频清晰度无效');return {width:p.width,scale:p.width/1080,bitrate:Math.round(Math.min(p.max,Math.max(p.min,p.width*p.width*fps*.24)))}}",
    "function videoLayout(width=1080){const margin=20,gap=10,legendWidth=210,plotSize=width-2*margin-gap-legendWidth;return {width,height:width,plotX:margin,plotY:(width-plotSize)/2,plotSize,legendX:margin+plotSize+gap,legendWidth}}",
    "function wrapLegendText(c,text,width){const words=String(text).split(/\\s+/),lines=[];let line='';for(const word of words){let next=line?line+' '+word:word;if(c.measureText(next).width<=width){line=next;continue}if(line){lines.push(line);line=''}for(const letter of word){if(c.measureText(line+letter).width>width&&line){lines.push(line);line=''}line+=letter}}if(line)lines.push(line);return lines.length?lines:['']}",
    "function drawRecordingLegend(c,layout){",
    " const visible=labels.map((label,group)=>({label,group})).filter(item=>enabled[item.group]);",
    " let font=17,rows,total;",
    " do{c.font=`${font}px system-ui`;rows=visible.map(item=>({...item,lines:wrapLegendText(c,item.label,layout.legendWidth-28)}));total=rows.reduce((sum,row)=>sum+row.lines.length*(font+4)+9,0);if(total<=layout.height-205||font<=10)break;font--}while(true);",
    " if(total>layout.height-205)throw new Error('当前图例过多，无法在正方形视频中清晰放下；请减少展示的细胞群');",
    " const start=Math.max(28,Math.min(layout.height-total-140,(layout.sceneCenterY||layout.height/2)-(total+115)/2));let y=start;",
    " c.save();c.font=`${font}px system-ui`;c.textBaseline='top';c.textAlign='left';",
    " for(const row of rows){c.fillStyle=colors[row.group];c.beginPath();c.arc(layout.legendX+6,y+font*.55,5,0,2*Math.PI);c.fill();c.fillStyle='#e6eefc';for(const line of row.lines){c.fillText(line,layout.legendX+20,y);y+=font+4}y+=9}",
    " c.restore();return {x:layout.legendX+75,y:Math.min(layout.height-75,y+65)};",
    "}",
    "function captureViewport(){const box=canvas.getBoundingClientRect();return {width:Math.max(1,box.width),height:Math.max(1,box.height),dpr:Math.min(devicePixelRatio||1,2)}}",
    "function exportContentBounds(reference,motion='fixed'){",
    " const w=reference.width,h=reference.height,rotation=viewMatrix(ax,ay,az),distance=1.35/zoom,k=Math.min(w,h)*.875*(w<h?.85:1)/.72;",
    " let left=Infinity,right=-Infinity,top=Infinity,bottom=-Infinity;",
    " const include=p=>{",
    "  const q=applyViewMatrix(p,rotation);",
    "  if(motion==='rotate'){",
    "   const radius=Math.hypot(q[0],q[2]);if(distance-radius<=.02){left=0;right=w;top=0;bottom=h;return}",
    "   const X=k*radius/Math.sqrt(distance*distance-radius*radius),Y1=h/2-k*q[1]/(distance-radius),Y2=h/2-k*q[1]/(distance+radius);",
    "   left=Math.min(left,w/2-X);right=Math.max(right,w/2+X);top=Math.min(top,Y1,Y2);bottom=Math.max(bottom,Y1,Y2);",
    "  }else{",
    "   const depth=distance-q[2];if(depth<=.02)return;",
    "   const X=w/2+k*q[0]/depth,Y=h/2-k*q[1]/depth;left=Math.min(left,X);right=Math.max(right,X);top=Math.min(top,Y);bottom=Math.max(bottom,Y);",
    "  }",
    " };",
    " points.forEach((p,i)=>{if(visiblePoint(i))include(p)});",
    " if(membrane.checked)for(const mesh of meshes)if(enabled[mesh.group])for(const triangle of mesh.triangles)for(const p of triangle.p){const fit=1.025-.05*Number(membraneFit.value);include(p.map((v,j)=>mesh.center[j]+(v-mesh.center[j])*fit))}",
    " if(!Number.isFinite(left)||right<0||left>w||bottom<0||top>h)return {x:0,y:0,width:w,height:h};",
    " const bloom=bloomSettings(),haloPad=bloom.level>0?bloom.blur*2:8;",
    " const labelPad=showGroupLabels.checked?26:0;",
    " const pointPad=Math.max(Number(size.value),geneActive()?Number(geneSize.value)*(geneStyle.value==='soft'?Number(geneSpread.value):1):0)*4,pad=Math.max(12,haloPad,labelPad,pointPad);",
    " left=Math.max(0,left-pad);right=Math.min(w,right+pad);top=Math.max(0,top-pad-(geneActive()?42:0));bottom=Math.min(h,bottom+pad);",
    " return {x:left,y:top,width:Math.max(1,right-left),height:Math.max(1,bottom-top)};",
    "}",
    "function sceneLayout(layout,reference,scale,crop={x:0,y:0,width:reference.width,height:reference.height}){",
    " const gap=22,margin=24,availableWidth=layout.width-2*margin-gap-layout.legendWidth,availableHeight=layout.height-2*margin;",
    " const factor=Math.min(availableWidth/crop.width,availableHeight/crop.height),pw=crop.width*factor,ph=crop.height*factor,totalWidth=pw+gap+layout.legendWidth;",
    " const x=(layout.width-totalWidth)/2,y=(layout.height-ph)/2,width=Math.max(1,Math.round(pw*scale)),height=Math.max(1,Math.round(ph*scale));",
    " layout.plotX=x;layout.plotY=y;layout.plotWidth=pw;layout.plotHeight=ph;layout.legendX=x+pw+gap;layout.sceneCenterY=y+ph/2;",
    " return {x:Math.round(x*scale),y:Math.round(y*scale),width,height,dpr:height/crop.height,screenDpr:reference.dpr,reference,crop,",
    "  transform:[reference.width/crop.width,reference.height/crop.height,(reference.width-2*crop.x-crop.width)/crop.width,(2*crop.y+crop.height-reference.height)/crop.height]};",
    "}",
    "function setExportViewport(target,reference,motion='fixed'){target.scene=sceneLayout(target.layout,reference,target.renderScale||1,exportContentBounds(reference,motion));exportViewport=target.scene}",
    "function bloomSettings(){const cellLevel=showPoints.checked&&(!geneActive()||geneBasePoints.checked)?Number(halo.value)/8:0,level=Math.min(.9,Math.max(0,cellLevel+(membrane.checked?Number(membraneBrightness.value)/14:0)));return {level,blur:(4+10*level)*(exportViewport?exportViewport.screenDpr:Math.min(devicePixelRatio||1,2))}}",
    "function drawRecordingFrame(target,targetCtx){",
    " const l=target.layout,s=target.renderScale||1,r=target.scene;",
    " targetCtx.save();try{",
    "  targetCtx.clearRect(0,0,target.width,target.height);targetCtx.fillStyle='#000000';targetCtx.fillRect(0,0,target.width,target.height);",
    "  // The GPU bitmap is already at final resolution; do not resize the particle layer.",
    "  targetCtx.drawImage(canvas,r.x,r.y);",
    "  const bloom=bloomSettings();if(bloom.level>0&&glState){targetCtx.save();targetCtx.globalCompositeOperation='screen';targetCtx.globalAlpha=bloom.level;targetCtx.filter=`blur(${bloom.blur*r.dpr}px)`;targetCtx.drawImage(bloomCanvas,r.x,r.y);targetCtx.restore()}",
    "  targetCtx.save();targetCtx.translate(r.x,r.y);drawGroupLabels(targetCtx,r.width,r.height,s);drawGeneScale(targetCtx,r.width,r.height,s);targetCtx.restore();",
    "  targetCtx.scale(s,s);const anchor=drawRecordingLegend(targetCtx,l);drawDirectionAxes(targetCtx,l.width,l.height,1,anchor,false);",
    " }finally{targetCtx.restore()}",
    "}",
    "function imageQualityPlan(ratio,quality){",
    " const aspect={'4:3':4/3,'1:1':1,'16:9':16/9}[ratio],width={standard:1600,high:2400,ultra:3600}[quality];",
    " if(!aspect||!width)throw new Error('图片比例或清晰度无效');",
    " const height=Math.round(width/aspect),baseWidth=1080*aspect,margin=20,gap=10,legendWidth=230;",
    " const plotSize=Math.min(1080-2*margin,baseWidth-2*margin-gap-legendWidth),plotX=margin+(baseWidth-2*margin-gap-legendWidth-plotSize)/2;",
    " return {width,height,scale:width/baseWidth,layout:{width:baseWidth,height:1080,plotX,plotY:(1080-plotSize)/2,plotSize,legendX:plotX+plotSize+gap,legendWidth}};",
    "}",
    "async function savePng(){",
    " if(recording)return;",
    " const button=document.querySelector('#save'),status=document.querySelector('#imageStatus'),oldView={ax,ay,az,zoom},wasAuto=autoRotate.checked,wasBreathing=geneBreath.checked,oldPulse=genePulse;",
    " const controls=[...document.querySelectorAll('header input,header select,header button,#legend input,#legend button')].map(el=>({el,disabled:el.disabled}));",
    " const target=document.createElement('canvas'),reference=captureViewport();",
    " recording=true;stopAutoRotate();stopGeneBreathing(false);controls.forEach(({el})=>el.disabled=true);button.textContent='正在导出…';",
    " try{",
    "  const plan=imageQualityPlan(document.querySelector('#imageRatio').value,document.querySelector('#imageQuality').value);",
    "  target.width=plan.width;target.height=plan.height;target.layout=plan.layout;target.renderScale=plan.scale;",
    "  setExportViewport(target,reference);",
    "  draw();drawRecordingFrame(target,target.getContext('2d'));",
    "  const blob=await new Promise((resolve,reject)=>target.toBlob(value=>value?resolve(value):reject(new Error('PNG编码失败')),'image/png'));",
    "  const link=document.createElement('a'),url=URL.createObjectURL(blob);link.href=url;link.download=`${activeEmbedding==='tsne'?'tSNE_3D':'UMAP_3D'}_${mode.value}_${plan.width}x${plan.height}.png`;link.click();setTimeout(()=>URL.revokeObjectURL(url),60000);",
    "  status.textContent=`已导出 PNG · ${plan.width} × ${plan.height} · 含细胞群图例`;",
    " }catch(e){status.textContent=`导出失败：${e.message}`}",
    " finally{target.width=target.height=0;exportViewport=null;recording=false;controls.forEach(({el,disabled})=>el.disabled=disabled);button.textContent='保存 PNG';setViewPose(oldView);zoom=oldView.zoom;genePulse=oldPulse;draw();if(wasBreathing)syncGeneBreathing(true);if(wasAuto)toggleAutoRotate()}",
    "}",
    "",
    "// WebCodecs + a small WebM muxer: timestamps come from frame indices, never wall time.",
    "function ebmlSize(n){for(let len=1;len<=8;len++){if(n<2**(7*len)-1){let out=new Uint8Array(len),v=n;for(let j=len-1;j>=0;j--){out[j]=v%256;v=Math.floor(v/256)}out[0]|=1<<(8-len);return out}}throw new Error('WebM文件过大')}",
    "function ebmlUInt(n){let len=1;while(n>=2**(8*len))len++;const out=new Uint8Array(len);for(let j=len-1;j>=0;j--){out[j]=n%256;n=Math.floor(n/256)}return out}",
    "function ebmlFloat(n){const out=new Uint8Array(8);new DataView(out.buffer).setFloat64(0,n);return out}",
    "function ebml(id,payload){const bytes=Uint8Array.from(id.match(/../g).map(x=>parseInt(x,16))),body=payload instanceof Blob?payload:new Blob([payload]);return new Blob([bytes,ebmlSize(body.size),body])}",
    "function webmVideo(chunks,width,height,fps,frameCount,codec){",
    " const uint=(id,n)=>ebml(id,ebmlUInt(n)),str=(id,s)=>ebml(id,new TextEncoder().encode(s));",
    " const header=ebml('1a45dfa3',new Blob([uint('4286',1),uint('42f7',1),uint('42f2',4),uint('42f3',8),str('4282','webm'),uint('4287',4),uint('4285',2)]));",
    " const info=ebml('1549a966',new Blob([uint('2ad7b1',1000),str('4d80','scToolbox'),str('5741','scToolbox'),ebml('4489',ebmlFloat(frameCount*1e6/fps))]));",
    " const video=ebml('e0',new Blob([uint('b0',width),uint('ba',height)]));",
    " const tracks=ebml('1654ae6b',ebml('ae',new Blob([uint('d7',1),uint('73c5',1),uint('83',1),uint('9c',0),str('86',codec==='vp8'?'V_VP8':'V_VP9'),uint('23e383',Math.round(1e9/fps)),video])));",
    " const clusters=chunks.map(c=>ebml('1f43b675',new Blob([uint('e7',c.timestamp),ebml('a3',new Blob([new Uint8Array([0x81,0,0,c.key?0x80:0]),c.data]))])));",
    " return new Blob([header,ebml('18538067',new Blob([info,tracks,...clusters]))],{type:'video/webm'});",
    "}",
    "// H.264/AVC samples in decode order, no B-frames (WebCodecs realtime mode).",
    "function mp4Video(chunks,width,height,fps,description){",
    " if(!description)throw new Error('H.264 编码器没有提供 MP4 配置，请切换 WebM');",
    " const cat=(...a)=>{const out=new Uint8Array(a.reduce((n,x)=>n+x.length,0));let p=0;for(const x of a){out.set(x,p);p+=x.length}return out};",
    " const u32=(...a)=>{const out=new Uint8Array(a.length*4),v=new DataView(out.buffer);a.forEach((x,i)=>v.setUint32(i*4,x));return out};",
    " const u16=(...a)=>{const out=new Uint8Array(a.length*2),v=new DataView(out.buffer);a.forEach((x,i)=>v.setUint16(i*2,x));return out};",
    " const str=s=>new TextEncoder().encode(s),zero=n=>new Uint8Array(n);",
    " const box=(name,...a)=>{const body=cat(...a);return cat(u32(body.length+8),str(name),body)};",
    " const full=(name,flags,...a)=>box(name,u32(flags),...a);",
    " const matrix=u32(65536,0,0,0,65536,0,0,0,1073741824);",
    " const count=chunks.length,scale=fps*1000,duration=count*1000;",
    " if(duration>0xffffffff)throw new Error('视频过长，请减少时长');",
    " for(let i=1;i<count;i++)if(chunks[i].timestamp<=chunks[i-1].timestamp)throw new Error('编码器产生了重排帧，请切换 WebM');",
    " const ftyp=box('ftyp',str('isom'),u32(512),str('isomiso2avc1mp41'));",
    " const avc1=box('avc1',zero(6),u16(1),zero(16),u16(width,height),u32(4718592,4718592,0),u16(1),zero(32),u16(24,65535),box('avcC',new Uint8Array(description)));",
    " const stsd=full('stsd',0,u32(1),avc1),stts=full('stts',0,u32(1,count,1000));",
    " const stsc=full('stsc',0,u32(1,1,count,1));",
    " const sizes=new Uint8Array(count*4),sv=new DataView(sizes.buffer);chunks.forEach((c,i)=>sv.setUint32(i*4,c.data.length));",
    " const stsz=full('stsz',0,u32(0,count),sizes);",
    " const keys=chunks.flatMap((c,i)=>c.key?[i+1]:[]),stss=full('stss',0,u32(keys.length),u32(...keys));",
    " const moov=offset=>box('moov',",
    "  full('mvhd',0,u32(0,0,scale,duration,65536),u16(256,0),zero(8),matrix,zero(24),u32(2)),",
    "  box('trak',full('tkhd',7,u32(0,0,1,0,duration),zero(8),u16(0,0,0,0),matrix,u32(width*65536,height*65536)),",
    "   box('mdia',full('mdhd',0,u32(0,0,scale,duration),u16(21956,0)),full('hdlr',0,u32(0),str('vide'),zero(12),str('scToolbox Video\\0')),",
    "    box('minf',full('vmhd',1,u16(0,0,0,0)),box('dinf',full('dref',0,u32(1),full('url ',1))),",
    "     box('stbl',stsd,stts,stsc,stsz,full('stco',0,u32(1,offset)),stss)))));",
    " const header=moov(0),offset=ftyp.length+header.length+8,total=chunks.reduce((n,c)=>n+c.data.length,8);",
    " if(offset+total>0xffffffff)throw new Error('MP4 超过4GB，请缩短视频');",
    " return new Blob([ftyp,moov(offset),u32(total),str('mdat'),...chunks.map(c=>c.data)],{type:'video/mp4'});",
    "}",
    "",
    "let exportCancelled=false;",
    "async function recordRotation(){",
    " const button=document.querySelector('#record');if(recording){exportCancelled=true;return}",
    " if(!window.VideoEncoder||!window.VideoFrame){window.alert('逐帧导出需要支持 WebCodecs 的新版 Chrome 或 Edge，请用它打开此网页。');return}",
    " let plan;try{plan=recordingPlan(Number(recordFps.value),recordStop.value,Number(recordTurns.value),Number(recordDuration.value),Number(rotateSpeed.value),recordMotion.value)}catch(e){window.alert(e.message);return}",
    " const frameCount=Math.max(1,Math.round(plan.duration*plan.fps));",
    " const oldView={ax,ay,az,zoom},wasAuto=autoRotate.checked,wasBreathing=geneBreath.checked,oldGenePulse=genePulse,oldPointer=canvas.style.pointerEvents;",
    " const controls=[...document.querySelectorAll('header input,header select,header button,#legend input,#legend button')].filter(el=>el!==button).map(el=>({el,disabled:el.disabled}));",
    " const target=document.createElement('canvas'),reference=captureViewport();draw();",
    " const quality=videoQualityPlan(document.querySelector('#videoQuality').value,plan.fps);target.layout=videoLayout();target.renderScale=quality.scale;target.width=target.height=quality.width;",
    " const targetCtx=target.getContext('2d'),chunks=[];let encoder=null,encodeError=null,bytes=0,avcDescription=null;const isMp4=videoFormat.value==='mp4';",
    " const embeddingName=activeEmbedding==='tsne'?'tSNE_3D':'UMAP_3D';",
    " recording=true;exportCancelled=false;autoRotate.checked=false;stopAutoRotate();stopGeneBreathing(false);controls.forEach(({el})=>el.disabled=true);canvas.style.pointerEvents='none';button.textContent='取消生成';",
    " try{",
    "  let config;",
    "  for(const codec of (isMp4?['avc1.420034','avc1.4d0034','avc1.640034']:['vp09.00.10.08','vp8'])){",
    "   const candidate={codec,width:target.width,height:target.height,framerate:plan.fps,bitrate:quality.bitrate,latencyMode:'realtime',...(isMp4?{avc:{format:'avc'}}:{})};",
    "   if((await VideoEncoder.isConfigSupported(candidate)).supported){config=candidate;break}",
    "  }",
    "  if(!config)throw new Error(isMp4?'当前浏览器不支持所选尺寸/帧率的 H.264 编码，请降低清晰度或帧率，或选择 WebM；建议使用新版 Chrome / Edge。':'当前浏览器没有可用的VP8/VP9编码器，请降低清晰度后重试');",
    "  encoder=new VideoEncoder({output(chunk,metadata){if(metadata?.decoderConfig?.description)avcDescription=new Uint8Array(metadata.decoderConfig.description).slice();const data=new Uint8Array(chunk.byteLength);chunk.copyTo(data);bytes+=data.byteLength;chunks.push({timestamp:chunk.timestamp,key:chunk.type==='key',data});if(bytes>1024**3)encodeError=new Error('视频已超过1GB，请缩短时长或减少圈数')},error(e){encodeError=e}});",
    "  encoder.configure(config);setExportViewport(target,reference,plan.motion);const origin={ax:oldView.ax,ay:oldView.ay,az:oldView.az};",
    "  for(let i=0;i<frameCount;i++){",
    "   if(exportCancelled)throw new Error('已取消生成');if(encodeError)throw encodeError;",
    "   const pose=plan.motion==='fixed'?origin:rotationPose(origin,plan.angle*i/frameCount);setViewPose(pose);genePulse=wasBreathing?geneBreathFromDarkAt(i/plan.fps,Number(geneBreathPeriod.value)):1;draw();drawRecordingFrame(target,targetCtx);",
    "   const timestamp=Math.round(i*1e6/plan.fps),duration=Math.round((i+1)*1e6/plan.fps)-timestamp;",
    "   const frame=new VideoFrame(target,{timestamp,duration});try{encoder.encode(frame,{keyFrame:i%(plan.fps*2)===0})}finally{frame.close()}",
    "   if(encoder.encodeQueueSize>=4)await encoder.flush();",
    "   recordProgress.textContent=`逐帧生成 ${i+1} / ${frameCount} 帧 · ${((i+1)/frameCount*100).toFixed(0)}% · 成片 ${plan.fps} FPS`;",
    "   await new Promise(resolve=>setTimeout(resolve,0));",
    "  }",
    "  await encoder.flush();if(encodeError)throw encodeError;if(exportCancelled)throw new Error('已取消生成');",
    "  if(chunks.length!==frameCount)throw new Error(`编码帧数不完整：${chunks.length}/${frameCount}，没有导出不完整视频`);",
    "  if(!isMp4)chunks.sort((a,b)=>a.timestamp-b.timestamp);",
    "  const blob=isMp4?mp4Video(chunks,target.width,target.height,plan.fps,avcDescription):webmVideo(chunks,target.width,target.height,plan.fps,frameCount,config.codec);const url=URL.createObjectURL(blob),a=document.createElement('a'),motionName=plan.motion==='fixed'?'fixed':'rotate';",
    "  a.download=`${embeddingName}_${motionName}_${quality.width}x${quality.width}_${plan.fps}fps_${frameCount}frames.${isMp4?'mp4':'webm'}`;a.href=url;a.click();setTimeout(()=>URL.revokeObjectURL(url),60000);",
    "  recordProgress.textContent=plan.motion==='fixed'?`已生成静止镜头视频 · ${quality.width} × ${quality.width} · ${frameCount} 帧 · ${plan.fps} FPS · ${(frameCount/plan.fps).toFixed(2)} 秒`:`已生成 ${quality.width} × ${quality.width} · ${frameCount} 帧 · ${plan.fps} FPS · ${(frameCount/plan.fps).toFixed(2)} 秒 · ${(plan.angle/(2*Math.PI)).toFixed(2)} 圈`;",
    " }catch(e){recordProgress.textContent=exportCancelled?'已取消生成':`导出失败：${e.message}`}",
    " finally{if(encoder&&encoder.state!=='closed')encoder.close();target.width=target.height=0;exportViewport=null;recording=false;controls.forEach(({el,disabled})=>el.disabled=disabled);canvas.style.pointerEvents=oldPointer;button.textContent='生成视频';setViewPose(oldView);zoom=oldView.zoom;genePulse=oldGenePulse;draw();if(wasBreathing)syncGeneBreathing(true);if(wasAuto){autoRotate.checked=true;toggleAutoRotate()}}",
    "}",
    "mode.onchange=halo.oninput=size.oninput=pointOpacity.oninput=pointBrightness.oninput=pointGrain.oninput=membraneFit.oninput=membraneOpacity.oninput=membraneBrightness.oninput=schedule;embedding.onchange=()=>prepareEmbedding(embedding.value);autoRotate.onchange=toggleAutoRotate;membrane.onchange=showPoints.onchange=showInside.onchange=syncVisibility;syncVisibility();document.querySelector('#view12').onclick=()=>{stopAutoRotate();autoRotate.checked=false;ax=0;ay=0;az=0;rotationAnchor={ax,ay,az};zoom=1;schedule()};document.querySelector('#bestView').onclick=bestView;document.querySelector('#reset').onclick=()=>{stopAutoRotate();autoRotate.checked=false;ax=.2;ay=-.35;az=0;rotationAnchor={ax,ay,az};zoom=1;enabled.fill(true);document.querySelectorAll('#legend input[type=\"checkbox\"]').forEach(c=>c.checked=true);schedule()};",
    "// WebGL renderer: the density surface is built once above, while dragging only",
    "// updates uniforms. This avoids reprojecting every membrane triangle on the CPU.",
    "function glShader(gl,type,source){const shader=gl.createShader(type);gl.shaderSource(shader,source);gl.compileShader(shader);if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS))throw new Error(gl.getShaderInfoLog(shader));return shader}",
    "function glProgram(gl,vertex,fragment){const program=gl.createProgram();gl.attachShader(program,glShader(gl,gl.VERTEX_SHADER,vertex));gl.attachShader(program,glShader(gl,gl.FRAGMENT_SHADER,fragment));gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(program));return program}",
    "function hexRgb(hex){const value=parseInt(String(hex).replace('#',''),16);return [((value>>16)&255)/255,((value>>8)&255)/255,(value&255)/255]}",
    "function rebuildWebGL(){if(!glState)return;const old=glState;for(const item of old.pointBuffers)old.gl.deleteBuffer(item.buffer);for(const item of old.membraneBuffers)old.gl.deleteBuffer(item.buffer);old.gl.deleteProgram(old.pointProgram);old.gl.deleteProgram(old.membraneProgram);glState=null;initWebGL()}",
    "function updateGeneBuffers(){",
    " if(!glState)return;const gl=glState.gl;",
    " for(const item of glState.pointBuffers){for(let k=0;k<item.ids.length;k++)item.vertices[k*5+4]=geneLevel(item.ids[k]);gl.bindBuffer(gl.ARRAY_BUFFER,item.buffer);gl.bufferSubData(gl.ARRAY_BUFFER,0,item.vertices)}",
    "}",
    "function initWebGL(){",
    " try{",
    "  const gl=canvas.getContext('webgl2',{antialias:true,alpha:false,preserveDrawingBuffer:true});if(!gl)return;",
    "  const common=`#version 300 es\\nprecision highp float;\\nuniform mat3 uView;vec3 rotatePoint(vec3 p){return uView*p;}`;",
    "  const pointProgram=glProgram(gl,common+`\\nin vec4 aPoint;in float aExpression;out float vExpression;uniform vec4 uFrame;uniform float uZoom,uAspect,uPointSize,uDpr,uPointFloor,uHideInterior,uExpression,uGeneStyle,uGeneSpread;flat out float vSeed;out float vDepth;void main(){vExpression=aExpression;if(uHideInterior>.5&&aPoint.w>.5){gl_Position=vec4(2.0);gl_PointSize=0.;return;}vec3 q=rotatePoint(aPoint.xyz);float depth=1.35/uZoom-q.z;vSeed=float(gl_VertexID);vDepth=q.z;gl_Position=vec4(q.x/(depth*.72*uAspect)*1.75*min(1.,uAspect)*(uAspect<1.?.85:1.),q.y/(depth*.72)*1.75*min(1.,uAspect)*(uAspect<1.?.85:1.),0.,1.);gl_Position.xy=gl_Position.xy*uFrame.xy+uFrame.zw;gl_PointSize=max(uPointFloor,uPointSize*uDpr*2.5/depth)*(uExpression>.5&&uGeneStyle>.5?uGeneSpread:1.);}`,",
    "   `#version 300 es\\nprecision highp float;flat in float vSeed;in float vDepth;in float vExpression;uniform float uDepthCue;uniform vec3 uColor;uniform float uOpacity,uBrightness,uHalo,uGlowMode,uGrain,uExpression,uGeneGain,uGeneStyle;out vec4 outColor;float hash12(vec2 p){return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}void main(){vec2 uv=gl_PointCoord;float d=distance(uv,vec2(.5))*2.;if(d>1.)discard;",
    "if(uExpression>.5&&uGeneStyle>.5){",
    " vec2 q=(uv-.5)*2.;float r2=dot(q,q);",
    " float taper=1.-smoothstep(.75,1.,d);",
    " float softProfile=(.65*exp(-28.*r2)+.18*exp(-5.*r2))*taper;",
    " float e=clamp(vExpression,0.,1.);",
    " float energy=(.008+e*uGeneGain*uBrightness)*softProfile;",
    " outColor=vec4(uColor,clamp(uOpacity*.3*energy,0.,1.));return;",
    "}",
    "float n=hash12(floor(uv*9.0)+vSeed*vec2(.071,.113));float facet=.96+.04*sin(atan(uv.y-.5,uv.x-.5)*5.0+vSeed*.17);float core=1.-smoothstep(.39*facet,.72*facet,d);float edge=pow(max(0.,1.-d),2.65);float glow=(.12+.78*uHalo)*edge*(uGlowMode>.5?1.28:1.);float fleck=step(.84,n)*(1.-smoothstep(.24,.96,d))*uGrain;float spec=pow(max(0.,1.-distance(uv,vec2(.32,.29))*3.8),9.0);float intensity=.22+uBrightness*(.82*core+.28*glow+.72*fleck+.42*spec);float alpha=uOpacity*(.95*core+.20*glow+.42*fleck+.30*spec);if(uExpression>.5){float e=clamp(vExpression,0.,1.);intensity=.015+uGeneGain*e*uBrightness*(.82*core+.28*glow+.72*fleck+.42*spec);alpha*=mix(.08,1.,sqrt(e));}vec3 lit=1.-exp(-uColor*intensity);float cue=mix(1.,mix(.28,1.18,smoothstep(-.65,.65,vDepth)),uDepthCue);outColor=vec4(lit*cue,alpha*mix(1.,cue,.6));}`);",
    "  const membraneProgram=glProgram(gl,common+`\\nin vec3 aPosition,aNormal;uniform vec4 uFrame;uniform float uZoom,uAspect,uFitScale;uniform vec3 uCenter;out vec3 vNormal;out float vDepth;void main(){vec3 p=uCenter+(aPosition-uCenter)*uFitScale;vec3 q=rotatePoint(p);float depth=1.35/uZoom-q.z;vNormal=normalize(rotatePoint(aNormal));vDepth=q.z;gl_Position=vec4(q.x/(depth*.72*uAspect)*1.75*min(1.,uAspect)*(uAspect<1.?.85:1.),q.y/(depth*.72)*1.75*min(1.,uAspect)*(uAspect<1.?.85:1.),q.z/depth,1.);gl_Position.xy=gl_Position.xy*uFrame.xy+uFrame.zw;}`,",
    "   `#version 300 es\\nprecision highp float;in vec3 vNormal;in float vDepth;uniform float uDepthCue;uniform vec3 uColor;uniform float uOpacity,uBrightness;out vec4 outColor;void main(){float facing=abs(dot(normalize(vNormal),vec3(0.,0.,1.)));float rim=pow(1.-facing,2.35);float fill=.035+.11*(1.-facing);float intensity=.20+uBrightness*(.45*fill+1.75*rim);float alpha=uOpacity*(fill+1.55*rim);if(alpha<.004)discard;vec3 n=normalize(vNormal);float near=mix(.38,1.15,smoothstep(-.65,.65,vDepth));float side=n.z < 0. ? .45 : 1.;float light=.7+.3*max(0.,dot(n,normalize(vec3(-.5,.7,1.))));float cue=mix(1.,near*side*light,uDepthCue);outColor=vec4((1.-exp(-uColor*intensity))*cue,alpha*mix(1.,cue,.65));}`);",
    "  const pointBuffers=[],membraneBuffers=[];",
    "  for(let group=0;group<labels.length;group++){",
    "   const pv=[],ids=members[group];for(const i of ids)pv.push(points[i][0],points[i][1],points[i][2],interior[i]?1:0,geneLevel(i));",
    "   const pb=gl.createBuffer(),vertices=new Float32Array(pv);gl.bindBuffer(gl.ARRAY_BUFFER,pb);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.DYNAMIC_DRAW);pointBuffers.push({buffer:pb,count:ids.length,color:hexRgb(colors[group]),ids,vertices});",
    "   const mesh=meshes.find(item=>item.group===group),mv=[];if(mesh)for(const tri of mesh.triangles)for(let k=0;k<3;k++)mv.push(...tri.p[k],...tri.normal);",
    "   const mb=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,mb);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(mv),gl.STATIC_DRAW);membraneBuffers.push({buffer:mb,count:mv.length/6,color:hexRgb(colors[group]),center:mesh?mesh.center:[0,0,0]});",
    "  }",
    "  glState={gl,pointProgram,membraneProgram,pointBuffers,membraneBuffers};draw=()=>{drawWebGL();drawAxisOverlay()};document.querySelector('#status').textContent+=` · GPU渲染已启用`;",
    " }catch(error){ctx=canvas.getContext('2d');document.querySelector('#status').textContent+=` · GPU不可用，使用兼容模式`;console.warn('WebGL unavailable; using 2D fallback',error)}",
    "}",
    "function glResize(gl){const box=canvas.getBoundingClientRect(),dpr=Math.min(devicePixelRatio||1,2),w=exportViewport?exportViewport.width:Math.max(1,Math.round(box.width*dpr)),h=exportViewport?exportViewport.height:Math.max(1,Math.round(box.height*dpr));if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;bloomCanvas.width=w;bloomCanvas.height=h;gl.viewport(0,0,w,h)}return exportViewport?exportViewport.dpr:dpr}",
    "function drawWebGL(){pending=false;const {gl,pointProgram,membraneProgram,pointBuffers,membraneBuffers}=glState,dpr=glResize(gl),rotation=viewMatrix(ax,ay,az),viewUniform=new Float32Array([rotation[0],rotation[3],rotation[6],rotation[1],rotation[4],rotation[7],rotation[2],rotation[5],rotation[8]]);gl.enable(gl.DEPTH_TEST);gl.depthMask(false);gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE);gl.clearColor(0,0,0,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);",
    " const aspect=exportViewport?exportViewport.reference.width/exportViewport.reference.height:canvas.width/Math.max(1,canvas.height),frame=exportViewport?exportViewport.transform:[1,1,0,0],fitScale=1.025-.05*Number(membraneFit.value);",
    " if(membrane.checked){gl.useProgram(membraneProgram);gl.uniform1f(gl.getUniformLocation(membraneProgram,'uDepthCue'),Number(depthCue.value));gl.uniformMatrix3fv(gl.getUniformLocation(membraneProgram,'uView'),false,viewUniform);gl.uniform1f(gl.getUniformLocation(membraneProgram,'uZoom'),zoom);gl.uniform1f(gl.getUniformLocation(membraneProgram,'uAspect'),aspect);gl.uniform4fv(gl.getUniformLocation(membraneProgram,'uFrame'),frame);gl.uniform1f(gl.getUniformLocation(membraneProgram,'uFitScale'),fitScale);gl.uniform1f(gl.getUniformLocation(membraneProgram,'uOpacity'),Number(membraneOpacity.value));gl.uniform1f(gl.getUniformLocation(membraneProgram,'uBrightness'),Number(membraneBrightness.value));const pos=gl.getAttribLocation(membraneProgram,'aPosition'),normal=gl.getAttribLocation(membraneProgram,'aNormal');for(const [group,item] of membraneBuffers.entries()){if(!enabled[group]||!item.count)continue;gl.bindBuffer(gl.ARRAY_BUFFER,item.buffer);gl.enableVertexAttribArray(pos);gl.vertexAttribPointer(pos,3,gl.FLOAT,false,24,0);gl.enableVertexAttribArray(normal);gl.vertexAttribPointer(normal,3,gl.FLOAT,false,24,12);gl.uniform3fv(gl.getUniformLocation(membraneProgram,'uColor'),item.color);gl.uniform3fv(gl.getUniformLocation(membraneProgram,'uCenter'),item.center);gl.drawArrays(gl.TRIANGLES,0,item.count)}}",
    " let drawPointPass=null;",
    " if(showPoints.checked){gl.useProgram(pointProgram);gl.uniform1f(gl.getUniformLocation(pointProgram,'uGeneGain'),Number(geneGain.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uGeneStyle'),geneStyle.value==='solid'?0:1);gl.uniform1f(gl.getUniformLocation(pointProgram,'uGeneSpread'),Number(geneSpread.value));gl.uniformMatrix3fv(gl.getUniformLocation(pointProgram,'uView'),false,viewUniform);gl.uniform1f(gl.getUniformLocation(pointProgram,'uZoom'),zoom);gl.uniform1f(gl.getUniformLocation(pointProgram,'uAspect'),aspect);gl.uniform4fv(gl.getUniformLocation(pointProgram,'uFrame'),frame);gl.uniform1f(gl.getUniformLocation(pointProgram,'uDpr'),dpr);gl.uniform1f(gl.getUniformLocation(pointProgram,'uPointFloor'),exportViewport?dpr/exportViewport.screenDpr:1);const expression=gl.getAttribLocation(pointProgram,'aExpression'),point=gl.getAttribLocation(pointProgram,'aPoint');drawPointPass=(isExpression)=>{gl.uniform1f(gl.getUniformLocation(pointProgram,'uExpression'),isExpression?1:0);gl.uniform1f(gl.getUniformLocation(pointProgram,'uDepthCue'),isExpression?0:Number(depthCue.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uHideInterior'),!isExpression&&membrane.checked&&!showInside.checked?1:0);gl.uniform1f(gl.getUniformLocation(pointProgram,'uPointSize'),Number(isExpression?geneSize.value:size.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uOpacity'),isExpression?genePulse*Number(geneOpacity.value):Number(pointOpacity.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uBrightness'),isExpression?1:Number(pointBrightness.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uHalo'),isExpression?0:Number(halo.value));gl.uniform1f(gl.getUniformLocation(pointProgram,'uGlowMode'),!isExpression&&mode.value==='glow'?1:0);gl.uniform1f(gl.getUniformLocation(pointProgram,'uGrain'),isExpression?.65:Number(pointGrain.value));for(const [group,item] of pointBuffers.entries()){if(!enabled[group]||!item.count)continue;gl.bindBuffer(gl.ARRAY_BUFFER,item.buffer);gl.enableVertexAttribArray(point);gl.vertexAttribPointer(point,4,gl.FLOAT,false,20,0);gl.enableVertexAttribArray(expression);gl.vertexAttribPointer(expression,1,gl.FLOAT,false,20,16);gl.uniform3fv(gl.getUniformLocation(pointProgram,'uColor'),isExpression?hexRgb(geneColor.value):item.color);gl.drawArrays(gl.POINTS,0,item.count)}};if(!geneActive()||geneBasePoints.checked)drawPointPass(false)}",
    " // Screen-space bloom is a cached bitmap pass: it adds halo without changing geometry.",
    " const bloom=bloomSettings();",
    " if(bloom.level>0){bloomCtx.clearRect(0,0,bloomCanvas.width,bloomCanvas.height);bloomCtx.drawImage(canvas,0,0);bloomCanvas.style.opacity=String(bloom.level);bloomCanvas.style.filter=`blur(${bloom.blur}px)`}else bloomCanvas.style.opacity='0';",
    " if(geneActive()&&drawPointPass)drawPointPass(true);",
    "}",
    "",
    "const defaultsKey='scToolbox.standalone.3d.defaults.v2';",
    "const preferenceIds=['embedding','geneGain','geneSize','geneOpacity','geneColor','geneStyle','geneBasePoints','geneSpread','geneBreath','geneBreathPeriod','halo','size','pointOpacity','pointBrightness','pointGrain','membrane','membraneFit','membraneOpacity','membraneBrightness','showPoints','showInside','autoRotate','rotateSpeed','showAxes','showGroupLabels','depthCue','recordMotion','recordFps','videoQuality','videoFormat','recordStop','recordTurns','recordDuration','imageRatio','imageQuality'];",
    "const defaultsStatus=document.querySelector('#defaultsStatus');",
    "function collectPreferences(){const controls={};for(const id of preferenceIds){const el=document.getElementById(id);controls[id]=el.type==='checkbox'?el.checked:(el.type==='range'||el.type==='number'?Number(el.value):el.value)}return {format:'scToolbox-3d-defaults',version:1,updatedAt:new Date().toISOString(),controls,view:{ax,ay,az,zoom},palette:Object.fromEntries([...Object.entries(paletteMemory),...labels.map((label,j)=>[label,colors[j]])])}}",
    "function validatePreferences(raw){",
    " if(!raw||raw.format!=='scToolbox-3d-defaults'||raw.version!==1||!raw.controls||typeof raw.controls!=='object'||!raw.view||!Number.isFinite(Date.parse(raw.updatedAt)))throw new Error('不是有效的三维默认参数文件');",
    " const controls={};for(const [id,v] of Object.entries(raw.controls)){",
    "  if(id==='mode'||id==='geneDiffraction')continue;const value=id==='geneStyle'&&v==='diffraction'?'soft':v;if(!preferenceIds.includes(id))throw new Error(`未知选项：${id}`);const el=document.getElementById(id);",
    "  if(el.type==='checkbox'){if(typeof value!=='boolean')throw new Error(`开关无效：${id}`)}",
    "  else if(el.type==='color'){if(typeof value!=='string'||!/^#[0-9a-f]{6}$/i.test(value))throw new Error(`颜色无效：${id}`)}",
    "  else if(el.tagName==='SELECT'){if(![...el.options].some(o=>o.value===String(value)))throw new Error(`选项无效：${id}`)}",
    "  else if(typeof value!=='number'||!Number.isFinite(value)||value<Number(el.min)||value>Number(el.max))throw new Error(`参数越界：${id}`);",
    "  controls[id]=value;",
    " }",
    " const {ax,ay,zoom}=raw.view,az=raw.view.az===undefined?0:raw.view.az;if(![ax,ay,az,zoom].every(v=>typeof v==='number'&&Number.isFinite(v))||zoom<.2||zoom>6)throw new Error('视角参数无效');",
    " const palette=raw.palette||{};if(typeof palette!=='object'||Array.isArray(palette)||Object.entries(palette).some(([label,color])=>!label||typeof color!=='string'||!/^#[0-9a-f]{6}$/i.test(color)))throw new Error('配色格式无效');",
    " return {format:raw.format,version:1,updatedAt:new Date(raw.updatedAt).toISOString(),controls,view:{ax,ay,az,zoom},palette:Object.fromEntries(Object.entries(palette))};",
    "}",
    "function setPreferenceControls(p){applyPalette(p.palette);for(const [id,value] of Object.entries(p.controls)){const el=document.getElementById(id);if(el.type==='checkbox')el.checked=value;else el.value=String(value)}syncGeneAppearance();syncGeneBreathing(true);if(embedding.selectedOptions[0]?.disabled)embedding.value=[...embedding.options].find(o=>!o.disabled).value}",
    "function applyPreferences(p,rebuild=true){",
    " stopAutoRotate();setPreferenceControls(p);",
    " if(rebuild&&activeEmbedding!==embedding.value)prepareEmbedding(embedding.value);",
    " setViewPose(p.view);zoom=p.view.zoom;rotationAnchor={ax,ay,az};document.querySelector('#bestViewStatus').textContent='';",
    " syncVisibility();syncRecordingOptions();schedule();",
    " if(autoRotate.checked){autoClock=performance.now();autoFrame=requestAnimationFrame(autoTick)}",
    "}",
    "const factoryPreferences=collectPreferences();",
    "const recommendedPreferences=validatePreferences(recommendedDefaults);",
    "document.querySelector('#recommendedDefaults').onclick=()=>{applyPreferences(recommendedPreferences);defaultsStatus.textContent='已应用推荐参数（2026-10-12保存的方案）'};",
    "let startupPreferences=null;",
    "try{if(embeddedDefaults)startupPreferences=validatePreferences(embeddedDefaults)}catch(e){defaultsStatus.textContent=`内置参数未加载：${e.message}`}",
    "try{const raw=(localStorage.getItem(defaultsKey)||localStorage.getItem('scToolbox.3d.defaults.v1')||localStorage.getItem('scToolbox.standalone.3d.defaults.v1'));if(raw){const saved=validatePreferences(JSON.parse(raw));if(!startupPreferences||saved.updatedAt>startupPreferences.updatedAt)startupPreferences=saved}}catch(e){defaultsStatus.textContent='浏览器本地保存不可用；仍可导出参数文件'}",
    "if(startupPreferences)setPreferenceControls(startupPreferences);",
    "document.querySelector('#saveDefaults').onclick=()=>{",
    " try{",
    "  const prefs=validatePreferences(collectPreferences());let localSaved=true;",
    "  try{localStorage.setItem(defaultsKey,JSON.stringify(prefs))}catch(e){localSaved=false}",
    "  const blob=new Blob([JSON.stringify(prefs,null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='UMAP3D_显示参数.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);",
    "  defaultsStatus.textContent=(localSaved?'浏览器默认参数已保存；':'')+'参数文件已请求下载；可在其他网页导入。';",
    " }catch(e){defaultsStatus.textContent=`保存失败：${e.message}`}",
    "};",
    "document.querySelector('#factoryDefaults').onclick=()=>{applyPreferences(factoryPreferences);defaultsStatus.textContent='当前页面已恢复初始参数；点击“保存为默认参数”可更新后续默认值'};",
    "document.querySelector('#importDefaults').onclick=()=>document.querySelector('#defaultsFile').click();",
    "document.querySelector('#defaultsFile').onchange=async event=>{try{const file=event.target.files[0];if(!file)return;if(file.size>100000)throw new Error('参数文件过大');const p=validatePreferences(JSON.parse(await file.text()));applyPreferences(p);defaultsStatus.textContent='参数已应用；点击“保存为默认参数”以用于后续网页'}catch(e){defaultsStatus.textContent=`导入失败：${e.message}`}finally{event.target.value=''}};",
    " document.querySelector('#save').onclick=savePng;document.querySelector('#record').onclick=recordRotation;prepareEmbedding(embedding.value);initWebGL();if(!glState){const cpuDraw=draw;draw=()=>{cpuDraw();drawAxisOverlay()}}new ResizeObserver(schedule).observe(canvas);if(geneData)selectGene();if(startupPreferences){applyPreferences(startupPreferences,false);defaultsStatus.textContent=\"已加载保存的默认参数\"}else schedule();",
    "if(headerChordSelectorData){",
    " const bar=document.querySelector('#headerChordSelector'),button=document.querySelector('#selectHeaderView'),backButton=document.querySelector('#headerChordBack'),selectorStatus=document.querySelector('#headerChordStatus');",
    " bar.hidden=false;",
    " document.querySelector('#headerChordProject').textContent=`项目 ${headerChordSelectorData.index}/${headerChordSelectorData.total}：${headerChordSelectorData.project}`;",
    " document.querySelector('#headerChordTarget').textContent=`表头群：${headerChordSelectorData.headerGroup}`;",
    " backButton.disabled=Number(headerChordSelectorData.index)<=1;",
    " backButton.onclick=()=>{",
    "  if(backButton.disabled)return;",
    "  backButton.disabled=true;button.disabled=true;selectorStatus.textContent='正在返回上一个项目……';",
    "  if(window.parent&&window.parent.Shiny&&typeof window.parent.Shiny.setInputValue==='function')window.parent.Shiny.setInputValue('header_chord_back',{index:Number(headerChordSelectorData.index),project:headerChordSelectorData.project},{priority:'event'});",
    "  else selectorStatus.textContent='当前页面未连接R选择器，无法返回。';",
    " };",
    " const initial=headerChordSelectorData.initialSelection;",
    " if(initial){",
    "  if(initial.palette)applyPalette(initial.palette);",
    "  if(initial.view&&[initial.view.ax,initial.view.ay,initial.view.zoom].every(Number.isFinite)){setViewPose(initial.view);zoom=initial.view.zoom;rotationAnchor={ax,ay,az}}",
    "  if(initial.embedding&&embeddingData[initial.embedding]){embedding.value=initial.embedding;prepareEmbedding(initial.embedding)}",
    "  selectorStatus.textContent='已载入上次保存的角度和配色，可直接修改后重新选定。';schedule();",
    " }",
    " button.onclick=()=>{",
    "  stopAutoRotate();autoRotate.checked=false;",
    "  const groupIndex=labels.indexOf(headerChordSelectorData.headerGroup);",
    "  if(groupIndex<0){selectorStatus.textContent='目标表头群不在当前标签中，无法保存。';return}",
    "  const selection={format:'scToolbox-header-chord-view',version:1,project:headerChordSelectorData.project,headerColumn:headerChordSelectorData.headerColumn,headerGroup:headerChordSelectorData.headerGroup,embedding:headerChordSelectorData.embedding,view:{ax,ay,az,zoom},color:colors[groupIndex],palette:Object.fromEntries(labels.map((label,j)=>[label,colors[j]])),selectedAt:new Date().toISOString()};",
    "  if(window.parent&&window.parent.Shiny&&typeof window.parent.Shiny.setInputValue==='function'){",
    "   button.disabled=true;selectorStatus.textContent='已选定，正在保存并进入下一项目……';",
    "   window.parent.Shiny.setInputValue('header_chord_selected',selection,{priority:'event'});",
    "  }else selectorStatus.textContent='当前页面未连接R选择器，请从实验性功能模块5打开本页。';",
    " };",
    "}",
    "</script></html>"
  ), collapse = "\n")
  resolve_expression_features <- function(features, available, assay, layer, notify = TRUE) {
    if (identical(features, 'all')) {
      if (!length(available)) stop('表达层没有基因：assay=', assay, '，layer=', layer, call. = FALSE)
      return(available)
    }
    index <- match(features, available); unresolved <- which(is.na(index))
    if (length(unresolved)) {
      folded <- tolower(available)
      for (i in unresolved) {
        candidates <- which(folded == tolower(features[i]))
        if (length(candidates) > 1L)
          stop('基因名大小写匹配不唯一：', features[i], ' -> ',
            paste(available[candidates], collapse = ', '),
            '。请填写完整名称（assay=', assay, '，layer=', layer, '）。', call. = FALSE)
        if (length(candidates) == 1L) index[i] <- candidates
      }
    }
    missing <- features[is.na(index)]
    if (length(missing))
      stop('表达层缺少基因：', paste(head(missing, 20L), collapse = ', '),
        if (length(missing) > 20L) ' 等' else '',
        '（assay=', assay, '，layer=', layer, '）。\n',
        '请核对该层的基因行名及assay；只看细胞类型可设gene_load="none"。', call. = FALSE)
    resolved <- available[index]; changed <- which(features != resolved)
    if (notify && length(changed))
      message('基因名称匹配：', paste(paste0(features[changed], ' -> ', resolved[changed]), collapse = ', '),
        '（assay=', assay, '，layer=', layer, '）。')
    unique(resolved)
  }
  expression_data <- function(d, features, assay = NULL, layer = 'data', validate_only = FALSE) {
    if (!length(features)) return(NULL)
    features <- unique(trimws(as.character(features)))
    if (anyNA(features) || any(!nzchar(features))) stop('基因名不能为空。')
    assay <- assay %||% SeuratObject::DefaultAssay(d$object)
    if (!assay %in% SeuratObject::Assays(d$object)) stop('不存在 assay：', assay)
    a <- d$object[[assay]]
    if (inherits(a, 'Assay5')) {
      available <- SeuratObject::Layers(a, search = NA)
      if (!layer %in% available) stop('找不到唯一表达层 ', layer, '。请先在对象副本 JoinLayers，或选择实际层。可用：', paste(available,collapse=', '))
      available_features <- SeuratObject::Features(a, layer = layer)
      features <- resolve_expression_features(features, available_features, assay, layer, !validate_only)
      if (validate_only) return(invisible(features))
      m <- SeuratObject::LayerData(a, layer = layer, features = features, cells = d$cells)
    } else {
      m <- if (utils::packageVersion('SeuratObject') >= '5.0.0')
        SeuratObject::GetAssayData(d$object, assay = assay, layer = layer) else
        SeuratObject::GetAssayData(d$object, assay = assay, slot = layer)
      features <- resolve_expression_features(features, rownames(m), assay, layer, !validate_only)
      if (validate_only) return(invisible(features))
      m <- m[features, , drop = FALSE]
    }
    if (is.null(m) || !length(features) || is.null(rownames(m)) || is.null(colnames(m)) || anyDuplicated(rownames(m)) ||
        anyDuplicated(colnames(m)) || !all(features %in% rownames(m)) || !all(d$cells %in% colnames(m)))
      stop('表达矩阵需要 genes × cells，并完整覆盖请求的基因和细胞。')
    m <- m[features, d$cells, drop = FALSE]
    if (as.double(nrow(m)) * ncol(m) > 2e6) {
      need('Matrix')
      mt <- methods::as(Matrix::t(m), 'dgCMatrix')
      if (any(!is.finite(mt@x)) || any(mt@x < 0)) stop('表达光效需要非负有限数值；请选择 data/counts，不要使用 scale.data。')
      maxima <- numeric(nrow(m))
      values <- lapply(seq_len(nrow(m)), function(i) {
        index <- if (mt@p[i] < mt@p[i+1L]) seq.int(mt@p[i]+1L, mt@p[i+1L]) else integer()
        maxima[i] <<- if (length(index)) max(mt@x[index]) else 0
        list(indices = as.list(unname(mt@i[index])), values = as.list(unname(mt@x[index])))
      })
      encoding <- 'sparse'; message('载入表达量：', nrow(m), '个基因 × ', ncol(m), '个细胞（稀疏存储）。')
    } else {
      m <- as.matrix(m)
      if (any(!is.finite(m)) || any(m < 0)) stop('表达光效需要非负有限数值；请选择 data/counts，不要使用 scale.data。')
      values <- unname(lapply(seq_len(nrow(m)), function(i) unname(m[i,])))
      maxima <- unname(apply(m, 1L, max)); encoding <- 'dense'
    }
    list(names = as.list(features), values = unname(values), maxima = as.list(maxima),
      layer = layer, assay = assay, encoding = encoding)
  }
  render <- function(obj, label_column = 'celltype', output_dir = 'UMAP3D',
      features = character(), assay = NULL, layer = 'data', colors = NULL,
      defaults_file = NULL, launch.browser = TRUE, overwrite = TRUE, methods = NULL) {
    d <- as_data(obj)
    expression_data(d, features, assay, layer, validate_only = TRUE)
    if (is.null(methods)) methods <- c('umap', 'tsne')[vapply(c('umap', 'tsne'), function(method)
      file.exists(file.path(output_dir, 'results', paste0(method, '3d_coordinates.rds'))), logical(1))]
    if (!length(methods)) stop('缺少3D坐标缓存，请先运行compute。')
    methods <- normalize_methods(methods)
    embeddings <- setNames(lapply(methods, function(method) load3d(output_dir, obj, method)), methods)
    cached <- embeddings[[1L]]; labels <- labels_for(d, label_column)
    levels <- sort(unique(labels)); cols <- palette(length(levels))
    if (!is.null(colors)) {
      if (is.null(names(colors)) || !all(levels %in% names(colors))) stop('colors 需要按群名命名，并覆盖全部群。')
      cols <- as.character(colors[levels]); invisible(grDevices::col2rgb(cols))
    }
    expression <- expression_data(d, features, assay, layer)
    payload <- list(xyz = unname((embeddings$umap %||% embeddings$tsne)$xyz),
      tsne = if (is.null(embeddings$tsne)) NULL else unname(embeddings$tsne$xyz),
      availableEmbeddings = as.list(methods),
      groups = as.list(unname(match(labels, levels)-1L)), labels = as.list(levels), colors = as.list(cols),
      stage = if (is.null(expression)) 'apply' else 'expression', expression = expression)
    path <- file.path(mkdir(file.path(output_dir, 'plots')), 'UMAP3D.html')
    path <- write_page(template3d, payload, path, overwrite, defaults_file)
    open_page(path, launch.browser); message('3D网页：', path)
    invisible(list(path = path, coordinate_file = cached$paths$coordinates,
      coordinate_files = lapply(embeddings, function(x) x$paths$coordinates), methods = methods,
      label_column = label_column, features = if (is.null(expression)) character() else unlist(expression$names)))
  }
  run <- function(obj, label_column = 'celltype', output_dir = 'UMAP3D', features = character(),
      assay = NULL, layer = 'data', colors = NULL, defaults_file = NULL, launch.browser = TRUE,
      overwrite = FALSE, methods = 'umap', ...) {
    compute3d(obj, output_dir, overwrite = overwrite, methods = methods, ...)
    render(obj, label_column, output_dir, features, assay, layer, colors, defaults_file, launch.browser,
      overwrite = TRUE, methods = methods)
  }
  execute <- function(input, options = list()) {
    if (!is.list(options) || is.data.frame(options) ||
        (length(options) && (is.null(names(options)) || anyNA(names(options)) ||
          any(!nzchar(names(options))) || anyDuplicated(names(options)))))
      stop('umap3d_options需要每项都有名称且名称唯一的list。', call. = FALSE)
    allowed <- c('action', 'gene_load', setdiff(union(names(formals(compute3d)), names(formals(render))), 'obj'))
    unknown <- setdiff(names(options), allowed)
    if (length(unknown)) stop('未知参数：', paste(unknown, collapse = ', '),
      '。请核对脚本开头的 umap3d_options。', call. = FALSE)
    action <- options$action %||% 'run'
    if (!is.character(action) || length(action) != 1L || is.na(action) ||
        !action %in% c('run', 'compute', 'render'))
      stop('action只能填写 run、compute 或 render。', call. = FALSE)
    # NULL使用默认值。
    arguments <- options[setdiff(names(options), c('action', 'gene_load'))]
    arguments <- arguments[!vapply(arguments, is.null, logical(1))]
    mode <- options$gene_load
    if (!is.null(mode)) {
      if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
          !mode %in% c('all', 'none', 'selected'))
        stop('gene_load 只能填写 all、none 或 selected。', call. = FALSE)
      if (action != 'compute') {
        if (mode == 'selected') {
          features <- arguments$features
          if (!is.character(features) || !length(features) || anyNA(features) ||
              any(!nzchar(trimws(features))) || identical(unique(trimws(features)), 'all'))
            stop('gene_load="selected"时，features必须填写非空基因名列表。', call. = FALSE)
          arguments$features <- unique(trimws(features))
        } else arguments$features <- if (mode == 'all') 'all' else character()
      }
    }
    output_dir <- arguments$output_dir %||% 'UMAP3D'
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
        !nzchar(trimws(output_dir))) stop('output_dir必须是一个非空目录路径。', call. = FALSE)
    obj <- resolve_input(input); d <- as_data(obj)
    if (action == 'run') {
      # 计算前检查标签列和基因。
      labels_for(d, arguments$label_column %||% 'celltype')
      expression_data(d, arguments$features %||% character(), arguments$assay,
        arguments$layer %||% 'data', validate_only = TRUE)
      return(do.call(run, c(list(obj = obj), arguments)))
    }
    fun <- if (action == 'compute') compute3d else render
    arguments <- arguments[intersect(names(arguments), names(formals(fun)))]
    arguments$output_dir <- output_dir
    if (action == 'render') arguments$overwrite <- TRUE
    do.call(fun, c(list(obj = obj), arguments))
  }
  demo <- function(output_dir = 'UMAP3D_demo', launch.browser = TRUE) {
    run(demo_data(), output_dir = output_dir, features = c('CD3D','MS4A1'),
      n_neighbors = 12L, n_epochs = 100L, verbose = FALSE, launch.browser = launch.browser, overwrite = TRUE)
  }
  list(execute = execute, inspect = inspect, compute = compute3d, render = render, run = run, demo = demo,
    example_data = demo_data, version = '1.7.4')

})

# -------------------------------------------------------------------------------
# 执行
# -------------------------------------------------------------------------------

umap3d_cli_usage <- function() {
  cat('
3D UMAP 通用绘图脚本
用法:
  Rscript 3dumap.R <rds> [基因1 基因2 ...] [选项]

位置参数:
  rds       必需。任意单细胞 RDS 对象（Seurat v3/v4/v5、SingleCellExperiment 等）
  基因      可选，0~N 个。空格 / 逗号 / 分号分隔；不传则只按细胞分组着色

选项:
  --label=列名     分组列，默认自动探测 CellType / annotation / cluster 等
  --out=目录       输出目录，默认 <rds所在目录>/<rds文件名>_UMAP3D
  --assay=名称     表达 assay，默认自动选择含目标基因的 data/counts 层
  --layer=层名     表达层，默认自动
  --methods=umap   umap | tsne | umap,tsne
  --force          忽略已有坐标缓存，强制重算 3D 坐标
  --browser        生成后尝试打开浏览器（服务器默认不开）
  --help           显示本帮助

示例:
  Rscript 3dumap.R /data/obj.rds
  Rscript 3dumap.R /data/obj.rds WRAP73 SEC22B YKT6
  Rscript 3dumap.R /data/obj.rds "WRAP73,SEC22B,YKT6" --label=CellType

输出:
  <out>/plots/UMAP3D.html   交互式 3D 网页
  <out>/results/           3D 坐标 CSV / RDS 缓存 / 参数 JSON
')
}

umap3d_cli_die <- function(...) stop(paste0(...), call. = FALSE)
umap3d_cli_say <- function(...) message('[3dumap] ', ...)
umap3d_cli_need <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) umap3d_cli_die('缺少 R 包：', p, '。请先安装。')
  invisible(TRUE)
}

umap3d_cli_parse <- function(argv) {
  opts <- list(); positional <- character(); i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c('--help', '-h')) { umap3d_cli_usage(); quit(save = 'no', status = 0L) }
    if (identical(a, '--force')) { opts$force <- TRUE; i <- i + 1L; next }
    if (identical(a, '--browser')) { opts$browser <- TRUE; i <- i + 1L; next }
    if (startsWith(a, '--')) {
      kv <- sub('^--', '', a)
      if (grepl('=', kv, fixed = TRUE)) {
        key <- sub('=.*$', '', kv); value <- sub('^[^=]*=', '', kv)
      } else {
        key <- kv
        if (i == length(argv)) umap3d_cli_die('选项 --', key, ' 缺少取值。')
        i <- i + 1L; value <- argv[i]
      }
      key <- tolower(key)
      if (!key %in% c('label', 'out', 'assay', 'layer', 'methods', 'genes'))
        umap3d_cli_die('未知选项：--', key, '；用 --help 查看用法。')
      opts[[key]] <- value; i <- i + 1L; next
    }
    positional <- c(positional, a); i <- i + 1L
  }
  if (!length(positional)) {
    umap3d_cli_usage()
    umap3d_cli_die('必须提供 RDS 路径作为第一个参数。')
  }
  opts$rds <- positional[1L]
  tokens <- c(positional[-1L], if (is.null(opts$genes)) character() else opts$genes)
  genes <- character()
  if (length(tokens)) {
    tokens <- gsub('和', ' ', tokens)
    genes <- unique(unlist(strsplit(tokens, '[,;\\s]+', perl = TRUE)))
    genes <- trimws(genes); genes <- genes[nzchar(genes)]
  }
  list(opts = opts, genes = genes)
}

umap3d_cli_to_seurat <- function(sce) {
  umap3d_cli_need('Seurat'); umap3d_cli_need('SummarizedExperiment')
  an <- SummarizedExperiment::assayNames(sce)
  if (!length(an)) umap3d_cli_die('SingleCellExperiment 对象没有任何 assay。')
  prefer <- c('logcounts', 'normcounts', 'counts', 'data')
  pick <- prefer[prefer %in% an]
  pick <- if (length(pick)) pick[1L] else an[1L]
  m <- SummarizedExperiment::assay(sce, pick)
  m <- tryCatch(methods::as(m, 'dgCMatrix'),
    error = function(e) Matrix::Matrix(as.matrix(m), sparse = TRUE))
  obj <- Seurat::CreateSeuratObject(counts = m, assay = 'RNA', min.cells = 0, min.features = 0)
  if (pick != 'counts') {
    obj <- tryCatch({
      SeuratObject::LayerData(obj, assay = 'RNA', layer = 'data') <- m
      obj
    }, error = function(e) obj)
  }
  cd <- tryCatch(as.data.frame(SummarizedExperiment::colData(sce)), error = function(e) NULL)
  if (!is.null(cd) && nrow(cd) == ncol(obj)) {
    rownames(cd) <- colnames(obj)
    keep <- vapply(cd, function(v) is.atomic(v) && is.null(dim(v)), logical(1))
    cd <- cd[, keep, drop = FALSE]
    cd <- cd[, setdiff(names(cd), names(obj[[]])), drop = FALSE]
    if (ncol(cd)) obj <- Seurat::AddMetaData(obj, cd)
  }
  if (requireNamespace('SingleCellExperiment', quietly = TRUE)) {
    for (r in SingleCellExperiment::reducedDimNames(sce)) {
      e <- tryCatch(as.matrix(SingleCellExperiment::reducedDim(sce, r)), error = function(e2) NULL)
      if (is.null(e) || !is.numeric(e) || nrow(e) != ncol(obj) || ncol(e) < 2L) next
      rownames(e) <- colnames(obj)
      obj <- tryCatch({
        obj[[r]] <- SeuratObject::CreateDimReducObject(embeddings = e, key = 'PC_', assay = 'RNA')
        obj
      }, error = function(e2) obj)
    }
  }
  umap3d_cli_say('SingleCellExperiment 已转换为 Seurat（assay=RNA，来源层=', pick, '）。')
  obj
}

umap3d_cli_generic <- function(x) {
  umap3d_cli_need('Seurat'); umap3d_cli_need('Matrix')
  nms <- names(x)
  find <- function(keys) {
    if (is.null(nms)) return(NULL)
    for (k in keys) {
      hit <- which(tolower(nms) == tolower(k))
      if (length(hit)) return(x[[hit[1L]]])
    }
    for (k in keys) {
      hit <- grep(tolower(k), tolower(nms), fixed = TRUE)
      if (length(hit)) return(x[[hit[1L]]])
    }
    NULL
  }
  emb <- find(c('embedding', 'embeddings', 'coordinates', 'coords', 'xyz', 'umap_3d', 'umap', 'pca', 'harmony'))
  if (is.data.frame(emb)) emb <- as.matrix(emb)
  if (is.null(emb) || !is.matrix(emb) || !is.numeric(emb) || ncol(emb) < 2L)
    umap3d_cli_die('list/data.frame 型 RDS 需要含数值坐标矩阵（字段 embedding/coordinates/umap 等，行为细胞）。')
  cells <- rownames(emb)
  if (is.null(cells) || anyDuplicated(cells) || any(!nzchar(cells)))
    cells <- sprintf('cell_%06d', seq_len(nrow(emb)))
  md <- find(c('metadata', 'meta.data', 'meta', 'coldata', 'colData'))
  md <- tryCatch(as.data.frame(md), error = function(e) data.frame(row.names = cells))
  if (nrow(md) != length(cells)) md <- data.frame(row.names = cells)
  rownames(md) <- cells
  keep <- vapply(md, function(v) is.atomic(v) && is.null(dim(v)), logical(1))
  md <- md[, keep, drop = FALSE]
  counts <- Matrix::Matrix(0, nrow = 1L, ncol = length(cells), sparse = TRUE,
    dimnames = list('placeholder', cells))
  args <- list(counts = counts, assay = 'RNA', min.cells = 0, min.features = 0)
  if (ncol(md)) args$meta.data <- md
  obj <- do.call(Seurat::CreateSeuratObject, args)
  rownames(emb) <- cells
  obj <- tryCatch({
    obj[['pca']] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = 'PC_', assay = 'RNA')
    obj
  }, error = function(e) obj)
  umap3d_cli_say('已从通用坐标矩阵构造 Seurat 对象（无表达矩阵，仅分组着色）。')
  obj
}

umap3d_cli_load <- function(path) {
  if (!file.exists(path) || dir.exists(path)) umap3d_cli_die('RDS 文件不存在：', path)
  obj <- tryCatch(readRDS(normalizePath(path, winslash = '/', mustWork = TRUE)),
    error = function(e) umap3d_cli_die('读取 RDS 失败：', conditionMessage(e)))
  if (inherits(obj, 'Seurat')) return(list(obj = obj, expressionless = FALSE))
  if (inherits(obj, 'SingleCellExperiment') || inherits(obj, 'SummarizedExperiment'))
    return(list(obj = umap3d_cli_to_seurat(obj), expressionless = FALSE))
  if (is.list(obj) || is.data.frame(obj))
    return(tryCatch(list(obj = umap3d_cli_generic(obj), expressionless = TRUE),
      error = function(e) umap3d_cli_die('无法解析该 RDS 对象：', conditionMessage(e))))
  umap3d_cli_die('不支持的 RDS 对象：', paste(class(obj), collapse = ', '),
    '。支持 Seurat(v3/v4/v5)、SingleCellExperiment/SummarizedExperiment、含坐标矩阵的 list/data.frame。')
}

umap3d_cli_join_layers <- function(obj) {
  for (a in SeuratObject::Assays(obj)) {
    x <- obj[[a]]
    if (!inherits(x, 'Assay5')) next
    lys <- tryCatch(SeuratObject::Layers(x, search = NA), error = function(e) character())
    if (anyDuplicated(lys) || any(grepl('^(counts|data|scale\\.data)\\.', lys))) {
      umap3d_cli_say('检测到分层的 Assay5（', a, '），合并表达层 JoinLayers。')
      obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = a), error = function(e) {
        message('[3dumap] JoinLayers 失败：', conditionMessage(e)); obj })
    }
  }
  obj
}

umap3d_cli_label <- function(obj, requested = NULL) {
  md <- obj[[]]; nms <- names(md)
  if (!is.null(requested)) {
    hit <- nms[tolower(nms) == tolower(requested)]
    if (length(hit)) return(list(obj = obj, column = hit[1L]))
    umap3d_cli_die('指定的分组列不存在：', requested, '。可用列：', paste(nms, collapse = ', '))
  }
  norm <- function(s) tolower(gsub('[^A-Za-z0-9]', '', s))
  score <- vapply(nms, function(n) {
    v <- md[[n]]
    if (!is.atomic(v) || !is.null(dim(v))) return(-Inf)
    u <- length(unique(as.character(v[!is.na(v)])))
    if (u < 2L || u > 100L) return(-Inf)
    k <- norm(n); s <- 0
    if (k %in% c('celltype', 'celltypes', 'celltypeannotation', 'cellannotation',
                 'annotation', 'annotations', 'cellclass', 'celltypefinal', 'majorcelltype'))
      s <- s + 1000
    if (grepl('celltype|annotation|cellclass|celllabel', k)) s <- s + 300
    if (grepl('cluster|subtype', k)) s <- s + 100
    if (grepl('phase|qc|batch|sample|ident|status|score|class|remove', k)) s <- s - 300
    s - abs(u - 12L) / 10
  }, numeric(1))
  if (any(is.finite(score))) {
    col <- nms[which.max(score)]
    umap3d_cli_say('自动探测分组列：', col)
    return(list(obj = obj, column = col))
  }
  umap3d_cli_say('未找到合适的分组列，使用单一分组 All。')
  obj <- Seurat::AddMetaData(obj,
    data.frame(CellType = rep('All', ncol(obj)), row.names = colnames(obj)))
  list(obj = obj, column = 'CellType')
}

umap3d_cli_layers <- function(obj) {
  out <- list()
  for (a in SeuratObject::Assays(obj)) {
    x <- obj[[a]]
    if (inherits(x, 'Assay5')) {
      lys <- tryCatch(unique(SeuratObject::Layers(x, search = NA)), error = function(e) character())
      lys <- lys[lys %in% c('data', 'counts', 'scale.data') | !grepl('\\.', lys)]
      for (l in lys) {
        fs <- tryCatch(as.character(SeuratObject::Features(x, layer = l)), error = function(e) character())
        if (length(fs)) out[[length(out) + 1L]] <- list(assay = a, layer = l, features = fs)
      }
    } else {
      for (l in c('data', 'counts', 'scale.data')) {
        m <- tryCatch(SeuratObject::GetAssayData(obj, assay = a, layer = l),
          error = function(e) tryCatch(SeuratObject::GetAssayData(obj, assay = a, slot = l),
            error = function(e2) NULL))
        fs <- if (is.null(m)) character() else rownames(m)
        if (length(fs)) out[[length(out) + 1L]] <- list(assay = a, layer = l, features = fs)
      }
    }
  }
  out
}

umap3d_cli_match_genes <- function(genes, features) {
  idx <- match(genes, features)
  miss <- which(is.na(idx))
  if (length(miss)) {
    folded <- tolower(features)
    for (i in miss) {
      hit <- which(folded == tolower(genes[i]))
      if (length(hit)) idx[i] <- hit[1L]
    }
  }
  ok <- !is.na(idx)
  list(features = unique(features[idx[ok]]), missing = genes[!ok])
}

umap3d_cli_expression <- function(obj, genes, assay_opt = NULL, layer_opt = NULL) {
  layers <- umap3d_cli_layers(obj)
  if (!length(layers)) {
    if (length(genes)) umap3d_cli_die('对象中没有可用的表达层，无法绘制基因表达。')
    return(list(assay = NULL, layer = NULL, features = character(), missing = character()))
  }
  default_assay <- SeuratObject::DefaultAssay(obj)
  rank <- vapply(layers, function(z) {
    lp <- if (z$layer == 'data') 0L else if (z$layer == 'counts') 1L else if (z$layer == 'scale.data') 2L else 3L
    ap <- if (identical(z$assay, default_assay)) 0L else 1L
    lp * 10L + ap
  }, integer(1))
  layers <- layers[order(rank)]
  if (!is.null(assay_opt) || !is.null(layer_opt)) {
    layers <- Filter(function(z)
      (is.null(assay_opt) || tolower(z$assay) == tolower(assay_opt)) &&
      (is.null(layer_opt) || tolower(z$layer) == tolower(layer_opt)), layers)
    if (!length(layers)) umap3d_cli_die('指定的 assay/layer 不存在。')
  }
  if (!length(genes))
    return(list(assay = layers[[1L]]$assay, layer = layers[[1L]]$layer,
      features = character(), missing = character()))
  hits <- lapply(layers, function(z) umap3d_cli_match_genes(genes, z$features))
  counts <- vapply(hits, function(h) length(h$features), integer(1))
  pick <- which.max(counts)
  z <- layers[[pick]]; h <- hits[[pick]]
  if (length(h$missing)) {
    umap3d_cli_say('以下基因在 ', z$assay, '/', z$layer, ' 中不存在，已跳过：',
      paste(h$missing, collapse = ', '))
    for (g in head(h$missing, 5L)) {
      near <- head(grep(substr(g, 1L, 3L), z$features, value = TRUE, ignore.case = TRUE), 5L)
      if (length(near)) umap3d_cli_say('  ', g, ' 近似基因：', paste(near, collapse = ', '))
    }
  }
  if (!length(h$features))
    umap3d_cli_die('所有指定基因都不在表达层中。检查基因名/物种/ID 类型；只看分组可不传基因。')
  umap3d_cli_say('表达来源：', z$assay, '/', z$layer, '；载入基因：',
    paste(h$features, collapse = ', '))
  list(assay = z$assay, layer = z$layer, features = h$features, missing = h$missing)
}

umap3d_cli_reduction <- function(obj) {
  reds <- SeuratObject::Reductions(obj)
  if (!length(reds)) return(NULL)
  nd <- vapply(reds, function(r)
    tryCatch(ncol(SeuratObject::Embeddings(obj, r)), error = function(e) 0L), integer(1))
  keys <- tolower(gsub('[^A-Za-z0-9]', '', reds))
  vis <- grepl('umap|tsne|visual|2d', keys)
  score <- ifelse(keys == 'harmony', 100,
           ifelse(keys == 'pca', 90,
           ifelse(grepl('harmony|integrated|mnn|scvi|fastmnn', keys), 80,
           ifelse(grepl('pca|lsi', keys), 70, 0))))
  score[vis] <- score[vis] - 50
  score[nd < 3L] <- score[nd < 3L] - 30
  if (any(score > 0L)) {
    r <- reds[which.max(score)]
    if (nd[match(r, reds)] >= 2L) return(r)
  }
  usable <- reds[nd >= 2L]
  if (length(usable)) return(usable[1L])
  NULL
}

umap3d_cli_pca <- function(obj, assay, layer, npcs = 30L) {
  umap3d_cli_need('Matrix'); umap3d_cli_need('irlba')
  m <- tryCatch(SeuratObject::LayerData(obj[[assay]], layer = layer), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m <- tryCatch(methods::as(m, 'dgCMatrix'), error = function(e) Matrix::Matrix(as.matrix(m), sparse = TRUE))
  m <- m[Matrix::rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) < 5L) return(NULL)
  if (nrow(m) > 2000L) {
    mu <- Matrix::rowMeans(m)
    v <- Matrix::rowMeans(m ^ 2) - mu ^ 2
    m <- m[order(v, decreasing = TRUE)[seq_len(2000L)], , drop = FALSE]
  }
  k <- min(npcs, nrow(m) - 1L, ncol(m) - 1L)
  if (k < 2L) return(NULL)
  X <- Matrix::t(m)
  fit <- irlba::irlba(X, nv = k, center = Matrix::colMeans(X))
  emb <- as.matrix(fit$u %*% diag(fit$d))
  rownames(emb) <- colnames(obj)
  colnames(emb) <- paste0('PC_', seq_len(ncol(emb)))
  obj <- tryCatch({
    obj[['pca']] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = 'PC_', assay = assay)
    obj
  }, error = function(e) {
    message('[3dumap] 注入 PCA 失败：', conditionMessage(e)); NULL
  })
  obj
}

umap3d_cli_main <- function(argv) {
  parsed <- umap3d_cli_parse(argv)
  opts <- parsed$opts; genes <- parsed$genes
  umap3d_cli_need('SeuratObject'); umap3d_cli_need('jsonlite'); umap3d_cli_need('Matrix')
  loaded <- umap3d_cli_load(opts$rds)
  obj <- loaded$obj
  if (!inherits(obj, 'Seurat')) umap3d_cli_die('对象转换后不是 Seurat 对象。')
  obj <- umap3d_cli_join_layers(obj)
  lab <- umap3d_cli_label(obj, opts$label)
  obj <- lab$obj
  if (loaded$expressionless && length(genes)) {
    umap3d_cli_say('该 RDS 无表达矩阵，忽略基因参数。')
    genes <- character()
  }
  ex <- umap3d_cli_expression(obj, genes, opts$assay, opts$layer)
  src <- umap3d_cli_reduction(obj)
  if (is.null(src)) {
    umap3d_cli_say('对象缺少可用降维，尝试从表达矩阵计算 PCA ...')
    obj <- tryCatch(umap3d_cli_pca(obj, ex$assay, ex$layer),
      error = function(e) { message('[3dumap] PCA 计算失败：', conditionMessage(e)); NULL })
    if (is.null(obj))
      umap3d_cli_die('无法获得 3D UMAP 的输入降维。请先在对象上运行 RunPCA 或 RunHarmony 后重新保存 RDS。')
    src <- 'pca'
  }
  umap3d_cli_say('3D 坐标输入降维：', src)
  rds_abs <- normalizePath(opts$rds, winslash = '/', mustWork = TRUE)
  out <- if (is.null(opts$out))
    file.path(dirname(rds_abs), paste0(tools::file_path_sans_ext(basename(rds_abs)), '_UMAP3D'))
  else opts$out
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE, showWarnings = FALSE))
    umap3d_cli_die('无法创建输出目录：', out)
  methods <- if (is.null(opts$methods)) 'umap' else
    unique(trimws(unlist(strsplit(opts$methods, '[,;\\s]+', perl = TRUE))))
  methods <- methods[nzchar(methods)]
  if (!length(methods) || any(!methods %in% c('umap', 'tsne')))
    umap3d_cli_die('--methods 只能是 umap、tsne 或 umap,tsne。')
  options <- list(
    action = 'run',
    gene_load = if (length(ex$features)) 'selected' else 'none',
    features = if (length(ex$features)) ex$features else NULL,
    label_column = lab$column,
    output_dir = out,
    assay = ex$assay,
    layer = ex$layer,
    source_reduction = src,
    methods = methods,
    launch.browser = isTRUE(opts$browser),
    overwrite = isTRUE(opts$force),
    reuse_embedding = TRUE,
    verbose = FALSE
  )
  umap3d_cli_say('输入对象：', rds_abs)
  umap3d_cli_say('细胞数：', ncol(obj), '；分组列：', lab$column,
    '；基因：', if (length(ex$features)) paste(ex$features, collapse = ', ') else '无')
  umap3d_cli_say('输出目录：', out)
  res <- tryCatch(umap3d$execute(obj, options), error = function(e) {
    msg <- conditionMessage(e)
    if (!grepl('缓存|cache|signature|已存在|不匹配', msg)) stop(e)
    umap3d_cli_say('已有 3D 坐标缓存与当前输入不匹配，自动重算（手动强制重算用 --force）。')
    options$overwrite <- TRUE
    umap3d$execute(obj, options)
  })
  umap3d_cli_say('完成，3D 网页：', res$path)
  invisible(res)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  umap3d_cli_usage()
  if (!interactive() && sys.nframe() <= 1L) quit(save = 'no', status = 0L)
} else {
  umap3d_cli_main(args)
}
