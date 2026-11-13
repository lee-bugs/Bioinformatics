#!/usr/bin/env Rscript
# Adapted from hdWGCNA.R for the WRAP73 glioma Seurat object.
# Key changes vs. the original template:
#   1. No Assay5 -> Assay coercion and NO CheckSeurat5 override.
#      SeuratObject >= 5.0 removed GetAssayData(slot=); hdWGCNA 0.4.12 runs
#      natively on Seurat 5 through LayerData(layer=).
#   2. All hard-coded `celltype` / `sample` references are parameterised
#      (this object uses CellType / orig.ident, 22 samples).
#   3. Assay defaults to `decontX` because `RNA` has no `data` layer here.
#   4. Extra CSV exports (hub genes, MEs, module scores, session info).
#   5. Per-figure tryCatch: a plotting failure never kills the whole run.
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Rcpp); library(ggplot2)
  library(ggrepel); library(igraph); library(hdWGCNA); library(tidyverse)
  library(cowplot); library(patchwork); library(optparse)
})
# ---------------------------------------------------------------------------
# WGCNA 1.74 bugfix (hdWGCNA 0.4.12 always calls blockwiseConsensusModules,
# even with consensus = FALSE). Its gene-reassignment step does
#   apply(scoreAll[modGenes[candidates], ], 1, which.max)
# WITHOUT drop = FALSE. When exactly one gene is a candidate the subset
# degenerates to a vector and apply() aborts with
# "dim(X) must have a positive length". WGCNA::blockwiseModules (the
# single-set twin) already has drop = FALSE. Patch the namespace copy.
# ---------------------------------------------------------------------------
patch_wgcna_consensus <- function() {
  ns <- asNamespace('WGCNA')
  txt <- deparse(get('blockwiseConsensusModules', envir = ns))
  i <- grep('whichBest = apply(scoreAll', txt, fixed = TRUE)
  if (length(i) != 1L) {
    warning('WGCNA patch target not found; continuing unpatched')
    return(invisible(FALSE))
  }
  if (!grepl('which.max', txt[i + 1L], fixed = TRUE))
    stop('unexpected WGCNA source layout; refusing to patch')
  txt[i] <- '                whichBest = apply(scoreAll[modGenes[candidates], , drop = FALSE], 1, which.max)'
  txt[i + 1L] <- ''
  f <- eval(parse(text = txt), envir = ns)
  assignInNamespace('blockwiseConsensusModules', f, ns = 'WGCNA')
  invisible(TRUE)
}
if (patch_wgcna_consensus()) cat('[WGCNA patch] blockwiseConsensusModules: drop = FALSE applied\n')
theme_set(theme_cowplot())
set.seed(12345)

option_list <- list(
  make_option(c('-r','--rds'), help='input Seurat RDS'),
  make_option(c('-i','--idents'), default='CellType', help='cell type column'),
  make_option(c('-g','--groupby'), default='orig.ident', help='sample column for metacells / ME harmonisation'),
  make_option(c('-u','--use_cluster'), default='all', help='comma separated cell types, default all'),
  make_option(c('-k','--kvalue'), default='25', help='k for metacells'),
  make_option(c('-m','--min_cells'), default='100', help='min cells per group'),
  make_option(c('--assays'), default='decontX', help='assay holding counts+data'),
  make_option(c('-s','--soft_power'), default='auto', help='soft power or auto'),
  make_option(c('-o','--outdir'), default='./result', help='output directory'),
  make_option(c('--threads'), type='integer', default=8, help='WGCNA threads'),
  make_option(c('--smoke'), action='store_true', default=FALSE, help='small subset, full pipeline test'),
  make_option(c('--smoke_cells'), type='integer', default=800, help='max cells per type in smoke mode'),
  make_option(c('--smoke_genes'), type='integer', default=300, help='n variable genes in smoke mode'),
  make_option(c('--skip_power_sweep'), action='store_true', default=FALSE, help='skip TestSoftPowers'),
  make_option(c('--resume'), action='store_true', default=FALSE, help='resume from _checkpoint_power.rds'),
  make_option(c('--no_checkpoint'), action='store_true', default=FALSE, help='do not write checkpoint files')
)
opt <- parse_args(OptionParser(option_list = option_list))

rds <- opt$rds; 
idents <- opt$idents; 
groupby <- opt$groupby
use_cluster <- opt$use_cluster; 
kvalue <- as.numeric(opt$kvalue)
assays <- opt$assays; 
min_cells <- as.numeric(opt$min_cells)
soft_power <- opt$soft_power; 
outdir <- opt$outdir

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# CRITICAL: all downstream paths (including the log file) must stay valid
# after ConstructNetwork() temporarily changes the working directory.
outdir <- normalizePath(outdir, mustWork = TRUE)
logf <- file.path(outdir, 'run.log')
step <- function(...) {
  msg <- sprintf('[%s] %s', format(Sys.time(), '%H:%M:%S'), paste0(..., collapse = ' '))
  cat(msg, '\n')
  tryCatch(cat(msg, '\n', file = logf, append = TRUE),
           error = function(e) warning('log write failed: ', conditionMessage(e)))
  flush.console()
}
safe <- function(label, expr) {
  step('START ', label)
  t0 <- Sys.time()
  out <- tryCatch(withCallingHandlers(force(expr), error = function(e) {
    cc <- tail(sys.calls(), 10)
    step('  TRACE ', paste(vapply(cc, function(z) substr(paste(deparse(z), collapse=' '), 1, 150), ''), collapse=' >> '))
  }), error = function(e) {
    step('  FAILED ', label, ' :: ', conditionMessage(e)); NULL
  })
  step('  DONE  ', label, sprintf('(%.1fs)', as.numeric(difftime(Sys.time(), t0, units='secs'))))
  out
}
cat('', file = logf)
step('=== hdWGCNA run ===')
step('rds      =', rds)
step('idents   =', idents, '| groupby =', groupby, '| assay =', assays)
step('outdir   =', outdir, '| smoke =', opt$smoke)
step('versions: R', R.version.string, '| Seurat', as.character(packageVersion('Seurat')),
     '| hdWGCNA', as.character(packageVersion('hdWGCNA')))
enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------- checkpoint ----
# Preparation (readRDS -> MetacellsByGroups -> TestSoftPowers) costs ~20 min
# on this 75k-cell object. Cache it right before network construction so a
# later failure can be retried with --resume in seconds.
cp_file <- file.path(outdir, '_checkpoint_power.rds')
pick_power <- function(obj) {
  pt <- tryCatch(GetPowerTable(obj), error = function(e) NULL)
  if (is.null(pt) || !nrow(pt) || !all(c('Power','SFT.R.sq') %in% colnames(pt))) return(NULL)
  pw <- suppressWarnings(as.numeric(pt$Power)); fit <- suppressWarnings(as.numeric(pt$SFT.R.sq))
  ok <- is.finite(pw) & is.finite(fit) & fit >= 0.8 & pw > 3
  if (any(ok)) return(min(pw[ok]))
  ok2 <- is.finite(pw) & is.finite(fit)
  if (any(ok2)) return(pw[ok2][which.max(fit[ok2])])
  NULL
}
sig <- paste(basename(as.character(rds)), idents, groupby, assays, kvalue, min_cells,
             opt$smoke, opt$smoke_genes, opt$skip_power_sweep, soft_power, sep='|')
resumed <- FALSE
if (opt$resume && file.exists(cp_file)) {
  ck <- tryCatch(readRDS(cp_file), error = function(e) NULL)
  if (!is.null(ck) && identical(ck$sig, sig)) {
    step('RESUME from checkpoint:', cp_file)
    seurat_obj <- ck$obj; use_cluster <- ck$use_cluster; new_name <- ck$new_name
    sp <- ck$sp; resumed <- TRUE
    step('checkpoint:', ncol(seurat_obj), 'cells |', length(use_cluster),
         'cell types | soft_power =', sp)
  } else step('checkpoint signature mismatch -> full run')
}

if (!resumed) {
# ---------------------------------------------------------------- load ----
seurat_obj <- safe('readRDS', readRDS(rds))
if (is.null(seurat_obj)) stop('readRDS failed')
step('object:', ncol(seurat_obj), 'cells x', nrow(seurat_obj), 'features |',
     'DefaultAssay =', DefaultAssay(seurat_obj))
if (!(idents %in% colnames(seurat_obj@meta.data)))
  stop(sprintf('idents column "%s" not found. Available: %s', idents, paste(colnames(seurat_obj@meta.data), collapse=', ')))
if (!(groupby %in% colnames(seurat_obj@meta.data)))
  stop(sprintf('groupby column "%s" not found. Available: %s', groupby, paste(colnames(seurat_obj@meta.data), collapse=', ')))
if (!(assays %in% names(seurat_obj@assays)))
  stop(sprintf('assay "%s" not found. Available: %s', assays, paste(names(seurat_obj@assays), collapse=', ')))
step('reductions:', paste(names(seurat_obj@reductions), collapse=', '))
step('layers of', assays, ':', paste(Layers(seurat_obj[[assays]]), collapse=', '))

DefaultAssay(seurat_obj) <- assays
# hdWGCNA needs a normalised `data` layer for SetDatExpr / ModuleExprScore.
d <- tryCatch(dim(SeuratObject::LayerData(seurat_obj, assay=assays, layer='data')), error=function(e) c(0L,0L))
if (prod(d) == 0) {
  step('no data layer in', assays, '-> NormalizeData')
  seurat_obj <- safe('NormalizeData', NormalizeData(seurat_obj, assay=assays, verbose=FALSE))
}

# --------------------------------------------------------------- subset ----
if (opt$smoke) {
  ct <- as.character(seurat_obj@meta.data[[idents]])
  keep_types <- names(sort(table(ct), decreasing=TRUE))[1:min(2, length(unique(ct)))]
  pick <- unlist(lapply(keep_types, function(x) {
    ix <- which(ct == x); sample(ix, min(length(ix), opt$smoke_cells))
  }), use.names=FALSE)
  step('SMOKE subset:', length(pick), 'cells from', paste(keep_types, collapse=', '))
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[pick])
  gc()
}

# ---------------------------------------------------------------- umap ----
safe('DimPlot', {
  p <- DimPlot(seurat_obj, group.by=idents, reduction='umap', label=FALSE)
  ggsave(file.path(outdir,'umap.pdf'), p, w=10, h=8)
  ggsave(file.path(outdir,'umap.png'), p, w=10, h=8, bg='white')
})

if (use_cluster == 'all') {
  use_cluster <- sort(unique(as.character(seurat_obj@meta.data[[idents]])))
  new_name <- 'M'
} else {
  use_cluster <- unlist(strsplit(use_cluster, split=','))
  new_name <- 'sub-M'
  seurat_obj <- seurat_obj[, seurat_obj@meta.data[[idents]] %in% use_cluster]
}
step('requested cell types:', paste(use_cluster, collapse=', '))
print(table(seurat_obj@meta.data[[idents]], seurat_obj@meta.data[[groupby]]))

# ------------------------------------------------------------- features ----
nfeat <- if (opt$smoke) opt$smoke_genes else 2000
seurat_obj <- safe('FindVariableFeatures',
  FindVariableFeatures(seurat_obj, selection.method='vst', nfeatures=nfeat, assay=assays, verbose=FALSE))
seurat_obj <- safe('ScaleData',
  ScaleData(seurat_obj, features=VariableFeatures(seurat_obj), assay=assays, verbose=FALSE))
seurat_obj <- safe('SetupForWGCNA',
  SetupForWGCNA(seurat_obj, fraction=0.05, wgcna_name='tutorial'))
step('network genes:', length(GetWGCNAGenes(seurat_obj, wgcna_name='tutorial')))

reduct <- if ('harmony' %in% names(seurat_obj@reductions)) 'harmony' else if ('pca' %in% names(seurat_obj@reductions)) 'pca' else names(seurat_obj@reductions)[1]
step('reduction for metacells =', reduct)

gb <- if (opt$smoke) c(idents) else c(idents, groupby)
step('metacell group.by:', paste(gb, collapse=' + '))
seurat_obj <- safe('MetacellsByGroups',
  MetacellsByGroups(seurat_obj, group.by=gb, reduction=reduct,
                    k=kvalue, max_shared=10, ident.group=idents, min_cells=min_cells,
                    assay=assays))
if (is.null(seurat_obj)) stop('MetacellsByGroups failed')
seurat_obj <- safe('NormalizeMetacells', NormalizeMetacells(seurat_obj))
if (is.null(seurat_obj)) stop('NormalizeMetacells failed')

meta_obj <- seurat_obj@misc$tutorial$wgcna_metacell_obj
use_cluster <- intersect(use_cluster, unique(as.character(meta_obj@meta.data[[idents]])))
step('cell types after min_cells filter:', paste(use_cluster, collapse=', '))
if (!length(use_cluster)) stop('No cell type survived the min_cells filter')

seurat_obj <- safe('SetDatExpr',
  SetDatExpr(seurat_obj, group_name=use_cluster, group.by=idents,
             assay=assays))
if (is.null(seurat_obj)) stop('SetDatExpr failed')
safe('datExpr diagnostics', {
  de <- GetDatExpr(seurat_obj, wgcna_name='tutorial')
  step('  GetDatExpr class:', class(de)[1], '| length:', length(de))
  step('   datExpr dim:', paste(dim(de), collapse=' x '), '| NA:', sum(is.na(de)),
       '| zero-variance cols:', sum(apply(as.matrix(de), 2, function(z) length(unique(z)) < 2)))
})

# ----------------------------------------------------------- soft power ----
if (!opt$skip_power_sweep) {
  seurat_obj <- safe('TestSoftPowers', TestSoftPowers(seurat_obj, networkType='signed'))
  safe('PlotSoftPowers', {
    pl <- PlotSoftPowers(seurat_obj)
    p <- wrap_plots(pl, ncol=2)
    ggsave(file.path(outdir,'parameter_sweep.pdf'), p, w=12, h=10)
    ggsave(file.path(outdir,'parameter_sweep.png'), p, w=12, h=10, bg='white')
    invisible(NULL)
  })
}


  # ------------------------------------------------------- soft power ----
  if (soft_power == 'auto' && !opt$skip_power_sweep) {
    sp <- pick_power(seurat_obj)
    if (is.null(sp) || !is.finite(sp)) { sp <- 8; step('auto power unavailable -> fallback soft_power =', sp) }
    else step('auto-selected soft_power =', sp)
  } else {
    sp <- if (soft_power == 'auto') 8 else as.numeric(soft_power)
    step('using fixed soft_power =', sp)
  }
  if (!opt$no_checkpoint) {
    safe('save checkpoint', saveRDS(
      list(obj=seurat_obj, use_cluster=use_cluster, new_name=new_name, sp=sp, sig=sig),
      cp_file))
    step('checkpoint written:', cp_file)
  }
}

# ------------------------------------------------------------ construct ----
# hdWGCNA 0.4.12 stores the TOM path as paste0(getwd(), '/', <tom_outdir>...).
# An absolute tom_outdir produces a doubled path and GetTOM() fails. Run
# ConstructNetwork from inside outdir with tom_outdir = '.', then restore.
wd0 <- getwd()
seurat_obj <- tryCatch({
  setwd(outdir)
  safe('ConstructNetwork',
    ConstructNetwork(seurat_obj, soft_power=sp, tom_outdir='.',
                     tom_name='hdWGCNA', overwrite_tom=TRUE))
}, finally = setwd(wd0))
if (is.null(seurat_obj)) stop('ConstructNetwork failed')

safe('PlotDendrogram', {
  png(file.path(outdir,'hdWGCNA_Dendrogram.png'), width=1800, height=1200, res=200)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram'); dev.off()
  pdf(file.path(outdir,'hdWGCNA_Dendrogram.pdf'), width=9, height=6)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram'); dev.off()
})
TOM <- safe('GetTOM', GetTOM(seurat_obj))

# NOTE: the original template overrode hdWGCNA:::CheckSeurat5 here to force the
# Seurat-v4 code path. With SeuratObject 5.4 that path calls the defunct
# GetAssayData(slot=) and fails. Keep the native v5 code path intact.
step('CheckSeurat5 (native) =', hdWGCNA::CheckSeurat5())

seurat_obj <- safe('ModuleEigengenes',
  ModuleEigengenes(seurat_obj, group.by.vars=groupby))
if (is.null(seurat_obj)) stop('ModuleEigengenes failed')
seurat_obj <- safe('ModuleConnectivity',
  ModuleConnectivity(seurat_obj, group.by=idents, group_name=use_cluster))
seurat_obj <- safe('ResetModuleNames', ResetModuleNames(seurat_obj, new_name=new_name))
if (is.null(seurat_obj)) stop('ModuleConnectivity/ResetModuleNames failed')

# --------------------------------------------------------------- output ----
modules <- GetModules(seurat_obj)
write.csv(modules, file.path(outdir,'wgcna_modules.csv'), row.names=TRUE)
step('modules table:', nrow(modules), 'genes |',
     paste(sort(unique(na.omit(modules$module))), collapse=', '))

safe('GetHubGenes', {
  hub <- GetHubGenes(seurat_obj, n_hubs=10)
  write.csv(hub, file.path(outdir,'hub_genes.csv'), row.names=FALSE)
})

safe('PlotKMEs', {
  p <- PlotKMEs(seurat_obj, ncol=4, n_hubs=10)
  ggsave(file.path(outdir,'kME.png'), plot=p, width=18, height=7, bg='white')
  ggsave(file.path(outdir,'kME.pdf'), p, width=18, height=7, bg='white')
})

seurat_obj <- safe('ModuleExprScore',
  ModuleExprScore(seurat_obj, n_genes=25, method='Seurat'))

safe('ModuleFeaturePlot', {
  pl <- ModuleFeaturePlot(seurat_obj, features='hMEs', order_points=TRUE)
  p <- wrap_plots(pl, ncol=4)
  ggsave(file.path(outdir,'ModuleFeaturePlot.pdf'), p, height=7, width=18)
  ggsave(file.path(outdir,'ModuleFeaturePlot.png'), p, height=7, width=18, bg='white')
})

safe('ModuleCorrelogram', {
  png(file.path(outdir,'ModuleCorrelogram.png'), width=1600, height=1400, res=200)
  ModuleCorrelogram(seurat_obj); dev.off()
  pdf(file.path(outdir,'ModuleCorrelogram.pdf'), width=8, height=7)
  ModuleCorrelogram(seurat_obj); dev.off()
})

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
write.csv(MEs, file.path(outdir,'module_eigengenes.csv'))
mods <- setdiff(colnames(MEs), 'grey')
step('module eigengenes:', length(mods), 'non-grey modules x', nrow(MEs), 'cells')
safe('DotPlot', {
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  p <- DotPlot(seurat_obj, features=mods, group.by=idents) + RotatedAxis() +
       scale_color_gradient2(high='red', mid='grey95', low='blue')
  ggsave(file.path(outdir,'Dotplot.pdf'), p, height=7, width=9)
  ggsave(file.path(outdir,'Dotplot.png'), p, height=7, width=9, bg='white')
})
safe('module scores export', {
  sc <- as.data.frame(seurat_obj@meta.data)
  cc <- grep('^(M|sub-M)[0-9]+$|^hME|^ME|Score', colnames(sc), value=TRUE)
  if (length(cc)) write.csv(sc[, cc, drop=FALSE], file.path(outdir,'module_scores.csv'))
})

safe('saveRDS', saveRDS(seurat_obj, file.path(outdir,'WGCNA.rds')))
writeLines(capture.output(sessionInfo()), file.path(outdir,'sessionInfo.txt'))
writeLines(sprintf('rds=%s\nidents=%s\ngroupby=%s\nassay=%s\nk=%s\nmin_cells=%s\nsoft_power_used=%s\nresumed=%s\nsmoke=%s\ncelltypes=%s',
  rds, idents, groupby, assays, kvalue, min_cells, sp, resumed, opt$smoke,
  paste(use_cluster, collapse=',')), file.path(outdir,'run_params.txt'))
step('=== RUN COMPLETE ===')
