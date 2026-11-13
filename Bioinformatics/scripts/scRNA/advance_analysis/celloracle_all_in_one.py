#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""celloracle_all_in_one.py - single-file, fully parameterised CellOracle pipeline.

One script replaces the whole 00/01/02/03/04 toolchain:

    Seurat RDS  --(embedded R exporter)-->  MatrixMarket intermediates
                -->  AnnData (raw counts)  -->  Oracle + base GRN
                -->  PCA / knn imputation  -->  cluster-wise GRN (Links)
                -->  ridge GRN fitting     -->  in-silico TF perturbation
                -->  transition probability / embedding shift / figures / tables

ALL important files are supplied on the command line, including the promoter
base-GRN database (parquet). Run with --help for the full reference and
--check-env to verify the environment before a long run.

Author: generated for the sw2 xiyoucloud project, 2026-11-07
Version: 1.1.0
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

VERSION = "1.1.0"

DEFAULT_RSCRIPT = "/home/shazhi11/miniforge3/envs/R_Seurat/bin/Rscript"
DEFAULT_BASE_GRN = (
    "/home/shazhi11/celloracle_data/promoter_base_GRN/"
    "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"
)

# ---------------------------------------------------------------------------
# Embedded R code. Written to <workdir>/_celloracle_export.R and executed with
# the Rscript given by --rscript. Raw string so R's backslash escapes survive.
# ---------------------------------------------------------------------------
R_EXPORT = r'''
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript _celloracle_export.R <rds> <outdir> [assay] [reduction] [cluster] [n_per_group] [which_cells] [seed]")
infile    <- args[1]
outdir    <- args[2]
assay     <- if (length(args) >= 3) args[3] else "RNA"
reduction <- if (length(args) >= 4) args[4] else "umap"
cluster_col <- if (length(args) >= 5) args[5] else ""
cells_per_group <- suppressWarnings(if (length(args) >= 6) as.integer(args[6]) else 0L)
which_cells <- if (length(args) >= 7) args[7] else ""
seed <- suppressWarnings(if (length(args) >= 8) as.integer(args[8]) else 0L)
if (is.na(cells_per_group)) cells_per_group <- 0L
if (is.na(seed)) seed <- 0L

if (!file.exists(infile)) stop("input RDS not found: ", infile)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("[R] reading ", infile, "\n", sep = "")
obj <- readRDS(infile)
DefaultAssay(obj) <- assay

# Optional in-R subsetting. Essential for large objects: exporting all
# 50k+ cells would produce a multi-GB MatrixMarket text file.
if (nchar(cluster_col) > 0 && cluster_col %in% colnames(obj@meta.data)) {
  grp <- as.character(obj@meta.data[[cluster_col]])
  if (nchar(which_cells) > 0) {
    wanted <- trimws(strsplit(which_cells, ",", fixed = TRUE)[[1]])
    keep <- grp %in% wanted
    cat("[R] which_cells filter: ", sum(keep), " of ", length(keep), " cells\n", sep = "")
    if (sum(keep) == 0) stop("no cell matches --rds-which-cells: ", which_cells)
    obj <- subset(obj, cells = colnames(obj)[keep])
    grp <- as.character(obj@meta.data[[cluster_col]])
  }
  if (cells_per_group > 0) {
    set.seed(seed)
    idx <- unlist(lapply(split(seq_len(ncol(obj)), grp), function(ii) {
      if (length(ii) > cells_per_group) sample(ii, cells_per_group) else ii
    }), use.names = FALSE)
    obj <- subset(obj, cells = colnames(obj)[sort(idx)])
    cat("[R] downsampled to ", ncol(obj), " cells (<= ", cells_per_group,
        " per ", cluster_col, ")\n", sep = "")
  }
}
cat("[R] dim: ", nrow(obj), " genes x ", ncol(obj), " cells\n", sep = "")

counts <- GetAssayData(obj, layer = "counts")
if (is.null(counts)) stop("assay '", assay, "' has no 'counts' layer")

# Matrix::writeMM accepts a FILENAME, not a connection. Write plain, then gzip.
mtx <- file.path(outdir, "counts.mtx")
cat("[R] writing MatrixMarket (cells x genes) ...\n")
writeMM(t(counts), mtx)
cat("[R] plain size: ", file.info(mtx)$size, " bytes\n", sep = "")
unlink(paste0(mtx, ".gz"))
status <- system2("gzip", c("-f", shQuote(mtx)))
if (status != 0 || !file.exists(paste0(mtx, ".gz"))) stop("gzip failed")
cat("[R] gz size: ", file.info(paste0(mtx, ".gz"))$size, " bytes\n", sep = "")

writeLines(colnames(obj), file.path(outdir, "barcodes.txt"))
writeLines(rownames(obj), file.path(outdir, "genes.txt"))
write.csv(obj@meta.data, file.path(outdir, "metadata.csv"), row.names = TRUE)

emb <- Embeddings(obj, reduction = reduction)
if (is.null(emb) || ncol(emb) < 2) stop("reduction '", reduction, "' has no 2D embedding")
emb <- emb[, 1:2, drop = FALSE]
colnames(emb) <- c("UMAP_1", "UMAP_2")
write.csv(emb, file.path(outdir, "embedding.csv"), row.names = TRUE)

vf <- try(VariableFeatures(obj), silent = TRUE)
if (!inherits(vf, "try-error") && length(vf) > 0) {
  writeLines(vf, file.path(outdir, "variable_genes.txt"))
  cat("[R] exported ", length(vf), " variable features\n", sep = "")
}
cat("[R] done. intermediates in ", outdir, "\n", sep = "")
'''

_LOG_PATH = None


def log(msg, tag="MAIN"):
    line = "[%s] %s" % (tag, msg)
    print(line, flush=True)
    if _LOG_PATH is not None:
        with open(_LOG_PATH, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
DESCRIPTION = """\
Single-file CellOracle pipeline for Seurat RDS / AnnData.

Stages
  1) RDS -> intermediates   embedded R exporter, called through --rscript
  2) intermediates -> h5ad  raw-count AnnData + X_umap + cluster column
  3) h5ad + base GRN -> Oracle -> PCA -> knn imputation -> cluster-wise GRN
  4) ridge GRN fitting -> in-silico TF knockout(s) -> embedding shift
  5) tables, figures, and serialised Oracle objects

Input modes (pick exactly one; evaluated in this priority order)
  --reuse-oracle FILE   load a fitted Oracle and go straight to perturbations
  --h5ad FILE           start from an existing AnnData (raw counts, X_umap)
  --intermediates DIR   start from pre-exported MatrixMarket intermediates
  --rds FILE            full pipeline starting from a Seurat RDS

The promoter base GRN database (--base-grn) is a parquet file with columns
'peak_id', 'gene_short_name', and one 0/1 column per TF. The official human
hg38 file shipped with CellOracle works out of the box.
"""

EPILOG = """\
EXAMPLES
  # 0) Show help / inspect the environment
  python celloracle_all_in_one.py --help
  python celloracle_all_in_one.py --check-env --base-grn GRN.parquet

  # 1) Full run from a Seurat RDS (all important files passed in)
  python celloracle_all_in_one.py \\
      --rds GSE182109.NPC_like_run2.MORC4_CDC45_4groups.rds \\
      --base-grn hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet \\
      --rscript /home/shazhi11/miniforge3/envs/R_Seurat/bin/Rscript \\
      --cluster-column MORC4_CDC45_group --assay RNA --reduction umap \\
      --outdir results --bagging-number 3 --n-jobs 8

  # 2) Only convert RDS -> h5ad (fast, no GRN)
  python celloracle_all_in_one.py --rds input.rds --base-grn GRN.parquet \\
      --cluster-column CellType --outdir results --stop-after h5ad

  # 3) Reuse an existing h5ad, knock out specific TFs
  python celloracle_all_in_one.py --h5ad celloracle_input.h5ad \\
      --base-grn GRN.parquet --perturb-genes KLF6,ATF3,TFAP2A \\
      --outdir results_ko --bagging-number 3

  # 4) Reuse a fitted Oracle: seconds per knockout, no GRN re-fitting
  python celloracle_all_in_one.py \\
      --reuse-oracle results/oracle_after_grn.celloracle.oracle \\
      --base-grn GRN.parquet --perturb-genes KLF6,ATF3 --outdir results_screen

  # 5) Non-TF gene virtual knockout (e.g. YKT6). The gene is NOT a TF in the
  #    base GRN, so inject an EXPLORATORY co-expression prior. --gene-list
  #    keeps it in the HVG set; --rds-cells-per-group avoids a multi-GB matrix.
  python celloracle_all_in_one.py --rds Glioma_3B_Fully_Annotated.rds \\
      --base-grn GRN.parquet --assay decontX --cluster-column CellType \\
      --rds-cells-per-group 400 --max-cells-per-group 0 --gene-list ykt6.txt \\
      --custom-regulator YKT6 --n-custom-targets 500 --perturb-genes YKT6 \\
      --outdir results_YKT6 --bagging-number 1

  # 6) Print resolved configuration without running anything
  DRY RUN: add --dry-run to any command.

OUTPUTS (all under --outdir)
  pipeline.log, run_config.json, run_summary.txt, base_grn_info.txt
  celloracle_input.h5ad                 (when built by this script)
  links_<group>.csv / links_<group>_filtered.csv
  perturbable_tfs.csv, active_regulatory_genes.txt
  oracle_after_grn.celloracle.oracle, oracle_object.celloracle.oracle
  embedding_shift_<TF>.csv, quiver_<TF>.png, delta_magnitude_<TF>.png
  ko_screen_summary.csv, ko_screen_matrix.csv
  custom_regulator_targets.csv, custom_regulator_correlation.csv (non-TF mode)

NOTES
  * Requires conda env celloracle_env (celloracle 0.20.0) and the R_Seurat
    R environment. Activate with: conda activate celloracle_env
  * --base-grn is REQUIRED, even in reuse mode (it is validated and logged).
  * bagging_number=1 produces NO p-values; the script detects this and filters
    by top coef_abs edges instead. Use >=3 for p-value based filtering.
  * --custom-regulator lets you knock out a gene that is NOT a TF in the base
    GRN. It adds a CO-EXPRESSION prior (top |Pearson r| targets by default),
    NOT motif/ChIP evidence. Treat those effect sizes as exploratory.
    Use --gene-list to force-keep that gene in the HVG set, otherwise
    simulate_shift will reject it as "not included in expression matrix".
  * For large RDS files use --rds-cells-per-group (subset INSIDE R before the
    MatrixMarket export), otherwise the .mtx can reach several GB.
  * --seed fixes subsampling and the main-process RNG, but CellOracle's
    parallel bagging workers keep their own random state: two GRN fits can
    still differ slightly. For strict reproducibility, fit ONCE (--stop-after
    grn or a normal run), then reuse the saved
    oracle_after_grn.celloracle.oracle with --reuse-oracle for every KO.
  * de novo motif scanning is not performed; only the precomputed base GRN.
"""


def build_parser():
    p = argparse.ArgumentParser(
        prog="celloracle_all_in_one.py",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version="celloracle_all_in_one.py %s" % VERSION)

    g = p.add_argument_group("INPUT FILES (choose one input mode)")
    g.add_argument("--rds", metavar="FILE",
                   help="Seurat .rds object to convert (full pipeline)")
    g.add_argument("--h5ad", metavar="FILE",
                   help="existing AnnData with RAW counts, obsm['X_umap'] and a cluster column")
    g.add_argument("--intermediates", metavar="DIR",
                   help="directory with counts.mtx(.gz), barcodes.txt, genes.txt, metadata.csv, embedding.csv")
    g.add_argument("--reuse-oracle", metavar="FILE",
                   help="fitted *.celloracle.oracle cache; skip conversion + GRN fitting")
    g.add_argument("--base-grn", metavar="FILE", required=True,
                   help="promoter base GRN database, parquet (REQUIRED). Default location: %s"
                        % DEFAULT_BASE_GRN)
    g.add_argument("--rscript", metavar="FILE", default=DEFAULT_RSCRIPT,
                   help="Rscript executable used for RDS export (default: %(default)s)")
    g.add_argument("--gene-list", metavar="FILE",
                   help="optional one-gene-per-line file; these genes are force-kept when selecting HVGs")

    g = p.add_argument_group("OUTPUT FILES / DIRECTORIES")
    g.add_argument("--outdir", metavar="DIR", default="results",
                   help="all result files go here (default: %(default)s)")
    g.add_argument("--workdir", metavar="DIR", default=None,
                   help="intermediates directory (default: <outdir>/intermediates)")
    g.add_argument("--h5ad-out", metavar="FILE", default=None,
                   help="where to write the generated h5ad (default: <outdir>/celloracle_input.h5ad)")

    g = p.add_argument_group("DATA / ANNDATA OPTIONS")
    g.add_argument("--cluster-column", metavar="NAME", default=None,
                   help="obs column defining GRN units; default: uns['cluster_column_name'] or the Oracle's own")
    g.add_argument("--assay", default="RNA", help="Seurat assay for RDS export (default: %(default)s)")
    g.add_argument("--reduction", default="umap", help="2D reduction for RDS export (default: %(default)s)")
    g.add_argument("--n-hvg", type=int, default=3000,
                   help="number of highly variable genes, 0 keeps all (default: %(default)s)")
    g.add_argument("--max-cells-per-group", type=int, default=1000,
                   help="subsample size per cluster in the Python stage, 0 keeps all (default: %(default)s)")
    g.add_argument("--rds-cells-per-group", type=int, default=0,
                   help="subsample INSIDE R before exporting; 0 exports all cells. Use this for large RDS (default: %(default)s)")
    g.add_argument("--rds-which-cells", default="",
                   help="comma separated labels of --cluster-column to keep in the R stage, e.g. 'AC-like,NPC-like' (default: all)")
    g.add_argument("--seed", type=int, default=0,
                   help="seed for subsampling and the main-process RNG; see NOTES on reproducibility (default: %(default)s)")

    g = p.add_argument_group("NON-TF / CUSTOM REGULATOR (optional, exploratory)")
    g.add_argument("--custom-regulator", metavar="GENE", default=None,
                   help="inject a gene that is NOT a TF in the base GRN (e.g. YKT6) so it can be knocked out. "
                        "A synthetic prior is added: the gene becomes a regulator of a data-driven target set")
    g.add_argument("--custom-targets", metavar="FILE", default=None,
                   help="one-gene-per-line target list for --custom-regulator; default: top correlated genes")
    g.add_argument("--n-custom-targets", type=int, default=500,
                   help="number of co-expression targets when --custom-targets is omitted (default: %(default)s)")

    g = p.add_argument_group("BASE GRN / NETWORK OPTIONS")
    g.add_argument("--links-alpha", type=float, default=10,
                   help="ridge alpha for cluster-wise GRN / Links (default: %(default)s)")
    g.add_argument("--fit-alpha", type=float, default=1,
                   help="ridge alpha for fit_GRN_for_simulation (default: %(default)s)")
    g.add_argument("--bagging-number", type=int, default=20,
                   help="Links bagging rounds; >=3 yields p-values, 1 is a fast smoke test (default: %(default)s)")
    g.add_argument("--p-threshold", type=float, default=0.001,
                   help="p-value cutoff for link filtering when p-values exist (default: %(default)s)")
    g.add_argument("--top-links", type=int, default=3000,
                   help="max edges kept per GRN unit by coef_abs (default: %(default)s)")

    g = p.add_argument_group("PERTURBATION OPTIONS")
    g.add_argument("--perturb-genes", metavar="G1,G2,...", default=None,
                   help="comma separated TF symbols to knock out; default: auto-pick highest out-degree active TFs")
    g.add_argument("--auto-top", type=int, default=1,
                   help="how many TFs to pick automatically when --perturb-genes is omitted (default: %(default)s)")
    g.add_argument("--perturb-value", type=float, default=0.0,
                   help="expression value after perturbation, 0.0 = knockout (default: %(default)s)")
    g.add_argument("--n-propagation", type=int, default=3,
                   help="signal propagation steps (default: %(default)s)")
    g.add_argument("--n-neighbors", type=int, default=None,
                   help="neighbours for transition probability (default: min(200, n_cells/5))")
    g.add_argument("--sigma-corr", type=float, default=0.05,
                   help="kernel scale for embedding shift (default: %(default)s)")

    g = p.add_argument_group("EXECUTION CONTROL")
    g.add_argument("--n-jobs", type=int, default=8, help="parallel workers (default: %(default)s)")
    g.add_argument("--stop-after", choices=["h5ad", "grn", "all"], default="all",
                   help="stop after building h5ad / fitting GRN / finish everything (default: %(default)s)")
    g.add_argument("--dry-run", action="store_true",
                   help="print the resolved configuration and exit without running")
    g.add_argument("--check-env", action="store_true",
                   help="verify python/R/dependencies/input files and exit")
    g.add_argument("--no-plots", action="store_true", help="skip figure generation")
    g.add_argument("--no-cache", action="store_true", help="do not serialise Oracle objects")
    return p


# ---------------------------------------------------------------------------
# Environment check / validation
# ---------------------------------------------------------------------------
def check_env(a):
    print("=" * 68)
    print("celloracle_all_in_one.py %s  -- environment check" % VERSION)
    print("=" * 68)
    print("python            : %s" % sys.version.split()[0])
    print("executable        : %s" % sys.executable)
    ok = True
    for name in ("numpy", "pandas", "scipy", "anndata", "scanpy",
                 "matplotlib", "pyarrow", "celloracle", "velocyto"):
        try:
            mod = __import__(name)
            print("  %-14s : %s" % (name, getattr(mod, "__version__", "OK")))
        except Exception as exc:
            ok = False
            print("  %-14s : MISSING (%s)" % (name, exc))
    try:
        import celloracle as co
        st = co.check_python_requirements(return_detail=True, print_warning=False)
        bad = st[st.requirement_satisfied == False]
        print("celloracle requirements unsatisfied : %d" % len(bad))
        if len(bad):
            print(bad.to_string(index=False))
    except Exception as exc:
        print("celloracle requirements check failed : %s" % exc)
    rscript = shutil.which(a.rscript) if not os.path.isabs(a.rscript) else a.rscript
    print("Rscript           : %s (%s)" %
          (a.rscript, "FOUND" if rscript and os.path.isfile(rscript) else "NOT FOUND"))
    print("base GRN          : %s (%s)" %
          (a.base_grn, "FOUND" if Path(a.base_grn).is_file() else "NOT FOUND"))
    if not Path(a.base_grn).is_file():
        ok = False
    for label, path in (("rds", a.rds), ("h5ad", a.h5ad),
                        ("intermediates", a.intermediates),
                        ("reuse-oracle", a.reuse_oracle),
                        ("gene-list", a.gene_list)):
        if path:
            good = Path(path).exists()
            print("%-18s: %s (%s)" % (label, path, "FOUND" if good else "NOT FOUND"))
            if not good:
                ok = False
    print("=" * 68)
    print("RESULT: %s" % ("OK" if ok else "PROBLEMS FOUND - see above"))
    return 0 if ok else 2


def validate(a, parser):
    modes = [m for m, v in (("--reuse-oracle", a.reuse_oracle),
                            ("--h5ad", a.h5ad),
                            ("--intermediates", a.intermediates),
                            ("--rds", a.rds)) if v]
    if not modes:
        parser.error("one input is required: --rds / --h5ad / --intermediates / --reuse-oracle")
    if len(modes) > 1:
        parser.error("only ONE input mode is allowed, got: %s" % ", ".join(modes))
    for label, value in (("--reuse-oracle", a.reuse_oracle), ("--h5ad", a.h5ad),
                         ("--rds", a.rds), ("--gene-list", a.gene_list),
                         ("--custom-targets", a.custom_targets)):
        if value and not Path(value).is_file():
            parser.error("%s file does not exist: %s" % (label, value))
    if a.custom_regulator and modes[0] == "--reuse-oracle":
        parser.error("--custom-regulator needs a fresh GRN fit; do not combine with --reuse-oracle")
    if a.rds_cells_per_group < 0 or a.n_custom_targets < 1:
        parser.error("--rds-cells-per-group must be >= 0 and --n-custom-targets must be >= 1")
    if a.intermediates and not Path(a.intermediates).is_dir():
        parser.error("--intermediates directory does not exist: %s" % a.intermediates)
    if not Path(a.base_grn).is_file():
        parser.error("--base-grn database file does not exist: %s" % a.base_grn)
    if a.perturb_value < 0:
        parser.error("--perturb-value must be >= 0")
    if a.auto_top < 1:
        parser.error("--auto-top must be >= 1")
    if a.n_hvg < 0 or a.max_cells_per_group < 0:
        parser.error("--n-hvg and --max-cells-per-group must be >= 0")
    if modes[0] == "--rds":
        rscript = a.rscript if os.path.isabs(a.rscript) else shutil.which(a.rscript)
        if not rscript or not Path(rscript).is_file():
            parser.error("--rscript not found: %s" % a.rscript)
    return modes[0]


# ---------------------------------------------------------------------------
# Stage 1: RDS -> intermediates (embedded R)
# ---------------------------------------------------------------------------
def export_rds(a, indir):
    import subprocess
    rscript = a.rscript if os.path.isabs(a.rscript) else (shutil.which(a.rscript) or a.rscript)
    script = indir / "_celloracle_export.R"
    script.write_text(R_EXPORT, encoding="utf-8")
    cmd = [str(rscript), str(script), str(Path(a.rds).resolve()),
           str(indir), a.assay, a.reduction,
           a.cluster_column or "",
           str(a.rds_cells_per_group), a.rds_which_cells, str(a.seed)]
    log("running embedded R exporter: %s" % " ".join(cmd), "R")
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    for stream, text in (("stdout", proc.stdout), ("stderr", proc.stderr)):
        for line in (text or "").splitlines():
            log("[R:%s] %s" % (stream, line), "R")
    if proc.returncode != 0:
        raise RuntimeError("R exporter failed with exit code %d" % proc.returncode)
    required = ["barcodes.txt", "genes.txt", "metadata.csv", "embedding.csv"]
    missing = [f for f in required if not (indir / f).is_file()]
    if not (indir / "counts.mtx.gz").is_file() and not (indir / "counts.mtx").is_file():
        missing.append("counts.mtx.gz")
    if missing:
        raise RuntimeError("R export incomplete, missing: %s" % ", ".join(missing))
    log("R export complete -> %s" % indir, "R")


# ---------------------------------------------------------------------------
# Stage 2: intermediates -> AnnData
# ---------------------------------------------------------------------------
def build_h5ad(a, indir, out):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.io as sio
    from anndata.utils import make_index_unique

    gz, plain = indir / "counts.mtx.gz", indir / "counts.mtx"
    log("reading MatrixMarket ...", "H5AD")
    if gz.is_file():
        with gzip.open(gz, "rt") as fh:
            X = sio.mmread(fh).tocsr().astype("float32")   # cells x genes
    elif plain.is_file():
        X = sio.mmread(str(plain)).tocsr().astype("float32")
    else:
        raise FileNotFoundError("no counts.mtx(.gz) in %s" % indir)

    barcodes = [str(x) for x in (indir / "barcodes.txt").read_text().splitlines()]
    genes = [str(x) for x in (indir / "genes.txt").read_text().splitlines()]
    meta = pd.read_csv(indir / "metadata.csv", index_col=0)
    emb = pd.read_csv(indir / "embedding.csv", index_col=0)
    meta.index = meta.index.astype(str)
    emb.index = emb.index.astype(str)
    if len(set(genes)) != len(genes):
        genes = make_index_unique(pd.Index(genes)).tolist()
    if len(set(barcodes)) != len(barcodes):
        barcodes = make_index_unique(pd.Index(barcodes)).tolist()
    if X.shape != (len(barcodes), len(genes)):
        raise ValueError("shape mismatch: matrix %s vs %d barcodes / %d genes" %
                         (X.shape, len(barcodes), len(genes)))

    meta = meta.loc[barcodes]
    emb = emb.loc[barcodes]
    adata = ad.AnnData(X=X.copy(), obs=meta.copy(),
                       var=pd.DataFrame(index=pd.Index(genes)))
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["X_umap"] = emb.values.astype("float32")

    col = a.cluster_column
    if col is None:
        raise ValueError("--cluster-column is required when building an h5ad")
    if col not in adata.obs.columns:
        raise ValueError("cluster column '%s' not found; available: %s" %
                         (col, list(adata.obs.columns)))
    adata.obs[col] = adata.obs[col].astype(str)
    adata = adata[~adata.obs[col].isin(["nan", "NA", "", "None"])].copy()
    log("cluster counts:\n%s" % adata.obs[col].value_counts().to_string(), "H5AD")

    if a.max_cells_per_group > 0:
        rng = np.random.default_rng(a.seed)
        idx = []
        for _, sub in adata.obs.groupby(col, observed=True):
            names = sub.index.to_numpy()
            if len(names) > a.max_cells_per_group:
                names = rng.choice(names, a.max_cells_per_group, replace=False)
            idx.extend(names.tolist())
        adata = adata[idx].copy()
        log("downsampled to %d cells" % adata.n_obs, "H5AD")

    if 0 < a.n_hvg < adata.n_vars:
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(tmp, n_top_genes=a.n_hvg, flavor="seurat")
        hvg = tmp.var_names[tmp.var.highly_variable].tolist()
        # --gene-list comes FIRST so force-kept genes survive the [:n_hvg] cut.
        forced = []
        if a.gene_list:
            forced += Path(a.gene_list).read_text().splitlines()
        vf = indir / "variable_genes.txt"
        if vf.is_file():
            forced += vf.read_text().splitlines()
        hvg = [g for g in list(dict.fromkeys(forced + hvg))
               if g in adata.var_names][:a.n_hvg]
        adata = adata[:, hvg].copy()
        log("kept %d genes (%d force-kept requested)" %
            (adata.n_vars, len(forced)), "H5AD")

    adata.uns["cluster_column_name"] = col
    adata.uns["embedding_name"] = "X_umap"
    adata.X = adata.X.astype("float32")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    log("wrote %s : %d cells x %d genes" % (out, adata.n_obs, adata.n_vars), "H5AD")
    return adata


# ---------------------------------------------------------------------------
# Base GRN
# ---------------------------------------------------------------------------
def load_base_grn(path, outdir):
    import pandas as pd
    log("reading base GRN database: %s" % path, "GRN")
    grn = pd.read_parquet(path)
    if "peak_id" not in grn.columns or "gene_short_name" not in grn.columns:
        raise ValueError("--base-grn must contain 'peak_id' and 'gene_short_name' columns; got %s"
                         % list(grn.columns[:6]))
    tf_cols = [c for c in grn.columns if c not in ("peak_id", "gene_short_name")]
    if len(tf_cols) < 10:
        raise ValueError("--base-grn has only %d TF columns; is this the right file?" % len(tf_cols))
    info = ("path: %s\nshape: %s\nTF columns: %d\ntarget genes: %d\n"
            "example TFs: %s\n" %
            (path, grn.shape, len(tf_cols), grn["gene_short_name"].nunique(),
             ", ".join(map(str, tf_cols[:12]))))
    (outdir / "base_grn_info.txt").write_text(info, encoding="utf-8")
    log("base GRN %s, %d TFs, %d target genes" %
        (str(grn.shape), len(tf_cols), grn["gene_short_name"].nunique()), "GRN")
    return grn, set(map(str, tf_cols))


def augment_base_grn(a, adata, grn, tf_cols, outdir):
    """Inject a NON-TF gene as a synthetic regulator into the base GRN.

    CellOracle can only knock out genes present as regulators (TF columns) in
    the base GRN. Genes such as YKT6 are not TFs and have no motif prior, so we
    add a clearly-labelled, EXPLORATORY co-expression prior: the gene becomes
    a candidate regulator of either a user-supplied target list, or of the top
    |Pearson r| genes computed from log-normalised counts.

    This is NOT motif/ChIP evidence. Downstream effect sizes must be
    interpreted as exploratory, not as transcriptional regulation.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    gene = a.custom_regulator
    genes = list(map(str, adata.var_names))
    if gene not in genes:
        raise ValueError("--custom-regulator '%s' is not in adata.var_names; "
                         "use --gene-list to force-keep it in the HVG set" % gene)
    if gene in tf_cols:
        log("'%s' is already a TF in the base GRN; no injection needed" % gene, "CUSTOM")
        return grn, [], {"mode": "already_tf", "n_targets": 0}

    if a.custom_targets:
        wanted = [x.strip() for x in Path(a.custom_targets).read_text().splitlines() if x.strip()]
        targets = [g for g in dict.fromkeys(wanted) if g in genes and g != gene]
        source = "user file %s" % a.custom_targets
        corr_table = None
    else:
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        M = tmp.X.toarray() if sp.issparse(tmp.X) else np.asarray(tmp.X, dtype=float)
        j = genes.index(gene)
        x = M[:, j].astype(float)
        xc = x - x.mean()
        Mc = M - M.mean(axis=0)
        denom = np.sqrt((Mc * Mc).sum(axis=0) * (xc * xc).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = (Mc * xc[:, None]).sum(axis=0) / denom
        corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
        order = np.argsort(-np.abs(corr))
        corr_table = pd.DataFrame({"gene": [genes[i] for i in order],
                                   "pearson_r": [float(corr[i]) for i in order]})
        corr_table.to_csv(outdir / "custom_regulator_correlation.csv", index=False)
        targets = [genes[i] for i in order if genes[i] != gene][:a.n_custom_targets]
        source = "top %d |Pearson r| on log-normalised counts" % a.n_custom_targets

    if not targets:
        raise RuntimeError("no usable targets for --custom-regulator")
    pd.DataFrame({"target_gene": targets}).to_csv(
        outdir / "custom_regulator_targets.csv", index=False)

    grn = grn.copy()
    grn[gene] = 0.0
    mask = grn["gene_short_name"].astype(str).isin(set(targets))
    grn.loc[mask, gene] = 1.0
    info = {"mode": "injected", "gene": gene, "n_targets": len(targets),
            "source": source, "n_peaks_marked": int(mask.sum())}
    log("INJECTED synthetic regulator %s -> %d target genes (%s), %d peak rows" %
        (gene, len(targets), source, int(mask.sum())), "CUSTOM")
    log("WARNING: this is a co-expression prior, NOT motif evidence. "
        "Treat the KO effect as exploratory.", "CUSTOM")
    return grn, targets, info


# ---------------------------------------------------------------------------
# Stages 3-5: Oracle / GRN / perturbation
# ---------------------------------------------------------------------------
def fit_grn(a, adata, grn, tf_cols, cluster_col, outdir,
            custom_targets=None, custom_regulator=None):
    import random
    import celloracle as co
    import numpy as np
    import pandas as pd

    # CellOracle's bagging / ridge fitting draws random subsets. Seeding the
    # main-process RNGs removes one source of variation, but parallel loky
    # workers keep their own state: use --reuse-oracle for exact reproducibility.
    random.seed(a.seed)
    np.random.seed(a.seed)
    log("random seed fixed at %d" % a.seed, "GRN")

    overlap = sorted(set(map(str, adata.var_names)) & tf_cols)
    log("TF symbols shared with h5ad: %d (e.g. %s)" %
        (len(overlap), ", ".join(overlap[:12])), "GRN")
    if not overlap:
        raise RuntimeError("no TF overlap between h5ad genes and --base-grn")

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata, cluster_column_name=cluster_col, embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=grn)
    log("TFdict size: %d target genes" % len(getattr(oracle, "TFdict", {})), "GRN")

    oracle.perform_PCA(n_components=min(30, adata.n_obs - 1, adata.n_vars - 1))
    oracle.knn_imputation(n_pca_dims=min(30, adata.n_obs - 1),
                          #k=min(50, max(5, adata.n_obs // 20)),
                          #balanced=False,
                          k=int(0.025*adata.n_obs),
                          balanced=True,
                          n_jobs=a.n_jobs)
    log("imputed matrix: %s" % str(oracle.adata.shape), "GRN")

    log("fitting cluster-wise GRN (bagging=%d, alpha=%g)" %
        (a.bagging_number, a.links_alpha), "GRN")
    links = oracle.get_links(cluster_name_for_GRN_unit=cluster_col,
                             alpha=a.links_alpha,
                             bagging_number=a.bagging_number,
                             verbose_level=0, n_jobs=a.n_jobs)
    log("raw links per group: %s" %
        {k: len(v) for k, v in links.links_dict.items()}, "GRN")

    first = next(iter(links.links_dict.values()))
    has_p = bool(first["p"].notna().any()) if "p" in first.columns else False
    p_filter_used = False
    if has_p:
        links.filter_links(p=a.p_threshold, weight="coef_abs",
                           threshold_number=a.top_links)
        n_try = sum(len(v) for v in links.filtered_links.values())
        if n_try < max(100, a.top_links // 10):
            log("p-value filter kept only %d edges; falling back to top %d coef_abs edges" %
                (n_try, a.top_links), "GRN")
            links.filter_links(p=None, weight="coef_abs",
                               threshold_number=a.top_links)
        else:
            p_filter_used = True
    else:
        log("p-values unavailable (bagging=%d); filtering by top %d coef_abs edges" %
            (a.bagging_number, a.top_links), "GRN")
        links.filter_links(p=None, weight="coef_abs",
                           threshold_number=a.top_links)
    n_filtered = {k: len(v) for k, v in links.filtered_links.items()}
    log("filtered links per group: %s" % n_filtered, "GRN")
    if sum(n_filtered.values()) == 0:
        raise RuntimeError("all links were filtered out; raise --top-links or use --bagging-number >= 3")

    for grp, df in links.links_dict.items():
        df.to_csv(outdir / ("links_%s.csv" % str(grp).replace("/", "_")), index=False)
    for grp, df in links.filtered_links.items():
        df.to_csv(outdir / ("links_%s_filtered.csv" % str(grp).replace("/", "_")), index=False)

    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    if custom_regulator and custom_targets:
        # Link filtering can delete every synthetic edge. Re-insert the custom
        # regulator into each cluster-specific TFdict so the ridge fit (and
        # therefore active_regulatory_genes) actually contains it.
        import numpy as np
        n_added = 0
        for _cluster, tf_dict in oracle.cluster_specific_TFdict.items():
            for tgt in custom_targets:
                if tgt in tf_dict:
                    vals = [str(x) for x in tf_dict[tgt]]
                    if custom_regulator not in vals:
                        vals.append(custom_regulator)
                        tf_dict[tgt] = np.array(vals, dtype=object)
                        n_added += 1
                else:
                    tf_dict[tgt] = np.array([custom_regulator], dtype=object)
                    n_added += 1
        log("re-inserted custom regulator %s into cluster-specific TFdict (%d entries)" %
            (custom_regulator, n_added), "CUSTOM")
    log("fitting GRN for simulation (alpha=%g)" % a.fit_alpha, "GRN")
    oracle.fit_GRN_for_simulation(alpha=a.fit_alpha,
                                  use_cluster_specific_TFdict=True,
                                  verbose_level=0)

    active = [str(g) for g in (getattr(oracle, "active_regulatory_genes", None) or [])]
    allreg = [str(g) for g in (getattr(oracle, "all_regulatory_genes_in_TFdict", None) or [])]
    var = set(map(str, adata.var_names))
    cand = sorted(set(active) & set(allreg) & var)
    src = pd.concat([d[["source"]] for d in links.filtered_links.values()],
                    ignore_index=True)["source"].astype(str)
    degree = src.value_counts()
    ranked = sorted(cand, key=lambda g: (-int(degree.get(g, 0)), g))
    pd.DataFrame({"TF": ranked,
                  "out_degree": [int(degree.get(g, 0)) for g in ranked]}) \
        .to_csv(outdir / "perturbable_tfs.csv", index=False)
    (outdir / "active_regulatory_genes.txt").write_text("\n".join(active) + "\n")
    log("active regulators: %d; perturbable TFs: %d; top: %s" %
        (len(active), len(cand), ", ".join(ranked[:12])), "GRN")
    if not cand:
        raise RuntimeError("no perturbable TF (empty active_regulatory_genes); "
                           "check --top-links / --bagging-number")

    meta = {"has_p": has_p, "p_filter_used": p_filter_used,
            "n_filtered": n_filtered, "n_active": len(active),
            "ranked": ranked, "degree": degree.to_dict()}
    return oracle, meta


def save_cache(oracle, path):
    import celloracle as co  # noqa: F401  (ensures package is loaded)
    oracle.to_hdf5(str(path))
    return path


def perturb_one(oracle, gene, a, cluster_col, outdir, make_plots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    log("simulating KO of %s (value=%g, n_propagation=%d)" %
        (gene, a.perturb_value, a.n_propagation), "SIM")
    oracle.simulate_shift(perturb_condition={gene: a.perturb_value},
                          n_propagation=a.n_propagation)
    n_neighbors = a.n_neighbors or min(200, max(20, oracle.adata.n_obs // 5))
    oracle.estimate_transition_prob(n_neighbors=n_neighbors,
                                    knn_random=True, 
                                    sampled_fraction=0.7,
                                    n_jobs=a.n_jobs, 
                                    calculate_randomized=True)
    oracle.calculate_embedding_shift(sigma_corr=a.sigma_corr)

    if hasattr(oracle, "delta_embedding"):
        delta = np.asarray(oracle.delta_embedding)
    elif "delta_embedding" in oracle.adata.obsm:
        delta = np.asarray(oracle.adata.obsm["delta_embedding"])
    else:
        raise RuntimeError("delta_embedding was not produced")
    emb = np.asarray(oracle.embedding)
    mag = np.sqrt((delta ** 2).sum(axis=1))

    frame = pd.DataFrame({
        "cell": oracle.adata.obs_names,
        cluster_col: oracle.adata.obs[cluster_col].values,
        "umap_1": emb[:, 0], "umap_2": emb[:, 1],
        "delta_1": delta[:, 0], "delta_2": delta[:, 1],
        "delta_magnitude": mag,
    })
    frame.to_csv(outdir / ("embedding_shift_%s.csv" % gene), index=False)
    rows = [{"TF": gene, "group": "ALL", "mean_delta_magnitude": float(mag.mean())}]
    for grp, value in frame.groupby(cluster_col)["delta_magnitude"].mean().items():
        rows.append({"TF": gene, "group": grp, "mean_delta_magnitude": float(value)})
    log("  %s mean |delta| = %.6f" % (gene, mag.mean()), "SIM")

    if make_plots:
        fig, ax = plt.subplots(figsize=(7, 6))
        groups = oracle.adata.obs[cluster_col].astype(str).values
        for g in pd.unique(groups):
            m = groups == g
            ax.scatter(emb[m, 0], emb[m, 1], s=6, alpha=0.5, label=str(g))
        idx = np.random.default_rng(0).choice(len(emb), min(400, len(emb)), replace=False)
        ax.quiver(emb[idx, 0], emb[idx, 1], delta[idx, 0], delta[idx, 1],
                  color="black", alpha=0.55, width=0.0025)
        ax.set_title("CellOracle in-silico KO: %s" % gene)
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        ax.legend(fontsize=7, markerscale=2)
        fig.tight_layout()
        fig.savefig(outdir / ("quiver_%s.pdf" % gene))
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(mag, bins=40, color="#3A4892")
        ax.set_xlabel("|delta_embedding|")
        ax.set_ylabel("cells")
        ax.set_title("Perturbation effect: %s" % gene)
        fig.tight_layout()
        fig.savefig(outdir / ("delta_magnitude_%s.pdf" % gene))
        plt.close(fig)
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    global _LOG_PATH
    parser = build_parser()
    a = parser.parse_args(argv)
    if a.check_env:
        return check_env(a)

    a.outdir = Path(a.outdir).resolve()
    a.workdir = Path(a.workdir).resolve() if a.workdir else a.outdir / "intermediates"
    a.h5ad_out = Path(a.h5ad_out).resolve() if a.h5ad_out else a.outdir / "celloracle_input.h5ad"
    a.cache_path = a.outdir / "oracle_after_grn.celloracle.oracle"

    if a.dry_run:
        printable = {k: (str(v) if isinstance(v, Path) else v)
                     for k, v in vars(a).items()}
        print(json.dumps(printable, indent=2, ensure_ascii=False, default=str))
        return 0

    mode = validate(a, parser)
    a.outdir.mkdir(parents=True, exist_ok=True)
    a.workdir.mkdir(parents=True, exist_ok=True)
    _LOG_PATH = a.outdir / "pipeline.log"
    _LOG_PATH.write_text("", encoding="utf-8")
    t0 = time.time()

    log("celloracle_all_in_one.py v%s | mode=%s | outdir=%s" %
        (VERSION, mode, a.outdir))
    (a.outdir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v)
                    for k, v in vars(a).items()},
                   indent=2, ensure_ascii=False, default=str), encoding="utf-8")

    grn, tf_cols = load_base_grn(a.base_grn, a.outdir)
    import pandas as pd

    # ---- resolve the Oracle object -------------------------------------
    used_genes = []
    meta = {"has_p": None, "p_filter_used": None,
            "n_filtered": {}, "n_active": 0, "ranked": [], "degree": {}}
    custom_info = None
    oracle = None
    cluster_col = a.cluster_column
    cache_ok = False

    if mode == "--reuse-oracle":
        import celloracle as co
        log("loading fitted Oracle: %s" % a.reuse_oracle, "MAIN")
        ref = co.load_hdf5(a.reuse_oracle)
        cluster_col = cluster_col or ref.cluster_column_name
        oracle = ref
        cache_ok = True
        a.cache_path = Path(a.reuse_oracle)
        active = [str(g) for g in (getattr(ref, "active_regulatory_genes", None) or [])]
        allreg = [str(g) for g in (getattr(ref, "all_regulatory_genes_in_TFdict", None) or [])]
        var = set(map(str, ref.adata.var_names))
        cand = sorted(set(active) & set(allreg) & var)
        ranked = []
        for table in (a.outdir / "perturbable_tfs.csv",
                      Path(a.reuse_oracle).parent / "perturbable_tfs.csv"):
            if table.is_file():
                ranked = [g for g in pd.read_csv(table)["TF"].astype(str).tolist()
                          if g in cand]
                if ranked:
                    log("loaded %d ranked TFs from %s" % (len(ranked), table), "MAIN")
                    break
        if not ranked:
            ranked = sorted(cand)
        meta.update(n_active=len(active), ranked=ranked)
        log("loaded Oracle: %d cells x %d genes, cluster=%s, %d active regulators" %
            (ref.adata.n_obs, ref.adata.n_vars, cluster_col, len(active)))
    else:
        # ---- stage 1 + 2 -------------------------------------------------
        if mode == "--rds":
            export_rds(a, a.workdir)
            adata = build_h5ad(a, a.workdir, a.h5ad_out)
        elif mode == "--intermediates":
            adata = build_h5ad(a, Path(a.intermediates).resolve(), a.h5ad_out)
        else:  # --h5ad
            import anndata as ad
            log("reading AnnData: %s" % a.h5ad, "H5AD")
            adata = ad.read_h5ad(a.h5ad)
            cluster_col = cluster_col or adata.uns.get("cluster_column_name")
            log("AnnData: %d cells x %d genes, cluster=%s" %
                (adata.n_obs, adata.n_vars, cluster_col), "H5AD")

        if cluster_col is None:
            raise ValueError("no cluster column; pass --cluster-column")
        if cluster_col not in adata.obs.columns:
            raise ValueError("cluster column '%s' not found; available: %s" %
                             (cluster_col, list(adata.obs.columns)))
        if "X_umap" not in adata.obsm:
            raise ValueError("adata.obsm['X_umap'] is required")
        if adata.X.min() < 0:
            raise ValueError("adata.X must be raw counts (found negative values)")

        if a.stop_after == "h5ad":
            log("--stop-after h5ad; stopping. h5ad: %s" %
                (a.h5ad_out if mode != "--h5ad" else a.h5ad), "MAIN")
            (a.outdir / "run_summary.txt").write_text(
                "mode %s\nstop_after h5ad\nelapsed_sec %.1f\n" %
                (mode, time.time() - t0), encoding="utf-8")
            log("DONE in %.1f s" % (time.time() - t0))
            return 0

        # ---- stage 3 ------------------------------------------------------
        custom_targets = []
        if a.custom_regulator:
            grn, custom_targets, custom_info = augment_base_grn(
                a, adata, grn, tf_cols, a.outdir)
            tf_cols = tf_cols | {a.custom_regulator}
        oracle, meta = fit_grn(a, adata, grn, tf_cols, cluster_col, a.outdir,
                               custom_targets, a.custom_regulator)
        try:
            if not a.no_cache:
                save_cache(oracle, a.cache_path)
                cache_ok = True
                log("cached %s" % a.cache_path, "GRN")
        except Exception as exc:
            log("cache skipped: %s: %s" % (type(exc).__name__, exc), "GRN")

    if a.stop_after == "grn":
        log("--stop-after grn; stopping.", "MAIN")
        (a.outdir / "run_summary.txt").write_text(
            "mode %s\nstop_after grn\nactive_regulators %d\nelapsed_sec %.1f\n" %
            (mode, meta["n_active"], time.time() - t0), encoding="utf-8")
        log("DONE in %.1f s" % (time.time() - t0))
        return 0

    # ---- stage 4: perturbation(s) ----------------------------------------
    ranked = meta["ranked"]
    if a.perturb_genes:
        used_genes = [g.strip() for g in a.perturb_genes.split(",") if g.strip()]
    else:
        used_genes = ranked[:a.auto_top]
    if not used_genes:
        raise RuntimeError("no TF to perturb; check perturbable_tfs.csv")
    bad = [g for g in used_genes if g not in set(ranked)]
    if bad:
        raise ValueError("not perturbable: %s; candidates: %s" %
                         (bad, ", ".join(ranked[:20])))
    if len(used_genes) > 1 and not cache_ok:
        raise RuntimeError("multiple --perturb-genes needs a working Oracle cache; "
                           "drop --no-cache or run one gene at a time")
    log("KO panel: %s" % ", ".join(used_genes), "SIM")

    import celloracle as co
    all_rows = []
    last_oracle = None
    for i, gene in enumerate(used_genes):
        if cache_ok and Path(a.cache_path).is_file():
            current = co.load_hdf5(str(a.cache_path))   # fresh object per gene
        elif i == 0:
            current = oracle
        else:
            raise RuntimeError("no Oracle cache available for gene %d" % i)
        all_rows += perturb_one(current, gene, a, cluster_col,
                                a.outdir, not a.no_plots)
        last_oracle = current

    summary = pd.DataFrame(all_rows)
    summary.to_csv(a.outdir / "ko_screen_summary.csv", index=False)
    pivot = summary.pivot(index="TF", columns="group",
                          values="mean_delta_magnitude")
    pivot.to_csv(a.outdir / "ko_screen_matrix.csv")
    log("KO summary:\n%s" % pivot.round(6).to_string(), "SIM")

    if last_oracle is not None and not a.no_cache:
        try:
            save_cache(last_oracle, a.outdir / "oracle_object.celloracle.oracle")
            log("saved oracle_object.celloracle.oracle", "SIM")
        except Exception as exc:
            log("final cache skipped: %s: %s" % (type(exc).__name__, exc), "SIM")

    report = (
        "celloracle_all_in_one.py %s\n"
        "mode %s\nbase_grn %s\ncluster_column %s\n"
        "bagging_number %d\nlinks_alpha %g\nfit_alpha %g\n"
        "p_values_available %s\np_value_filter_used %s\n"
        "filtered_links %s\nactive_regulators %d\n"
        "custom_regulator %s\ncustom_regulator_info %s\n"
        "perturb_genes %s\nelapsed_sec %.1f\n" %
        (VERSION, mode, a.base_grn, cluster_col,
         a.bagging_number, a.links_alpha, a.fit_alpha,
         meta["has_p"], meta["p_filter_used"], str(meta["n_filtered"]),
         meta["n_active"], a.custom_regulator or "none",
         json.dumps(custom_info, ensure_ascii=False) if custom_info else "none",
         ",".join(used_genes), time.time() - t0))
    (a.outdir / "run_summary.txt").write_text(report, encoding="utf-8")
    log("DONE in %.1f s -> %s" % (time.time() - t0, a.outdir))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
