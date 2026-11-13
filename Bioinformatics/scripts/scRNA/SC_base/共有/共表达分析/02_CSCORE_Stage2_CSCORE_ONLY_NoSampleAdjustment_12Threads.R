# ==============================================================================
# STAGE 2 | 正式 target-centered CS-CORE 共表达
# ==============================================================================
# 使用方法：
#   1) 必须先运行 Stage 1；
#   2) Stage 2 自动读取 Stage1_QC_Result.rds；
#   3) 若目标基因 Detection >=5%，按常规基因集运行；
#   4) 若目标基因 0<Detection<5%，仅强制加入目标基因；其他基因仍要求 Detection >=5%；
#   5) 不做 SampleID/批次校正（covariates = NULL），只运行原始 CS-CORE；
#   6) 不做 ORA/GSEA，只输出共表达基因。
# ==============================================================================

rm(list=ls()); gc(); options(stringsAsFactors=FALSE, warn=1)

# ============================== 只改这里 ======================================
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE271959(scRNA)Object1/共表达基因分析"
target_gene <- "UBE2J2"
cell_type <- "NPC-like"
# ============================================================================

safe_tag <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
analysis_dir <- file.path(work_dir,"CSCORE_TwoStage",safe_tag(cell_type),safe_tag(target_gene))
manifest_file <- file.path(analysis_dir,"Stage1_QC","Stage1_QC_Result.rds")
if (!file.exists(manifest_file)) stop("Stage 1 result not found:\n",manifest_file,"\nPlease run Stage 1 first.")
manifest <- readRDS(manifest_file)
if (!identical(manifest$target_gene,target_gene)) stop("Target gene does not match Stage 1: ",manifest$target_gene)
if (!identical(manifest$cell_type,cell_type)) stop("Cell type does not match Stage 1: ",manifest$cell_type)
if (!isTRUE(manifest$stage2_allowed)) stop("Stage 1 does not recommend Stage 2: ",manifest$advice)

# Stage 1 is authoritative from here onward.
rds_file <- manifest$rds_file
out_dir <- manifest$stage2_dir; plot_dir <- file.path(out_dir,"Plots")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
celltype_col <- manifest$celltype_col; sample_col <- manifest$sample_col
detection_cutoff_pct <- manifest$detection_cutoff_pct
final_genes <- manifest$selected_genes; validation_genes <- manifest$validation_genes; gene_selection <- manifest$gene_qc
validation_n_top <- manifest$validation_n_top
validation_est_tol <- manifest$validation_est_tol; validation_test_tol <- manifest$validation_test_tol; validation_p_tol <- manifest$validation_p_tol
irls_n_iter <- manifest$irls_n_iter; irls_eps <- manifest$irls_eps; irls_verbose <- manifest$irls_verbose
auto_install_cscore <- manifest$auto_install_cscore; fdr_cutoff <- manifest$fdr_cutoff; strong_positive_estimate <- manifest$strong_positive_estimate
target_omp_threads <- 12L; blas_threads <- 12L

logical_cores <- parallel::detectCores(logical=TRUE); physical_cores <- parallel::detectCores(logical=FALSE)
if (is.na(logical_cores)) logical_cores <- 8L
if (is.na(physical_cores)) physical_cores <- max(1L,floor(logical_cores/2L))
Sys.setenv(OMP_NUM_THREADS=target_omp_threads,OPENBLAS_NUM_THREADS=blas_threads,MKL_NUM_THREADS=blas_threads,
  VECLIB_MAXIMUM_THREADS=blas_threads,NUMEXPR_NUM_THREADS=blas_threads)

# 1. Packages
# ==============================================================================

required_pkgs <- c("Seurat", "Matrix", "ggplot2", "Rcpp", "RcppArmadillo")
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop( "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before running this script." )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})
if (!requireNamespace("CSCORE", quietly = TRUE)) {
  if (!auto_install_cscore) {
    stop( "CSCORE is not installed. Install with:\n", "  install.packages('remotes')\n", "  remotes::install_github('ChangSuBiostats/CS-CORE')" )
  }
  message("[SETUP] CSCORE not found. Installing from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github( "ChangSuBiostats/CS-CORE", upgrade = "never", dependencies = NA )
}
cscore_version <- packageVersion("CSCORE")
if (cscore_version < package_version("1.0.1")) {
  stop( "CSCORE version ", as.character(cscore_version), " is too old.\n", "Please update with:\n", "  remotes::install_github('ChangSuBiostats/CS-CORE')"
  )
}
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(blas_threads), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(target_omp_threads), silent = TRUE)
}

# ==============================================================================
# 2. Logging and helpers
# ==============================================================================

log_file <- file.path(out_dir, "00_CSCORE_TargetOnly_Run_Log.txt")
sink(log_file, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)
cat("======================================================================\n")
cat(target_gene," | ",cell_type," | FAST target-only CS-CORE co-expression only\n",sep="")
cat("======================================================================\n")
cat("Start time             : ", format(Sys.time()), "\n", sep = "")
cat("CSCORE version         : ", as.character(cscore_version), "\n", sep = "")
cat("Logical CPUs           : ", logical_cores, "\n", sep = "")
cat("Physical CPUs          : ", physical_cores, "\n", sep = "")
cat("Target OpenMP threads  : ", target_omp_threads, "\n", sep = "")
cat("BLAS thread cap        : ", blas_threads, "\n", sep = "")
cat("Sample adjustment      : NONE (original CS-CORE; covariates = NULL)\n")
cat("Primary gene filter    : detection >= ",detection_cutoff_pct,"% of ",cell_type," cells\n",sep="")
cat("Target forced into set : ",manifest$target_forced,"\n",sep="")
cat("Validation subset      : top ", validation_n_top, " + target\n", sep = "")
cat("======================================================================\n\n")
save_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
safe_fdr_log <- function(x) {
  -log10(pmax(x, .Machine$double.xmin))
}
sparse_mean_relative_expression <- function(counts_mat, seq_depth) {
  inv_depth <- 1 / pmax(seq_depth, 1)
  scaled <- counts_mat %*% Matrix::Diagonal(x = inv_depth)
  out <- Matrix::rowMeans(scaled)
  rm(scaled)
  gc(verbose = FALSE)
  out
}

# ==============================================================================
# 3. Compile the target-only WLS covariance helper
# ==============================================================================

compile_target_wls_cpp <- function() {
  cpp_code <- '
  #include <RcppArmadillo.h>
  #include <algorithm>
  #ifdef _OPENMP
  #include <omp.h>
  #endif

  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::plugins(cpp14)]]
  // [[Rcpp::plugins(openmp)]]

  // Target-only version of the k=1 WLS covariance calculation used by
  // CS-CORE. "denom" is the matrix whose reciprocal is the gene-specific WLS
  // weight matrix (w_cov for covariance estimation; Sigma for testing).
  //
  // [[Rcpp::export]]
  Rcpp::List target_wls_cov_inverse_denom_cpp(
      const arma::vec& d,
      const arma::mat& X_centered,
      const arma::mat& denom,
      const int target_idx,
      const int n_threads) {

    const arma::uword n = X_centered.n_rows;
    const arma::uword p = X_centered.n_cols;

    if (d.n_elem != n) Rcpp::stop("d length does not match cell number");
    if (denom.n_rows != n || denom.n_cols != p) {
      Rcpp::stop("denom dimension mismatch");
    }
    if (target_idx < 0 || target_idx >= static_cast<int>(p)) {
      Rcpp::stop("target_idx out of range");
    }

    arma::vec coef(p, arma::fill::zeros);
    arma::vec test_stat(p, arma::fill::zeros);

    #ifdef _OPENMP
    omp_set_num_threads(std::max(1, n_threads));
    #pragma omp parallel for schedule(static)
    #endif
    for (int jj = 0; jj < static_cast<int>(p); ++jj) {
      double numerator = 0.0;
      double denominator = 0.0;

      for (arma::uword i = 0; i < n; ++i) {
        const double pair_weight =
          1.0 / (denom(i, target_idx) * denom(i, jj));
        const double di = d(i);
        const double xprod =
          X_centered(i, target_idx) * X_centered(i, jj);

        numerator   += di * pair_weight * xprod;
        denominator += di * di * pair_weight;
      }

      coef(jj) = numerator / denominator;
      test_stat(jj) = coef(jj) * std::sqrt(denominator);
    }

    return Rcpp::List::create(
      Rcpp::Named("coef") = coef,
      Rcpp::Named("test_stat") = test_stat
    );
  }

  // [[Rcpp::export]]
  bool target_cscore_openmp_enabled() {
    #ifdef _OPENMP
      return true;
    #else
      return false;
    #endif
  }
  '
  tryCatch( {
      Rcpp::sourceCpp( code = cpp_code, env = globalenv(), rebuild = TRUE, verbose = FALSE )
      TRUE
    }, error = function(e) {
      cat("\n[ERROR] Runtime C++ compilation failed:\n", e$message, "\n", sep = "")
      FALSE
    } )
}
cat("[0] Compiling optimized target-only OpenMP helper (C++14 + OpenMP)...\n")
compile_ok <- compile_target_wls_cpp()
if (!compile_ok) {
  stop( "Could not compile the target-only C++ helper.\n", "Because the full 5001-gene official covariance calculation is too slow, ",
    "this script will NOT silently fall back to the full-network method.\n",
    "On WSL, confirm that build-essential, Rcpp and RcppArmadillo are installed, and that the compiler supports C++14." )
}
openmp_enabled <- isTRUE(target_cscore_openmp_enabled())
cat("OpenMP enabled         : ", openmp_enabled, "\n", sep = "")
if (!openmp_enabled) {
  warning( "OpenMP was not enabled by the compiler. The target-only implementation ",
    "will still be much faster than the full network, but it will use one CPU ", "thread for the covariance loops." )
}

# ==============================================================================
# 4. Exact target-only CS-CORE implementation
# ==============================================================================

run_target_cscore_exact <- function(sparse_gene_by_cell, genes, seq_depth, target_gene, n_threads, label = "TargetOnly") {
  if (!identical(rownames(sparse_gene_by_cell), genes)) {
    stop("Gene order mismatch in ", label)
  }
  cat("\n------------------------------------------------------------------\n")
  cat("[", label, "] Preparing dense selected matrix...\n", sep = "")
  cat("Cells: ", ncol(sparse_gene_by_cell), "\n", sep = "")
  cat("Genes: ", nrow(sparse_gene_by_cell), "\n", sep = "")
  dense_gib <- as.numeric(ncol(sparse_gene_by_cell)) * as.numeric(nrow(sparse_gene_by_cell)) * 8 / 1024^3
  cat("Approx. dense X size: ", sprintf("%.2f GiB", dense_gib), "\n", sep = "")

  # CS-CORE expects cells x genes.
  X <- t(as.matrix(sparse_gene_by_cell))
  storage.mode(X) <- "double"
  colnames(X) <- genes
  n <- nrow(X)
  p <- ncol(X)
  target_idx <- match(target_gene, genes)
  if (is.na(target_idx)) stop("Target gene missing in ", label)
  if (length(seq_depth) != n) {
    stop("seq_depth length does not match cell number in ", label)
  }

  # Original CS-CORE: no SampleID/batch covariates.
  # Equivalent to CSCORE_IRLS_cpp(..., covariates = NULL).
  seq_depth_sq <- seq_depth^2
  D_mu <- matrix(seq_depth, ncol = 1)
  D_sigma2 <- matrix(seq_depth_sq, ncol = 1)
  D_sigma <- matrix(seq_depth_sq, ncol = 1)
  irls_start <- proc.time()[["elapsed"]]
  cat("[", label, "] Starting exact CS-CORE IRLS mean/variance stage at ", format(Sys.time()), "\n", sep = "")

  # ---- Initial mean / variance estimates: mirrors official source ----
  W <- matrix(1, nrow = n, ncol = p)
  mu_beta <- CSCORE::WLS_mean(D_mu, X, W)
  rm(W)
  gc(verbose = FALSE)
  M <- D_mu %*% mu_beta
  mu <- as.numeric(mu_beta[1, ])
  X_centered <- X - M
  response_var <- X_centered^2 - M
  W <- matrix(1, nrow = n, ncol = p)
  sigma2_beta <- CSCORE::WLS_mean(D_sigma2, response_var, W)
  rm(W, response_var)
  gc(verbose = FALSE)
  sigma2 <- as.numeric(sigma2_beta[1, ])
  theta <- mu^2 / sigma2
  iter <- 0L
  delta <- Inf
  while (delta > irls_eps && iter <= irls_n_iter) {
    theta_prev <- theta
    theta_pos <- theta[is.finite(theta) & theta > 0]
    if (length(theta_pos) == 0) {
      stop("No positive finite theta values during IRLS in ", label)
    }
    theta_median <- median(theta_pos)
    mu_sq_over_theta <- mu^2 / theta_median

    # weight = M + seq_depth_sq outer mu_sq_over_theta
    weight <- M + outer(seq_depth_sq, mu_sq_over_theta)
    weight[!is.finite(weight) | weight <= 0] <- 1
    W_mu <- 1 / weight
    rm(weight)
    mu_beta <- CSCORE::WLS_mean(D_mu, X, W_mu)
    rm(W_mu)
    gc(verbose = FALSE)
    mu <- as.numeric(mu_beta[1, ])
    M <- D_mu %*% mu_beta
    X_centered <- X - M
    h <- (M^2 / theta_median + M)^2
    h[!is.finite(h) | h <= 0] <- 1
    W_var <- 1 / h
    rm(h)
    response_var <- X_centered^2 - M
    sigma2_beta <- CSCORE::WLS_mean(D_sigma2, response_var, W_var)
    rm(W_var, response_var)
    gc(verbose = FALSE)
    sigma2 <- as.numeric(sigma2_beta[1, ])
    theta <- mu^2 / sigma2
    valid <- is.finite(theta) & is.finite(theta_prev) & theta > 0 & theta_prev > 0
    if (!any(valid)) {
      stop("No valid theta ratios during IRLS in ", label)
    }
    ratio <- log(theta[valid] / theta_prev[valid])
    delta_vec <- abs(ratio)
    delta <- max(delta_vec)
    iter <- iter + 1L
    if (irls_verbose) {
      cat("[", label, "] IRLS iter ", iter, " delta=", signif(delta, 5), "\n", sep = "")
    }
  }
  if (iter > irls_n_iter && delta > irls_eps) {
    warning( "[", label, "] IRLS did not converge after ", irls_n_iter, " iterations. Final delta=", delta )
  } else {
    cat("[", label, "] IRLS converged after ", iter, " iterations; delta=", signif(delta, 5), "\n", sep = "")
  }
  irls_elapsed <- proc.time()[["elapsed"]] - irls_start
  cat("[", label, "] IRLS elapsed: ", sprintf("%.1f", irls_elapsed), " sec\n", sep = "")

  # ---- Target-only covariance estimation: same equations, only target-gene pairs ----
  theta_pos <- theta[is.finite(theta) & theta > 0]
  if (length(theta_pos) == 0) stop("No positive theta values after IRLS.")
  theta_median <- median(theta_pos)
  mu_sq_over_theta <- mu^2 / theta_median
  w_cov_denom <- M + outer(seq_depth_sq, mu_sq_over_theta)
  w_cov_denom[!is.finite(w_cov_denom) | w_cov_denom <= 0] <- 1
  X_centered <- X - M
  d_sigma <- as.numeric(D_sigma[, 1])
  cat("[", label, "] Target covariance pass 1/2 using ", n_threads, " OpenMP threads...\n", sep = "")
  cov_start <- proc.time()[["elapsed"]]
  cov_fit <- target_wls_cov_inverse_denom_cpp( d = d_sigma, X_centered = X_centered, denom = w_cov_denom, target_idx = target_idx - 1L,
    n_threads = n_threads )
  covar_target <- as.numeric(cov_fit$coef)
  rm(cov_fit, w_cov_denom)
  gc(verbose = FALSE)

  # Statistical-testing weights from the official source:
  # Sigma = M + seq_depth_sq outer sigma2; ele_inv_Sigma = 1 / Sigma.
  Sigma <- M + outer(seq_depth_sq, sigma2)
  cat("[", label, "] Target covariance pass 2/2 (test statistic)...\n", sep = "")
  test_fit <- target_wls_cov_inverse_denom_cpp( d = d_sigma, X_centered = X_centered, denom = Sigma, target_idx = target_idx - 1L, n_threads = n_threads )
  test_stat <- as.numeric(test_fit$test_stat)
  rm(test_fit, Sigma, X_centered, M, X)
  gc(verbose = FALSE)
  cov_elapsed <- proc.time()[["elapsed"]] - cov_start
  cat("[", label, "] Target-only covariance/testing elapsed: ", sprintf("%.1f", cov_elapsed), " sec\n", sep = "")

  # Same normalization used by official CS-CORE.
  sigma <- sqrt(pmax(sigma2, 0))
  raw_est <- covar_target / (sigma[target_idx] * sigma)
  p_value <- 2 * pnorm(-abs(test_stat))

  # Approximate the package post-processing for the target row. Invalid or
  # non-positive biological variances imply insufficient biological variation.
  invalid_variance <- !is.finite(sigma2) | sigma2 <= 0
  processed_est <- raw_est
  if (invalid_variance[target_idx]) {
    processed_est[] <- 0
  } else {
    processed_est[invalid_variance] <- 0
  }
  processed_est[processed_est > 1] <- 1
  processed_est[processed_est < -1] <- -1
  processed_est[target_idx] <- 1
  names(raw_est) <- genes
  names(processed_est) <- genes
  names(test_stat) <- genes
  names(p_value) <- genes
  names(sigma2) <- genes
  target_df <- data.frame( Gene = genes, CSCORE_Estimate = processed_est, Raw_CSCORE_Estimate = raw_est, Test_Statistic = test_stat, P_value = p_value,
    Biological_Variance = sigma2, Invalid_Biological_Variance = invalid_variance, stringsAsFactors = FALSE )
  target_df <- target_df[target_df$Gene != target_gene, , drop = FALSE]
  target_df$FDR <- p.adjust(target_df$P_value, method = "BH")
  target_df$Direction <- ifelse( target_df$CSCORE_Estimate > 0, "Positive", ifelse(target_df$CSCORE_Estimate < 0, "Negative", "Zero") )
  target_df$Significant_FDR005 <- !is.na(target_df$FDR) & target_df$FDR < fdr_cutoff
  total_elapsed <- irls_elapsed + cov_elapsed
  list( target = target_df, raw_estimate = raw_est, test_stat = test_stat, p_value = p_value, sigma2 = sigma2,
    target_variance_valid = !invalid_variance[target_idx], invalid_variance_genes = genes[invalid_variance], irls_iterations = iter, irls_delta = delta,
    irls_seconds = irls_elapsed, covariance_seconds = cov_elapsed, total_seconds = total_elapsed )
}

# ==============================================================================
# 5. Validation against official CSCORE_IRLS_cpp
# ==============================================================================

validate_against_official <- function(validation_sparse, genes, seq_depth, target_gene) {
  cat("\n======================================================================\n")
  cat("[VALIDATION] Comparing FAST target-only vs official full CS-CORE\n")
  cat("======================================================================\n")
  fast_fit <- run_target_cscore_exact(sparse_gene_by_cell = validation_sparse, genes = genes, seq_depth = seq_depth,
    target_gene = target_gene, n_threads = target_omp_threads, label = "Validation_FAST")
  Xv <- t(as.matrix(validation_sparse))
  storage.mode(Xv) <- "double"
  colnames(Xv) <- genes
  irls_par <- list(n_iter = irls_n_iter, eps = irls_eps, verbose = FALSE, conv = "max")
  cat("\n[Validation_OFFICIAL] Starting official full CSCORE_IRLS_cpp (covariates=NULL)...\n")
  official_start <- proc.time()[["elapsed"]]
  official_fit <- CSCORE::CSCORE_IRLS_cpp(X = Xv, seq_depth = seq_depth, covariates = NULL,
    post_process = FALSE, IRLS_par = irls_par, return_all = FALSE)
  official_seconds <- proc.time()[["elapsed"]] - official_start
  cat("[Validation_OFFICIAL] Finished in ", sprintf("%.1f", official_seconds), " sec\n", sep = "")
  target_idx <- match(target_gene, genes)
  official_est <- as.numeric(official_fit$est[target_idx, ])
  official_test <- as.numeric(official_fit$test_stat[target_idx, ])
  official_p <- as.numeric(official_fit$p_value[target_idx, ])
  fast_est <- as.numeric(fast_fit$raw_estimate[genes])
  fast_test <- as.numeric(fast_fit$test_stat[genes])
  fast_p <- as.numeric(fast_fit$p_value[genes])
  keep <- genes != target_gene
  validation_table <- data.frame( Gene = genes[keep], FAST_Estimate = fast_est[keep], Official_Estimate = official_est[keep],
    Delta_Estimate = fast_est[keep] - official_est[keep], FAST_TestStat = fast_test[keep], Official_TestStat = official_test[keep],
    Delta_TestStat = fast_test[keep] - official_test[keep], FAST_P = fast_p[keep], Official_P = official_p[keep],
    Delta_P = fast_p[keep] - official_p[keep], stringsAsFactors = FALSE )
  shared_finite_est <- is.finite(validation_table$FAST_Estimate) & is.finite(validation_table$Official_Estimate)
  shared_finite_test <- is.finite(validation_table$FAST_TestStat) & is.finite(validation_table$Official_TestStat)
  shared_finite_p <- is.finite(validation_table$FAST_P) & is.finite(validation_table$Official_P)
  max_est_delta <- if (any(shared_finite_est)) {
    max(abs(validation_table$Delta_Estimate[shared_finite_est]))
  } else Inf
  max_test_delta <- if (any(shared_finite_test)) {
    max(abs(validation_table$Delta_TestStat[shared_finite_test]))
  } else Inf
  max_p_delta <- if (any(shared_finite_p)) {
    max(abs(validation_table$Delta_P[shared_finite_p]))
  } else Inf
  finite_pattern_est_ok <- identical( is.finite(fast_est[keep]), is.finite(official_est[keep]) )
  finite_pattern_test_ok <- identical( is.finite(fast_test[keep]), is.finite(official_test[keep]) )
  validation_pass <- finite_pattern_est_ok && finite_pattern_test_ok && is.finite(max_est_delta) && max_est_delta <= validation_est_tol &&
    is.finite(max_test_delta) && max_test_delta <= validation_test_tol && is.finite(max_p_delta) && max_p_delta <= validation_p_tol
  summary_df <- data.frame( Metric = c( "ValidationGenes", "FAST_TotalSeconds", "Official_TotalSeconds", "FinitePatternEstimateMatch",
      "FinitePatternTestStatMatch", "MaxAbsDeltaEstimate", "MaxAbsDeltaTestStat", "MaxAbsDeltaP", "EstimateTolerance", "TestStatTolerance", "PTolerance",
      "ValidationPass" ), Value = c( sum(keep), fast_fit$total_seconds, official_seconds, finite_pattern_est_ok, finite_pattern_test_ok, max_est_delta,
      max_test_delta, max_p_delta, validation_est_tol, validation_test_tol, validation_p_tol, validation_pass ), stringsAsFactors = FALSE )
  save_csv( validation_table, file.path(out_dir, "Validation_FAST_vs_Official_PerGene.csv") )
  save_csv( summary_df, file.path(out_dir, "Validation_FAST_vs_Official_Summary.csv") )
  cat("\nValidation max |delta estimate| : ", format(max_est_delta, scientific = TRUE), "\n", sep = "")
  cat("Validation max |delta test|     : ", format(max_test_delta, scientific = TRUE), "\n", sep = "")
  cat("Validation max |delta p|        : ", format(max_p_delta, scientific = TRUE), "\n", sep = "")
  cat("Validation pass                 : ", validation_pass, "\n", sep = "")
  rm(Xv, official_fit)
  gc(verbose = FALSE)
  if (!validation_pass) {
    stop( "FAST target-only implementation did NOT match official CS-CORE within ", "the predefined tolerances.\n",
      "Final Detection>=5% analysis was NOT started.\n", "Inspect Validation_FAST_vs_Official_Summary.csv and PerGene.csv." )
  }
  invisible(list( pass = validation_pass, fast_fit = fast_fit, summary = summary_df ))
}

# ==============================================================================

# ==============================================================================
# 6. Reconstruct the exact Stage 1 analysis set
# ==============================================================================
cat("\n[Stage2-1] Loading the Stage 1-locked cells and genes...\n")
if (!file.exists(rds_file)) stop("RDS not found: ", rds_file)
obj <- readRDS(rds_file)
if (!celltype_col %in% colnames(obj@meta.data)) stop("Metadata column not found: ", celltype_col)
if (!sample_col %in% colnames(obj@meta.data)) stop("Metadata column not found: ", sample_col)
if (!"RNA" %in% Assays(obj)) stop("RNA assay not found.")
Idents(obj) <- celltype_col; ac <- subset(obj, idents=cell_type); rm(obj); gc(verbose=FALSE)
missing_cells <- setdiff(manifest$included_cells, colnames(ac))
if (length(missing_cells)) stop(length(missing_cells), " Stage 1 cells are missing from the current RDS. Do not continue.")
ac <- subset(ac, cells=manifest$included_cells); DefaultAssay(ac) <- "RNA"
if (sum(grepl("^counts($|\\.)", Layers(ac[["RNA"]]))) > 1) ac <- JoinLayers(ac, assay="RNA")
if (!"counts" %in% Layers(ac[["RNA"]])) stop("RNA counts layer not found. CS-CORE requires raw UMI counts.")
counts <- LayerData(ac, assay="RNA", layer="counts")
if (!inherits(counts,"sparseMatrix")) counts <- Matrix::Matrix(counts,sparse=TRUE)
if (!all(final_genes %in% rownames(counts))) stop("Some Stage 1-selected genes are missing from the current RDS.")
if (!identical(colnames(counts), manifest$included_cells)) counts <- counts[, manifest$included_cells, drop=FALSE]
ac@meta.data[[sample_col]] <- droplevels(factor(ac@meta.data[[sample_col]]))
cat("Locked cells : ",ncol(counts),"\nLocked genes : ",length(final_genes),"\nSamples      : ",nlevels(ac@meta.data[[sample_col]]),"\n",sep="")

# ==============================================================================
# 7. Sequencing depth and target QC (sample labels only for QC; NOT adjusted in CS-CORE)
# ==============================================================================
seq_depth <- as.numeric(Matrix::colSums(counts)); names(seq_depth) <- colnames(counts)
if (any(!is.finite(seq_depth)) || any(seq_depth<=0)) stop("Invalid/zero sequencing depth detected.")
target_raw <- as.numeric(counts[target_gene,,drop=TRUE]); sample_names <- levels(ac@meta.data[[sample_col]])
sample_qc <- do.call(rbind,lapply(sample_names,function(s){ idx<-ac@meta.data[[sample_col]]==s; x<-target_raw[idx]
  data.frame(Sample=s,Cells=sum(idx),Target_Positive_Cells=sum(x>0),Target_Detection_Percent=mean(x>0)*100,Target_Mean_UMI=mean(x),Target_Total_UMI=sum(x),stringsAsFactors=FALSE)}))
save_csv(sample_qc,file.path(out_dir,"00_Stage2_Cell_Sample_QC.csv"))
parameters <- data.frame(Parameter=c("Stage1Manifest","TargetGene","CellType","LockedCells","LockedGenesIncludingTarget","TargetDetectionPercent_Stage1",
  "GeneDetectionCutoffPercent","SampleAdjustmentMode","ValidationGenesIncludingTarget","CSCORE_Version","TargetOpenMPThreads","BLASThreadCap",
  "FDR_Method","FDR_Cutoff","StrongPositiveEstimateCutoff","Stage1_MinCellsPerSample","Stage1_ExcludedSamples","TargetForcedBecauseBelowDetectionCutoff"),
  Value=c(manifest_file,target_gene,cell_type,ncol(counts),length(final_genes),manifest$target_detection_pct,detection_cutoff_pct,"none (covariates=NULL)",
  length(validation_genes),as.character(cscore_version),target_omp_threads,blas_threads,"BH",fdr_cutoff,strong_positive_estimate,manifest$min_cells_per_sample,paste(manifest$excluded_samples,collapse=";"),manifest$target_forced),stringsAsFactors=FALSE)
save_csv(parameters,file.path(out_dir,"00_Stage2_Analysis_Parameters.csv"))

# ==============================================================================
# 8. Mandatory validation against official CS-CORE
# ==============================================================================
cat("\n[Stage2-2] Validating FAST target-only implementation...\n")
validation_sparse <- counts[validation_genes,,drop=FALSE]
validation_result <- validate_against_official(validation_sparse=validation_sparse,genes=validation_genes,
  seq_depth=seq_depth,target_gene=target_gene)
rm(validation_sparse,validation_result); gc(verbose=FALSE)
cat("[VALIDATION] PASS. Starting final target-only analysis.\n")

# ==============================================================================
# 9. Final Stage 2 target-centered CS-CORE
# ==============================================================================
final_sparse <- counts[final_genes,,drop=FALSE]; rm(counts,ac); gc(verbose=FALSE)
cat("\n[Stage2-3] Running final ",target_gene,"-centered CS-CORE...\n",sep="")
final_fit <- run_target_cscore_exact(sparse_gene_by_cell=final_sparse,genes=final_genes,seq_depth=seq_depth,
  target_gene=target_gene,n_threads=target_omp_threads,label="Stage2_Final_FAST_TargetCentered")
rm(final_sparse); gc(verbose=FALSE)
if (!final_fit$target_variance_valid) stop("Invalid/non-positive biological variance for ",target_gene,"; results should not be interpreted.")
invalid_df <- data.frame(Gene=final_fit$invalid_variance_genes,IsTarget=final_fit$invalid_variance_genes==target_gene,stringsAsFactors=FALSE)
save_csv(invalid_df,file.path(out_dir,"01_Stage2_Invalid_Biological_Variance_Genes.csv"))

# ==============================================================================
# 10. Export co-expression genes
# ==============================================================================
target_res <- final_fit$target
gene_meta <- gene_selection[match(target_res$Gene,gene_selection$Gene),c("Gene","Mean_Relative_Expression","Mean_UMI","Detection_Percent","Expression_Rank","Pass_Detection"),drop=FALSE]
target_res <- merge(target_res,gene_meta,by="Gene",all.x=TRUE,sort=FALSE)
target_res <- target_res[order(target_res$FDR,-target_res$CSCORE_Estimate,na.last=TRUE),,drop=FALSE]
save_csv(target_res,file.path(out_dir,paste0("02_",safe_tag(target_gene),"_CSCORE_All_Genes.csv")))
positive_sig <- target_res[!is.na(target_res$FDR)&target_res$FDR<fdr_cutoff&target_res$CSCORE_Estimate>0&!target_res$Invalid_Biological_Variance,,drop=FALSE]
positive_sig <- positive_sig[order(-positive_sig$CSCORE_Estimate,positive_sig$FDR),,drop=FALSE]
save_csv(positive_sig,file.path(out_dir,paste0("03_",safe_tag(target_gene),"_CSCORE_Positive_FDR_lt_0.05.csv")))
strong_positive <- positive_sig[positive_sig$CSCORE_Estimate>=strong_positive_estimate,,drop=FALSE]
save_csv(strong_positive,file.path(out_dir,paste0("04_",safe_tag(target_gene),"_CSCORE_Positive_Est_ge_",strong_positive_estimate,"_FDR_lt_",fdr_cutoff,".csv")))
negative_sig <- target_res[!is.na(target_res$FDR)&target_res$FDR<fdr_cutoff&target_res$CSCORE_Estimate<0&!target_res$Invalid_Biological_Variance,,drop=FALSE]
negative_sig <- negative_sig[order(negative_sig$CSCORE_Estimate,negative_sig$FDR),,drop=FALSE]
save_csv(negative_sig,file.path(out_dir,paste0("05_",safe_tag(target_gene),"_CSCORE_Negative_FDR_lt_0.05.csv")))

# ==============================================================================
# 11. Global QC + plots
# ==============================================================================
finite_est <- target_res$CSCORE_Estimate[is.finite(target_res$CSCORE_Estimate)&!target_res$Invalid_Biological_Variance]
qc_summary <- data.frame(Metric=c("TargetVarianceValid","GenesTested","GenesWithValidBiologicalVariance","GenesWithInvalidBiologicalVariance","FiniteValidEstimates",
  "PositiveEstimates","NegativeEstimates","ZeroEstimates","PositivePercent","NegativePercent","MedianEstimate","MeanEstimate","MaxPositiveEstimate",
  "MinNegativeEstimate","FDR_lt_0.05_Positive","FDR_lt_0.05_Negative","Estimate_ge_0.10_and_FDR_lt_0.05","IRLS_Iterations","IRLS_FinalDelta",
  "IRLS_Seconds","TargetCovarianceTesting_Seconds","FinalTotalSeconds"), Value=c(final_fit$target_variance_valid,nrow(target_res),
  sum(!target_res$Invalid_Biological_Variance),sum(target_res$Invalid_Biological_Variance),length(finite_est),sum(finite_est>0),sum(finite_est<0),sum(finite_est==0),
  mean(finite_est>0)*100,mean(finite_est<0)*100,median(finite_est),mean(finite_est),ifelse(any(finite_est>0),max(finite_est[finite_est>0]),NA_real_),
  ifelse(any(finite_est<0),min(finite_est[finite_est<0]),NA_real_),nrow(positive_sig),nrow(negative_sig),nrow(strong_positive),final_fit$irls_iterations,
  final_fit$irls_delta,final_fit$irls_seconds,final_fit$covariance_seconds,final_fit$total_seconds),stringsAsFactors=FALSE)
save_csv(qc_summary,file.path(out_dir,paste0("06_",safe_tag(target_gene),"_CSCORE_Global_QC.csv")))

plot_df <- target_res[is.finite(target_res$CSCORE_Estimate)&!is.na(target_res$FDR)&!target_res$Invalid_Biological_Variance,,drop=FALSE]
p_hist <- ggplot(plot_df,aes(CSCORE_Estimate))+geom_histogram(bins=80,fill="grey75",color="white",linewidth=.2)+geom_vline(xintercept=0,linetype="dashed",linewidth=.5)+
  labs(title=paste0("Stage 2: ",target_gene," CS-CORE Co-expression Distribution"),subtitle=paste0(cell_type," | Background genes Detection >= ",detection_cutoff_pct,"% | target forced: ",manifest$target_forced," | covariates=NULL"),x="CS-CORE estimate",y="Gene count")+
  theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"),plot.subtitle=element_text(hjust=.5))
ggsave(file.path(plot_dir,paste0("01_",safe_tag(target_gene),"_CSCORE_Estimate_Distribution.pdf")),p_hist,width=7,height=5)
plot_df$PlotClass <- "Not significant"; plot_df$PlotClass[plot_df$FDR<fdr_cutoff&plot_df$CSCORE_Estimate>0] <- "Positive FDR < 0.05"
plot_df$PlotClass[plot_df$FDR<fdr_cutoff&plot_df$CSCORE_Estimate<0] <- "Negative FDR < 0.05"; plot_df$MinusLog10FDR <- safe_fdr_log(plot_df$FDR)
label_df <- head(positive_sig,12)
p_vol <- ggplot(plot_df,aes(CSCORE_Estimate,MinusLog10FDR,color=PlotClass))+geom_point(alpha=.55,size=1.2)+geom_vline(xintercept=0,linetype="dashed",linewidth=.4)+
  geom_hline(yintercept=-log10(fdr_cutoff),linetype="dashed",linewidth=.4)+scale_color_manual(values=c("Positive FDR < 0.05"="#D2352C","Negative FDR < 0.05"="#3376B2","Not significant"="grey80"))+
  labs(title=paste0("Stage 2: ",target_gene,"-centric CS-CORE"),x="CS-CORE estimate",y="-Log10(FDR)",color=NULL)+theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"),legend.position="top")
if(nrow(label_df)>0&&requireNamespace("ggrepel",quietly=TRUE)) p_vol <- p_vol+ggrepel::geom_text_repel(data=label_df,aes(CSCORE_Estimate,safe_fdr_log(FDR),label=Gene),inherit.aes=FALSE,size=3,color="black",max.overlaps=Inf)
ggsave(file.path(plot_dir,paste0("02_",safe_tag(target_gene),"_CSCORE_Volcano.pdf")),p_vol,width=7.2,height=5.6)
top_pos <- head(positive_sig,20)
if(nrow(top_pos)>0){ top_pos$Gene<-factor(top_pos$Gene,levels=rev(top_pos$Gene)); p_top<-ggplot(top_pos,aes(CSCORE_Estimate,Gene))+geom_col(width=.72)+
  labs(title=paste0("Top Positive ",target_gene," Co-expression Genes"),subtitle="CS-CORE | FDR < 0.05",x="CS-CORE estimate",y=NULL)+theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"),plot.subtitle=element_text(hjust=.5));
  ggsave(file.path(plot_dir,paste0("03_",safe_tag(target_gene),"_CSCORE_Top_Positive.pdf")),p_top,width=6.5,height=6) }

cat("\n======================================================================\n")
cat("STAGE 2 completed | ",target_gene," CS-CORE co-expression\n",sep="")
cat("======================================================================\n")
cat("Locked cells                : ",length(manifest$included_cells),"\n",sep="")
cat("Locked genes incl. target   : ",length(final_genes),"\n",sep="")
cat("Target forced (<5%)         : ",manifest$target_forced,"\n",sep="")
cat("Sample/batch adjustment     : NONE (covariates=NULL)\n")
cat("Positive FDR < 0.05         : ",nrow(positive_sig),"\n",sep="")
cat("Negative FDR < 0.05         : ",nrow(negative_sig),"\n",sep="")
cat("Positive est>=0.10 FDR<0.05 : ",nrow(strong_positive),"\n",sep="")
cat("OpenMP / BLAS threads       : ",target_omp_threads," / ",blas_threads,"\n",sep="")
cat("Output                      : ",out_dir,"\n",sep="")
cat("======================================================================\n")
rm(final_fit); gc(verbose=FALSE); sink()
