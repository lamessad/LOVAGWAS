#' Compute direct genetic effect (MVMR-LOVA–adjusted) GWAS effects with optional Manhattan plot
#'
#' Adjust GWAS SNP effect sizes to estimate the direct genetic effect using
#' the Multivariable Mendelian Randomization Latent Outcome Variable Approach (MVMR-LOVA).
#' Optionally, generate Manhattan plots of the original and adjusted results via **qqman**.
#'
#' @param betaY Numeric vector of outcome SNP effect sizes (length = n_snps).
#' @param betaX Numeric matrix of exposure SNP effect sizes (n_snps x n_exposures).
#' @param betaYse Numeric vector of standard errors for \code{betaY} (length = n_snps).
#' @param tau_t Numeric vector of causal effect estimates for exposures (length = n_exposures).
#' @param pcor Optional numeric correlation matrix with outcome in the first row/column
#'   and exposures following. If \code{NULL}, correlations are computed from \code{cbind(betaY, betaX)}.
#' @param ny Numeric scalar, sample size of the outcome GWAS (must be > 0).
#' @param snp_info Data frame with SNP information, aligned with \code{betaY}.
#'   Must contain columns \code{SNP}, \code{CHR}, \code{BP}, \code{A1}, \code{A2}
#'   if \code{plot_manhattan = TRUE}.
#' @param plot_manhattan Logical; if \code{TRUE}, generate Manhattan plots of original
#'   and LOVA-adjusted p-values (default = \code{FALSE}).
#' @param save_path Optional file path to save the Manhattan plots as PNG. If \code{NULL},
#'   plots are displayed to the active device.
#'
#' @return A data frame with, for each SNP:
#'   \itemize{
#'     \item \code{b_c} – MVMR-LOVA–adjusted effect size (direct genetic effect)
#'     \item \code{se_c} – adjusted standard error
#'     \item \code{z_c} – adjusted z-score
#'     \item \code{pval_c} – adjusted p-value
#'     \item \code{pval} – original p-value
#'     \item all columns from \code{snp_info}
#'   }
#'
#' @importFrom stats cor pchisq
#' @importFrom qqman manhattan
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#' @export
lova_gwas <- function(betaY, betaX, betaYse, tau_t, pcor = NULL, ny, snp_info,
                      plot_manhattan = FALSE, save_path = NULL) {
  # --- Basic checks ---
  if (!is.numeric(betaY) || !is.numeric(betaYse)) stop("betaY and betaYse must be numeric.")
  if (!is.matrix(betaX) || !is.numeric(betaX)) stop("betaX must be a numeric matrix.")
  if (!is.numeric(tau_t)) stop("tau_t must be numeric.")
  if (!is.data.frame(snp_info)) stop("snp_info must be a data frame.")
  if (length(betaY) != nrow(betaX)) stop("betaY length must match number of rows in betaX.")
  if (length(betaY) != nrow(snp_info)) stop("betaY length must match number of SNPs in snp_info.")
  if (length(tau_t) != ncol(betaX)) stop("tau_t length must match number of exposures (columns in betaX).")
  if (!is.numeric(ny) || length(ny) != 1L || ny <= 0) stop("ny must be a positive numeric scalar.")

  # Extra check only if plotting
  if (plot_manhattan) {
    required_cols <- c("SNP","CHR","BP","A1","A2")
    if (!all(required_cols %in% colnames(snp_info))) {
      stop("For Manhattan plotting, snp_info must contain: SNP, CHR, BP, A1, A2.")
    }
  }

  nexp <- length(tau_t)

  # --- Correlation matrix ---
  if (is.null(pcor)) {
    data_mat <- cbind(betaY, betaX)
    pcor <- stats::cor(data_mat)
  }

  # --- Covariance adjustment (loop, base R) ---
  cov_tmp <- 0
  for (i in 1:nexp) {
    cov_tmp <- cov_tmp - 2 * tau_t[i] * pcor[1, 1 + i]
    for (j in 1:nexp) {
      if (i != j) cov_tmp <- cov_tmp + 2 * tau_t[i] * tau_t[j] * pcor[1 + i, 1 + j]
    }
  }

  # --- Adjusted effect sizes & SEs ---
  b_c <- as.numeric(betaY - betaX %*% tau_t)
  se_c <- sqrt(betaYse^2 + sum(tau_t^2) / ny - cov_tmp / ny)

  # --- Assemble output ---
  z_c <- b_c / se_c
  pval_c <- stats::pchisq(z_c^2, df = 1, lower.tail = FALSE)
  pval   <- stats::pchisq((betaY / betaYse)^2, df = 1, lower.tail = FALSE)

  gwas <- data.frame(
    snp_info,
    b_c    = b_c,
    se_c   = se_c,
    z_c    = z_c,
    pval_c = pval_c,
    pval   = pval,
    stringsAsFactors = FALSE
  )

  # --- Optional Manhattan plots ---
  if (plot_manhattan) {
    if (!is.null(save_path)) {
      dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
      png(save_path, height = 7, width = 12, units = "in", res = 300)
    }

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1, 2))

    # Original p-values
    qqman::manhattan(gwas, chr = "CHR", bp = "BP", snp = "SNP", p = "pval",
                     genomewideline = -log10(5e-8),
                     suggestiveline = -log10(1e-5),
                     main = "Manhattan (Original)\nhighlighted are significicant after correction", highlight = c(gwas$SNP[gwas$pval_c<5*10**-8]),
                     annotateTop = TRUE)

    # Adjusted p-values
    qqman::manhattan(gwas, chr = "CHR", bp = "BP", snp = "SNP", p = "pval_c",
                     genomewideline = -log10(5e-8),
                     suggestiveline = -log10(1e-5),
                     main = "Manhattan (adjusted)\nhighlighted are significant in the before correction", highlight = c(gwas$SNP[gwas$pval<5*10**-8]),
                     annotateTop = TRUE)

    if (!is.null(save_path)) dev.off()
  }

  return(gwas)
}
