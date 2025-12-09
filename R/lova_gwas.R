#' @title LOVAGWAS: a locus direct effect estimation
#'
#' @description
#' Adjust GWAS SNP effects to estimate the direct genetic effect
#' using the Multivariable Mendelian Randomization Latent Outcome Variable
#' Approach (MVMR-LOVA).
#' Automatically detects column names such as `b`, `beta`, `BETA`, `effect`,
#' `logodds`, and their SE variants, as well as `SNP`, `CHR`, `BP`, `A1`, `A2`.
#' Optionally generates Manhattan plots (original vs adjusted) using qqman.
#'
#' @param gwasY Data frame of outcome GWAS summary statistics.
#'   Must contain SNP IDs and effect estimates (any variant of b/beta/logodds etc.).
#'   For Manhattan plotting, must also contain columns for chromosome (`CHR`)
#'   and base position (`BP`).
#' @param betaX Data frame of exposure SNP effects: first column \code{SNP},
#'   remaining columns are exposures.
#' @param tau_t Numeric vector of causal estimates for exposures
#'   (length must match the number of exposure columns in \code{betaX} minus one).
#' @param pcor Optional correlation matrix with outcome in the first row/column
#'   and exposures following. If \code{NULL}, it is computed from
#'   \code{cbind(b, betaX_exposures)}.
#' @param ny Numeric; outcome GWAS sample size (must be > 0).
#' @param plot_manhattan Logical; if \code{TRUE}, generates Manhattan plots.
#' @param save_path Optional file path to save the plots as PNG; if \code{NULL},
#'   plots are displayed to the active graphics device.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{b_c}{LOVAGWAS SNP effect (direct effect)}
#'     \item{se_c}{Adjusted standard error}
#'     \item{z_c}{Adjusted z-score}
#'     \item{pval_c}{Adjusted p-value}
#'     \item{missing_count}{Number of missing exposure effects replaced by 0}
#'     \item{missing_index}{Indices of missing exposures}
#'   }
#'
#' @details
#' Direct effects are computed as \eqn{b_c = b_Y - B_X \tau_t}.
#' Standard errors use \eqn{se_c^2 = se_Y^2 + ||\tau_t||^2 / n_y - \mathrm{cov\_tmp}/n_y},
#' where \code{cov_tmp} aggregates outcome–exposure and exposure–exposure
#' correlation terms.
#'
#' @importFrom stats cor pchisq na.omit
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#' @import qqman
#'
#' @examples
#' # Simulated example with 100 SNPs and 2 exposures
#' set.seed(1)
#' n_snps <- 100
#' gwasY <- data.frame(
#'   SNP = paste0("rs", 1:n_snps),
#'   CHR = sample(1:22, n_snps, TRUE),
#'   BP  = sample(1e6:2e6, n_snps),
#'   A1  = sample(c("A", "C", "G", "T"), n_snps, TRUE),
#'   A2  = sample(c("A", "C", "G", "T"), n_snps, TRUE),
#'   BETA = rnorm(n_snps, 0, 0.05),
#'   SE   = runif(n_snps, 0.01, 0.03)
#' )
#'
#' betaX <- data.frame(
#'   SNP = gwasY$SNP,
#'   exp1 = rnorm(n_snps, 0, 0.05),
#'   exp2 = rnorm(n_snps, 0, 0.05)
#' )
#'
#' tau_t <- c(0.1, 0.2)
#'
#' res <- lova_gwas(gwasY, betaX, tau_t, ny = 50000)
#' head(res)
#'
#' @export

lova_gwas <- function(gwasY, betaX, tau_t, pcor = NULL, ny,
                      plot_manhattan = FALSE, save_path = NULL) {

  ## ---- Helper: auto-detect column names ------------------------------------
  detect_columns <- function(df) {
    nms <- tolower(names(df))
    snp_candidates <- c("snp", "rsid", "rs_id", "markername", "id")
    chr_candidates <- c("chr", "chrom", "chromosome")
    bp_candidates  <- c("bp", "pos", "position", "basepair")
    a1_candidates  <- c("a1", "allele1", "effect_allele", "ea")
    a2_candidates  <- c("a2", "allele2", "other_allele", "nea")
    b_candidates   <- c("b", "beta", "effect", "logodds", "estimate", "betahat")
    se_candidates  <- c("se", "stderr", "sebeta", "standard_error", "se_", "sigma")
    pick <- function(cands) names(df)[match(TRUE, nms %in% cands, nomatch = 0)]
    list(
      SNP  = pick(snp_candidates),
      CHR  = pick(chr_candidates),
      BP   = pick(bp_candidates),
      A1   = pick(a1_candidates),
      A2   = pick(a2_candidates),
      b    = pick(b_candidates),
      se   = pick(se_candidates)
    )
  }

  ## ---- Input validation -----------------------------------------------------
  if (!is.data.frame(gwasY)) stop("gwasY must be a data frame.")
  if (!is.data.frame(betaX)) stop("betaX must be a data frame.")
  if (!is.numeric(tau_t)) stop("tau_t must be numeric.")
  if (!is.numeric(ny) || ny <= 0) stop("ny must be a positive scalar.")

  # detect all columns
  cols <- detect_columns(gwasY)
  if (is.na(cols$SNP) || cols$SNP == "") stop("No SNP column found (SNP/rsid/id/etc.).")
  if (is.na(cols$b)   || cols$b == "")   stop("No effect column found (b/beta/effect/logodds/etc.).")
  if (is.na(cols$se)  || cols$se == "")  stop("No SE column found (se/SE/stderr/etc.).")

  # optional (for plotting)
  if (plot_manhattan) {
    for (v in c("CHR", "BP", "A1", "A2")) {
      if (is.na(cols[[v]]) || cols[[v]] == "")
        stop(sprintf("Column for %s missing: required for Manhattan plot.", v))
    }
  }

  ## ---- Standardize GWAS columns --------------------------------------------
  cols_to_keep <- unique(c(
    cols$SNP, cols$CHR, cols$BP, cols$A1, cols$A2, cols$b, cols$se
  ))
  cols_to_keep <- cols_to_keep[cols_to_keep %in% names(gwasY)]
  if ("pval" %in% names(gwasY)) cols_to_keep <- c(cols_to_keep, "pval")

  gwasY <- gwasY[, cols_to_keep, drop = FALSE]

  names(gwasY)[names(gwasY) == cols$SNP] <- "SNP"
  if (!is.na(cols$CHR)) names(gwasY)[names(gwasY) == cols$CHR] <- "CHR"
  if (!is.na(cols$BP))  names(gwasY)[names(gwasY) == cols$BP]  <- "BP"
  if (!is.na(cols$A1))  names(gwasY)[names(gwasY) == cols$A1]  <- "A1"
  if (!is.na(cols$A2))  names(gwasY)[names(gwasY) == cols$A2]  <- "A2"
  names(gwasY)[names(gwasY) == cols$b]  <- "b"
  names(gwasY)[names(gwasY) == cols$se] <- "se"


  ## ---- Basic consistency checks --------------------------------------------
  if (!"SNP" %in% names(betaX)) stop("betaX must contain a 'SNP' column.")
  exp_cols <- setdiff(names(betaX), "SNP")
  if (length(tau_t) != length(exp_cols))
    stop("tau_t length must equal number of exposure columns in betaX (excluding SNP).")
  if (!identical(gwasY$SNP, betaX$SNP))
    stop("SNP order mismatch: gwasY$SNP and betaX$SNP must be identical and aligned.")

  ## ---- Handle missing exposures --------------------------------------------
  betaX$missing_count <- rowSums(is.na(betaX[, exp_cols, drop = FALSE]))
  has_miss <- betaX$missing_count > 0L
  betaX$missing_index <- NA_character_
  if (any(has_miss)) {
    idx_mat <- is.na(betaX[has_miss, exp_cols, drop = FALSE])
    betaX$missing_index[has_miss] <- apply(idx_mat, 1L, function(z) paste(which(z), collapse = ", "))
  }
  for (nm in exp_cols) {
    x <- betaX[[nm]]; x[is.na(x)] <- 0; betaX[[nm]] <- x
  }

  ## ---- Extract matrices/vectors --------------------------------------------
  by  <- as.vector(gwasY$b)
  sey <- as.vector(gwasY$se)
  bX  <- as.matrix(betaX[, exp_cols, drop = FALSE])
  nexp <- length(tau_t)

  ## ---- Compute correlation matrix ------------------------------------------
  if (is.null(pcor)) {
    data_mat <- cbind(by, bX)
    pcor <- stats::cor(data_mat)
  }

  ## ---- Covariance term -----------------------------------------------------
  cov_tmp <- 0
  for (i in 1:nexp) {
    cov_tmp <- cov_tmp - 2 * tau_t[i] * pcor[1, 1 + i]
    for (j in 1:nexp) {
      if (i != j) cov_tmp <- cov_tmp + 2 * tau_t[i] * tau_t[j] * pcor[1 + i, 1 + j]
    }
  }

  ## ---- Adjusted effects ----------------------------------------------------
  b_c  <- as.numeric(by - bX %*% tau_t)
  se_c <- sqrt(sey^2 + sum(tau_t^2) / ny - cov_tmp / ny)
  z_c  <- b_c / se_c
  pval_c <- stats::pchisq(z_c^2, df = 1, lower.tail = FALSE)
  if (!"pval" %in% names(gwasY))
    gwasY$pval <- stats::pchisq((by / sey)^2, df = 1, lower.tail = FALSE)

  ## ---- Assemble output -----------------------------------------------------
  out <- data.frame(
    gwasY,
    b_c  = b_c,
    se_c = se_c,
    z_c  = z_c,
    pval_c = pval_c,
    missing_count = betaX$missing_count,
    missing_index = betaX$missing_index,
    stringsAsFactors = FALSE
  )

  ## ---- Optional Manhattan plots --------------------------------------------
  if (plot_manhattan) {
    if (!is.null(save_path)) {
      dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
      png(save_path, height = 7, width = 12, units = "in", res = 300)
    }
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1, 2))

    qqman::manhattan(
      out, chr = "CHR", bp = "BP", snp = "SNP", p = "pval",
      genomewideline = -log10(5e-8),
      suggestiveline = -log10(1e-5),
      main = "Manhattan (Original)\nHighlight: significant after adjustment",
      highlight = na.omit(out$SNP[out$pval_c < 5e-8]),
      annotateTop = TRUE,
      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
    )
    qqman::manhattan(
      out, chr = "CHR", bp = "BP", snp = "SNP", p = "pval_c",
      genomewideline = -log10(5e-8),
      suggestiveline = -log10(1e-5),
      main = "Manhattan (Adjusted)\nHighlight: significant before adjustment",
      highlight = na.omit(out$SNP[out$pval < 5e-8]),
      annotateTop = TRUE,
      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
    )
    if (!is.null(save_path)) dev.off()
  }

  return(out)
}
