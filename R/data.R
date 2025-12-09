#' Simulated GWAS summary data for outcome Y
#'
#' @format A data frame with 100 rows and 7 columns:
#' \describe{
#'   \item{SNP}{SNP identifier}
#'   \item{CHR}{Chromosome}
#'   \item{BP}{Base-pair position}
#'   \item{A1, A2}{Effect and other alleles}
#'   \item{BETA}{Estimated effect size}
#'   \item{SE}{Standard error}
#' }
"gwasY"

#' Simulated exposure effect estimates
#'
#' @format A data frame with 100 rows and 3 columns:
#' \describe{
#'   \item{SNP}{SNP identifier}
#'   \item{exp1}{Effect on exposure 1}
#'   \item{exp2}{Effect on exposure 2}
#' }
"betaX"
