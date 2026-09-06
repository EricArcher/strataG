#' @title Split Alleles For Diploid Data
#' @description Split loci stored in one column to two columns for each allele 
#'   in a matrix of diploid data.
#'   
#' @param x a matrix or data.frame containing diploid data. Every column 
#'   represents one locus with two alleles.
#' @param sep separator used between alleles of a locus. If \code{NULL}, then 
#'  alleles should be of equal length (e.g., 145095 = 145 and 095, or 
#'  AG = A and G).
#' 
#' @details
#' Separators are interpreted literally. Non-missing separated genotypes must
#' contain exactly two non-empty alleles. Unseparated genotypes must have even
#' width after removing spaces. In this legacy unseparated encoding, zero-valued
#' allele codes and the string NA denote missing alleles; separated zero alleles
#' are retained. Actual NA genotypes return two missing alleles.
#'
#' @return matrix with alleles for each locus in one column split into 
#'   separate columns.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' # A sample SNP data set with no separators between nucleotides in a genotype
#' snps <- do.call(cbind, lapply(1:3, function(i) {
#'   a1 <- sample(c("A", "G"), 10, rep = TRUE)
#'   a2 <- sample(c("A", "G"), 10, rep = TRUE)
#'   paste(a1, a2, sep = "")
#' }))
#' colnames(snps) <- paste("Loc", LETTERS[1:3], sep = "_")
#' snps
#' alleleSplit(snps)
#' 
#' # A sample microsatellie data set with alleles separated by "/"
#' alleles <- seq(100, 150, 2)
#' msats <- do.call(cbind, lapply(1:3, function(i) {
#'   a1 <- sample(alleles, 10, rep = TRUE)
#'   a2 <- sample(alleles, 10, rep = TRUE)
#'   paste(a1, "/", a2, sep = "")
#' }))
#' colnames(msats) <- paste("Loc", LETTERS[1:3], sep = "_")
#' msats
#' alleleSplit(msats, sep = "/")
#' 
#' @export
#' 
alleleSplit <- function(x, sep = NULL) {
  if(!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data.frame.", call. = FALSE)
  }
  if(!is.null(sep)) {
    if(!is.character(sep) || length(sep) != 1L || is.na(sep)) {
      stop("'sep' must be NULL or a single non-missing string.", call. = FALSE)
    }
    if(sep == "") sep <- NULL
  }
  locus.names <- colnames(x)
  if(is.null(locus.names)) locus.names <- if(ncol(x)) paste0("Locus", seq_len(ncol(x))) else character()
  result <- matrix(NA_character_, nrow(x), 2L * ncol(x),
    dimnames = list(rownames(x),
      if(ncol(x)) paste(rep(locus.names, each = 2L), rep(1:2, ncol(x)), sep = ".") else character()))
  missing.code <- function(a) {
    numeric <- suppressWarnings(as.numeric(a))
    if(a == "NA" || a == "" || (!is.na(numeric) && numeric == 0)) NA_character_ else a
  }
  for(i in seq_len(ncol(x))) {
    values <- as.character(x[, i])
    for(j in which(!is.na(values))) {
      value <- values[j]
      if(is.null(sep)) {
        value <- gsub(" ", "", value, fixed = TRUE)
        if(value == "") next
        if(nchar(value) %% 2L != 0L) {
          stop("Unseparated genotypes must have even width (row ", j,
               ", locus ", locus.names[i], ").", call. = FALSE)
        }
        half <- nchar(value) / 2L
        alleles <- c(missing.code(substr(value, 1L, half)),
                     missing.code(substr(value, half + 1L, nchar(value))))
      } else {
        alleles <- strsplit(value, sep, fixed = TRUE)[[1]]
        if(length(alleles) != 2L || any(!nzchar(alleles)) ||
           startsWith(value, sep) || endsWith(value, sep)) {
          stop("Expected two non-empty alleles (row ", j,
               ", locus ", locus.names[i], ").", call. = FALSE)
        }
      }
      result[j, c(2L * i - 1L, 2L * i)] <- alleles
    }
  }
  result
}
