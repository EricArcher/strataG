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
  if(!is.null(sep)) if(sep == "") sep <- NULL
  
  locus.names <- if(is.null(colnames(x))) {
    paste("Locus", 1:ncol(x), sep = "") 
  } else {
    colnames(x)
  }
  locus.names <- paste(rep(locus.names, each = 2), c(1, 2), sep = ".")
  
  x <- do.call(cbind, lapply(1:ncol(x), function(col) as.character(x[, col])))

  split.alleles <- lapply(1:ncol(x), function(i) {
    if(!is.null(sep)) {
      do.call(rbind, strsplit(x[, i], split = sep))
    } else {
      t(sapply(x[, i], function(a) {
        a <- sub(" ", "", a)
        if (is.na(a)) return(c(NA, NA))
        end <- nchar(a)
        half <- end / 2
        
        a1 <- substr(a, 1, half)
        a1.num <- suppressWarnings(as.numeric(a1))
        a1 <- if(is.na(a1.num)) a1 else if(a1.num == 0) NA else a1
        a1 <- if((a1 == "NA") | (a1 == "")) NA else a1
        
        a2 <- substr(a, half + 1, end)
        a2.num <- suppressWarnings(as.numeric(a2))
        a2 <- if(is.na(a2.num)) a2 else if(a2.num == 0) NA else a2
        a2 <- if((a2 == "NA") | (a2 == "")) NA else a2
        
        c(a1, a2)
      }))
    }
  })
  split.alleles <- do.call(cbind, split.alleles)
  colnames(split.alleles) <- locus.names
  rownames(split.alleles) <- NULL
  return(split.alleles)
}