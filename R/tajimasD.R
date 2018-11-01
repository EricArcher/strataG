#' @title Tajima's D
#' @description Calculate Tajima's D for a set of sequences to test 
#'   for selection.
#' 
#' @param x set of DNA sequences or a haploid \linkS4class{gtypes} 
#'   object with sequences.
#' @param CI desired central confidence interval.
#' 
#' @return A named vector with the estimate for \code{D} and 
#'   the \code{p.value} that it is different from 0.
#' 
#' @references Tajima, F. 1989. Statistical method for testing the neutral 
#'   mutation hypothesis by DNA polymorphism. Genetics 123:585-595.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' 
#' tajimasD(dolph.seqs)
#' 
#' @export
#' 
tajimasD <- function(x, CI = 0.95) {
  # making life easier:
  # D.to.x <- function(D, Dmin, Dmax) (D - Dmin) / (Dmax - Dmin) # conversion 1
  x.to.D <- function(x, Dmin, Dmax) x * (Dmax - Dmin) + Dmin   # conversion 2
  
  # distribution function from Tajima paper:
  beta.D <- function(tajima.D, alpha, beta, Dmin, Dmax) { 
    # return P(D | alpha, beta, Dmin, Dmax)
    gamma(alpha + beta) * (Dmax - tajima.D) ^ (alpha - 1) * 
      (tajima.D - Dmin) ^ (beta - 1) / 
      (gamma(alpha) * gamma(beta) * (Dmax - Dmin) ^ (alpha + beta - 1))
  }
  
  x <- if(inherits(x, "gtypes")) {
    getSequences(x, as.haplotypes = FALSE, as.multidna = TRUE)
  } else {
    as.multidna(x)
  }
  x <- apex::getSequences(x, simplify = FALSE)
  
  purrr::map(x, function(dna) {
    dna <- as.matrix(dna)
    num.sites <- ncol(dna)
    
    pws.diff <- ape::dist.dna(
      dna, 
      model = "N", 
      pairwise.deletion = TRUE, 
      as.matrix = TRUE
    )
    pi <- mean(pws.diff[lower.tri(pws.diff)])
    S <- ncol(variableSites(dna)$site.freqs) # number segregating sites
    n <- nrow(dna) # number individuals
    
    # now for all the pieces...
    n.vec <- 1:(n-1) # common vector used
    a1 <- sum(1 / n.vec)
    a2 <- sum(1 / n.vec ^ 2)
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n ^ 2 + n + 3) / (9 * n * (n - 1))
    c1 <- b1 - 1 / a1
    c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1 ^ 2
    e1 <- c1 / a1
    e2 <- c2 / (a1 ^ 2 + a2)
    
    # which allows you to calculate D!
    D_obs <- (pi - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))
    
    if(is.nan(D_obs)) {
      warning("D cannot be computed (division by zero in final equation)")
      return(c(D = NA, p.value = NA))
    }
    
    # compare with expected value from beta distribution:
    DMin <- (2 / n - 1 / a1) / e2 ^ (1 / 2) #S approaches infinity
    DMax <- if(n %% 2 == 0) {
      (n / (2 * (n - 1)) - 1 / a1) / e2 ^ (1 / 2)
    } else {
      ((n + 1)/(2 * n) - 1 / a1) / e2 ^ (1 / 2)
    }
    Alpha <- -(1 + DMin * DMax) * DMax / (DMax - DMin)
    Beta <- (1 + DMin * DMax) * DMin / (DMax - DMin)
    
    # 95% confidence interval:
    lci.p <- (1 - CI) / 2
    LCI <- x.to.D(stats::qbeta(lci.p, Beta, Alpha), DMin, DMax)
    UCI <- x.to.D(stats::qbeta(1 - lci.p, Beta, Alpha), DMin, DMax)
    
    # important probabilities 
    # D negative: intregrate from Dmin to D, 
    # D positive: integrate from D to Dmax
    tryCatch({
      prob <- stats::integrate(
        beta.D, 
        lower = ifelse(D_obs < 0, DMin, D_obs),
        upper = ifelse(D_obs < 0, D_obs, DMax), 
        alpha = Alpha, 
        beta = Beta, 
        Dmin = DMin, 
        Dmax = DMax
      )
      tibble::tibble(D = D_obs, p.value = prob$value, LCI = LCI, UCI = UCI)
    }, error = function(e) {
      warning("error in Tajima's D integration, NA returned")
      tibble::tibble(D = NA, p.value = NA, LCI = NA, UCI = NA)
    })
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(locus = names(x)) %>% 
    as.data.frame()
}