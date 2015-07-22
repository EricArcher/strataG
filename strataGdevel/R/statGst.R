#' @rdname popStructStat
#' @export
#' 
statGst <- function(g, strata = NULL, ...) {
  if(ploidy(g) < 2) return(c(Gst = NA))
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  est <- statGst_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata(g)) - 1,
    ploidy(g)
  )
  
#   hets <- Hstats(g, strata)
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   est <- 1 - (Hs / Ht) 
#   if(is.nan(est)) est <- NA
  
  names(est) <- "Gst"
  est
}


#' @rdname popStructStat
#' @export
#' 
statGstPrime <- function(g, strata = NULL, prime.type = c("nei", "hedrick"), ...) { 
  if(ploidy(g) < 2) return(c('G\'st' = NA))
  prime.type <- match.arg(prime.type)
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  est <- statGstPrime_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata(g)) - 1,
    ploidy(g),
    switch(prime.type, nei = 0, hedrick = 1)
  )
  
#   hets <- Hstats(g, strata)
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
#   est <- if(prime.type == "nei") {
#     (n * (Ht - Hs)) / ((n * Ht) - Hs) 
#   } else {
#     gst.max <- ((n - 1) * (1 - Hs)) / (n - 1 + Hs)
#     (1 - (Hs / Ht)) / gst.max
#   }
#   if(is.nan(est)) est <- NA
  
  names(est) <- "G'st"
  est
}


#' @rdname popStructStat
#' @export
#' 
statGstDblPrime <- function(g, strata = NULL, ...) {
  if(ploidy(g) < 2) return(c('G\'\'st' = NA))
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  est <- statGstDblPrime_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata(g)) - 1,
    ploidy(g)
  )
  
#   hets <- Hstats(g, strata)  
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
#   est <- (n * (Ht - Hs)) / ((n * Ht) - Hs) / (1 - Hs)
#   if(is.nan(est)) est <- NA
  
  names(est) <- "G''st"
  est
}