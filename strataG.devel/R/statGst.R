#' @rdname popStructStat
#' @export
#' 
statGst <- function(g, strata = NULL, ...) {
  if(ploidy(g) == 1) return(c(Gst = NA))
  hets <- Hstats(g, strata)
  Hs <- mean(hets["Hs", ], na.rm = TRUE)
  Ht <- mean(hets["Ht", ], na.rm = TRUE)
  est <- 1 - (Hs / Ht) 
  if(is.nan(est)) est <- NA
  names(est) <- "Gst"
  est
}


#' @rdname popStructStat
#' @export
#' 
statGstPrime <- function(g, strata = NULL, prime.type = "nei", ...) { 
  if(ploidy(g) == 1) return(c('G\'st' = NA))
  if(!prime.type %in% c("nei", "hedrick")) {
    stop("'prime.type' must be either 'nei' or 'hedrick'")
  }
  hets <- Hstats(g, strata)
  Hs <- mean(hets["Hs", ], na.rm = TRUE)
  Ht <- mean(hets["Ht", ], na.rm = TRUE)
  n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
  est <- if(prime.type == "nei") {
    (n * (Ht - Hs)) / ((n * Ht) - Hs) 
  } else {
    gst.max <- ((n - 1) * (1 - Hs)) / (n - 1 + Hs)
    1 - (Hs / Ht) / gst.max
  }
  if(is.nan(est)) est - NA
  names(est) <- "G'st"
  est
}


#' @rdname popStructStat
#' @export
#' 
statGstDblPrime <- function(g, strata = NULL, ...) {
  if(ploidy(g) == 1) return('G\'\'st' = NA)
  hets <- Hstats(g, strata)  
  Hs <- mean(hets["Hs", ], na.rm = TRUE)
  Ht <- mean(hets["Ht", ], na.rm = TRUE)
  n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
  est <- (n * (Ht - Hs)) / ((n * Ht) - Hs) / (1 - Hs)
  if(is.nan(est)) est - NA
  names(est) <- "G''st"
  est
}