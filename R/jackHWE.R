#' @title Hardy-Weinberg Equlibrium Jackknife
#' @description Test influence of samples on Hardy-Weinberg equilibrium via 
#'   jackknife.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param exclude.num Number of samples to exclude at a time.
#' @param min.hwe.samples minimum samples needed to calculate HWE.
#' @param show.progress logical. Show progress of jackknife?
#' @param ... other arguments to be passed to \code{\link{hweTest}}.
#' @param jack.result result from run of \code{jackHWE}.
#' @param alpha critical value to determine if exclusion is "influential".
#' @param x result from a call to \code{jackInfluential}.
#' @param main main title for influential sample plots from 
#'   \code{plot.jack.influential}.
#' 
#' @details \describe{
#'   \item{\code{jackHWE}}{performs a HWE jackknife where all combinations 
#'     of \code{exclude.num} samples are left out and HWE is recalculated}
#'   \item{\code{jackInfluential}}{calculates odds.ratios between jackknife 
#'     HWE and observed HWE and identifies "influential" samples. Samples 
#'     are "influential" if the observed HWE p-value is < \code{alpha}, but 
#'     is > \code{alpha} when the samples are not present}
#'   \item{\code{plot.jack.influential}}{creates a cumulative frequency plot 
#'     of all odds-ratios from \code{jack.influential}. A vertical dashed 
#'     line marks the smallest influential exclusion}
#' }
#' 
#' @return \code{jackHWE} returns a list with:
#' \item{obs}{a named vector of HWE p-values for each locus.}
#' \item{jack}{a \code{data.frame} of HWE p-values where each row is an 
#'   exclusion and columns are loci.}
#' \item{gtypes}{the original \code{gtypes} object.}\cr
#' \code{jackInfluential} returns a list with:
#' \item{influential}{a \code{data.frame} of influential exclusions.}
#' \item{allele.freqs}{a \code{data.frame} listing the allele frequencies of 
#'   influential exclusions.}
#' \item{odds.ratio}{a \code{matrix} of odds ratios between exclusions (rows) 
#'   and loci (columns).} 
#' 
#' @references Morin, P.A., R.G. LeDuc, F.I. Archer, K.K. Martien, 
#'   R. Huebinger, J.W. Bickham, and B.L. Taylor. 2009. Significant deviations 
#'   from Hardy-Weinberg equilibirum caused by low levels of microsatellite 
#'   genotyping errors. Molecular Ecology Resources 9:498-504.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{hweTest}}
#' 
#' @export
#' 
jackHWE <- function(g, exclude.num = 1, min.hwe.samples = 5, 
                    show.progress = TRUE, ...) {  
  
  if((getNumInd(g) - exclude.num) < min.hwe.samples){
    stop(paste("'exclude.num' or 'min.HWE.samples' is too large to analyze this data"))
  }
  
  # Setup array of samples to exclude
  exclude.arr <- utils::combn(getIndNames(g), exclude.num)
  
  start.time <- Sys.time()
  nsteps <- ncol(exclude.arr) + 1 
  if(show.progress) {
    cat(format(start.time, "%Y-%m-%d %H:%M:%S"), "|", 1, "/", nsteps, "\n")
  }
  
  # calculate observed HWE
  obs <- hweTest(g, ...)
  
  # jackknife HWE
  jack <- sapply(1:ncol(exclude.arr), function(i) {           
    if(show.progress) {
      now <- Sys.time()    
      avg.time <- difftime(now, start.time, units = "secs") / (i + 1)
      est.complete.time <- now + (avg.time * (nsteps - i))
      cat(format(now, "%Y-%m-%d %H:%M:%S"), "|", i + 1, "/", nsteps, "| ETA =", 
          format(est.complete.time, "%Y-%m-%d %H:%M:%S\n")
      )
    }
    to.keep <- setdiff(getIndNames(g), exclude.arr[, i])
    jack.gtypes <- g[to.keep, , ]
    hweTest(jack.gtypes, ...)
  })

  result <- list(
    obs = obs, 
    jack = data.frame(
      excluded = apply(exclude.arr, 2, paste, collapse = ", "), 
      t(jack), 
      stringsAsFactors = FALSE
    ), 
    gtypes = g
  )
  class(result) <- c("jack.result", class(result))
  result
}


#' @rdname jackHWE
#' @export
#' 
jackInfluential <- function(jack.result, alpha = 0.05) {  
  if(!inherits(jack.result, "jack.result")) {
    stop("'jack.result' is not from 'jack.hwe'")
  }
  
  obs <- jack.result$obs
  exclude.vec <- jack.result$jack$excluded
  jack <- t(as.matrix(jack.result$jack[, -1]))
  
  # Find influential samples and calculate odds and odds ratios
  obs.arr <- matrix(obs, length(obs), length(exclude.vec))
  which.infl <- which((obs.arr <= alpha) & (jack > alpha), arr.ind = TRUE)
  influential <- if(nrow(which.infl) > 0) {
    tibble::tibble(
      excluded = exclude.vec[which.infl[, "col"]],
      locus = names(obs)[which.infl[, "row"]],
      obs.pval = obs[which.infl[, "row"]],
      jack.pval = apply(which.infl, 1, function(i) jack[i[1], i[2]])
    ) %>% 
      dplyr::mutate(
        obs.odds = swfscMisc::odds(.data$obs.pval),
        jack.odds = swfscMisc::odds(.data$jack.pval),
        odds.ratio = .data$jack.odds / .data$obs.odds
      ) %>% 
      dplyr::arrange(dplyr::desc(.data$odds.ratio)) %>% 
      as.data.frame()
  } else NULL
  
  # Calculate allele frequencies of influential samples 
  allele.freqs <- if(!is.null(influential)) {
    purrr::map(1:nrow(influential), function(i) {
      ids <- unlist(strsplit(influential$excluded[i], ", "))
      tibble::tibble(id = ids, locus = influential$locus[i])
    }) %>% 
      dplyr::bind_rows() %>% 
      unique() %>% 
      .alleleFreqFormat(jack.result$gtypes)
  } else NULL
  
  # Calculate full odds ratio matrix
  odds.ratio <- t(swfscMisc::odds(jack) / swfscMisc::odds(obs.arr))
  rownames(odds.ratio) <- exclude.vec
  colnames(odds.ratio) <- names(obs)
  
  result <- list(
    influential = influential, 
    allele.freqs = allele.freqs, 
    odds.ratio = odds.ratio
  )
  class(result) <- c("jack.influential", class(result))
  result
}


#' @rdname jackHWE
#' @method plot jack.influential
#' @export
#' 
plot.jack.influential <- function(x, main = NULL, ...) {
  or.tbl <- tibble::tibble(or = as.vector(x$odds.ratio)) %>% 
    dplyr::filter(
      !is.na(.data$or) & 
      !is.nan(.data$or) & 
      !is.infinite(.data$or)
    ) %>% 
    dplyr::arrange(.data$or) %>% 
    dplyr::mutate(freq = 1:dplyr::n() / dplyr::n()) 
    
  p <- ggplot2::ggplot(or.tbl, ggplot2::aes_string(x = "or", y = "freq")) + 
    ggplot2::geom_line() +
    ggplot2::labs(x = "Odds ratio", y = "Cumulative frequency")
  
  if(!is.null(main)) p <- p + ggplot2::ggtitle(main)
  
  if(!is.null(x$influential)) {
    p <- p + ggplot2::geom_vline(
      xintercept = min(x$influential$odds.ratio, na.rm = TRUE), 
      linetype = 2, 
      color = "red"
    ) 
  }
  
  print(p)
}

#' @rdname jackHWE 
#' @param x a matrix or data.frame where first column is sample id and 
#'   second colum is locus name.
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.alleleFreqFormat <- function(x, g) {
  x <- as.matrix(x)
  
  freqs <- alleleFreqs(g, by.strata = FALSE, type = "prop")
  fmtd <- rep(as.character(NA), nrow(x))
  for(i in 1:length(fmtd)) {
    # get genotype of this id at this locus
    gt <- g@data %>% 
      dplyr::filter(.data$id ==  x[i, 1] & .data$locus == x[i, 2]) %>% 
      dplyr::pull(.data$allele) 
    # if the genotype is NA skip and leave format as NA
    if(any(is.na(gt))) next
    # get frequency and round
    f <- round(freqs[[x[i, 2]]][gt], 3)
    f <- paste(gt, " (", format(f, nsmall = 3), ")", sep = "")
    # return a single frequency if homozygote
    fmtd[i] <- if(length(unique(gt)) == 1) f[1] else paste(f, collapse = " / ")
  }
  cbind(x, allele.freqs = fmtd)
}