#' @title gtype Random Forest
#' @description Conduct Random Forest on a gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} slot.
#' @param pairwise do analysis on all pairwise combinations of strata?
#' @param conf.level confidence level for the \code{\link{binom.test}} confidence interval
#' @param unk vector of strata to be treated as "unknowns" for prediction with 
#' @param ... arguments passed to \code{\link{sequenceRF}} and \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}). If \code{pairwise} is \code{TRUE} then 
#'   the \code{rf} element is a list of \code{randomForest} results for each row in
#'   \code{smry}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' library(strataG)
#' data(dloop.g)
#' 
#' rf.result <- gtypesRF(dloop.g, pairwise = TRUE)
#' 
#' rf.result
#' 
#' @export

gtypesRF <- function(g, gene = 1, pairwise = FALSE, conf.level = 0.95, 
                     unk = NULL, ...) {
  # check for and set names to sampsize
  args <- list(...)
  ss.match <- which(pmatch(names(args), "sampsize") == 1)
  if(length(ss.match) == 1) {
    temp.ss <- args[[ss.match]]
    if(is.null(names(temp.ss))) {
      names(temp.ss) <- getStrataNames(g)
      args[[ss.match]] <- temp.ss
      names(args)[ss.match] <- "sampsize"
    }
  }
  
  if(pairwise) {
    sp <- .strataPairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- g[, , unlist(sp[i, ])]
      do.call(gtypesRF, c(list(g = pair.g, pairwise = FALSE), args))
    })
    list(
      smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
      rp = lapply(result, function(x) x$rp)
    )
  } else {
    seq.df <- gtypes2rfDF(g, gene = gene)
    rf.df <- if(!is.null(unk)) {
      df <- seq.df[!seq.df$stratum %in% unk, ]
      df$stratum <- droplevels(df$stratum)
      df
    } else seq.df
    args <- c(list(x = rf.df), args)
    result <- do.call(sequenceRF, args)
    result$seq.df <- seq.df
    result$unk <- if(!is.null(unk) & !is.null(result$rp)) {
      unk.df <- seq.df[seq.df$stratum %in% unk, ]
      pred <- data.frame(stats::predict(result$rp, unk.df, type = "prob"))
      pred$predicted <- stats::predict(result$rp, unk.df, type = "response")
      pred
    } else NULL
    result
  }
}
