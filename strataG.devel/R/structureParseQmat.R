#' @title Parse Q-Matrix Text
#' 
#' @param q.mat.txt character vector of Q matrix from STRUCTURE output file.
#' @param pops vector of population labels to be used in place of numbers in 
#'   STRUCTURE file.
#' 
#' @details Function is used internally by \code{\link{structure}} and 
#'   \code{\link{clumpp}} to parse output files.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
structureParseQmat <- function(q.mat.txt, pops) {
  q.mat.txt <- sub("[*]+", "", q.mat.txt)
  q.mat.txt <- sub("[(]", "", q.mat.txt)
  q.mat.txt <- sub("[)]", "", q.mat.txt)
  q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
  
  # Parse population assignment portion of table to create single-line 
  #   data.frame
  do.call(rbind, lapply(q.mat.txt, function(x) {
    # Split on spaces and remove empty spaces and colons
    x <- strsplit(x, " ")[[1]]
    x <- x[!x %in% c("", ":")]
    p <- as.numeric(x[4])
  
    df <- data.frame(row = as.numeric(x[1]), id = x[2], 
      pct.miss = as.numeric(x[3]), 
      orig.pop = if(is.null(pops)) p else pops[p], 
      stringsAsFactors = FALSE
    )
    pop.prob <- as.data.frame(rbind(as.numeric(x[-(1:4)])))
    colnames(pop.prob) <- paste("prob", 1:ncol(pop.prob), sep = ".")
    cbind(df, pop.prob) 
  }))
}