#' @title Plot STRUCTURE Results
#' @description Plot Q-matrix from a call to \code{\link{structure}} or 
#'   \code{\link{clumpp}}.
#' 
#' @param q.mat matrix or data.frame of assignment probabilities.
#' @param pop.col column number identifying original population designations.
#' @param prob.col column number of first assignment probabilities to first 
#'  group. It is assumed that the remainder of columns 
#'  (\code{prob.col:ncol(q.mat)}) contain all assignment probabilities, 
#'  and thus \emph{k} = \code{ncol(q.mat) - prob.col + 1}.
#' @param sort.probs logical. Sort individuals by probabilities within 
#'   populations? 
#' If \code{FALSE} individuals will be plotted as in \code{q.mat}.
#' @param label.pops logical. Label the populations on the plot? If 
#'   \code{FALSE}, then population labels are omitted so the user can 
#'   customize their format.
#' @param col colors to use for each group.
#' @param horiz logical. Plot horizontal bars?
#' @param legend logical. Include a legend?
#' @param ... optional arguments to be passed to 
#'   \code{\link[graphics]{barplot}}.
#' 
#' @return invisibly, a list containing:
#' \tabular{ll}{
#'   \code{q.mat} \tab the sorted matrix or data.frame of assignment 
#'     probabilities used in the plot.\cr
#'   \code{bar.centers} \tab a vector of the centers of bars for each 
#'     individual.\cr
#'   \code{pop.ticks} \tab a vector of the tick marks separating populations.\cr
#'  }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{structure}}, \code{\link{clumpp}}
#'
#' @importFrom graphics par strwidth barplot axis mtext
#' @importFrom grDevices rainbow
#' @export
#' 
structurePlot <- function(q.mat, pop.col = 3, prob.col = 4, sort.probs = TRUE, 
                          label.pops = TRUE, col = NULL, horiz = TRUE,
                          legend = TRUE, ...) {    
  # sort q.mat within strata by probability
  prob.cols <- prob.col:ncol(q.mat)
  q.mat <- q.mat[order(q.mat[, pop.col]), ]
  if(sort.probs) {      
    q.mat <- do.call(rbind, by(q.mat, list(q.mat[, pop.col]), function(x) {
      prob.list <- lapply(prob.cols, function(i) x[, i])
      x[do.call(order, prob.list), ]
    }, simplify = FALSE))
  }
  
  # make sure probs sum to 1 and get matrix
  q.mat[, prob.cols] <- prop.table(as.matrix(q.mat[, prob.cols]), 1)
  assign.mat <- t(q.mat[, prob.cols])  
  
  # create barplot
  if(is.null(col)) col <- rainbow(length(prob.cols))
  mai <- par("mai")
  mai[1] = 0.7
  mai[3] = 0.2
  if(horiz) {
    mai[2] <- max(strwidth(unique(q.mat[, pop.col]), "inches")) + 0.4
  }
  mai[4] <- if(legend) {
    text.width <- max(strwidth(colnames(q.mat)[prob.cols], "inches"))
    text.width + 0.8
  } else 0.5
  op <- par(mai = mai, las = 1)
  if(legend) par(xpd = TRUE)
  bp <- barplot(assign.mat, axes = FALSE, axisnames = FALSE, col = col, horiz = horiz, ...)
  prob.side <- if(horiz) 1 else 2
  pop.side <- if(horiz) 2 else 1
  tx <- tapply(bp, q.mat[, pop.col], min) + 0.1
  tx <- c(tx - tx[1], max(bp) + tx[1])
  axis(prob.side, pos = min(tx) - 1.5)
  axis(pop.side, at = tx, pos = -0.005, labels = FALSE)
  if(label.pops) {
    lbl.x <- sapply(1:(length(tx) - 1),
                    function(i) tx[i] + (tx[i + 1] - tx[i]) / 2
    )
    names(lbl.x) <- names(tx)[1:(length(tx) - 1)]
    mtext(names(lbl.x), side = pop.side, at = lbl.x, line = 1)
  }
  if(legend) {
    x <- if(horiz) par("usr")[2] else max(tx)
    y <- if(horiz) max(tx) else par("usr")[4]
    legend(x, y, colnames(q.mat)[prob.cols],
           fill = col)
  }
  par(op)
  invisible(list(q.mat = q.mat, bar.centers = bp, pop.ticks = tx))
}