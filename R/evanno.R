#' @title Run Evanno Method on STRUCTURE Results
#' @description Calculate first and second order rates of changes of LnPr(K) 
#'   from STRUCTURE results based on Evanno et al. 2005.
#' 
#' @param sr output from a call to \code{\link{structure}}.
#' @param plot logical. Generate a plot of Evanno metrics?
#' 
#' @return a list with:
#' \describe{
#'   \item{\code{df}}{data.frame with Evanno log-likelihood metrics for each value of K.}
#'   \item{\code{plots}}{list of four ggplot objects for later plotting.}
#' }
#' 
#' @references Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the 
#'   number of clusters of individuals using the software STRUCTURE: a 
#'   simulation study. Molecular Ecology 14:2611-2620.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{structure}} \code{\link{clumpp}}
#' 
#' @examples \dontrun{
#' data(msats.g)
#' 
#' # Run STRUCTURE
#' sr <- structureRun(msats, k.range = 1:4, num.k.rep = 10)
#' 
#' # Calculate Evanno metrics
#' evno <- evanno(sr)
#' evno
#' }
#' 
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot aes_string geom_line geom_segment geom_point xlim ylab theme element_blank theme_void geom_text ggplot_gtable ggplot_build
#' @importFrom grid unit.pmax
#' @importFrom gridExtra grid.arrange
#' @export
#' 
evanno <- function(sr, plot = TRUE) {
  if(!"structure.result" %in% class(sr)) {
    stop("'sr' is not a result from 'structure.run'.")
  }
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  if(length(k.tbl) < 3) stop("must have at least two values of k.")
  
  # collect summary statistics
  sr.smry <- t(sapply(sr, function(x) x$summary))
  
  # calculate mean and sd of LnPr(K)
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], sd)
  
  # Ln'(K)
  ln.pk <- diff(ln.k)
  # Ln''(K)
  ln.ppk <- abs(diff(ln.pk))
  # Delta-K
  delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
    abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1]) / sd.ln.k[i]
  })
  
  df <- data.frame(
    k = as.numeric(names(ln.k)),
    reps = as.numeric(table(sr.smry[, "k"])),
    mean.ln.k = as.numeric(ln.k),
    sd.ln.k = as.numeric(sd.ln.k),
    ln.pk = c(NA, ln.pk),
    ln.ppk = c(NA, ln.ppk, NA),
    delta.k = c(NA, delta.k,  NA)
  )
  rownames(df) <- NULL
  
  # Build plots
  df$sd.min <- df$mean.ln.k - df$sd.ln.k
  df$sd.max <- df$mean.ln.k + df$sd.ln.k
  
  plot.list <- list(
    mean.ln.k = ggplot(df, aes_string(x = "k", y = "mean.ln.k")) +
      ylab("mean LnP(K)") +
      geom_segment(aes_string(x = "k", xend = "k", y = "sd.min", yend = "sd.max")),
    ln.pk = ggplot(df[!is.na(df$ln.pk), ], aes_string(x = "k", y = "ln.pk")) +
      ylab("LnP'(K)"),
    ln.ppk = ggplot(df[!is.na(df$ln.ppk), ], aes_string(x = "k", y = "ln.ppk")) +
      ylab("LnP''(K)")
  )
  if(!all(is.na(df$delta.k))) {
    plot.list$delta.k <- ggplot(df[!is.na(df$delta.k), ], aes_string(x = "k", y = "delta.k")) +
        ylab(expression(Delta(K)))
  }
  
  for(i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + 
      geom_line() +
      geom_point(fill = "white", shape = 21, size = 3) +
      xlim(c(1, max(df$k))) +
      theme(axis.title.x = element_blank())
  }

  if(plot) {
    p <- lapply(plot.list, function(x) ggplot_gtable(ggplot_build(x)))
    maxWidth <- do.call(unit.pmax, lapply(p, function(x) x$widths[2:3]))
    for(i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
    p$bottom <- "K"
    p$ncol <- 2
    do.call(grid.arrange, p)
  } 
  
  df$sd.min <- df$sd.max <- NULL
  print(df)
  invisible(list(df = df, plots = plot.list))
}