#' @title Run Evanno Method on STRUCTURE Results
#' @description Calculate first and second order rates of changes of LnPr(K) 
#'   from STRUCTURE results based on Evanno et al. 2005.
#' 
#' @param sr output from a call to \code{\link{structure}}.
#' @param plot logical. Generate a plot of Evanno metrics?
#' 
#' @return a list with: \tabular{ll}{
#'   \code{df} \tab data.frame with Evanno log-likelihood metrics for each
#'     value of K.\cr
#'   \code{plots} \tab list of four ggplot objects for later plotting.\cr
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
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], stats::sd)
  
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
    mean.ln.k = ggplot2::ggplot(
        df, 
        ggplot2::aes_string(x = "k", y = "mean.ln.k")
      ) +
        ggplot2::ylab("mean LnP(K)") +
        ggplot2::geom_segment(
          ggplot2::aes_string(
            x = "k", 
            xend = "k", 
            y = "sd.min", 
            yend = "sd.max"
          )
        ),
    ln.pk = ggplot2::ggplot(
      df[!is.na(df$ln.pk), ], 
      ggplot2::aes_string(x = "k", y = "ln.pk")
    ) +
      ggplot2::ylab("LnP'(K)"),
    ln.ppk = ggplot2::ggplot(
      df[!is.na(df$ln.ppk), ], 
      ggplot2::aes_string(x = "k", y = "ln.ppk")
    ) +
      ggplot2::ylab("LnP''(K)")
  )
  if(!all(is.na(df$delta.k))) {
    plot.list$delta.k <- ggplot2::ggplot(
      df[!is.na(df$delta.k), ], 
      ggplot2::aes_string(x = "k", y = "delta.k")
    ) +
      ggplot2::ylab(expression(Delta(K)))
  }
  
  for(i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + 
      ggplot2::geom_line() +
      ggplot2::geom_point(fill = "white", shape = 21, size = 3) +
      ggplot2::xlim(c(1, max(df$k))) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }

  if(plot) {
    p <- plot.list %>% 
      purrr::map(function(x) {
        ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
      })
    maxWidth <- do.call(
      grid::unit.pmax, 
      purrr::map(p, function(x) x$widths[2:3])
    )
    for(i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
    p$bottom <- "K"
    p$ncol <- 2
    do.call(gridExtra::grid.arrange, p)
  } 
  
  df$sd.min <- df$sd.max <- NULL
  invisible(list(df = df, plots = plot.list))
}