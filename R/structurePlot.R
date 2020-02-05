#' @title Plot STRUCTURE Results
#' @description Plot Q-matrix from a call to \code{\link{structure}} or 
#'   \code{\link{clumpp}}.
#' 
#' @param q.mat matrix or data.frame of assignment probabilities.
#' @param pop.col column number identifying original population designations.
#' @param prob.col column number of first assignment probabilities to first 
#'  group. It is assumed that the remainder of columns 
#'  (\code{prob.col:ncol(q.mat)}) contain all assignment probabilities.
#' @param sort.probs logical. Sort individuals by probabilities within
#'   populations? If \code{FALSE} individuals will be plotted as in
#'   \code{q.mat}.
#' @param label.pops logical. Label the populations on the plot?
#' @param col colors to use for each group.
#' @param horiz logical. Plot bars horizontally.
#' @param type either \code{"area"} for stacked continuous area plot or
#'   \code{"bar"} for discrete stacked bar chart. The latter is prefered for
#'   small numbers of samples. If not specified, a bar chart will be used if
#'   there are <= 100 samples.
#' @param legend.position the position of the legend (\code{"top", "left", 
#'   "right", "bottom"}, or two-element numeric vector).
#' @param plot display plot?
#' 
#' @return invisibly, the ggplot object
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{structure}}, \code{\link{clumpp}}
#'
#' @export
#' 
structurePlot <- function(
  q.mat, 
  pop.col = 3, 
  prob.col = 4, 
  sort.probs = TRUE,
  label.pops = TRUE, 
  col = NULL, 
  horiz = TRUE, 
  type = NULL,
  legend.position = c("top", "left", "right", "bottom", "none"),
  plot = TRUE
) {
  
  legend.position <- match.arg(legend.position)
  
  # convert q.mat to sorted data.table
  prob.cols <- prob.col:ncol(q.mat)
  qm <- as.data.frame(q.mat)[, c(pop.col, prob.cols), drop = FALSE]
  qm[, 1] <- factor(
    qm[, 1], 
    levels = sort(unique(qm[, 1]), decreasing = !horiz)
  )
  sort.cols <- c(1, if(sort.probs) 2:ncol(qm) else NULL)
  i <- do.call(
    order, 
    c(as.list(qm[, sort.cols, drop = FALSE]), decreasing = TRUE)
  )
  qm <- qm[i, ]
  qm$x <- 1:nrow(qm)
  
  # Get population frequencies, centers and dividing points
  pop.freq <- table(qm[, 1])
  levels(qm[, 1]) <- paste(
    levels(qm[, 1]), "\n(n = ", pop.freq, ")", sep = ""
  )
  pop.cntr <- tapply(qm$x, qm[, 1], mean)
  pop.div <- rev(tapply(qm$x, qm[, 1], min))[-1] - 0.5
  
  # Create data.frame for plotting
  df <- melt(qm[, c("x", colnames(qm)[-ncol(qm)])], id.vars = c(1, 2),
             variable.name = "Group", value.name = "probability")
  colnames(df)[1:2] <- c("x", "population")
  df <- df[order(-as.numeric(df$Group), df$probability), ]
  
  type <- if(is.null(type)) {
    if(nrow(df) <= 100) "bar" else "area"
  } else {
    match.arg(type, c("bar", "area"))
  }
  
  # Plot stacked bar graphs
  g <- ggplot2::ggplot(df, ggplot2::aes_string("x", "probability")) +  
    switch(
      type,
      area = ggplot2::geom_area(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      ),
      bar = ggplot2::geom_bar(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      )
    ) +
    ggplot2::ylab("Pr(Group Membership)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  if(label.pops) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div) +
      ggplot2::scale_x_continuous(
        name = "", 
        breaks = pop.cntr, 
        labels = names(pop.cntr),
        expand = c(0, 0)
      )
  } else {
    g <- g + 
      ggplot2::xlab("") + 
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  if(horiz) g <- g + ggplot2::coord_flip()
  if(!is.null(col)) g <- g + ggplot2::scale_fill_manual(values = col)
  
  if(plot) print(g)
  invisible(g)
}