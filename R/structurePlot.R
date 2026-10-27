#' @title Plot STRUCTURE Results
#' @description Plot Q-matrix from a call to \code{\link{structure}} or 
#'   \code{\link{clumpp}}.
#' 
#' @param q.mat matrix or data.frame of assignment probabilities.
#' @param type the type of plot (\code{"bar"} or \code{"area"}) to create. 
#'   see \code{\link[swfscMisc]{plotAssignments}} for more information.
#' @param plot display plot?
#' 
#' @return invisibly, the ggplot object
#' 
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#' 
#' @seealso \code{\link{structure}}, \code{\link{clumpp}}
#'
#' @export
#' 
structurePlot <- function(
  q.mat, 
  type = NULL,
  plot = TRUE
) {
  p <- swfscMisc::plotAssignments(
    q.mat |> 
      dplyr::select(id, dplyr::starts_with('Group.')) |> 
      tibble::column_to_rownames('id'),
    q.mat$orig.pop,
    type = 'bar', 
    plot = FALSE
  ) + 
    ggplot2::theme(legend.position = 'none')
  
  if(plot) print(p)
  invisible(p)
}
