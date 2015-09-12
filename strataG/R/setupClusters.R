#' @importFrom parallel detectCores makePSOCKcluster makeForkCluster
#' 
.setupClusters <- function(num.cores) {
  # setup clusters
  num.cores <- max(1, num.cores)
  num.cores <- min(num.cores, detectCores() - 1)
  if(num.cores > 1) {
    is.windows <- .Platform$OS.type == "windows"
    cl.func <- ifelse(is.windows, makePSOCKcluster, makeForkCluster)
    cl.func(num.cores)
  } else NULL
}