#' @title GELATo - Group ExcLusion and Assignment Test
#' @description Run a GELATo test to evaluate assignment likelihoods of 
#'   groups of samples.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param unknown.strata a character vector listing to assign. Strata must 
#'   occur in \code{g}.
#' @param nrep number of permutation replicates for Fst distribution.
#' @param min.sample.size minimum number of samples to use to characterize 
#'   knowns. If the known sample size would be smaller than this after drawing 
#'   an equivalent number of unknowns for self-assignment, then the comparison 
#'   is not done.
#' @param gelato.result the result of a call to \code{gelato}.
#' @param unknown the names of the unknown strata in the \code{x$likelihoods}
#'   element to create plots. If \code{NULL} one plot for each stratum is
#'   created.
#' @param main main label for top of plot.#' 
#' @param num.cores The number of cores to use to distribute replicates over.
#'   If set to \code{NULL}, the value will be what is reported 
#'   by \code{\link[parallel]{detectCores} - 1}.
#' 
#' @return A list with the following elements:
#' \tabular{ll}{
#'   \code{assign.prob} \tab a data.frame of assignment probabilities.\cr
#'   \code{likelihoods} \tab a list of likelihoods.\cr
#' }
#' 
#' @references O'Corry-Crowe, G., W. Lucey, F.I. Archer, and B. Mahoney. 2015. 
#'   The genetic ecology and population origins of the beluga whales of 
#'   Yakutat Bay. Marine Fisheries Review 77(1):47-58
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples \dontrun{
#' data(msats.g)
#' 
#' # Run GELATo analysis
#' gelato.fine <- gelato(msats.g, unk = "Offshore.South", nrep = 20, num.cores = 2)
#' gelato.fine
#' 
#' # Plot results
#' gelatoPlot(gelato.fine, "Offshore.South")
#' }
#' 
#' @name gelato
#' @export
#' 
gelato <- function(g, unknown.strata, nrep = 1000, 
                   min.sample.size = 5, num.cores = 1) {
  # Check unknown strata
  all.strata <- getStrata(g)
  unknown.strata <- unique(as.character(unknown.strata))
  if(!all(unknown.strata %in% all.strata)) {
    stop("Some 'unknown.strata' could not be found in 'g'")
  }
  
  # identify known strata
  knowns <- sort(setdiff(all.strata, unknown.strata))
  
  # get number of cores
  num.cores <- .getNumCores(num.cores)
  
  strata.freq <- table(all.strata)
  # loop through every unknown stratum
  result <- sapply(unknown.strata, function(unknown) {
    # individual ids in this unknown
    unknown.ids <- names(all.strata)[all.strata == unknown]
    
    # loop through each known population and calculate distribution
    #   of Fst and log-likelihood of membership
    unknown.result <- sapply(knowns, function(known) {
      # check that enough samples exist in this known
      can.run <- (strata.freq[known] - length(unknown.ids)) >= min.sample.size
      if(!can.run) return(NULL)
      # individual ids in this known
      known.ids <- names(all.strata)[all.strata == known]
      known.g <- g[known.ids, , , drop = TRUE]
      
      # run permutations to collect Fst distributions 
      # fst.dist <- do.call(rbind, lapply(
      #   1:nrep, .gelatoPermFunc, known.ids = known.ids, 
      #   unknown.ids = unknown.ids, g = g, known.g = known.g
      # ))
      
      cl <- swfscMisc::setupClusters(num.cores)
      fst.dist <- tryCatch({
        if(is.null(cl)) {
          lapply(
            1:nrep, .gelatoPermFunc, known.ids = known.ids, 
            unknown.ids = unknown.ids, g = g, known.g = known.g
          )
        } else {
          parallel::clusterEvalQ(cl, require(strataG))
          parallel::clusterExport(
            cl, 
            c("known.ids", "unknown.ids", "g", "known.g"), 
            environment()
          )
          parallel::parLapplyLB(
            cl, 1:nrep, .gelatoPermFunc, known.ids = known.ids, 
            unknown.ids = unknown.ids, g = g, known.g = known.g
          )
        }
      }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)
      fst.dist <- do.call(rbind, fst.dist)
      fst.dist <- fst.dist[apply(fst.dist, 1, function(x) all(!is.na(x))), ]
      
      if(nrow(fst.dist) < 2) {
        NULL
      } else {
        # summarize Fst distribution
        null.mean <- mean(fst.dist[, "null"])
        null.sd <- stats::sd(fst.dist[, "null"])
        null.lik <- stats::dnorm(fst.dist[, "obs"], null.mean, null.sd)
        log.Lik <- sum(log(null.lik), na.rm = T)     
        obs.median <- stats::median(fst.dist[, "obs"], na.rm = T)
        obs.mean <- mean(fst.dist[, "obs"], na.rm = T)
        list(
          fst.dist = fst.dist, 
          log.Lik.smry = c(
            log.Lik = log.Lik, 
            mean.nreps = log.Lik / length(log.Lik),
            median = log(stats::dnorm(obs.median, null.mean, null.sd)), 
            mean = log(stats::dnorm(obs.mean, null.mean, null.sd))
          ),
          norm.coefs = c(mean = null.mean, sd = null.sd)
        )
      }
    }, simplify = FALSE)
    
    # calculate median logLikehood of assignment to each known
    log.Lik <- sapply(unknown.result, function(x) {
      if(is.null(x)) NA else x$log.Lik.smry["median"]
    })
    
    lik <- exp(log.Lik - max(log.Lik, na.rm = T))
    assign.prob <- lik / sum(lik, na.rm = T) 
    names(assign.prob) <- knowns
    
    list(assign.prob = assign.prob, likelihoods = unknown.result)
  }, simplify = F)
  
  assign.prob <- as.data.frame(t(sapply(result, function(x) x$assign.prob)))
  assign.prob$assignment <- apply(assign.prob, 1, function(x) {
    colnames(assign.prob)[which.max(x)]
  })
  result <- lapply(result, function(x) x$likelihoods)
  
  list(assign.prob = assign.prob, likelihoods = result)
}


#' @keywords internal
#'
.gelatoPermFunc <- function(i, known.ids, unknown.ids, g, known.g) {
  # select samples to self assign
  ran.knowns <- sample(known.ids, length(unknown.ids))
  
  # extract gtypes of base known strata
  known.to.keep <- setdiff(known.ids, ran.knowns)
  
  # gtypes for observed Fst
  obs.g <- g[c(unknown.ids, known.to.keep), , , drop = TRUE]
  
  # gtypes for null Fst
  st <- getStrata(known.g)
  st[ran.knowns] <- "<gelato.unknown>"
  setStrata(known.g) <- st
  
  c(obs = unname(statFst(obs.g, nrep = 0)$result["estimate"]), 
    null = unname(statFst(known.g, nrep = 0)$result["estimate"])
  )
}


#' @rdname gelato
#' @export
#' 
gelatoPlot <- function(gelato.result, unknown = NULL, main = NULL) { 
  if(is.null(unknown)) unknown <- names(gelato.result$likelihoods)
  for(unk in unknown) {
    lik <- gelato.result$likelihoods[[unk]]
    lik <- lik[!sapply(lik, is.null)]
    if(length(lik) == 0) {
      stop(paste("No likelihood distributions available for '", 
                 unknown, "'", sep = ""))
    }
    xticks <- pretty(unlist(sapply(lik, function(x) x$fst.dist)))
    xlim <- range(xticks)
    op <- graphics::par(mar = c(3, 3, 3, 2) + 0.1, 
              oma = c(2, 2, 0.1, 0.1), 
              mfrow = c(length(lik), 1)
    )
    high.prob <- gelato.result$assign.prob[unknown, "assignment"]
    for(known in names(lik)) {
      known.lik <- lik[[known]]
      null.max <- max(graphics::hist(known.lik$fst.dist[, "null"], plot = F)$density)
      obs.max <- max(graphics::hist(known.lik$fst.dist[, "obs"], plot = F)$density)
      lik.mean <- known.lik$norm.coefs["mean"]
      lik.sd <- known.lik$norm.coefs["sd"]
      norm.lik <- stats::dnorm(lik.mean, lik.mean, lik.sd)
      ylim <- range(pretty(c(0, null.max, obs.max, norm.lik)))
      graphics::hist(known.lik$fst.dist[, "null"], breaks = 10, freq = FALSE, 
           xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "", 
           col = "red", xaxt = "n")
      x <- NULL # To avoid R CMD CHECK warning about no global binding for 'x'
      graphics::curve(stats::dnorm(x, lik.mean, lik.sd), from = xlim[1], to = xlim[2],
            add = TRUE, col = "black", lwd = 3, ylim = ylim)
      graphics::par(new = TRUE)
      graphics::hist(known.lik$fst.dist[, "obs"], breaks = 10, freq = FALSE, 
           xlim = xlim, ylim = ylim, xlab = "", ylab = "", col = "darkgreen", 
           main = "", xaxt = "n", yaxt = "n") 
      graphics::axis(1, pretty(xlim))
      ll.median <- known.lik$log.Lik.smry["median"]
      log.lik <- if(!is.infinite(ll.median)) {
        format(ll.median, digits = 4) 
      } else "Inf"
      p.val <- format(gelato.result$assign.prob[unknown, known], digits = 2)
      pop <- paste(known, " (lnL = ", log.lik, ", p = ", p.val, ")", sep = "")
      graphics::mtext(pop, side = 3, line = 1, adj = 1, 
            font = ifelse(known == high.prob, 2, 1))
    }
    graphics::mtext("Fst", side = 1, outer = T, cex = 1.2)
    graphics::mtext("Density", side = 2, outer = T, cex = 1.2)
    graphics::par(op)
    if(!is.null(main)) graphics::mtext(main, side = 3, line = 3, adj = 0, font = 3)
  }
}