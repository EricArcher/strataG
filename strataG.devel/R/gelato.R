#' @name gelato
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
#' @param num.cores number of CPU cores to use. Value is passed to 
#'   \code{\link[parallel]{mclapply}}.
#' @param gelato.result the result of a call to \code{gelato}.
#' @param unknown the name of an unknown stratum in the \code{x$likelihoods} 
#'   element.
#' @param main main label for top of plot.
#' 
#' @return A list with the following elements:
#' \tabular{ll}{
#'   \code{assign.prob} \tab a data.frame of assignment probabilities.\cr
#'   \code{likelihoods} \tab a list of likelihoods.\cr
#' }
#' 
#' @references O'Corry-Crowe et. al. XXXX
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom parallel mclapply
#' @export
#' 
gelato <- function(g, unknown.strata, nrep = 1000, min.sample.size = 5, 
                   num.cores = 1) { 
#   # Check unknown strata
#   if(!is.character(unknown.strata) & !is.vector(unknown.strata)) {
#     stop("'unknown.strata' must be a character vector")
#   }
#   all.strata <- attr(g, "strata.names")[-1]
#   unknown.strata <- unique(unknown.strata)
#   if(!all(unknown.strata %in% all.strata)) {
#     stop("Some 'unknown.strata' could not be found in 'g'")
#   }
#   
#   opt <- options(mc.cores = num.cores)
#   
#   knowns <- sort(setdiff(all.strata, unknown.strata))
#   # loop through every unknown strata
#   result <- sapply(unknown.strata, function(unknown) {
#     unknown.gtypes <- subset(g, strata = unknown)
#     unknown.n <- nrow(unknown.gtypes$genotypes)
#     
#     # loop through each known population and calculate distribution
#     #   of Fst and log-likelihood of membership
#     unknown.result <- sapply(knowns, function(known) {
#       known.gtypes <- subset.gtypes(g, strata = known)
#       if((nrow(known.gtypes$genotypes) - unknown.n) >= min.sample.size) {
#         fst.dist <- do.call(rbind, mclapply(1:nrep, function(i) {
#           # select samples to self assign
#           ran.sample <- sample(rownames(known.gtypes$genotypes), unknown.n)
#           # extract gtypes of base known strata
#           known.to.keep <- setdiff(rownames(known.gtypes$genotypes), ran.sample)
#           known.sample <- subset(known.gtypes, ids = known.to.keep)
#           # gtypes for observed Fst
#           obs.gtypes <- merge(known.sample, unknown.gtypes)
#           # gtypes for null Fst
#           null.gtypes <- decode(known.gtypes)
#           null.gtypes$genotypes[ran.sample, "strata"] <- rep("gelato.unknown", unknown.n)
#           gen.mat <- cbind(id = rownames(null.gtypes$genotypes), null.gtypes$genotypes)
#           null.gtypes <- gtypes(gen.mat, dna.seq = null.gtypes$sequences)
#           c(obs = stat.fst(obs.gtypes)$estimate, null = stat.fst(null.gtypes)$estimate)
#         }))
#         fst.dist <- fst.dist[apply(fst.dist, 1, function(x) all(!is.na(x))), ]
#         
#         if(nrow(fst.dist) < 2) {
#           NULL
#         } else {
#           # summarize Fst distribution
#           null.mean <- mean(fst.dist[, "null"])
#           null.sd <- sd(fst.dist[, "null"])
#           log.Lik <- sum(log(dnorm(fst.dist[, "obs"], null.mean, null.sd)), na.rm = T)        
#           list(fst.dist = fst.dist, 
#                log.Lik.smry = c(log.Lik = log.Lik, mean.nreps = log.Lik / length(log.Lik),
#                                 median = log(dnorm(median(fst.dist[, "obs"], na.rm = T), null.mean, null.sd)), 
#                                 mean = log(dnorm(mean(fst.dist[, "obs"], na.rm = T), null.mean, null.sd))
#                ),
#                norm.coefs = c(mean = null.mean, sd = null.sd)
#           )
#         }
#       } else NULL
#     }, simplify = F)
#     
#     # calculate median logLikehood of assignment to each known
#     log.Lik <- sapply(unknown.result, function(x) {
#       if(is.null(x)) NA else x$log.Lik.smry["median"]
#     })
#     
#     print(log.Lik)
#     lik <- exp(log.Lik - max(log.Lik, na.rm = T))
#     assign.prob <- lik / sum(lik, na.rm = T) 
#     names(assign.prob) <- knowns
#     
#     list(assign.prob = assign.prob, likelihoods = unknown.result)
#   }, simplify = F)
#   
#   assign.prob <- as.data.frame(t(sapply(result, function(x) x$assign.prob)))
#   assign.prob$assignment <- apply(assign.prob, 1, function(x) colnames(assign.prob)[which.max(x)])
#   result <- lapply(result, function(x) x$likelihoods)
#   options(opt)
#   list(assign.prob = assign.prob, likelihoods = result)
}

#' @rdname gelato
#' @export
#' 
gelatoPlot <- function(gelato.result, unknown, main = NULL) { 
#   lik <- gelato.result$likelihoods[[unknown]]
#   lik <- lik[!sapply(lik, is.null)]
#   if(length(lik) == 0) stop(paste("No likelihood distributions available for '", unknown, "'", sep = ""))
#   xticks <- pretty(unlist(sapply(lik, function(x) x$fst.dist)))
#   xlim <- range(xticks)
#   op <- par(mar = c(3, 3, 3, 2) + 0.1, oma = c(2, 2, 0.1, 0.1), mfrow = c(length(lik), 1))
#   high.prob <- gelato.result$assign.prob[unknown, "assignment"]
#   for(known in names(lik)) {
#     known.lik <- lik[[known]]
#     null.max <- max(hist(known.lik$fst.dist[, "null"], plot = F)$density)
#     obs.max <- max(hist(known.lik$fst.dist[, "obs"], plot = F)$density)
#     lik.mean <- known.lik$norm.coefs["mean"]
#     lik.sd <- known.lik$norm.coefs["sd"]
#     ylim <- range(pretty(c(0, null.max, obs.max, dnorm(lik.mean, lik.mean, lik.sd))))
#     hist(known.lik$fst.dist[, "null"], breaks = 10, freq = FALSE, xlim = xlim, ylim = ylim, 
#          xlab = "", ylab = "", main = "", col = "red", xaxt = "n")
#     x <- NULL # To avoid R CMD CHECK warning about no global binding for 'x'
#     curve(dnorm(x, lik.mean, lik.sd), from = xlim[1], to = xlim[2],
#           add = TRUE, col = "black", lwd = 3, ylim = ylim)
#     par(new = TRUE)
#     hist(known.lik$fst.dist[, "obs"], breaks = 10, freq = FALSE, xlim = xlim, ylim = ylim, 
#          xlab = "", ylab = "", col = "darkgreen", main = "", xaxt = "n", yaxt = "n") 
#     axis(1, pretty(xlim))
#     ll.median <- known.lik$log.Lik.smry["median"]
#     log.lik <- if(!is.infinite(ll.median)) format(ll.median, digits = 4) else "Inf"
#     p.val <- format(gelato.result$assign.prob[unknown, known], digits = 2)
#     pop <- paste(known, " (lnL = ", log.lik, ", p = ", p.val, ")", sep = "")
#     mtext(pop, side = 3, line = 1, adj = 1, font = ifelse(known == high.prob, 2, 1))
#   }
#   mtext("Fst", side = 1, outer = T, cex = 1.2)
#   mtext("Density", side = 2, outer = T, cex = 1.2)
#   par(op)
#   if(!is.null(main)) mtext(main, side = 3, line = 3, adj = 0, font = 3)
}