#' @name sharedLoci
#' @title Shared Loci
#' @description Calculate proportion of alleles and number of loci shared 
#'   between pairs of individuals or strata.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param type a character vector determining type of pairwise comparsion. Can 
#'   be "strata" for strata or "ids" for individuals.
#' @param smry a character vector determining type of summary for 
#'   \code{sharedAlleles}. "which" returns the names of the alleles shared. 
#'   "num" returns the number of alleles shared.
#' @param ids character vector of two sample ids to compare.
#' @param num.cores number of CPU cores to use.
#' 
#' @return data.frame summary of pairwise shared loci.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' propSharedLoci(msats.g)
#' 
#' sharedAlleles(msats.g)
#' 
#' @importFrom utils combn
#' @importFrom stats na.omit
#' @importFrom parallel parLapply stopCluster
#' @export
#' 
propSharedLoci <- function(g, type = c("strata", "ids"), num.cores = 1) {
  type <- match.arg(type)
  type.pairs <- if(type == "strata") {
    as.matrix(.strataPairs(g))
  } else {
    id.vec <- indNames(g)
    if(length(id.vec) < 2) stop("'g' must have at least two individuals")
    t(combn(id.vec, 2))
  }
  cl <- .setupClusters(num.cores)

  shared <- tryCatch({
    prop.func <- if(type == "strata") {
      function(f, type.pair) {
        pair.f <- f[, "prop", type.pair]
        sum(apply(pair.f, 1, function(x) !all(x == 0))) / nrow(pair.f)
      }
    } else {
      function(i, type.pairs, g) propSharedIds(type.pairs[i, ], g)
    }
    
    if(type == "strata") {  
      freqs <- alleleFreqs(g, TRUE)
      sharedStrataList <- lapply(1:nrow(type.pairs), function(i) {
        prop.shared.loci <- if(num.cores > 1) {
          parLapply(cl, freqs, prop.func, type.pair = type.pairs[i, ])
        } else {
          lapply(freqs, prop.func, type.pair = type.pairs[i, ])
        }
        prop.shared.loci <- unlist(prop.shared.loci)
        names(prop.shared.loci) <- names(freqs)
        prop.shared.loci
      })
      do.call(rbind, sharedStrataList)
    } else {
      sharedIDlist <- if(num.cores > 1) {
        parLapply(cl, 1:nrow(type.pairs), prop.func, type.pairs = type.pairs, g = g)
      } else {
        lapply(1:nrow(type.pairs), prop.func, type.pairs = type.pairs, g = g)
      }
      do.call(rbind, sharedIDlist)
    }
  }, finally = if(num.cores > 1) stopCluster(cl))
  
  shared.summary <- do.call(rbind, lapply(1:nrow(shared), function(i) {  
    num.same <- sum(na.omit(shared[i, ]) == 1) 
    num.not.missing <- sum(!is.na(shared[i, ]))
    prop.same <- num.same / num.not.missing
    c(num.same = num.same, num.not.missing = num.not.missing, 
      prop.same = prop.same)
  }))
  shared.summary <- cbind(as.data.frame(type.pairs, stringsAsFactors = FALSE), 
                          shared.summary, shared)
  colnames(shared.summary)[1:2] <- paste(type, 1:2, sep = ".")
  shared.summary
}


#' @rdname sharedLoci
#' @importFrom utils combn
#' @export
#'  
sharedAlleles <- function(g, smry = "num") {
  st.vec <- levels(strata(g))
  if(length(st.vec) < 2) stop("'g' must have at least two strata")
  st.pairs <- t(combn(st.vec, 2))
  colnames(st.pairs) <- c("strata.1", "strata.2")
  freqs <- alleleFreqs(g, TRUE)
  
  shared <- do.call(rbind, lapply(1:nrow(st.pairs), function(i) {
    which.shared <- sapply(freqs, function(f) {
      pair.f <- f[, "prop", st.pairs[i, ]] 
      if(any(is.na(pair.f))) return(NA)
      is.shared <- apply(pair.f, 1, function(x) all(x != 0))
      if(sum(is.shared) > 0) {
        if(smry == "which") {
          paste(rownames(pair.f)[is.shared], collapse = ", ")
        } else sum(is.shared)
      } else NA
    }, USE.NAMES = TRUE)
    names(which.shared) <- names(freqs)
    which.shared
  }))
  
  st.pairs <- data.frame(st.pairs, stringsAsFactors = FALSE)
  cbind(st.pairs, shared)
}


#' @rdname sharedLoci
#' @export
#' 
propSharedIds <- function(ids, g) {
  sapply(locNames(g), function(x) {
    g1 <- unlist(loci(g, ids[1], x))
    g2 <- unlist(loci(g, ids[2], x))
    sum(g1 %in% g2) + sum(g2 %in% g1)
  }) / (2 * ploidy(g))
}
