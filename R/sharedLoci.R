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
#' @param num.cores number of CPU cores to use. Defaults to the number reported 
#'  by \code{\link[parallel]{detectCores}} - 1.
#' 
#' @return data.frame summary of pairwise shared loci.
#' 
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' sharedAlleles(msats.g)
#' 
#' \dontrun{
#' propSharedLoci(msats.g, num.cores = 2)
#' }
#' 
#' @name sharedLoci
#' 
NULL

#' @rdname sharedLoci
#' @export
#' 
propSharedLoci <- function(g, type = c("strata", "ids"), num.cores = NULL) {
  type <- match.arg(type)
  if(type == "strata" & getNumStrata(g) == 1) {
    stop("'type' cannot be 'strata' if only one stratum is present.")
  }
  type.pairs <- if(type == "strata") {
    as.matrix(.strataPairs(g))
  } else {
    id.vec <- getIndNames(g)
    if(length(id.vec) < 2) stop("'g' must have at least two individuals")
    t(utils::combn(id.vec, 2))
  }
  
  cl <- .setupClusters(num.cores)
  
  propFunc <- if(type == "strata") {
    function(f, type.pair, ...) {
      mean(apply(f[, type.pair], 1, function(x) all(x > 0)))
    }
  } else {
    function(i, type.pairs, g) {
      prop.shared <- .propSharedIds(type.pairs[i, ], g)
      stats::setNames(prop.shared$prop.shared, prop.shared$locus)
    }
  }
  
  shared <- tryCatch({
    pair.i <- 1:nrow(type.pairs)
    if(type == "strata") { 
      freqs <- alleleFreqs(g, by.strata = TRUE)
      lapply(pair.i, function(i) {
        prop.shared.loci <- if(!is.null(cl)) {
          parallel::parLapply(cl, freqs, propFunc, type.pair = type.pairs[i, ])
        } else {
          lapply(freqs, propFunc, type.pair = type.pairs[i, ])
        }
        stats::setNames(unlist(prop.shared.loci), names(freqs))
      })
    } else {
      cat("Comparing", nrow(type.pairs), "pairs of individuals...\n")
      if(!is.null(cl)) {
        parallel::parLapply(cl, pair.i, propFunc, type.pairs = type.pairs, g = g)
      } else {
        lapply(pair.i, propFunc, type.pairs = type.pairs, g = g)
      }
    }
  }, finally = if(!is.null(cl)) parallel::stopCluster(cl))
  shared <- do.call(rbind, shared)
  
  shared.summary <- do.call(rbind, lapply(1:nrow(shared), function(i) {  
    num.same <- sum(stats::na.omit(shared[i, ]) == 1) 
    num.not.missing <- sum(!is.na(shared[i, ]))
    prop.same <- num.same / num.not.missing
    c(
      num.loci.shared = num.same, 
      num.loci.genotyped = num.not.missing, 
      prop.loci.shared = prop.same
    )
  }))
  
  type.pairs %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    stats::setNames(paste(type, 1:2, sep = ".")) %>% 
    cbind(shared.summary, shared)
}


#' @rdname sharedLoci
#' @export
#'  
sharedAlleles <- function(g, smry = c("num", "which")) {
  st.pairs <- .strataPairs(g)
  if(is.null(st.pairs)) stop("'g' must have at least two strata")
  smry <- match.arg(smry)
  st.pairs <- as.matrix(st.pairs)
  freqs <- alleleFreqs(g, TRUE)
  
  shared <- do.call(rbind, lapply(1:nrow(st.pairs), function(i) {
    result <- sapply(freqs, function(f) {
      pair.f <- f[, st.pairs[i, ], drop = FALSE] 
      if(any(is.na(pair.f))) return(NA)
      is.shared <- apply(pair.f, 1, function(x) all(x != 0))
      num.shared <- sum(is.shared)
      if(smry == "which") {
        if(num.shared > 0) {
          paste(rownames(pair.f)[is.shared], collapse = ", ")
        } else NA
      } else num.shared
    }, USE.NAMES = TRUE)
    names(result) <- names(freqs)
    result
  }))
  
  return(data.frame(st.pairs, shared, stringsAsFactors = FALSE))
}


#' @rdname sharedLoci
#' @param ids character vector of two sample ids to compare.
#' @param g a \code{gtypes} object
#' @keywords internal
#' 
.propSharedIds <- function(ids, g) {
  g@data %>% 
    dplyr::group_by(.data$locus) %>% 
    dplyr::do({
      gt1 <- dplyr::filter(.data, .data$id == ids[1]) %>% 
        dplyr::pull(.data$allele)
      gt2 <- dplyr::filter(.data, .data$id == ids[2]) %>% 
        dplyr::pull(.data$allele)
      if(any(is.na(c(gt1, gt2)))) {
        data.frame(prop.shared = NA)
      } else {
        data.frame(
          prop.shared = sum(gt1 %in% gt2, gt2 %in% gt1) / (2 * ploidy(g))
        )
      }
    }) %>% 
    dplyr::ungroup() %>% 
    as.data.frame()
}
