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
#' 
#' @return data.frame summary of pairwise shared loci.
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
#' propSharedLoci(msats.g)
#' }
#' 
#' @name sharedLoci
#' 
NULL

#' @rdname sharedLoci
#' @export
#' 
propSharedLoci <- function(g, type = c("strata", "ids")) {
  g <- .checkHapsLabelled(g)
  
  type <- match.arg(type)
  
  .formatResult <- function(type.pairs, prop.shared) {
    as.data.frame(type.pairs, stringsAsFactors = FALSE) %>% 
      stats::setNames(paste(type, 1:2, sep = ".")) %>% 
      dplyr::mutate(
        num.loci.genotyped = rowSums(!is.na(prop.shared)),
        num.loci.shared = rowSums(prop.shared == 1, na.rm = TRUE),
        prop.loci.shared = .data$num.loci.shared / .data$num.loci.genotyped
      ) %>% 
      cbind(prop.shared) %>% 
      as.data.frame()
  }
  
  smry <- if(type == "ids") {
    if(getNumInd(g) == 1) stop("'g' must have at least two individuals")
    af <- table(g@data[, c("id", "locus", "allele")])
    id.pairs <- t(utils::combn(dimnames(af)[[1]], 2))
    
    # proportion of alleles at each locus present in both individuals
    prop.shared <- t(sapply(1:nrow(id.pairs), function(i) {
      x <- af[id.pairs[i, ], , ]
      is.shared <- x[1, ,] > 0 & x[2, , ] > 0
      prop <- (rowSums(is.shared * x[1, , ]) + rowSums(is.shared * x[2, , ]))
      not.genotyped <- apply(apply(x, c(1, 2), sum), 2, function(z) any(z == 0))
      prop[not.genotyped] <- NA
      prop
    })) / (2 * getPloidy(g))
    
    .formatResult(id.pairs, prop.shared)
  } else {
    if(getNumStrata(g) == 1) {
      stop("'type' cannot be 'strata' if only one stratum is present.")
    }
    na <- numAlleles(g)
    na <- stats::setNames(na$num.alleles, na$locus)
    af <- table(g@data[, c("stratum", "locus", "allele")])
    strata.freq <- .strataFreq(g)
    strata.pairs <- t(utils::combn(dimnames(af)[[1]], 2))
    
    # proportion of alleles at each locus present in both strata
    prop.shared <- t(sapply(1:nrow(strata.pairs), function(i) {
      is.shared <- af[strata.pairs[i, 1], ,] > 0 & af[strata.pairs[i, 2], , ] > 0
      prop <- rowSums(is.shared) / na
      not.genotyped <- apply(
        apply(af[strata.pairs[i, ], , ], c(1, 2), sum), 
        2, 
        function(z) any(z == 0)
      )
      prop[not.genotyped] <- NA
      prop
    }))
    
    .formatResult(strata.pairs, prop.shared)
  }
  
  unassigned <- getOther(g, "haps.unassigned")
  if(!is.null(unassigned)) attr(smry, "gtypes") <- g
  
  smry
}


#' @rdname sharedLoci
#' @export
#'  
sharedAlleles <- function(g, smry = c("num", "which")) {
  g <- .checkHapsLabelled(g)
  
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
  
  shared <- data.frame(st.pairs, shared, stringsAsFactors = FALSE)
  
  unassigned <- getOther(g, "haps.unassigned")
  if(!is.null(unassigned)) attr(shared, "gtypes") <- g
  
  shared
}