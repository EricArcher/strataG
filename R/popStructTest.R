#' @name popStructTest 
#' @title Population Structure Tests
#' @description Conduct multiple tests of population structure / differentiation. 
#'   Overall tests can be conducted for the current stratification scheme 
#'   (\code{overallTest()}), or can be conducted for all unique 
#'   pairs of strata (\code{pairwiseTest()}). All statistics appropriate to 
#'   the ploidy of the data are estimated at once. See \code{Note} for a 
#'   description of each statistic.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param by.locus return by-locus values of statistics? If \code{TRUE} the overall 
#'   value will be contained in the first row, labelled \code{"All"}. Only applies 
#'   if the ploidy of \code{g} is > 1 (non-haploid).
#' @param hap.locus which locus to use if \code{g} is haploid. Can be specified 
#'   by number or name. 
#' @param max.cores the maximum number of cores to use to distribute replicates
#'   for permutation tests over. If set to \code{NULL}, the value will be what is 
#'   reported by \code{\link[parallel]{detectCores} - 1}. If \code{detectCores} 
#'   reports \code{NA}, \code{max.cores} will be set to 1 and parallel 
#'   processing will not be done. 
#' @param quietly logical. print progress to screen?
#' @param ... parameters passed to \code{\link[ape]{dist.dna}} for 
#'  computation of pairwise distance matrix for AMOVA PHIst statistic.
#' @param pws a list returned from a call to \code{pairwiseTest()}.
#' @param stat the name of a statistic in the \code{$result} element for pairwise
#'   comparisons returned by \code{pairwiseTest()}.
#' @param locus the name of a single locus. If \code{"All"}, the overall result
#'   from all loci is returned. See \code{by.locus}.
#' 
#' @return \describe{
#'  \item{\code{overallTest()}}{a list containing: \describe{
#'    \item{\code{strata.freq}}{a table of the sample sizes for each stratum}
#'    \item{\code{result}}{an array with the statistic estimate and p-value 
#'      for each statistic. If \code{by.locus = FALSE} or \code{g} is a haploid dataset, 
#'      this is a two-dimensional array, with one row per statistic, 
#'      statistic estimate in the first column and permutation test p-value 
#'      in the second column. If \code{by.locus = TRUE} and \code{g} has ploidy > 1,
#'      then this is a three-dimensional array where the first dimension 
#'      is loci, second dimension is statistics, and third dimension is 
#'      statistic estimate and p-value.}
#'  }}
#'  \item{\code{pairwiseTest()}}{a list containing a list of results as described above 
#'    for \code{overallTest()} for each pairwise comparison.}
#'  \item{\code{pairwiseMatrix()}}{a matrix summarizing a chosen statistic 
#'    (\code{stat}) for a chosen locus (\code{locus}) between pairs of strata 
#'    with the statistic estimate in the lower left and the p-value in the upper right.}
#'  \item{\code{pairwiseSummary()}}{a data frame summarizing all pairwise 
#'    statistics and p-values along with strata sample sizes.}  
#' }
#' 
#' @note The computed statistics are:\tabular{ll}{
#'   \code{CHIsq} \tab chi-squared estimate measuring random allele frequency
#'     distribution distributed across strata (haploid and diploid) \cr
#'   \code{Ho, Hs, Ht} \tab Nei and Kumar 2002 : 
#'     observed heterozygosity (\code{Ho}), 
#'     within population diversity (\code{Hs}), overall diversity (\code{Ht}) 
#'     \cr
#'   \code{Ht_prime} \tab description \cr
#'   \code{Dst} \tab description \cr
#'   \code{Dst_prime} \tab description \cr
#'   \code{Fst} \tab For haploid data, equivalent to PHIst with pairwise 
#'     distances set to 1. For diploid data,  \cr
#'   \code{Fst_prime} \tab description \cr
#'   \code{Fis} \tab description \cr
#'   \code{Gst_prime} \tab description \cr
#'   \code{Gst_dbl_prime} \tab description \cr
#'   \code{Dest, Dest_Chao} \tab population differentiation (Jost 2008) \cr 
#'   \code{wcFit, wcFst, wcFit} \tab (Weir and Cockerham 1984) \cr
#'   \code{PHIst} \tab Haploid AMOVA estimate of differentiation derived 
#'     from matrix of pairwise distances between sequences. 
#'     See \code{\link[ape]{dist.dna}} for details on distance computation. 
#'     (Excoffier et al 1992) \cr
#' }
#'   
#' @seealso \code{\link[hierfstat]{basic.stats}}, \code{\link[pegas]{Fst}},
#'   \code{\link[pegas]{amova}}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references
#' Excoffier, L., Smouse, P.E. and Quattro, J.M. 1992. Analysis of 
#'   molecular variance inferred from metric distances among DNA haplotypes: 
#'   application to human mitochondrial DNA restriction data. 
#'   Genetics 131:479–491.  
#' Jost, L. 2008. GST and its relatives do not measure differentiation. 
#'   Molecular Ecology 17:4015-4026.  
#' Nei M. and Chesser R. 1983. Estimation of fixation indexes and gene 
#'   diversities. Annals of Human Genetics 47:253-259.  
#' Nei M. 1987. Molecular Evolutionary Genetics. Columbia University Press
#' Weir, B.S. and Cockerham, C.C. 1984. Estimating F-statistics for 
#'   the analysis of population structure. Evolution 38:1358–1370.  
#' Weir, B.S. and Hill, W.G. 2002. Estimating F-statistics. 
#'   Annual Review of Genetics 36:721–750.  
#' 
#' @examples 
#' # An overall test with microsatellite data
#' data(msats.g)
#' ovl <- overallTest(msats.g, nrep = 100)
#' ovl
#' 
#' #' A pairwise test with control region sequences
#' data(dloop.g)
#' pws <- pairwiseTest(dloop.g, nrep = 100)
#' pws
#'
NULL

#' @noRd
#' 
.dataPrep <- function(g, nrep = 0, hap.locus = 1, ...) {
  # format data for C input
  dt <- g@data %>% 
    dplyr::filter(!is.na(allele) & !is.na(stratum)) %>% 
    dplyr::transmute(
      id = as.numeric(factor(id)),
      stratum = factor(stratum),
      locus = factor(locus),
      allele = factor(paste0(locus, "_", allele))
    )
  ploidy <- strataG::getPloidy(g)
  
  # check data
  if(nrow(dt) == 0) stop(
    "no individuals left after removing missing genotypes ",
    "and unstratified individuals."
  )
  st.freq <- table(dt$stratum) / ploidy
  if(any(st.freq == 1)) {
    has.one <- names(st.freq)[st.freq == 1]
    stop(
      "the following strata have only one individual: ", 
      paste(has.one, collapse = ", ")
    )
  }
  
  st <- makeLookupVec(dt$id, as.numeric(dt$stratum) - 1)
  dt.list <- list(
    n_ind = dplyr::n_distinct(dt$id),
    n_allele = nlevels(dt$allele),
    n_loc = nlevels(dt$locus),
    n_st = dplyr::n_distinct(st),
    ploidy = ploidy,
    ind_allele_freq = table2D_int(dt$id, as.numeric(dt$allele) - 1),
    strata = cbind(
      st, 
      if(nrep > 0) replicate(nrep, sample(st)) else NULL, deparse.level = 0
    ),
    strata_names = levels(dt$stratum)
  )
  
  if(ploidy == 1) {
    if(is.numeric(hap.locus)) hap.locus <- getLociNames(g)[hap.locus]
    haps <- getAlleleNames(g)[[hap.locus]]
    hap.names <- paste0(hap.locus, "_", haps)
    
    # format distances for Fst (all 1s, and 0s on the diagonal)
    dt.list$unit_dist <- matrix(
      1, nrow = length(haps), ncol = length(haps), 
      dimnames = list(hap.names, hap.names)
    )
    diag(dt.list$unit_dist) <- 0    
    
    # get haplotype distances for PHIst
    dna.seq <- getSequences(g, as.haplotypes = TRUE, seqName = hap.locus)
    dt.list$hap_dist <- if(!is.null(dna.seq)) {
      hd <- ape::dist.dna(dna.seq, as.matrix = TRUE, ...)
      hd <- hd[haps, haps, drop = FALSE] ^ 2
      dimnames(hd) <- list(hap.names, hap.names)
      hd
    } else NULL
  } else {
    dt.list$locus <- makeLookupVec(
      as.numeric(dt$allele) - 1, 
      as.numeric(dt$locus) - 1
    )
    dt.list$loci_names <- levels(dt$locus)
  }
  
  dt.list
}

#' @noRd
#' 
.statChisq <- function(dt, s) {
  obsvd <- do.call(
    rbind,
    tapply(1:nrow(dt$strata), dt$strata[, s], function(i) {
      colSums(dt$ind_allele_freq[i, ])
    })
  )
  exptd <- tcrossprod(rowSums(obsvd), colSums(obsvd)) / sum(obsvd)
  chisq <- (obsvd - exptd) ^ 2 / exptd
  chisq <- if(dt$n_loc > 1) {
    by_locus <- tapply(
      1:length(dt$locus), 
      dt$locus, 
      function(i) sum(chisq[, i])
    )
    stats::setNames(
      c(sum(by_locus), by_locus),
      c("All", dt$loci_names)
    )
  } else sum(chisq)
  
  cbind(CHIsq = chisq)
}

#' @noRd
#' 
.statFstStats <- function(x, s) {
  result <- calc_FstStats(
    x$n_ind, x$n_allele, x$n_loc, x$n_st, x$ploidy, 
    x$strata[, s], x$locus, x$ind_allele_freq
  )
  rownames(result) <- c("All", x$loci_names)
  result
}

#' @noRd
#' 
.statPhist <- function(x, s) {
  calc_Phist(
    x$n_ind, x$n_allele, x$n_st, 
    x$strata[, s], x$ind_allele_freq, x$unit_dist, x$hap_dist
  )
}

#' @noRd
#' 
.runPopStructFuncs <- function(s, data, by.locus = TRUE) {
  mat <- if(data$ploidy == 1) {
    rbind(c(.statChisq(data, s)[, 1], .statPhist(data, s)))
  } else {
    result <- cbind(.statChisq(data, s), .statFstStats(data, s))
    if(by.locus) result else result[1, , drop = FALSE]
  }
  mat[is.nan(mat)] <- NA
  mat
}

#' @rdname popStructTest
#' @export
#' 
overallTest <- function(g, nrep = 1000, by.locus = FALSE, hap.locus = 1,
                        quietly = FALSE, max.cores = 1,  ...) {
  # check replicates
  if(!is.numeric(nrep) & length(nrep) != 1) {
    stop("'nrep' must be a single-element numeric vector")
  }
  if(nrep < 0) stop("'nrep' must be a positive number")
  
  # create input data (overwrite 'g')
  desc <- getDescription(g)
  g <- .dataPrep(g, nrep = nrep, hap.locus = hap.locus, ...)
  
  # compute observed statistics
  obs <- .runPopStructFuncs(1, g, by.locus)    
  
  # create permutation null distribution
  p.val <- if(nrep == 0) {
    x <- obs
    x[] <- NA
    x
  } else {
    # run permutation tests in parallel
    cl <- swfscMisc::setupClusters(nrep, max.cores)
    
    if(!quietly) {
      message("\n<<< ", desc, " >>>")
      message(format(Sys.time()), " : Overall tests (", nrep, " permutations)")
    }
    
    null.dist <- tryCatch({
      if(is.null(cl)) {
        lapply(2:nrep, .runPopStructFuncs, data = g, by.locus = by.locus)
      } else {
        parallel::clusterEvalQ(cl, require(strataG))
        parallel::clusterExport(cl, "g", environment())
        parallel::parLapplyLB(
          cl, 2:nrep, .runPopStructFuncs, data = g, by.locus = by.locus
        )
      }
    }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)
    if(!quietly) message(format(Sys.time()), " : Done")
    null.dist <- simplify2array(c(list(obs), null.dist))
    apply(null.dist, c(1, 2), function(x) mean(x >= x[1], na.rm = TRUE))
  }
  p.val[is.nan(p.val)] <- NA
  result <- simplify2array(list(estimate = obs, p.val = p.val))
  # simplify result if only one locus is present
  if(g$ploidy == 1 | !by.locus) result <- result[1, , ]
  
  list(
    strata.freq = table(g$strata_names[g$strata[, 1] + 1]), 
    result = result
  )
}

#' @rdname popStructTest
#' @export
#' 
pairwiseTest <- function(g, nrep = 1000, by.locus = FALSE, hap.locus = 1,
                         quietly = FALSE, max.cores = 1, ...) { 
  
  if(getNumStrata(g) == 1) stop("'g' must have more than one stratum defined.")
  
  if(!quietly) {
    message("\n<<< ", getDescription(g), " >>>")
    message(format(Sys.time()), " : Pairwise tests (", nrep, " permutations)")
  }
  
  # create strata pairs
  strata.pairs <- .strataPairs(g)
  
  # run permutation test on all pairwise gtypes subsets
  pws <- vector("list", length = nrow(strata.pairs))
  for(i in 1:nrow(strata.pairs)) {
    pair <- unlist(strata.pairs[i, ])
    if(!quietly) {
      message(format(Sys.time()), " : ", paste(pair, collapse = " v. "))
    }
    pws[[i]] <- overallTest(
      g[, , pair], nrep = nrep, by.locus = by.locus, hap.locus = hap.locus,
      quietly = TRUE, max.cores = max.cores, ...
    )
  }
  if(!quietly) message(format(Sys.time()), " : Done")
  
  pws
}

#' @rdname popStructTest
#' @export
#' 
pairwiseMatrix <- function(pws, stat, locus = "All") {
  res <- pws[[1]]$result
  
  # check that stat is present
  stat.dim <- length(dim(res)) - 1
  if(!stat %in% dimnames(res)[[stat.dim]]) {
    stop("statistic '", stat, "' not present.")
  }
  
  # check that locus is present
  locus.dim <- if(length(dim(res)) == 2) NULL else 1
  if(!is.null(locus.dim)) {
    if(!locus %in% dimnames(res)[[locus.dim]]) {
      stop("locus '", locus, "' not present.")
    }
  }
  
  # extract strata names and create empty matrix                          
  st.freqs <- do.call(
    rbind,
    lapply(pws, function(res) as.data.frame(res$strata.freq))
  )
  st.freqs <- unique(st.freqs)
  mat <- matrix(NA, nrow = nrow(st.freqs), ncol = nrow(st.freqs))
  rownames(mat) <- colnames(mat) <- st.freqs[, 1]
  
  # create and return matix
  for(res in pws) {
    pair <- names(res$strata.freq)
    res <- if(stat.dim == 1) res$result else res$result[locus, , ]
    mat[pair[2], pair[1]] <- res[stat, "estimate"]
    mat[pair[1], pair[2]] <- res[stat, "p.val"]
  }
  mat
}

#' @rdname popStructTest
#' @export
#' 
pairwiseSummary <- function(pws, locus = "All") {
  res <- pws[[1]]$result
  stat.dim <- length(dim(res)) - 1
  stats <- dimnames(res)[[stat.dim]]
  
  # check that locus is present
  locus.dim <- if(length(dim(res)) == 2) NULL else 1
  if(!is.null(locus.dim)) {
    if(!locus %in% dimnames(res)[[locus.dim]]) {
      stop("locus '", locus, "' not present.")
    }
  }
  
  result <- lapply(pws, function(pws) {
    df <- data.frame(
      strata.1 = names(pws$strata.freq)[1],
      strata.2 = names(pws$strata.freq)[2],
      n.1 = pws$strata.freq[1],
      n.2 = pws$strata.freq[2]
    )
    rownames(df) <- NULL
    mat <- if(is.null(locus.dim)) {
      t(pws$result)
    } else {
      t(pws$result[locus, , ])
    }
    cbind(
      label = paste0(
        df$strata.1, " (", df$n.1, ") v. ",
        df$strata.2, " (", df$n.2, ")"
      ),
      df,
      rbind(stats::setNames(
        c(mat),
        c(rbind(stats, paste0(stats, "_p.val")))
      ))
    )
  })
  
  do.call(rbind, result)
}