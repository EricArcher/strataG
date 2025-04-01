#' @title Diagnosability
#' @description Conduct Random Forest on a gtypes object to compute the
#'   diagnosability of each stratum (PD from Archer et al 2017).
#' 
#' @param g haploid \code{\link{gtypes}} object with aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} 
#'   slot. Defaults to the first gene in the object.
#' @param pairwise do analysis on all pairwise combinations of strata?
#' @param conf.level confidence level for the \code{\link{binom.test}} 
#'   confidence interval.
#' @param replace sample with replacement in Random Forest trees? 
#'    (see \code{\link[randomForest]{randomForest}}).
#' @param sampsize sample size for each Random Forest tree? 
#'    (see \code{\link[randomForest]{randomForest}}). If \code{NULL} a 
#'    balanced sample size is chosen 
#'    (see \code{\link[rfPermute]{balancedSampsize}}).
#' @param train.pct if \code{sampsize} is \code{NULL}, the percent of the 
#'    minimum strata size to use for \code{sampsize}.
#' @param min.n minimum sample size across all strata.
#' @param min.votes.pct numeric vector giving the minimum percent of votes for  
#'    the assigned strata for a sample to be considered correctly assigned.
#' @param rp.nrep number of replicates for \code{rfPermute} computation of 
#'    significance of site importance scores.
#' @param unk vector of strata to be treated as "unknowns" for prediction with 
#'   Random Forest model.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), 
#'   and the \code{randomForest} object (\code{rf}). If \code{pairwise} 
#'   is \code{TRUE} then the \code{rf} element is a list of 
#'   \code{randomForest} results for each row in \code{smry}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' \dontrun{
#' library(strataG)
#' data(dloop.g)
#' 
#' pd <- diagnosability(dloop.g, pairwise = TRUE)
#' 
#' lapply(pd, function(x) x$rf.confusion.mat)
#' }
#' 
#' @export
#' 
#' 
diagnosability <- function(g, gene = 1, pairwise = FALSE, conf.level = 0.95, 
                           replace = FALSE, sampsize = NULL, train.pct = 0.5,
                           min.n = 2, min.votes.pct = c(0.8, 0.9, 0.95), rp.nrep = 0,
                           unk = NULL) {

  arg.list <- as.list(environment())
  
  if(pairwise) {
    sp <- .strataPairs(g)
    arg.list$pairwise <- FALSE
    lapply(1:nrow(sp), function(i) {
      arg.list$g <- g[, , unlist(sp[i, ])]
      do.call(diagnosability, arg.list)
    })
  } else {
    seq.df <- .gtypes2rfDF(g, gene = gene)
    rf.df <- if(!is.null(unk)) {
      df <- seq.df[!seq.df$stratum %in% unk, ]
      df$stratum <- droplevels(df$stratum)
      df
    } else seq.df
    arg.list[c('g', 'gene', 'pairwise', 'unk')] <- NULL
    result <- do.call(.sequenceRF, c(list(x = rf.df), arg.list))
    result$rf.df <- rf.df
    result$unk <- if(!is.null(unk) & !is.null(result$rp)) {
      unk.df <- seq.df[seq.df$stratum %in% unk, ]
      pred <- data.frame(stats::predict(result$rp, unk.df, type = "prob"))
      pred$predicted <- stats::predict(result$rp, unk.df, type = "response")
      pred
    } else NULL
    result
  }
}


#' @keywords internal
#' @noRd
#'
.gtypes2rfDF <- function(g, gene = 1) {
  rf.df <- if(getPloidy(g) == 1) {
    if(is.null(getSequences(g))) stop("'g' must have aligned sequences")
    
    # extract stratified sequences
    if(is.numeric(gene)) gene <- getLociNames(g)[gene]
    df <- stats::na.omit(as.data.frame(g))[, c("id", "stratum", gene)]
    dna.seqs <- as.matrix(getSequences(g, seqName = gene))
    
    # extract variable sites for these sequences and create sequence matrix
    var.sites <- variableSites(dna.seqs)
    var.seq.mat <- tolower(
      do.call(rbind, as.character(as.list(var.sites$sites)))
    )
    sites <- paste("site", colnames(var.sites$site.freqs), sep = ".")
    
    # only use sites with valid bases
    to.keep <- apply(var.seq.mat, 2, function(x) {
      all(x %in% c("a", "c", "g", "t", "-", "."))
    })
    if(sum(to.keep) == 0) return(NULL)
    var.seq.mat <- cbind(var.seq.mat[, to.keep])
    colnames(var.seq.mat) <- sites[to.keep]
    rownames(var.seq.mat) <- df$id
    
    # create factors of variable site columns, add strata and ids, and return data frame
    var.seq.mat |> 
      as.data.frame() |> 
      dplyr::mutate(
        dplyr::across(dplyr::everything(), factor),
        stratum = df$stratum
      ) |>
      dplyr::select(stratum, dplyr::everything())
  } else {
    snp.df <- as.data.frame(g, one.col = TRUE, ids = TRUE, strata = TRUE) |> 
      tibble::column_to_rownames('id')
    all.biallelic <- all(
      sapply(snp.df[, -1], function(x) length(unique(x)) <= 3)
    )
    if(!all.biallelic) warning("some loci in 'g' may not be biallelic")
    snp.df
  }
  
  # remove any rows with missing data
  rf.df <- rf.df |> 
    stats::na.omit() |> 
    dplyr::mutate(dplyr::across(dplyr::everything(), factor))
  
  # remove columns where substitutions are represented by only one individual
  preds <- rf.df[, -1, drop = FALSE]
  to.keep <- apply(preds, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) {
    warning("all predictors are either constant or variable in only one individual. NULL returned")
    return(NULL)
  }
  
  # format and return data frame
  preds[, to.keep, drop = FALSE] |> 
    dplyr::mutate(stratum = rf.df$stratum) |>  
    dplyr::select('stratum', dplyr::everything()) 
}


#' @keywords internal
#' @noRd
#'
.sequenceRF <- function(x, replace = FALSE, sampsize = NULL, 
                        train.pct = 0.5, min.n = 2, rp.nrep = 0, 
                        conf.level = 0.95, min.votes.pct = 0.95) {
  if(is.null(x)) return(NULL)
  strata <- x$stratum
  if(length(unique(strata)) < 2) return(NULL)
  
  if(is.null(sampsize)) sampsize <- rfPermute::balancedSampsize(strata, train.pct)
  sampsize <- ifelse(sampsize < min.n, min.n, sampsize)
  rp <- rfPermute::rfPermute(
    stratum ~ ., data = x, replace = replace, 
    sampsize = sampsize, nrep = rp.nrep,
    num.cores = ifelse(rp.nrep == 0, 1, parallel::detectCores() - 1)
  )
  
  # Random Forest classifier
  rf <- rfPermute::as.randomForest(rp)
  predicted <- apply(rf$votes, 1, function(x) colnames(rf$votes)[which.max(x)])
  min.votes.pct <- min.votes.pct[min.votes.pct >= 1 / nlevels(rf$y)]
  if(length(min.votes.pct) == 0) min.votes.pct <- 0.95
  min.votes.pct <- stats::setNames(min.votes.pct, paste("pd",  min.votes.pct, sep = "."))
  strata.pd <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    is.correct = predicted[y.i] == rf$y[y.i]
    sapply(min.votes.pct, function(p) {
      100 * mean(is.correct & votes >= p)
    })
  })
  
  # Haplotype frequency classifier
  dna.seq <- apply(x[, -1, drop = FALSE], 1, paste, collapse = "")
  most.freq <- sort(table(dna.seq), decreasing = T)
  haps <- as.numeric(factor(dna.seq, levels = names(most.freq)))
  hap.freq <- prop.table(table(haps, x$stratum), 1)
  pred <- sapply(haps, function(i) {
    colnames(hap.freq)[which.max(hap.freq[as.character(i), ])]
  })
  hap.class.tbl <- table(
    factor(strata), 
    factor(pred, levels = sort(unique(strata)))
  )
  overall.diag <- 100 * sum(diag(hap.class.tbl)) / sum(hap.class.tbl)
  hap.class.tbl <- cbind(
    hap.class.tbl, 
    pct.correct = 100 * diag(hap.class.tbl) / rowSums(hap.class.tbl)
  )
  
  list(
    rf.confusion.mat = rfPermute::confusionMatrix(rp, conf.level = conf.level),
    strata.pd = strata.pd,
    hap.freq.confusion.mat = hap.class.tbl,
    num.vs = ncol(x) - 1,
    rp = rp
  )
}