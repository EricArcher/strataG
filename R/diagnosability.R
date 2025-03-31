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
#' @param unk vector of strata to be treated as "unknowns" for prediction with 
#'   Random Forest model.
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
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
#' diagnosability(dloop.g, pairwise = TRUE)
#' }
#' 
#' @export
#' 
diagnosability <- function(g, gene = 1, pairwise = FALSE, conf.level = 0.95, 
                     unk = NULL, ...) {
  # check for and set names to sampsize
  args <- list(...)
  ss.match <- which(pmatch(names(args), "sampsize") == 1)
  if(length(ss.match) == 1) {
    temp.ss <- args[[ss.match]]
    if(is.null(names(temp.ss))) {
      names(temp.ss) <- getStrataNames(g)
      args[[ss.match]] <- temp.ss
      names(args)[ss.match] <- "sampsize"
    }
  }
  
  if(pairwise) {
    sp <- .strataPairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- g[, , unlist(sp[i, ])]
      do.call(diagnosability, c(list(g = pair.g, pairwise = FALSE), args))
    })
    list(
      smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
      rp = lapply(result, function(x) x$rp)
    )
  } else {
    seq.df <- .gtypes2rfDF(g, gene = gene)
    rf.df <- if(!is.null(unk)) {
      df <- seq.df[!seq.df$stratum %in% unk, ]
      df$stratum <- droplevels(df$stratum)
      df
    } else seq.df
    args <- c(list(x = rf.df), args)
    result <- do.call(.sequenceRF, args)
    result$seq.df <- seq.df
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
#'
.gtypes2rfDF <- function(g, gene = 1, label = NULL) {
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
    var.seq.mat <- stats::setNames(
      cbind(var.seq.mat[, to.keep]),
      sites[to.keep]
    )
    
    # create factors of variable site columns
    seq.df <- do.call(
      data.frame, 
      lapply(colnames(var.seq.mat), function(x) factor(var.seq.mat[, x]))
    ) |> 
      stats::setNames(var.seq.mat)
    
    # add strata and ids
    cbind(
      df[, c('id', 'stratum')], 
      seq.df[df$id, ]
    ) |>
      tibble::column_to_rownames('id')
  } else {
    snp.df <- as.data.frame(g, one.col = TRUE, ids = TRUE, strata = TRUE) |> 
      tibble::column_to_rownames('id')
    all.biallelic <- all(
      sapply(snp.df[, -1], function(x) length(unique(x)) <= 3)
    )
    if(!all.biallelic) warning("some loci in 'g' may not be biallelic")
    snp.df
  }
  
  # add strata and remove any rows with missing data
  rf.df <- stats::na.omit(rf.df)
  
  # remove columns where substitutions are represented by only one individual
  preds <- rf.df[, -1, drop = FALSE]
  to.keep <- apply(preds, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) {
    warning("all predictors are either constant or variable in only one individual. NULL returned")
    return(NULL)
  }
  
  # format and return data frame
  st <- if(is.null(label)) rf.df$stratum else paste(label, rf.df$stratum)
  do.call(
    data.frame, 
    lapply(preds[, to.keep, drop = FALSE], factor)
  ) |> 
    dplyr::mutate(
      id = rownames(rf.df),
      stratum = factor(st)
    ) |> 
    tibble::column_to_rownames('id') |> 
    dplyr::select('stratum', dplyr::everything()) 
}


#' @keywords internal
#'
.sequenceRF <- function(x, replace = FALSE, sampsize = NULL, 
                       train.pct = 0.5, min.n = 2, nrep = 0, 
                       conf.level = 0.95, ...) {
  if(is.null(x)) return(NULL)
  strata <- x$stratum
  if(length(unique(strata)) < 2) return(NULL)
  
  if(is.null(sampsize)) sampsize <- rfPermute::balancedSampsize(strata)
  sampsize <- ifelse(sampsize < min.n, min.n, sampsize)
  rp <- rfPermute::rfPermute(
    stratum ~ ., data = x, replace = replace, 
    sampsize = sampsize, nrep = nrep, ...
  )
  
  rf <- rfPermute::as.randomForest(rp)
  overall.accuracy <- 1 - as.vector(rf$err.rate[nrow(rf$err.rate), "OOB"])
  ci <- rfPermute::confusionMatrix(rp, conf.level = conf.level)
  ci <- ci[, (length(rf$classes) + 1):(length(rf$classes) + 3)]
  min.diag <- which.min(ci[-nrow(ci), 1])
  diag.strata <- rownames(ci)[min.diag]
  
  dna.seq <- apply(x[, -1, drop = FALSE], 1, paste, collapse = "")
  most.freq <- sort(table(dna.seq), decreasing = T)
  haps <- as.numeric(factor(dna.seq, levels = names(most.freq)))
  hap.freq <- prop.table(table(haps, x$stratum), 1)
  pred <- sapply(haps, function(i) {
    colnames(hap.freq)[which.max(hap.freq[as.character(i), ])]
  })
  
  class.tbl <- table(
    factor(strata), 
    factor(pred, levels = unique(strata))
  )
  overall.diag <- sum(diag(class.tbl)) / sum(class.tbl)
  class.tbl <- cbind(
    class.tbl, 
    diagnosability = diag(class.tbl) / rowSums(class.tbl)
  )
  
  predicted <- apply(rf$votes, 1, function(x) colnames(rf$votes)[which.max(x)])
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = predicted[y.i] == rf$y[y.i])
  }, simplify = F)
  min.p <-  1 / nlevels(rf$y)
  pd.vec <- if(is.null(pd.vec)) min.p else c(min.p, pd.vec)
  pd.vec <- pd.vec[pd.vec >= min.p]
  pd95 <- sapply(pd.vec, function(p) {
    strata.pd <- sapply(vote.dfs, function(df) {
      mean(df$is.correct & df$votes >= p)
    })
    unname(strata.pd[which.min(strata.pd)])
  }) |> 
    stats::setNames(paste("pd",  pd.vec, sep = "."))
  
  smry <- data.frame(
    overall.accuracy = overall.accuracy * 100,
    diag.strata = diag.strata,
    diagnosability = ci[min.diag, 1],
    diag.lci = ci[min.diag, 2],
    diag.uci = ci[min.diag, 3],
    shared.hap.diag = class.tbl[diag.strata, "diagnosability"],
    pd95 = pd95[2],
    num.vs = ncol(x) - 1,
    num.haps = length(unique(haps)),
    hap.div = sprex::diversity(haps, type = "unb.gini"),
    eff.num.haps = sprex::diversity(haps, type = "effective.number")
  )
  rownames(smry) <- NULL
  
  list(smry = smry, rp = rp)
}