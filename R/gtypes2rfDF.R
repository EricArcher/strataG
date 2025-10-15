#' @title Create Sequence Data.Frame
#' @description Create data.frame of variable sites from gtypes object.
#' 
#' @param g a \code{\link[strataG]{gtypes}} object. If haploid, it must have 
#'   aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} slot.
#' @param label label to add to beginning of each stratum name.
#' 
#' @return a \code{data.frame} where the first column lists the (\code{strata}) 
#'   and every column afterwards is a variable site. All columns are factors.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
gtypes2rfDF <- function(g, gene = 1, label = NULL) {
  rf.df <- if(getPloidy(g) == 1) {
    if(is.null(getSequences(g))) stop("'g' must have aligned sequences")
    
    # extract stratified sequences
    if(is.numeric(gene)) gene <- getLociNames(g)[gene]
    df <- g |> 
      as.data.frame() |> 
      stats::na.omit() |> 
      dplyr::select(dplyr::all_of(c("id", "stratum", gene)))
    dna.seqs <- g |> 
      getSequences(seqName = gene) |> 
      as.matrix()
    
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
    
    # create factors of variable site columns
    seq.df <- do.call(
      data.frame, 
      lapply(colnames(var.seq.mat), function(x) factor(var.seq.mat[, x]))
    )
    colnames(seq.df) <- colnames(var.seq.mat)
    
    # add strata and Ids
    seq.df <- data.frame(cbind(stratum = df$stratum, seq.df[df$id, ]))
    rownames(seq.df) <- df$id
    seq.df
  } else {
    snp.df <- as.data.frame(g, one.col = TRUE, ids = TRUE, strata = TRUE)
    rownames(snp.df) <- snp.df$id
    snp.df$id <- NULL    
    all.biallelic <- all(
      sapply(snp.df[, -1], function(x) length(unique(x)) <= 3)
    )
    if (!all.biallelic) warning("some loci in 'g' may not be biallelic")
    snp.df
  }
  
  # add strata and remove any rows with missing data
  rf.df <- stats::na.omit(rf.df)
  
  # remove columns where substitutions are only represented by one individual
  preds <- rf.df[, -1, drop = FALSE]
  to.keep <- apply(preds, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) {
    warning("all predictors are either constant or variable in only one individual. NULL returned")
    return(NULL)
  }
  preds <- preds[, to.keep, drop = FALSE]
  
  # convert to factor
  preds <- do.call(data.frame, lapply(preds, factor))
  
  st <- if(is.null(label)) rf.df$stratum else paste(label, rf.df$stratum)
  result <- cbind(stratum = factor(st), preds)
  rownames(result) <- rownames(rf.df)
  result
}
