#' @title Find and label haplotypes
#' @description Identify and group sequences that share the same haplotype.
#'
#' @param x sequences in a \code{character matrix}, \code{list},
#'   or \code{\link{DNAbin}} object, or a haploid \linkS4class{gtypes} object
#'   with sequences.
#' @param prefix a character string giving prefix to be applied to numbered
#'   haplotypes. If NULL, haplotypes will be labeled with the first label
#'   from original sequences.
#' @param use.indels logical. Use indels when comparing sequences?
#' @param ... arguments to be passed to \code{labelHaplotypes.default}.
#'
#' @details If any sequences contain ambiguous bases (N's) they are first
#'   removed. Then haplotypes are assigned based on the remaining
#'   sequences. The sequences with N's that were removed are then assigned to
#'   the new haplotypes if it can be done unambiguously (they match only one
#'   haplotype with 0 differences once the N's have been removed). If this
#'   can't be done they are assigned NAs and listed in the
#'   \code{unassigned} element.
#'
#' @return
#'  For \code{character}, \code{list}, or \code{DNAbin}, a list with the following elements:
#'  \describe{
#'    \item{haps}{named vector (\code{DNAbin}) or list of named vectors
#'      (\code{multidina}) of haplotypes for each sequence in \code{x}.}
#'    \item{hap.seqs}{\code{DNAbin} or \code{multidna} object containing
#'      sequences for each haplotype.}
#'    \item{unassigned}{\code{data.frame} listing closest matching haplotypes
#'      for unassignable sequences with N's and the minimum number of
#'      substitutions between the two. Will be \code{NULL} if no sequences
#'      remain unassigned.}
#'  }
#'
#'  For \code{gtypes}, a new \code{gtypes} object with unassigned individuals
#'  stored in the @@other slot in an element named \code{'haps.unassigned'} (see
#'  \code{\link{getOther}}).
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{expandHaplotypes}}
#'
#' @examples
#' # create 5 example short haplotypes
#' haps <- c(
#'   H1 = "ggctagct",
#'   H2 = "agttagct",
#'   H3 = "agctggct",
#'   H4 = "agctggct",
#'   H5 = "ggttagct"
#' )
#
#' # draw and label 100 samples
#' sample.seqs <- sample(names(haps), 100, rep = TRUE)
#' ids <- paste(sample.seqs, 1:length(sample.seqs), sep = "_")
#' sample.seqs <- lapply(sample.seqs, function(x) strsplit(haps[x], "")[[1]])
#' names(sample.seqs) <- ids
#'
#' # add 1-2 random ambiguities
#' with.error <- sample(1:length(sample.seqs), 10)
#' for(i in with.error) {
#'   num.errors <- sample(1:2, 1)
#'   sites <- sample(1:length(sample.seqs[[i]]), num.errors)
#'   sample.seqs[[i]][sites] <- "n"
#' }
#'
#' hap.assign <- labelHaplotypes(sample.seqs, prefix = "Hap.")
#' hap.assign
#'
#' @name labelHaplotypes
#' @export
#'
labelHaplotypes <- function(x, prefix = NULL, use.indels = TRUE) {
  UseMethod("labelHaplotypes")
}


#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.default  <- function(x, prefix = NULL, use.indels = TRUE) {
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object.")
  x <- as.matrix(x)
  
  # return same data if only one sequence exists
  if(nrow(x) == 1) {
    haps <- rownames(x)
    names(haps) <- haps
    return(list(haps = haps, hap.seqs = x, unassigned = NULL))
  }
  
  # throw error if any sequence names are duplicated
  if(any(duplicated(rownames(x)))) {
    warning(
      "'x' cannot have duplicate sequence names. NULL returned.",
      call. = FALSE, 
      immediate. = TRUE
    )
    return(NULL)
  }
  
  # find sequences without ambiguities  
  good.bases <- c("a", "c", "g", "t", "-")
  has.amb <- apply(as.character(x), 1, function(bases) {
    !all(tolower(bases) %in% good.bases)
  })
  if(sum(!has.amb) <= 1) {
    warning(
      "There are fewer than two sequences without ambiguities (N's). Can't assign haplotypes. NULL returned.",
      call. = FALSE, 
      immediate. = TRUE
    )
    return(NULL)
  }
  
  # get pairwise distances and set all non-0 distances to 1
  x.no.amb <- x[!has.amb, , drop = FALSE]
  hap.dist <- ape::dist.dna(x.no.amb, model = "N", pairwise.deletion = TRUE)
  if(use.indels) {
    hap.dist <- hap.dist + ape::dist.dna(x.no.amb, model = "indelblock")
  }
  hap.dist <- as.matrix(hap.dist)
  hap.dist[hap.dist > 0] <- 1
  
  # create haplotype code out of 0s and 1s
  hap.code <- hap.dist %>% 
    apply(1, paste, collapse = "") %>% 
    factor() %>% 
    as.numeric() %>% 
    stats::setNames(rownames(hap.dist))
  
  # rename haplotypes
  hap.labels <- if(!is.null(prefix)) {
    # use prefix+number if prefix given
    # sort based on frequency first
    hap.order <- as.numeric(names(sort(table(hap.code), decreasing = TRUE)))
    hap.nums <- swfscMisc::zero.pad(1:length(hap.order))
    names(hap.order) <- paste(prefix, hap.nums, sep = "")
    names(sort(hap.order))
  } else {
    # if no prefix, use first sequence name for each haplotype
    names(hap.code)[!duplicated(hap.code)]
  }
  hap.code <- hap.labels[hap.code]
  names(hap.code) <- rownames(hap.dist)
  
  # get sequences for each haplotype
  unique.codes <- hap.code[!duplicated(hap.code)]
  hap.seqs <- x[names(unique.codes), , drop = FALSE]
  rownames(hap.seqs) <- unique.codes
  hap.seqs <- hap.seqs[order(rownames(hap.seqs)), , drop = FALSE] %>% 
    as.matrix() %>% 
    as.character()
  
  # get distance of all sequences with n's to other sequences (possible matching sequences)
  x.has.amb <- as.character(x)[has.amb, , drop = FALSE]
  unk.dist <- sapply(rownames(x.has.amb), function(i) {
    # get ambiguity sites in this sequence
    amb.sites <- which(!tolower(x.has.amb[i, ]) %in% good.bases)
    if(length(amb.sites) == 0) return(NULL)
    # create matrix of ambiguity sites and remove sequences with ambiguities
    possible.sites <- unique(hap.seqs[, amb.sites, drop = FALSE])
    # create sequences to test
    if(nrow(possible.sites) == 0) return(NULL)
    test.seqs <- do.call(rbind, lapply(1:nrow(possible.sites), function(j) {
      seq.j <- x.has.amb[i, ]
      seq.j[amb.sites] <- possible.sites[j, ]
      seq.j
    }))
    test.names <- paste("test.seq.", 1:nrow(test.seqs), sep = "")
    rownames(test.seqs) <- test.names
    # get distances between test sequences and sequences without n's
    new.mat <- ape::as.DNAbin(rbind(hap.seqs, test.seqs))
    hap.dist <- ape::dist.dna(
      new.mat, 
      model = "N", 
      pairwise.deletion = TRUE, 
      as.matrix = TRUE
    )
    if(use.indels) hap.dist <- hap.dist + 
      ape::dist.dna(
        new.mat, 
        model = "indelblock", 
        as.matrix = TRUE
      )
    hap.dist[test.names, !colnames(hap.dist) %in% test.names, drop = FALSE]
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  min.dist <- NULL
  hap.assign <- if(length(unk.dist) > 0) {
    # get minimum distance of unknowns to other sequences
    min.dist <- do.call(rbind, sapply(
      unk.dist, 
      function(mat) apply(mat, 2, min), 
      simplify = FALSE
    ))
    # assign sequences if only one of the minimum sequences is 0
    apply(min.dist, 1, function(counts) {
      num.matches <- sum(counts == 0)
      if(num.matches == 1) colnames(min.dist)[counts == 0] else NA
    })
  } else NULL
  
  # compile vector of haplotype assignments
  hap.vec <- c(hap.code, hap.assign)[rownames(x)]
  
  # create data.frame of unassigned sequences
  unassigned.df <- if(any(is.na(hap.vec))) {
    mat <- min.dist[names(hap.vec)[is.na(hap.vec)], , drop = FALSE]
    df <- do.call(rbind, lapply(1:nrow(mat), function(r) {
      min.dist <- min(mat[r, ])
      i <- which(mat[r, ] == min.dist)
      haps <- paste(colnames(mat)[i], collapse = ", ")
      data.frame(haplotype = haps, min.substitutions = min.dist)
    }))
    rownames(df) <- rownames(mat)
    df
  } else NULL
  
  list(
    haps = hap.vec, 
    hap.seqs = ape::as.DNAbin(hap.seqs), 
    unassigned = unassigned.df
  )
}

#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.list <- function(x, ...) labelHaplotypes(ape::as.DNAbin(x),...)


#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.character <- function(x, ...) labelHaplotypes(ape::as.DNAbin(x), ...)


#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.gtypes <- function(x, ...) {
  # check that sequences are present
  if(getPloidy(x) > 1 | is.null(getSequences(x))) {
    stop("'x' is not haploid or does not have any sequences")
  }

  # label haplotypes for each gene
  new.haps <- purrr::map(getSequences(x), labelHaplotypes, ...)
  has.errors <- sapply(new.haps, is.null)
  if(sum(has.errors) > 0) {
    has.errors <- paste(names(new.haps)[has.errors], collapse = ", ")
    stop("haplotypes could not be assigned for: ", has.errors)
  }

  # create haplotype data.frame
  hap.df <- as.data.frame(x)
  for(gene in names(new.haps)) {
    old.haps <- hap.df[, gene]
    hap.df[, gene] <- new.haps[[gene]]$haps[old.haps]
  }

  # collect sequences
  hap.seqs <- lapply(new.haps, function(x) x$hap.seqs)
  names(hap.seqs) <- names(new.haps)

  # collect unassigned
  unassigned <- lapply(new.haps, function(x) x$unassigned)

  # create new gtypes
  df2gtypes(
    hap.df, 
    ploidy = 1, 
    id.col = 1, 
    strata.col = 2, 
    loc.col = 3,
    sequences = hap.seqs, 
    description = getDescription(x), 
    schemes = getSchemes(x),
    other = c(getOther(x), list(haps.unassigned = unassigned))
  )
}
