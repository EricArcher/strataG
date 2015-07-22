
###############
# '[' operator
###############
## genind
setMethod("[", signature(x="genind", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., pop=NULL, loc=NULL, treatOther=TRUE, quiet=TRUE, drop=FALSE) {
  
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  
  ## HANDLE 'POP'
  if(!is.null(pop) && !is.null(pop(x))){
    if(is.factor(pop)) pop <- as.character(pop)
    if(!is.character(pop)) pop <- popNames(x)[pop]
    temp <- !pop %in% pop(x)
    if (any(temp)) { # if wrong population specified
      warning(paste("the following specified populations do not exist:", pop[temp]))
    }
    i <- pop(x) %in% pop
  }
  
  ## handle population factor
  if(!is.null(pop(x))) {
    pop <- factor(pop(x)[i])
  } else {
    pop <- NULL
  }
  
  tab       <- tab(x)
  old.other <- other(x)
  hier      <- x@strata
  prevcall  <- match.call()
  
  if (x@type == "codom"){
    ## handle loc argument
    if(!is.null(loc)){
      if(is.factor(loc)) loc <- as.character(loc)
      if(!is.character(loc)) loc <- locNames(x)[loc]
      temp <- !loc %in% locFac(x)
      if (any(temp)) { # if wrong loci specified
        warning(paste("the following specified loci do not exist:", loc[temp]))
      }
      j <- x$loc.fac %in% loc
    } # end loc argument
    if (drop){
      tab    <- tab[i, , ..., drop = FALSE]
      allNb  <- colSums(tab, na.rm=TRUE) # allele absolute frequencies
      toKeep <- (allNb > 1e-10)
      j      <- j & toKeep
      tab    <- tab[, j, ..., drop=FALSE]
    } else {
      tab <- tab[i, j, ..., drop=FALSE]
    }
  } else { # PA case
    tab <- tab[i, j, ..., drop = FALSE]
  }
  
  
  ## handle 'other' slot
  nOther <- length(other(x))
  namesOther <- names(other(x))
  counter <- 0
  if(treatOther){
    f1 <- function(obj,n=nrow(tab(x))){
      counter <<- counter+1
      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
        obj <- obj[i,,drop=FALSE]
      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
        obj <- obj[i]
        if(is.factor(obj)) {obj <- factor(obj)}
      } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}
      
      return(obj)
    } # end f1
    
    x@other <- lapply(other(x), f1) # treat all elements
    
  } else {
    other(x) <- old.other
  } # end treatOther
  
  x@tab    <- tab
  x@pop    <- pop
  x@call   <- prevcall
  x@type   <- x@type
  
  # Treat sample and strata
  x@ploidy    <- ploidy(x)[i]
  x@hierarchy <- x@hierarchy
  x@strata    <- hier[i, , drop = FALSE]
  
  if (x@type == "codom"){
    # Treat locus items
    x <- .drop_alleles(x, j)
  }
  return(x)
})
