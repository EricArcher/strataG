#' @title STRUCTURE
#' @description Run STRUCTURE to assess group membership of samples.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param k.range vector of values to use for \code{maxpops} (the number of
#'   clusters, \emph{K}) in multiple runs. Each value is run
#'   \code{num.k.rep} times. If \code{NULL} (the default), \code{maxpops} is
#'   used if it is supplied in \code{...}, otherwise \code{1:getNumStrata(g)}.
#'   If both \code{k.range} and \code{maxpops} are supplied, \code{k.range} is
#'   used and \code{maxpops} is ignored.
#' @param num.k.rep number of replicates for each value in \code{k.range}.
#' @param label label to use for input and output files
#' @param delete.files logical. Delete all files when STRUCTURE is finished?
#' @param exec name or full path of executable for STRUCTURE. Defaults to "structure".
#' @param ... arguments to be passed to \code{structureWrite}.
#' @param maxpops number of clusters (\emph{K}) for a single run.
#' @param burnin number of iterations for MCMC burnin.
#' @param numreps number of MCMC replicates.
#' @param noadmix logical. Use the no-admixture model, in which each individual
#'   is assumed to come entirely from one cluster? Note that STRUCTURE's own
#'   default is the admixture model (\code{noadmix = FALSE}).
#' @param freqscorr logical. Correlated frequencies?
#' @param randomize logical. If \code{TRUE}, STRUCTURE seeds its random number
#'   generator from the system clock and ignores \code{seed}. Results are then
#'   not reproducible, and replicate runs started within the same second can be
#'   identical, so \code{structureRun} issues a warning. The default
#'   (\code{FALSE}) uses \code{seed}, which is recommended.
#' @param seed integer seed for STRUCTURE's random number generator, used when
#'   \code{randomize = FALSE}. In \code{structureWrite}, the default is drawn
#'   from R's random number generator, so results can be reproduced with
#'   \code{set.seed()}. In \code{structureRun}, each run gets its own seed: if
#'   \code{seed} is supplied in \code{...}, runs use \code{seed},
#'   \code{seed + 1}, \code{seed + 2}, ..., in the order replicates are run
#'   (all replicates of the first value of \code{k.range}, then the next);
#'   otherwise distinct seeds are drawn at random. The seed used for each run
#'   is stored in its \code{seed} element, so any run can be repeated exactly
#'   with \code{structureWrite(..., seed = sr[[i]]$seed)}.
#' @param pop.prior a character specifying which population prior model to use:
#'   "locprior" or "usepopinfo".
#' @param locpriorinit parameterizes locprior parameter \emph{r} - how
#'   informative the populations are. Only used when \code{pop.prior} =
#'   "locprior".
#' @param maxlocprior specifies range of locprior parameter \emph{r}. Only used
#'   when \code{pop.prior} = "locprior".
#' @param gensback integer defining the number of generations back to test for
#'   immigrant ancestry. Only used when \code{pop.prior} = "usepopinfo".
#' @param migrprior numeric between 0 and 1 listing migration prior. Only used
#'   when \code{pop.prior} = "usepopinfo".
#' @param pfrompopflagonly logical. update allele frequencies from individuals
#'   specified by \code{popflag}. Only used when \code{pop.prior} =
#'   "usepopinfo".
#' @param popflag a vector of integers (0, 1) or logicals identifying whether
#'   or not to use strata information. Only used when \code{pop.prior} =
#'   "usepopinfo".
#' @param inferalpha logical. Infer the value of the admixture parameter
#'   \eqn{\alpha} from the data? If \code{FALSE}, \eqn{\alpha} is fixed at the
#'   value of \code{alpha}. Ignored under the no-admixture model
#'   (\code{noadmix = TRUE}). Small \eqn{\alpha} implies that most individuals
#'   are essentially from one population or another, while \eqn{\alpha} > 1
#'   implies that most individuals are admixed.
#' @param alpha Dirichlet parameter for degree of admixture. This is the initial
#'   value if \code{inferalpha = TRUE}.
#' @param unifprioralpha logical. Assume a uniform prior for \code{alpha} which
#'   runs between 0 and \code{alphamax}. This model seems to work fine; the
#'   alternative model (when \code{unifprioralpha = FALSE}) is to take \code{alpha}
#'   as having a Gamma prior, with mean \code{alphapriora × alphapriorb}, and
#'   variance \code{alphapriora × alphapriorb^2}.
#' @param alphamax maximum for uniform prior on \code{alpha} when
#'   \code{unifprioralpha = TRUE}.
#' @param alphapriora,alphapriorb parameters of Gamma prior on \code{alpha} when
#'   \code{unifprioralpha = FALSE}.
#' @param file name of the output file from STRUCTURE.
#' @param pops vector of population labels to be used in place of numbers in
#'   STRUCTURE file.
#' @param ids character vector of individual ids, used in place of the row
#'   numbers written as labels in the STRUCTURE input file.
#'
#' @return
#' \describe{
#'   \item{\code{structureRun}}{a list of class \code{structure.result} with
#'     one element per run, named by run label. Each element is the list
#'     returned by \code{structureRead}, plus:
#'     \describe{
#'       \item{\code{files}}{named vector of the input and output files used by
#'         STRUCTURE, or \code{NULL} if \code{delete.files = TRUE}}
#'       \item{\code{label}}{label for the run}
#'       \item{\code{seed}}{seed used by STRUCTURE for the run}
#'     }}
#'   \item{\code{structureWrite}}{invisibly, a list containing:
#'     \describe{
#'       \item{\code{files}}{named vector of the data, mainparams, extraparams,
#'         and output file names}
#'       \item{\code{pops}}{population labels in the order of the numeric codes
#'         written to the data file, with \code{"Unknown"} appended if any
#'         samples are unstratified}
#'       \item{\code{ids}}{individual ids in the order of the row numbers
#'         written as labels}
#'     }}
#'   \item{\code{structureRead}}{a list containing:
#'     \describe{
#'       \item{\code{summary}}{named numeric vector with \code{k} (number of
#'         clusters), \code{est.ln.prob} (estimated Ln Prob of Data), and
#'         \code{mean.lnL} and \code{var.lnL} (mean and variance of the log
#'         likelihood)}
#'       \item{\code{q.mat}}{data.frame with columns \code{id},
#'         \code{pct.miss}, \code{orig.pop}, and \code{Group.1} ...
#'         \code{Group.k} holding each individual's assignment probabilities}
#'       \item{\code{prior.anc}}{list of prior ancestry estimates for each
#'         individual where population priors were used, or \code{NULL}}
#'     }}
#' }
#'
#' @note STRUCTURE is not included with \code{strataG} and must be downloaded
#'   separately. Additionally, it must be installed such that it can be run from
#'   the command line in the current working directory. See the vignette for
#'   \code{external.programs} for installation instructions.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @references Pritchard, J.K., M. Stephens, P. Donnelly. 2000. Inference of
#'   population structure using multilocus genotype data. Genetics
#'   155:945-959.\cr
#'   \url{http://web.stanford.edu/group/pritchardlab/structure.html}
#'
#' @seealso \code{\link{structurePlot}}, \code{\link{evanno}}, \code{\link{clumpp}}
#'
#' @examples
#' \dontrun{
#' data(msats.g)
#'
#' # Run STRUCTURE
#' sr <- structureRun(msats.g, k.range = 1:4, num.k.rep = 10)
#'
#' # Calculate Evanno metrics
#' evno <- evanno(sr)
#' evno
#'
#' # Run CLUMPP to combine runs for K = 2
#' q.mat <- clumpp(sr, k = 3)
#' q.mat
#'
#' # Plot CLUMPP results
#' structurePlot(q.mat)
#' }
#'
#' @name structure
#'   
NULL

#' @rdname structure
#' @export
#' 
structureRun <- function(
  g, 
  k.range = NULL, 
  num.k.rep = 1, 
  label = NULL, 
  delete.files = TRUE, 
  exec = "structure", 
  ...
) {
  
  if(!getPloidy(g) > 1) stop("'g' must have a ploidy > 1")
  
  # setup folder
  label <- .getFileLabel(g, label)
  label <- paste(label, "structureRun", sep = ".")
  label <- gsub("[[:space:]]", ".", label)
  unlink(label, recursive = TRUE, force = TRUE)
  dir.create(label)
  if(!utils::file_test("-d", label)) {
    stop(paste("'", label, "' is not a valid folder.", sep = ""))
  }
  label <- file.path(label, label)
  
  # setup k and replicate data.frame to cycle through
  dots <- list(...)
  if(is.null(k.range)) {
    k.range <- if(is.null(dots$maxpops)) 1:getNumStrata(g) else dots$maxpops
  }
  rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
  rownames(rep.df) <- paste(label, ".k", rep.df$k, ".r", rep.df$rep, sep = "")
  
  # one distinct seed per run: independent replicates, reproducible with set.seed()
  randomize <- isTRUE(dots$randomize)
  if(randomize) {
    warning("randomize = TRUE: STRUCTURE seeds from the clock, so replicates ",
            "started in the same second can be identical.")
  }
  rep.df$seed <- if(!is.null(dots$seed)) {
    as.integer(dots$seed[1]) + seq_len(nrow(rep.df)) - 1L
  } else {
    sample.int(1e8, nrow(rep.df))
  }
  dots[c("maxpops", "seed", "randomize")] <- NULL
  
  out.files <- lapply(rownames(rep.df), function(x) {
    sw.out <- do.call(structureWrite, c(
      list(
        g, 
        label = x, 
        maxpops = rep.df[x, "k"],
        seed = rep.df[x, "seed"], 
        randomize = randomize
      ),
      dots
    ))
    files <- sw.out$files
    cmd <- paste0(
      exec, 
      " -m ", files["mainparams"], 
      " -e ", files["extraparams"], 
      " -i ", files["data"], 
      " -o ", files["out"]
    )
    
    err.code <- system(cmd)
    if(err.code == 127) {
      stop("You do not have STRUCTURE installed.")
    } else if(!err.code == 0) {
      stop(paste("Error running STRUCTURE. Error code", err.code, "returned."))
    }
    
    files["out"] <- paste(files["out"], "_f", sep = "")
    result <- structureRead(files["out"], sw.out$pops, sw.out$ids)
    
    if(file.exists("seed.txt")) file.remove("seed.txt")
    files <- if(delete.files) NULL else files
    
    result <- c(
      result, 
      list(files = files, label = basename(x), seed = rep.df[x, "seed"])
    )
    fname <- paste(x, ".ws.rdata", sep = "")
    save(result, file = fname)
    fname
  })
  
  run.result <- lapply(out.files, function(f) {
    result <- NULL
    load(f)
    result
  })
  names(run.result) <- sapply(run.result, function(x) x$label)
  class(run.result) <- c("structure.result", class(run.result))
  
  if(delete.files) unlink(dirname(label), recursive = TRUE, force = TRUE)
  run.result
}


#' @rdname structure
#' @export
#' 
structureWrite <- function(
  g, 
  label = NULL, 
  maxpops = getNumStrata(g), 
  burnin = 1000, 
  numreps = 1000, 
  noadmix = TRUE, 
  freqscorr = FALSE, 
  randomize = FALSE, 
  seed = sample.int(1e8, 1),
  pop.prior = NULL, 
  locpriorinit = 1, 
  maxlocprior = 20, 
  gensback = 2,
  migrprior = 0.05,
  pfrompopflagonly = TRUE, 
  popflag = NULL, 
  inferalpha = TRUE,
  alpha = 1.0,
  unifprioralpha = TRUE,
  alphamax = 20,
  alphapriora = 0.05,
  alphapriorb = 0.001,
  ...
) {
  
  if(getPloidy(g) != 2) stop("'g' must be diploid")
  
  # check parameters
  if(!is.null(pop.prior)) {
    if(!pop.prior %in% c("locprior", "usepopinfo")) {
      stop("'pop.prior' must be 'locprior' or 'usepopinfo'.")
    }
  }
  if(is.null(popflag)) popflag <- rep(1, getNumInd(g))
  if(length(popflag) != getNumInd(g)) {
    stop("'popflag' should be the same length as the number of individuals in 'g'.")
  }
  if(!all(popflag %in% c(0, 1))) {
    stop("all values in 'popflag' must be 0 or 1.")
  }
  popflag <- stats::setNames(as.integer(popflag), names(popflag))
  if(is.null(names(popflag))) names(popflag) <- getIndNames(g)
  if(!all(names(popflag) %in% getIndNames(g))) {
    stop("names of 'popflag' must match the individual ids in 'g'.")
  }
  popflag[is.na(getStrata(g)[names(popflag)])] <- 0L
  
  in.file <- ifelse(is.null(label), "data", paste(label, "data", sep = "_"))
  out.file <- ifelse(is.null(label), "out", paste(label, "out", sep = "_"))
  main.file <- ifelse(
    is.null(label), 
    "mainparams", 
    paste(label, "mainparams", sep = "_")
  )
  extra.file <- ifelse(
    is.null(label), 
    "extraparams", 
    paste(label, "extraparams", sep = "_")
  )
  
  # write data
  loci <- getLociNames(g)
  strata.names <- getStrataNames(g)
  unk.code <- length(strata.names) + 1   # code for unstratified samples
  mat <- .stackedAlleles(g, alleles2integer = TRUE, na.val = -9) |> 
    dplyr::select(-"allele") |> 
    dplyr::mutate(
      popflag = popflag[.data$id],
      stratum = as.numeric(factor(.data$stratum, levels = strata.names)),
      stratum = ifelse(is.na(.data$stratum), unk.code, .data$stratum),
      dplyr::across(dplyr::all_of(loci), \(x) replace(x, is.na(x), -9))
    ) |> 
    dplyr::select("id", "stratum", "popflag", dplyr::all_of(loci)) |> 
    as.matrix()
  has.unk <- any(mat[, "stratum"] == unk.code)
  
  # use row numbers as labels: STRUCTURE truncates long labels
  ind.ids <- unique(mat[, "id"])
  mat[, "id"] <- match(mat[, "id"], ind.ids)
  
  write(paste(loci, collapse = " "), file = in.file)
  for(i in 1:nrow(mat)) {
    write(paste(mat[i, ], collapse = " "), file = in.file, append = TRUE)
  }
  
  # write mainparams
  main.params <- c(
    paste("MAXPOPS", as.integer(maxpops)),
    paste("BURNIN", as.integer(burnin)),
    paste("NUMREPS", as.integer(numreps)),
    paste("INFILE", in.file),
    paste("OUTFILE", out.file),
    paste("NUMINDS", getNumInd(g)),
    paste("NUMLOCI", getNumLoci(g)),
    "MISSING -9",
    #"ONEROWPERIND 1",
    "LABEL 1",
    "POPDATA 1",
    "POPFLAG 1",
    "LOCDATA 0",
    "PHENOTYPE 0",
    "EXTRACOLS 0",
    "MARKERNAMES 1"
  )
  main.params <- paste("#define", main.params)
  write(main.params, file = main.file)
  
  # write extraparams
  extra.params <- c(
    paste("NOADMIX", as.integer(noadmix)),
    paste("FREQSCORR", as.integer(freqscorr)),
    paste("INFERALPHA", as.integer(inferalpha)),
    paste("ALPHA", as.numeric(alpha)),
    "FPRIORMEAN 0.01",
    "FPRIORSD 0.05",
    "LAMBDA 1.0",
    paste("UNIFPRIORALPHA", as.integer(unifprioralpha)),
    paste("ALPHAMAX", as.numeric(alphamax)),
    paste("ALPHAPRIORA", as.numeric(alphapriora)),
    paste("ALPHAPRIORB", as.numeric(alphapriorb)),
    "COMPUTEPROB 1",
    paste("ADMBURNIN", max(0, as.integer(burnin / 2))),
    "ALPHAPROPSD 0.025",
    "STARTATPOPINFO 0",
    paste("RANDOMIZE", as.integer(randomize)),
    paste("SEED", as.integer(seed)),
    "METROFREQ 10",
    "REPORTHITRATE 0" 
  )
  
  if(!is.null(pop.prior)) {
    pop.prior <- tolower(pop.prior)
    prior.params <- if(pop.prior == "locprior") {
      c("LOCPRIOR 1",
        "LOCISPOP 1",
        paste("LOCPRIORINIT", locpriorinit),
        paste("MAXLOCPRIOR", maxlocprior)
      )        
    } else if(pop.prior == "usepopinfo") {
      c("USEPOPINFO 1",
        paste("GENSBACK", trunc(gensback)),
        paste("MIGRPRIOR", migrprior),
        paste("PFROMPOPFLAGONLY", as.integer(pfrompopflagonly))
      )
    }
    extra.params <- c(extra.params, prior.params)
  }
  
  extra.params <- extra.params[!is.na(extra.params)]
  extra.params <- paste("#define", extra.params)
  write(extra.params, file = extra.file)
  
  invisible(list(
    files = c(
      data = in.file, 
      mainparams = main.file, 
      extraparams = extra.file, 
      out = out.file
    ),
    pops = if(has.unk) c(strata.names, "Unknown") else strata.names,
    ids = ind.ids
  ))
}


#' @rdname structure
#' @export
#' 
structureRead <- function(file, pops = NULL, ids = NULL) { 
  if(!file.exists(file)) {
    stop(paste("the file '", file, "' can't be found.", sep = ""))
  }
  
  # Read file to get single results and parameter values
  result <- scan(file, "character", quiet = TRUE)
  
  loc <- grep("Estimated", result, ignore.case = FALSE, value = FALSE)
  est.ln.prob <- as.numeric(result[loc[1] + 6])
  
  loc <- grep("likelihood", result, ignore.case = FALSE, value = FALSE)
  mean.lnL <- as.numeric(result[loc[1] + 2])
  var.lnL <- as.numeric(result[loc[2] + 2])
  
  loc <- grep("MAXPOPS", result, value = F)
  maxpops <- result[loc]
  maxpops <- sub("MAXPOPS=", "", maxpops)
  maxpops <- as.integer(sub(",", "", maxpops))
  
  loc <- grep("GENSBACK", result, value = F)
  gensback <- result[loc]
  gensback <- sub("GENSBACK=", "", gensback)
  gensback <- as.integer(sub(",", "", gensback))
  
  smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, var.lnL = var.lnL)
  
  # Read file to get population assignment probability table
  result <- scan(file, "character", sep = "\n", quiet = TRUE)
  first <- grep("(%Miss)", result, value = FALSE) + 1
  last <- grep("Estimated Allele", result, value = FALSE) - 1
  tbl.txt <- result[first:last]
  
  # Remove special characters from table and whitespace from end of lines
  tbl.txt <- sub("[*]+", "", tbl.txt)
  tbl.txt <- sub("[(]", "", tbl.txt)
  tbl.txt <- sub("[)]", "", tbl.txt)
  tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
  
  # Find which lines have had population priors
  prior.lines <- grep("[|]", tbl.txt)
  
  # Create table of population assignments for lines without population priors
  no.prior <- if(length(prior.lines) < length(tbl.txt)) {
    no.prior.q.txt <- if(length(prior.lines) == 0) tbl.txt else tbl.txt[-prior.lines]
    .structureParseQmat(no.prior.q.txt, pops, ids) 
  } else NULL
  
  # Return just this base table if MAXPOPS = 1
  if(maxpops == 1) {
    no.prior$row <- NULL
    return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
  }
  
  # Create table of population assignments for lines with population priors
  has.prior <- if(length(prior.lines) > 0) {
    prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
    
    # Get base table
    prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
    df <- .structureParseQmat(prior.q.txt, pops, ids) 
    
    # Parse ancestry assignments into matrix
    prior.anc <- purrr::map(prior.txt, function(x) {
      anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
      rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat), sep = ".")
      colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
      # Split on whitespace and colons
      x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
        y <- y[y != ""] # remove empty strings
        y[-1] # return vector with first element ("Pop") removed - vector has population # and ancestry assignments
      })
      # Populate ancestry matrix with probabilities
      for(i in 1:ncol(x)) {
        pop <- as.numeric(x[1, i])
        anc.mat[pop, ] <- as.numeric(x[-1, i])
      }
      anc.mat
    }) |> 
      stats::setNames(df$id)
    
    # Create population probability matrix for samples with priors and add to end of base table
    prob.mat <- t(sapply(1:nrow(df), function(i) {
      pop.probs <- rowSums(prior.anc[[i]])
      pop.probs[is.na(pop.probs)] <- df$Group.1[i]
      pop.probs
    }))
    colnames(prob.mat) <- paste("Group", 1:ncol(prob.mat), sep = ".")
    df$Group.1 <- NULL
    df <- cbind(df, prob.mat)
    
    list(df = df, prior.anc = prior.anc)
  } else NULL
  
  # Combine assignment probability matrices
  has.prior.df <- if(is.null(has.prior)) NULL else has.prior$df
  q.mat <- rbind(no.prior, has.prior.df)
  q.mat <- q.mat[order(q.mat$row), ]
  q.mat$row <- NULL
  rownames(q.mat) <- NULL
  # Make sure all probs sum to 1
  q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i / sum(i)))
  
  prior.anc <- if(is.null(has.prior)) NULL else has.prior$prior.anc
  
  list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
}


# Internal file used by 'stucture' to parse output files
# q.mat.txt: character vector of Q matrix from STRUCTURE output file.
# pops: vector of population labels to be used in place of numbers in STRUCTURE file.
#   
.structureParseQmat <- function(q.mat.txt, pops, ids = NULL) {
  q.mat.txt <- sub("[*]+", "", q.mat.txt)
  q.mat.txt <- sub("[(]", "", q.mat.txt)
  q.mat.txt <- sub("[)]", "", q.mat.txt)
  q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
  
  cols1to4 <- c("row", "id", "pct.miss", "orig.pop")
  
  strsplit(q.mat.txt, " ") |> 
    purrr::map(function(q) {
      q <- q[! q %in% c("", " ", ":")] |> 
        as.character() |> 
        rbind() |> 
        as.data.frame(stringsAsFactors = FALSE) 
      stats::setNames(q, c(cols1to4, paste("Group", 1:(ncol(q) - 4), sep = ".")))
    }) |> 
    dplyr::bind_rows() |> 
    dplyr::mutate_at(
      dplyr::vars("row", "pct.miss", "orig.pop", dplyr::starts_with("Group.")), 
      as.numeric
    ) |> 
    dplyr::mutate(orig.pop = if(!is.null(pops)) {
      pops[.data$orig.pop] 
    } else {
      .data$orig.pop
    }) |> 
    dplyr::mutate(id = if(!is.null(ids)) ids[as.integer(.data$id)] else .data$id)
}
