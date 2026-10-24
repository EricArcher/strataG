#' @title Run CLUMPP
#' @description Native R implementation of CLUMPP (Jakobsson & Rosenberg 2007)
#'   to align and aggregate multiple STRUCTURE runs from \code{structureRun()}.
#'
#' @param sr result from \code{structureRun()} (class
#'   \code{structure.result}).
#' @param k choice of \emph{k} in \code{sr} to combine.
#' @param align.algorithm algorithm to be used for aligning the runs. Can be
#'   "full.search", "greedy", or "large.k".
#' @param sim.stat pairwise matrix similarity statistic to be used. Can be
#'   "g" or "g.prime".
#' @param greedy.option input order of runs to be tested. Required if
#'   \code{align.algorithm} is "greedy" or "large.k". Valid choices are
#'   "all" (all R! input orders) or "ran.order" (\code{repeats} random
#'   input orders).
#' @param repeats the number of input orders of runs to be tested. Only used
#'   if \code{align.algorithm} is "greedy" or "large.k", and
#'   \code{greedy.option} is "ran.order".
#' @param order.by.run permute the clusters according to the cluster order of
#'   a specific run. Set this parameter to a number from 1 to the number of
#'   runs at \code{k} in \code{sr}. 0 keeps the cluster order of the first
#'   run in the input order that maximized H.
#'
#' @return a data.frame with columns \code{id}, \code{pct.miss},
#'   \code{orig.pop}, \code{Group.1} ... \code{Group.k} holding the mean of
#'   the aligned Q-matrices, as returned by \code{strataG::clumpp()}. The
#'   attribute \code{"clumpp"} additionally holds H, the column permutation
#'   applied to each run, and the settings used.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#' 
#' @references Jakobsson, M. and N. A. Rosenberg. 2007. CLUMPP: a cluster
#'   matching and permutation program for dealing with label switching and
#'   multimodality in analysis of population structure.
#'   Bioinformatics 23(14):1801-1806.
#'   
#' @seealso \code{\link{structure}}
#'
#' @export
#' 
clumpp <- function(sr, k, align.algorithm = "greedy", sim.stat = "g",
                   greedy.option = "ran.order", repeats = 100,
                   order.by.run = 0) {
  
  if(k < 2) stop("k must be greater than 1.")
  align.algorithm <- tolower(align.algorithm)
  sim.stat <- tolower(sim.stat)
  greedy.option <- tolower(greedy.option)
  if(!align.algorithm %in% c("full.search", "greedy", "large.k")) {
    stop("'align.algorithm' must be either 'full.search', 'greedy', or 'large.k'.")
  }
  if(!sim.stat %in% c("g", "g.prime")) {
    stop("'sim.stat' must be either 'g' or 'g.prime'.")
  }
  if(!"structure.result" %in% class(sr)) {
    stop("'sr' must be a result from 'structureRun'.")
  }
  
  eq.k <- vapply(sr, function(x) isTRUE(x$summary[["k"]] == k), logical(1))
  if(sum(eq.k) == 0) stop(paste("no entries for k =", k, "found in 'sr'."))
  sr <- sr[eq.k]
  num.runs <- length(sr)
  
  if(!is.numeric(order.by.run) || length(order.by.run) != 1 ||
     order.by.run < 0 || order.by.run > num.runs ||
     order.by.run != round(order.by.run)) {
    stop("'order.by.run' must be 0 or an integer from 1 to ", num.runs, ".")
  }
  
  # Q-matrices with rows in the id order of the first run
  id.order <- sr[[1]]$q.mat$id
  q.list <- lapply(sr, .clumppGetQ, id.order = id.order, k = k)
  
  # Find column permutations for each run
  result <- if(num.runs == 1) {
    list(perms = list(seq_len(k)), H = NA_real_)
  } else if(align.algorithm == "full.search") {
    .clumppFullSearch(q.list, sim.stat)
  } else {
    run.orders <- .clumppRunOrders(num.runs, greedy.option, repeats)
    .clumppGreedy(q.list, run.orders, sim.stat,
                  large.k = align.algorithm == "large.k")
  }
  perms <- result$perms
  
  # Reorder clusters to match a specific run
  if(order.by.run > 0) {
    ref <- order(perms[[order.by.run]])
    perms <- lapply(perms, function(p) p[ref])
  }
  
  # Mean of aligned Q-matrices
  aligned <- mapply(function(q, p) q[, p, drop = FALSE], q.list, perms,
                    SIMPLIFY = FALSE)
  mean.q <- Reduce(`+`, aligned) / num.runs
  colnames(mean.q) <- paste("Group", seq_len(k), sep = ".")
  
  q1 <- sr[[1]]$q.mat
  q.mat <- data.frame(
    id = as.character(q1$id),
    pct.miss = as.numeric(q1$pct.miss),
    orig.pop = as.character(q1$orig.pop),
    mean.q,
    stringsAsFactors = FALSE
  )
  rownames(q.mat) <- NULL
  
  names(perms) <- names(sr)
  attr(q.mat, "clumpp") <- list(
    H = result$H,
    perms = perms,
    align.algorithm = align.algorithm,
    sim.stat = sim.stat,
    num.runs = num.runs
  )
  q.mat
}


# ---- Ancillary functions ---------------------------------------------------

# Extract a run's Q-matrix (Group.* columns) with rows ordered by id.order
.clumppGetQ <- function(x, id.order, k) {
  q.mat <- x$q.mat
  i <- match(id.order, q.mat$id)
  if(anyNA(i)) {
    stop("Run '", x$label, "' does not contain the same individuals as the first run.")
  }
  q <- as.matrix(q.mat[i, grep("^Group\\.", names(q.mat)), drop = FALSE])
  if(ncol(q) != k) stop("Run '", x$label, "' does not have ", k, " groups.")
  storage.mode(q) <- "double"
  dimnames(q) <- NULL
  q
}

# All permutations of 1:n as rows of an integer matrix; row 1 is 1:n
.clumppPerms <- function(n) {
  if(n == 1) return(matrix(1L, 1, 1))
  sub <- .clumppPerms(n - 1)
  out <- do.call(
    rbind, 
    lapply(seq_len(n), function(i) {
      cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow(sub)))
    })
  )
  storage.mode(out) <- "integer"
  dimnames(out) <- NULL
  out
}

# Input orders of runs to test for the greedy algorithms
.clumppRunOrders <- function(num.runs, greedy.option, repeats) {
  if(!greedy.option %in% c("all", "ran.order")) {
    stop("'greedy.option' must be either 'all' or 'ran.order'.")
  }
  if(greedy.option == "all") {
    if(factorial(num.runs) > 40320) {
      warning("testing all ", factorial(num.runs), " input orders of ",
              num.runs, " runs; consider greedy.option = 'ran.order'.")
    }
    .clumppPerms(num.runs)
  } else {
    if(is.null(repeats)) {
      stop("'repeats' must be specified if 'greedy.option' is 'ran.order'.")
    }
    t(vapply(seq_len(repeats), function(i) sample.int(num.runs),
             integer(num.runs)))
  }
}

# Similarity of A to B[, perm] for every row of perm.mat (vectorised).
# Uses ||A - BP||^2 = ||A||^2 + ||B||^2 - 2 * sum_k (A'B)[k, perm[k]],
# and ||BP - W|| = ||B - W|| since W is constant.
.clumppSimAll <- function(A, B, perm.mat, sim.stat) {
  k <- ncol(A)
  n <- nrow(perm.mat)
  M <- crossprod(A, B)
  inner <- rowSums(matrix(M[cbind(rep(seq_len(k), each = n), as.vector(perm.mat))], n))
  d <- sqrt(pmax(sum(A^2) + sum(B^2) - 2 * inner, 0))
  denom <- if(sim.stat == "g") {
    w <- 1 / k
    sqrt(sqrt(sum((A - w)^2)) * sqrt(sum((B - w)^2)))
  } else {
    sqrt(2 * nrow(A))
  }
  if(denom == 0) return(ifelse(d == 0, 1, 0))
  1 - d / denom
}

# Similarity of two already-aligned Q-matrices
.clumppSim <- function(A, B, sim.stat) {
  .clumppSimAll(A, B, matrix(seq_len(ncol(B)), 1), sim.stat)
}

# H: mean pairwise similarity among aligned Q-matrices
.clumppH <- function(aligned, sim.stat) {
  pairs <- utils::combn(length(aligned), 2)
  mean(apply(pairs, 2, function(ij) {
    .clumppSim(aligned[[ij[1]]], aligned[[ij[2]]], sim.stat)
  }))
}

# FullSearch: first run fixed, every combination of the K! permutations of
# the remaining R - 1 runs, keep the combination with the highest H.
.clumppFullSearch <- function(q.list, sim.stat) {
  k <- ncol(q.list[[1]])
  num.runs <- length(q.list)
  perm.mat <- .clumppPerms(k)
  num.perms <- nrow(perm.mat)
  num.comb <- num.perms ^ (num.runs - 1)
  if(num.comb > 1e7) {
    stop("FullSearch would test ", format(num.comb, big.mark = ","),
         " combinations of permutations. Use 'greedy' or 'large.k'.")
  }
  
  # Pairwise similarity S[[i, j]][a, b] = sim(Q_i[, P_a], Q_j[, P_b])
  pairs <- utils::combn(num.runs, 2)
  sim.tbl <- lapply(seq_len(ncol(pairs)), function(m) {
    i <- pairs[1, m]
    j <- pairs[2, m]
    rows <- if(i == 1) 1L else seq_len(num.perms)
    out <- matrix(NA_real_, num.perms, num.perms)
    for(a in rows) {
      out[a, ] <- .clumppSimAll(q.list[[i]][, perm.mat[a, ], drop = FALSE],
                                q.list[[j]], perm.mat, sim.stat)
    }
    out
  })
  
  comb <- as.matrix(expand.grid(rep(list(seq_len(num.perms)), num.runs - 1)))
  comb <- cbind(1L, comb)
  total <- numeric(nrow(comb))
  for(m in seq_len(ncol(pairs))) {
    total <- total + sim.tbl[[m]][cbind(comb[, pairs[1, m]], comb[, pairs[2, m]])]
  }
  best <- which.max(total)
  
  list(
    perms = lapply(seq_len(num.runs), function(r) perm.mat[comb[best, r], ]),
    H = total[best] / ncol(pairs)
  )
}

# Greedy / LargeKGreedy: for each input order, fix the first run, then align
# each following run to all runs already fixed. Keep the order with the
# highest H.
.clumppGreedy <- function(q.list, run.orders, sim.stat, large.k = FALSE) {
  k <- ncol(q.list[[1]])
  num.runs <- length(q.list)
  perm.mat <- if(large.k) NULL else .clumppPerms(k)
  if(!large.k && k > 8) {
    warning("Greedy tests all ", factorial(k), " permutations per run; ",
            "consider align.algorithm = 'large.k'.")
  }
  
  best <- list(perms = NULL, H = -Inf)
  for(o in seq_len(nrow(run.orders))) {
    run.order <- run.orders[o, ]
    perms <- vector("list", num.runs)
    aligned <- vector("list", num.runs)
    perms[[run.order[1]]] <- seq_len(k)
    aligned[[run.order[1]]] <- q.list[[run.order[1]]]
    fixed <- run.order[1]
    
    for(x in run.order[-1]) {
      p <- if(large.k) {
        .clumppLargeKPerm(aligned[fixed], q.list[[x]])
      } else {
        s <- 0
        for(y in fixed) {
          s <- s + .clumppSimAll(aligned[[y]], q.list[[x]], perm.mat, sim.stat)
        }
        perm.mat[which.max(s), ]
      }
      perms[[x]] <- p
      aligned[[x]] <- q.list[[x]][, p, drop = FALSE]
      fixed <- c(fixed, x)
    }
    
    H <- .clumppH(aligned, sim.stat)
    if(H > best$H) best <- list(perms = perms, H = H)
  }
  best
}

# LargeKGreedy step: match columns of q to the fixed runs one pair at a time,
# taking the closest remaining (fixed column, q column) pair each time.
# Closeness is the squared distance between columns summed over fixed runs.
.clumppLargeKPerm <- function(fixed.list, q) {
  k <- ncol(q)
  cost <- Reduce(`+`, lapply(fixed.list, function(f) {
    outer(colSums(f^2), colSums(q^2), "+") - 2 * crossprod(f, q)
  }))
  p <- integer(k)
  for(i in seq_len(k)) {
    ab <- which(cost == min(cost), arr.ind = TRUE)[1, ]
    p[ab[1]] <- ab[2]
    cost[ab[1], ] <- Inf
    cost[, ab[2]] <- Inf
  }
  p
}