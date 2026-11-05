#' @title Theta from Expected Heterozygosity
#' @description Estimate the population-scaled mutation rate (theta) for each
#'   locus in a \code{gtypes} object by inverting the expected
#'   heterozygosity under the infinite alleles model (IAM) or stepwise
#'   mutation model (SMM), with or without Zouros's (1979) correction for the
#'   drift variance of homozygosity.
#'
#' @param g a gtypes object.
#' @param model one or more of \code{"IAM"}, \code{"SMM"}, \code{"ZourosIAM"},
#'   or \code{"ZourosSMM"}.
#' @param by.strata logical - calculate separately for each stratum?
#'
#' @return a data.frame with one row per locus (and stratum if \code{by.strata}),
#'   number of samples genotyped (\code{n}), expected heterozygosity (\code{H}),
#'   and for each model a \code{theta.<model>} column,
#'   and a \code{Ne.mu.<model>} column (theta / (2 * ploidy)).
#'
#' @details Expected heterozygosity (gene diversity) is calculated with Nei's (1987) 
#'   unbiased estimator. Homozygosity is taken as J = 1 - H. Theta is then:
#'   \describe{
#'     \item{IAM}{H / J, from E[J] = 1 / (1 + theta) (Kimura & Crow 1964).}
#'     \item{SMM}{0.5 * (1 / J^2 - 1), from
#'       E[J] = 1 / sqrt(1 + 2 * theta) (Ohta & Kimura 1973).}
#'     \item{ZourosIAM}{solves a = theta * (1 + 2 * (1 + theta) /
#'       ((2 + theta) * (3 + theta))) for theta, where a = 1 / J - 1
#'       (Zouros 1979, eqs. 10-12).}
#'     \item{ZourosSMM}{solves a = 2 * theta * (1 + 3 * (1 + 2 * theta)^2 /
#'       (2 * (3 + 11.25 * theta + 13 * theta^2 + 1.7 * theta^3))) for theta,
#'       where a = 1 / J^2 - 1 (Zouros 1979, eqs. 8, 13-14; Moran 1976).}
#'   }
#'
#' @references
#'   Kimura, M. and Crow, J.F. 1964. The number of alleles that can be
#'     maintained in a finite population. Genetics 49:725-738.
#'   Moran, P.A.P. 1976. Wandering distributions and the electrophoretic
#'     profile. II. Theoretical Population Biology 10:145-149.
#'   Nei, M. 1987. Molecular Evolutionary Genetics. Columbia University Press.
#'   Ohta, T. and Kimura, M. 1973. A model of mutation appropriate to estimate
#'     the number of electrophoretically detectable alleles in a finite
#'     population. Genetical Research 22:201-204.
#'   Zouros, E. 1979. Mutation rates, population sizes and amounts of
#'     electrophoretic variation at enzyme loci in natural populations.
#'     Genetics 92:623-646.
#'
#' @examples
#' data(msats.g)
#' thetaH(msats.g)
#' thetaH(msats.g, model = c("SMM", "ZourosSMM"), by.strata = TRUE)
#'
#' @export
thetaH <- function(g, model = c("IAM", "SMM", "ZourosIAM", "ZourosSMM"),
                   by.strata = FALSE) {
  model <- match.arg(model, several.ok = TRUE)

  ploidy <- getPloidy(g)
  if (length(ploidy) != 1 || is.na(ploidy) || ploidy < 1) {
    stop("could not read a valid ploidy from 'g'")
  }
  if(ploidy == 1) g <- labelHaplotypes(g)

  # Nei's unbiased expected heterozygosity
  .H <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2) return(c(n = n, H = NA_real_))
    p <- table(x) / n
    c(n = n, H = n * (1 - sum(p^2)) / (n - 1))
  }

  # Zouros (1979) approximations of E[a] as a function of theta (M)
  .EaIAM <- function(M) M * (1 + 2 * (1 + M) / ((2 + M) * (3 + M)))
  .EaSMM <- function(M) {
    2 * M * (1 + 3 * (1 + 2 * M)^2 / (2 * (3 + 11.25 * M + 13 * M^2 + 1.7 * M^3)))
  }

  # solve Ea(M) = a for M. Ea(M) >= scale * M, so the root lies in
  # [0, a / scale], which always brackets it.
  .invert <- function(a, Ea, scale) {
    vapply(a, function(a.i) {
      if (is.na(a.i)) return(NA_real_)
      if (a.i == 0) return(0)
      if (is.infinite(a.i)) return(Inf)
      stats::uniroot(
        function(M) Ea(M) - a.i,
        interval = c(0, a.i / scale),
        tol = .Machine$double.eps^0.75
      )$root
    }, numeric(1))
  }

  .theta <- list(
    IAM = function(H) H / (1 - H),
    SMM = function(H) 0.5 * (1 / (1 - H)^2 - 1),
    ZourosIAM = function(H) .invert(1 / (1 - H) - 1, .EaIAM, 1),
    ZourosSMM = function(H) .invert(1 / (1 - H)^2 - 1, .EaSMM, 2)
  )

  grp <- if(by.strata) c("stratum", "locus") else "locus"

  g@data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarize(
      res = list(.H(.data$allele)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ploidy = ploidy,
      n = purrr::map_int(.data$res, function(r) as.integer(r[["n"]])),
      H = purrr::map_dbl(.data$res, function(r) r[["H"]])
    ) |>
    dplyr::select(-"res") |>
    dplyr::mutate(
      dplyr::across(
        "H",
        .fns = .theta[model],
        .names = "theta.{.fn}"
      )
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("theta."),
        .fns = function(th) th / (2 * ploidy),
        .names = "{sub('theta.', 'Ne.mu.', .col, fixed = TRUE)}"
      )
    ) |> 
    dplyr::select(-"ploidy") |> 
    as.data.frame()
}
