permuteStrata <- function(g) {
  g@strata <- sample(g@strata)
  g
}