context("ldNe confidence level validation")

test_that("invalid confidence levels fail before genotype processing", {
  invalid <- list(NULL, numeric(), c(0.9, 0.95), NA_real_, NaN,
                  Inf, -Inf, -0.5, 0, 1, 1.2, "0.95", TRUE, 0.95 + 0i)
  for (ci in invalid) {
    expect_error(ldNe(NULL, ci = ci),
                 "'ci' must be a single finite number strictly between 0 and 1.",
                 fixed = TRUE)
  }
})

test_that("valid confidence levels give ordered bounds", {
  set.seed(126)
  d <- matrix(stats::rbinom(80 * 5, 2, 0.3), 80, 5)
  x <- data.frame(id = paste0("s", 1:80), pop = "P")
  for (i in 1:5) {
    x[[paste0("L", i, ".1")]] <- 1L + (d[, i] > 0)
    x[[paste0("L", i, ".2")]] <- 1L + (d[, i] == 2)
  }
  g <- df2gtypes(x, ploidy = 2)
  default <- ldNe(g)
  expect_equal(default, ldNe(g, ci = 0.95))
  for (ci in c(0.8, 0.9, 0.99)) {
    result <- ldNe(g, ci = ci)
    expect_equal(result$Ne, default$Ne)
    expect_true(result$param.lci <= result$param.uci)
  }
})
