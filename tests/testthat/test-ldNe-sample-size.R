context("ldNe sample size")

make_ldne_sample <- function(n) {
  # Deterministic, polymorphic diploid loci with non-identical dosages.
  x <- data.frame(id = paste0("s", seq_len(n)), pop = "P")
  patterns <- list(c(0, 1, 2), c(0, 2, 1, 1, 0), c(2, 0, 1, 0, 2, 1, 1))
  for (i in seq_along(patterns)) {
    dosage <- rep(patterns[[i]], length.out = n)
    x[[paste0("L", i, ".1")]] <- 1L + as.integer(dosage > 0)
    x[[paste0("L", i, ".2")]] <- 1L + as.integer(dosage == 2)
  }
  x
}

test_that("complete-data sample size does not depend on locus-pair count", {
  for (n in c(20L, 40L)) {
    x <- make_ldne_sample(n)
    for (k in 2:3) {
      g <- df2gtypes(x[, seq_len(2 + 2 * k)], ploidy = 2)
      result <- ldNe(g, num.cores = 1)
      expect_equal(result$S, as.numeric(n))
      expect_equal(result$num.comp, choose(k, 2))
    }
  }
})

test_that("dropping an extra incomplete locus preserves the result", {
  for (n in c(20L, 40L)) {
    x <- make_ldne_sample(n)
    complete <- ldNe(df2gtypes(x, ploidy = 2), num.cores = 1)
    x$L4.1 <- x$L1.1
    x$L4.2 <- x$L1.2
    x[1, c("L4.1", "L4.2")] <- NA
    dropped <- ldNe(df2gtypes(x, ploidy = 2),
                    drop.missing = TRUE, num.cores = 1)
    expect_equal(dropped$S, as.numeric(n))
    expect_equal(complete, dropped, tolerance = 1e-10)
  }
})
