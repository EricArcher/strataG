context("ldNe minor allele frequency filtering")

ldne_maf_gtypes <- function(dosages, strata = rep("P", nrow(dosages))) {
  x <- data.frame(id = paste0("s", seq_len(nrow(dosages))), pop = strata)
  for (i in seq_len(ncol(dosages))) {
    x[[paste0("L", i, ".1")]] <- 1L + (dosages[, i] > 0)
    x[[paste0("L", i, ".2")]] <- 1L + (dosages[, i] == 2)
  }
  df2gtypes(x, ploidy = 2)
}

test_that("rare alleles are removed irrespective of dosage orientation", {
  d <- cbind(c(1, rep(0, 39)), rep(0:2, length.out = 40),
             rep(c(0, 2, 1, 0, 1), 8))
  for (x in list(d, 2 - d)) {
    actual <- ldNe(ldne_maf_gtypes(x), maf.threshold = 0.05)
    expected <- ldNe(ldne_maf_gtypes(x[, 2:3]), maf.threshold = 0)
    expect_equal(actual$num.comp, 1)
    expect_equal(actual, expected)
  }
})

test_that("shared filtering works and computes MAF separately in each stratum", {
  # L1 is rare in A but common in B. L2 reverses its common allele
  # between strata and has MAF 0.10 in both. L3 is common in both.
  a <- cbind(c(1, rep(0, 19)), c(rep(1, 4), rep(2, 16)),
             rep(c(0, 1, 2, 1), 5))
  b <- cbind(rep(c(0, 1, 2, 1), 5), c(rep(1, 4), rep(0, 16)),
             rep(c(2, 1, 0, 1), 5))
  d <- rbind(a, b)
  st <- rep(c("A", "B"), each = 20)
  g <- ldne_maf_gtypes(d, st)
  shared <- ldNe(g, maf.threshold = 0.05, by.strata = TRUE)
  expected <- ldNe(ldne_maf_gtypes(d[, 2:3], st), maf.threshold = 0)
  expect_equal(shared, expected)
  expect_equal(shared$num.comp, c(1, 1))
  separate <- ldNe(g, maf.threshold = 0.05, by.strata = FALSE)
  expect_equal(separate$num.comp, c(1, 3))
})

test_that("shared MAF filtering retains matrix dimensions for one stratum", {
  d <- cbind(rep(0:2, length.out = 40), rep(c(0, 2, 1, 0, 1), 8))
  g <- ldne_maf_gtypes(d)
  expect_equal(ldNe(g, maf.threshold = 0.05, by.strata = TRUE),
               ldNe(g, maf.threshold = 0.05, by.strata = FALSE))
})
