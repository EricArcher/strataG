context("locus summaries")

test_that("allele frequencies are computed and formed correctly", {
  data(msats.g)
  af <- alleleFreqs(msats.g)
  af2 <- alleleFreqs(msats.g, by.strata = TRUE, type = "prop")
  
  expect_that(af, is_a("list"))
  expect_equal(names(af), getLociNames(msats.g))
  expect_equal(unique(sapply(af, class)), "table")
  expect_true(nrow(sapply(af2, dim)) == 2)
  expect_true(all(sapply(af2, dim)[2, ] == 3))
})