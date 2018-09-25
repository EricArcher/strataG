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

test_that("the number of alleles are computed and formed correctly", {
  data(msats.g)
  na <- numAlleles(msats.g)
  na2 <- numAlleles(msats.g, TRUE)
  
  expect_that(na, is_a("data.frame"))
  expect_equal(ncol(na), 2)
  expect_equal(ncol(na2), 3)
  expect_equal(nrow(na), getNumLoci(msats.g))
  expect_equal(sum(na$num.alleles), 68)
  expect_equal(nrow(na2), getNumLoci(msats.g) * getNumStrata(msats.g))
})

test_that("the number of individuals genotyped are computed and formed correctly", {  data(msats.g)
  ng <- numGenotyped(msats.g)
  ng2 <- numGenotyped(msats.g, TRUE)
  
  expect_that(ng, is_a("data.frame"))
  expect_equal(ncol(ng), 2)
  expect_equal(ncol(ng2), 3)
  expect_equal(nrow(ng), getNumLoci(msats.g))
  expect_true(all(ng$num.genotyped <= getNumInd(msats.g)))
  expect_equal(nrow(ng2), getNumLoci(msats.g) * getNumStrata(msats.g))
})