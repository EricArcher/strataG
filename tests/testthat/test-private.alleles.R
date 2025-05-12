context("private alleles")

test_that("private alleles are computed correctly", {
  data(msats.g)
  pa <- privateAlleles(msats.g)
  
  expect_equal(nrow(pa), getNumLoci(msats.g))
  expect_equal(ncol(pa), getNumStrata(msats.g))
  expect_true(all(pa[, "Coastal"] == 0))
  expect_true(all(pa[, "Offshore.South"] == 1))
})