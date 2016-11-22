context("Population structure metrics")

test_that("statFst equals fstat", {
  library(hierfstat)
  data(msats.g)
  fstat_fst <- fstat(gtypes2genind(msats.g))["Total", "pop"]
  strataG_fst <- statFst(msats.g)$result["estimate"]
  fst_diff <- abs(fstat_fst - strataG_fst)
  expect_true(fst_diff < 1e-10)
})