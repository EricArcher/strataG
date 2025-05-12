context("gtypes conversions")

test_that("gtypes can be converted to data.frames and matrices", {
  data(msats.g)
  
  df <- as.data.frame(msats.g)
  df2 <- as.data.frame(msats.g, one.col = TRUE)
  df3 <- as.data.frame(msats.g, one.col = TRUE, ids = FALSE)
  df4 <- as.data.frame(msats.g, strata = FALSE)
  mat <- as.matrix(msats.g)
  
  expect_that(df, is_a("data.frame"))
  expect_equal(nrow(df), getNumInd(msats.g))
  expect_equal(ncol(df2), getNumLoci(msats.g) + 2)
  expect_equal(colnames(df3)[1], "stratum")
  expect_equal(ncol(df4), (getNumLoci(msats.g) * getPloidy(msats.g)) + 1)
  expect_that(mat, is_a("matrix"))
})