context("population structure statistics")

test_that("haploid values are correct", {
  load("validation_results.rdata")
  data(dloop.g)
  dloop.ovl <- overallTest(dloop.g, 0, by.locus = TRUE)
  expect_equal(dloop.ovl$result["CHIsq", "estimate"], chisq.dloop.true)
  expect_equal(dloop.ovl$result["PHIst", "estimate"], phist.true)
})

test_that("diploid values are correct", {
  load("validation_results.rdata")
  data(msats.g)
  msats.ovl <- overallTest(msats.g, 0, by.locus = TRUE)
  x <- msats.ovl$result[, , "estimate"]
  expect_equal(x["All", "CHIsq"], chisq.msats.true)
  expect_lt(abs(x["All", "Ho"] - hstats.true["overall", "Ho"]), 1e-3)
  expect_lt(abs(x["All", "Hs"] - hstats.true["overall", "Hs"]), 1e-3)
  expect_lt(abs(x["All", "Ht"] - hstats.true["overall", "Ht"]), 1e-3)
  expect_lt(abs(x["All", "Dest"] - x["All", "Dest_Chao"]), 0.05)
  expect_lt(abs(x["D11t", "wcFit"] - wcFst.true["D11t", "Fit"]), 1e-2)
  expect_lt(abs(x["EV37", "wcFst"] - wcFst.true["EV37", "Fst"]), 1e-2)
  expect_lt(abs(x["Ttr34", "wcFis"] - wcFst.true["Ttr34", "Fis"]), 1e-2)
})