context("sequence summaries")

test_that("IUPAC codes are computed correctly", {  
  i1 <- iupacCode(c("a", "a", "g"))
  i2 <- iupacCode(c("t", "c", "g"))
  
  x1 <- validIupacCodes(c("c", "t", "c", "c"))
  x2 <- validIupacCodes(c("c", "y", "c", "c"))
  x3 <- validIupacCodes(c("a", "g", "t", "a"))
  
  expect_equal(i1, "r")
  expect_equal(i2, "b")
  
  expect_equal(x1[1], "y")
  expect_equal(x2[1], "y")
  expect_equal(x3[1], "d")
})

test_that("base frequencies are computed correctly", {
  data(dloop.g)
  
  bf <- baseFreqs(dloop.g)
  
  expect_true(bf$site.freqs["a", 1] == 0)
  expect_true(bf$site.freqs["g", 1] == 33)
  expect_true(bf$base.freqs["a"] == 3977)
  expect_true(bf$base.freqs["-"] == 61)
  expect_true(bf$ind.freqs["Hap.01", "a"] == 120)
  expect_true(bf$ind.freqs["Hap.33", "g"] == 52)
})

test_that("fixed sites are computed correctly", {
  data(dloop.g)
  
  fs <- fixedSites(dloop.g)
  
  expect_true(fs["30"] == "c")
  expect_true(!"20" %in% names(fs))
  expect_true(fs["389"] == "t")
})

test_that("variable sites are computed correctly", {
  data(dloop.g)
  
  vs <- variableSites(dloop.g)
  
  expect_equal(ncol(vs$site.freqs), 43)
  expect_true("20" %in% colnames(vs$site.freqs))
  expect_equal(vs$site.freqs["a", "20"], 2)
  expect_equal(vs$site.freqs["g", "20"], 31)
})