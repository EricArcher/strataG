context("sequence summaries")

test_that("IUPAC codes are computed correctly", {  
  i1 <- iupacCodes(c("a", "a", "g"))
  i2 <- iupacCodes(c("t", "c", "g"))
  
  x1 <- iupacCodes(c("c", "t", "c", "c"), all = TRUE)
  x2 <- iupacCodes(c("c", "y", "c", "c"), all = TRUE)
  x3 <- iupacCodes(c("a", "g", "t", "a"), all = TRUE)
  
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
  expect_true(bf$site.freqs["g", 1] == 126)
  expect_true(bf$base.freqs["t"] == 17166)
  expect_true(bf$base.freqs["-"] == 245)
  expect_true(bf$ind.freqs["23792", "a"] == 120)
  expect_true(bf$ind.freqs["4498", "g"] == 51)
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
  expect_equal(vs$site.freqs["g", "20"], 124)
})