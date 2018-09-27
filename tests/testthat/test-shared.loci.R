context("shared loci")

test_that("shared alleles are properly computed", {
  data(msats.g)
  
  sa.num <- sharedAlleles(msats.g, "num")
  sa.wch <- sharedAlleles(msats.g, "which")
  
  expect_equal(ncol(sa.num), getNumLoci(msats.g) + 2)
  expect_equal(unique(sapply(sa.wch, class)), "character")
  expect_equal(sa.num[1, "D11t"], 3)
  expect_equal(sa.num[3, "Ttr11"], 7)
  expect_equal(nchar(sa.wch[1, "D11t"]), 13)
  expect_equal(nchar(sa.wch[3, "Ttr34"]), 38)
})