context("HWE locus alignment")

hwe_fixture <- function() {
  df2gtypes(data.frame(id = paste0("s", 1:20), pop = rep(c("A", "B"), 10),
    L1.1 = rep(c(1,1,2,2),5), L1.2 = rep(c(1,2,1,2),5),
    L2.1 = rep(c(1,2,1,1,2),4), L2.2 = rep(c(2,2,1,1,2),4)), ploidy=2)
}

test_that("single locus retains its name", {
  g <- hwe_fixture()[, "L1", ]
  expect_named(hweTest(g, num.rep=20), "L1")
})

test_that("missing loci retain their positions", {
  g <- hwe_fixture()
  g@data$allele[g@data$locus == "L1"] <- NA
  result <- hweTest(g, num.rep=20)
  expect_named(result, getLociNames(g))
  expect_length(result, 2L)
  expect_true(is.na(result["L1"]))
  expect_true(is.finite(result["L2"]))
  g@data$allele[] <- NA
  expect_identical(hweTest(g, num.rep=20), stats::setNames(c(NA_real_,NA_real_), getLociNames(g)))
})

test_that("ordinary results preserve the pegas calculation", {
  g <- hwe_fixture()
  set.seed(12)
  expected <- pegas::hw.test(gtypes2genind(g), B=20)
  set.seed(12)
  actual <- hweTest(g, num.rep=20)
  expect_equal(actual, expected[,ncol(expected)])
})
