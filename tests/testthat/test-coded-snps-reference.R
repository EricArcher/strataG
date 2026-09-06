context("coded SNP reference allele orientation")

reference_test_gtypes <- function() {
  # L2 keeps samples with missing L1 genotypes in the object.
  x <- data.frame(id = letters[1:6], pop = "P",
                  L1.1 = c("A", "A", "G", "A", NA, "A"),
                  L1.2 = c("A", "G", "G", "A", NA, NA),
                  L2.1 = rep("C", 6), L2.2 = rep(c("C", "T"), 3))
  df2gtypes(x, ploidy = 2)
}

test_that("coded SNPs count alternate rather than reference alleles", {
  g <- reference_test_gtypes()
  result <- as.data.frame(g, coded.snps = TRUE,
                          ref.allele = c(L1 = "A", L2 = "C"))
  expect_equal(result$id, letters[1:6])
  expect_equal(result$L1, c(0, 1, 2, 0, NA, NA))
  expect_equal(result$L2, rep(c(0, 1), 3))
})

test_that("changing the reference reverses nonmissing dosages", {
  g <- reference_test_gtypes()
  a <- as.data.frame(g, coded.snps = TRUE,
                     ref.allele = c(L1 = "A", L2 = "C"))
  b <- as.data.frame(g, coded.snps = TRUE,
                     ref.allele = c(L2 = "T", L1 = "G"))
  expect_equal(b$L1, 2 - a$L1)
  expect_equal(b$L2, 2 - a$L2)
})

test_that("default major references and unnamed references agree", {
  g <- reference_test_gtypes()
  expected <- as.data.frame(g, coded.snps = TRUE,
                            ref.allele = c(L1 = "A", L2 = "C"))
  expect_equal(as.data.frame(g, coded.snps = TRUE), expected)
  expect_equal(as.data.frame(g, coded.snps = TRUE,
                             ref.allele = c("A", "C")), expected)
  expect_equal(as.data.frame(g, coded.snps = TRUE, ids = FALSE,
                             strata = FALSE), expected[, c("L1", "L2")])
})
