context("gtypes creation")

test_that("a haploid data frame can be converted to gtypes", {
  data(dolph.strata)
  data(dolph.haps)
  g <- df2gtypes(
    dolph.strata[, c("id", "fine", "dLoop")], 
    ploidy = 1, 
    sequences = dolph.haps
  )
  
  expect_that(g, is_a("gtypes"))
  expect_equal(getPloidy(g), 1)
  expect_equal(length(getSequences(g)), 1)
})

test_that("a diploid data frame can be converted to gtypes", {
  data(dolph.msats)
  data(dolph.strata)
  g <- df2gtypes(
    dolph.msats, 
    ploidy = 2, 
    strata.col = NULL, 
    loc.col = 2, 
    schemes = dolph.strata[, c("id", "broad")]
  )
  g <- stratify(g, "broad")
  
  expect_that(g, is_a("gtypes"))
  expect_equal(getPloidy(g), 2)
})

test_that("sequences can be converted to gtypes", {
  data(dolph.seqs)
  data(dolph.strata)
  g <- sequence2gtypes(
    dolph.seqs,
    schemes = dolph.strata[, c("id", "broad", "fine")]
  )
  g <- stratify(g, "fine")
  
  expect_that(g, is_a("gtypes"))
  expect_equal(getPloidy(g), 1)
  expect_equal(length(getSequences(g)), 1)
})