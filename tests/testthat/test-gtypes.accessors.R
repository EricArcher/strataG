context("gtypes accessors")

test_that("numerical accessors return right numbers", {
  data(dolph.strata)
  data(dolph.haps)
  data(dolph.msats)
  
  cr.g <- df2gtypes(
    dolph.strata[, c("id", "fine", "dLoop")], 
    ploidy = 1, 
    sequences = dolph.haps
  )

  ms.g <- df2gtypes(
    dolph.msats, 
    ploidy = 2, 
    strata.col = NULL, 
    loc.col = 2, 
    schemes = dolph.strata[, c("id", "broad")]
  )
  ms.g <- stratify(ms.g, "broad")
  
  expect_equal(getNumInd(cr.g), 126)
  expect_equal(getNumInd(ms.g), 126)
  expect_equal(getNumLoci(cr.g), 1)
  expect_equal(getNumLoci(ms.g), 5)
  expect_equal(getNumStrata(cr.g), 3)
  expect_equal(getNumStrata(ms.g), 2)
})

test_that("gtypes individual, strata, and locus names are correctly returned", {
  data(dolph.strata)
  data(dolph.msats)
  
  ms.g <- df2gtypes(
    dolph.msats, 
    ploidy = 2, 
    strata.col = NULL, 
    loc.col = 2, 
    schemes = dolph.strata[, c("id", "broad")]
  )
  ms.g <- stratify(ms.g, "broad")
  
  expect_true(all(c("18650", "41759", "78061") %in% getIndNames(ms.g)))
  expect_equal(length(getIndNames(ms.g)), getNumInd(ms.g))
  expect_true("Coastal" %in% getStrataNames(ms.g))
  expect_equal(length(getStrataNames(ms.g)), 2)
  expect_true(all(c("D11t", "EV94", "Ttr34") %in% getLociNames(msats.g)))
  expect_equal(length(getLociNames(msats.g)), getNumLoci(msats.g))
})

test_that("strata slot components are correctly returned and assigned", {
  data(msats.g)
  g <- msats.g
  setStrata(g) <- setNames(
    sample(letters[1:2], getNumInd(g), rep = TRUE), 
    getIndNames(g)
  )
  
  expect_equal(length(getStrata(g)), getNumInd(g))
  expect_true(all(names(getStrata(g)) %in% getIndNames(g)))
  expect_true(all(getStrata(g) %in% getStrataNames(g)))
  expect_error(setStrata(g) <- c("a", "b"))
  expect_error(setStrata(g) <- sample(letters, getNumInd(g), rep = TRUE))
})

test_that("indexing a gtypes object works for individuals, loci, and strata", {
  data(msats.g)
  data(dloop.g)
  
  expect_equal(getNumInd(msats.g[c("18650", "78055"), , ]), 2)
  expect_equal(getNumLoci(msats.g[, c("EV94", "D11t"), ]), 2)
  expect_equal(getStrataNames(msats.g[, , "Coastal"]), "Coastal")
})
