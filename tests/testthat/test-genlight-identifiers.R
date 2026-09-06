context("genlight sample identifiers")

test_that("genlight identifiers remain aligned with genotypes and strata", {
  d <- matrix(c(0,1,2, 2,0,1), nrow = 3)
  gl <- adegenet::as.genlight(d)
  adegenet::indNames(gl) <- c("fishZ", "fishA", "fishM")
  adegenet::locNames(gl) <- c("loc1", "loc2")
  adegenet::pop(gl) <- c("west", "east", "west")
  g <- genlight2gtypes(gl)
  result <- as.data.frame(g, one.col = TRUE)
  result <- result[match(adegenet::indNames(gl), result$id), ]
  expect_setequal(getIndNames(g), adegenet::indNames(gl))
  expect_equal(result$id, adegenet::indNames(gl))
  expect_equal(result$stratum, c("west", "east", "west"))
  expect_equal(result$loc1, c("A/A", "A/G", "G/G"))
  expect_equal(result$loc2, c("G/G", "A/A", "A/G"))
})

test_that("unnamed genlight samples retain numeric fallback identifiers", {
  gl <- adegenet::as.genlight(matrix(c(0,1,2,2,0,1), nrow = 3))
  gl@ind.names <- NULL
  g <- genlight2gtypes(gl)
  expect_setequal(getIndNames(g), as.character(1:3))
  result <- as.data.frame(g, one.col = TRUE)
  expect_true(all(result$stratum == "Default"))
})
