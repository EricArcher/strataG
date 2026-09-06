context("gtypes to genlight alignment")

test_that("genotypes and metadata stay aligned regardless of internal row order", {
  x <- data.frame(id = c("z", "a", "m", "b"),
                  pop = c("W", "E", "W", "E"),
                  zloc.1 = c("A", "A", "G", "A"),
                  zloc.2 = c("A", "G", "G", "A"),
                  aloc.1 = c("G", "G", "A", "A"),
                  aloc.2 = c("G", "G", "G", "A"))
  g <- df2gtypes(x, ploidy = 2)
  for (reverse in c(FALSE, TRUE)) {
    current <- g
    if (reverse) {
      current@data <- current@data[rev(seq_len(nrow(current@data))), ]
    }
    expect_true(methods::validObject(current))
    expected <- as.data.frame(current, coded.snps = TRUE)
    result <- gtypes2genlight(current)
    expect_equal(adegenet::indNames(result), expected$id)
    expect_equal(as.character(adegenet::pop(result)), expected$stratum)
    loci <- setdiff(names(expected), c("id", "stratum"))
    expect_equal(adegenet::locNames(result), loci)
    dosages <- as.matrix(expected[, loci, drop = FALSE])
    rownames(dosages) <- expected$id
    expect_equal(as.matrix(result), dosages)
  }
})
