context("fastsimcoal single-row population blocks")

# These fixtures exercise the ARP reader without running fastsimcoal.
parse_test_arp <- function(blocks, preamble = character()) {
  path <- tempfile(fileext = ".arp")
  on.exit(unlink(path))
  lines <- unlist(lapply(seq_along(blocks), function(i) {
    c(paste0('SampleName="deme', i, '"'),
      paste0("SampleSize=", length(blocks[[i]])),
      "SampleData= {", blocks[[i]], "}")
  }))
  writeLines(c("[Profile]", 'Title="Synthetic test"', "DataType=DNA",
               "GenotypicData=0", "[Data]", "[[Samples]]", preamble, lines), path)
  .fscParseArpFile(path)
}

test_that("a single sequence row remains a matrix", {
  result <- parse_test_arp(list("sample1 1 ACTG"))
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(1L, 3L))
  expect_equal(colnames(result), c("id", "deme", "col3"))
  expect_equal(unname(result[1, ]), c("sample1", "1", "ACTG"))
})

test_that("multiple singleton populations are retained", {
  result <- parse_test_arp(list("sample1 1 ACTG", "sample2 1 ACCG"))
  expect_equal(dim(result), c(2L, 3L))
  expect_equal(unname(result[, "deme"]), c("1", "2"))
  expect_equal(unname(result[, "col3"]), c("ACTG", "ACCG"))
})

test_that("mixed population sizes retain every row and deme", {
  result <- parse_test_arp(list(c("sample1 1 ACTG", "sample2 1 ACCG"),
                                "sample3 1 ATCG"))
  expect_equal(dim(result), c(3L, 3L))
  expect_equal(unname(result[, "id"]), paste0("sample", 1:3))
  expect_equal(unname(result[, "deme"]), c("1", "1", "2"))
})

test_that("multiple marker columns survive singleton extraction", {
  result <- parse_test_arp(list("sample1 1 0 1 0"))
  expect_equal(dim(result), c(1L, 5L))
  expect_equal(colnames(result), c("id", "deme", "col3", "col4", "col5"))
  expect_equal(unname(result[1, ]), c("sample1", "1", "0", "1", "0"))
})

test_that("ordinary multi-row input and position metadata are preserved", {
  result <- parse_test_arp(list(c("sample1 1 ACTG", "sample2 1 ACCG")),
                           c("# 2 polymorphic positions on chromosome 1", "1 3"))
  expect_equal(dim(result), c(2L, 3L))
  expect_equal(unname(result[, "col3"]), c("ACTG", "ACCG"))
  expect_equal(attr(result, "poly.pos"),
               cbind(chromosome = c(1L, 1L), position = c(1, 3)))
})

test_that("empty blocks are skipped without renumbering populations", {
  result <- parse_test_arp(list(character(), "sample1 1 ACTG"))
  expect_equal(unname(result[1, ]), c("sample1", "2", "ACTG"))
  expect_null(parse_test_arp(list(character())))
})
