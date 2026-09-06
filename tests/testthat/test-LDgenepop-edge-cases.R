context("GENEPOP LD edge cases")

ld_genepop_fixture <- function() {
  df2gtypes(data.frame(id = paste0("s", 1:20), pop = rep(c("A", "B"), 10),
                      L1.1 = 1, L1.2 = 1,
                      L2.1 = rep(c(1,1,2,2),5),
                      L2.2 = rep(c(1,2,1,2),5),
                      L3.1 = rep(c(1,2,1,1,2),4),
                      L3.2 = rep(c(2,2,1,1,2),4)), ploidy = 2)
}

test_that("one locus is rejected before any files are written", {
  folder <- tempfile("ld-test-")
  dir.create(folder)
  old <- setwd(folder)
  on.exit({setwd(old); unlink(folder, recursive = TRUE)})
  g <- ld_genepop_fixture()[, "L1", ]
  expect_error(LDgenepop(g), "At least two loci are required", fixed = TRUE)
  expect_length(list.files(), 0L)
})

test_that("real GENEPOP untestable pairs are parsed without warnings", {
  skip_if_not_installed("genepop")
  folder <- tempfile("ld-test-")
  dir.create(folder)
  old <- setwd(folder)
  on.exit({setwd(old); unlink(folder, recursive = TRUE)})
  warnings <- character()
  result <- withCallingHandlers(
    LDgenepop(ld_genepop_fixture(), dememorization = 100,
              batches = 10, iterations = 400, delete.files = FALSE, label = "test"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(warnings, 0L)
  expect_equal(nrow(result), 3L)
  expect_equal(result$Locus.1, c("L1", "L1", "L2"))
  expect_equal(result$Locus.2, c("L2", "L3", "L3"))
  expect_true(all(is.na(as.matrix(result[1:2, 3:5]))))
  expect_true(is.finite(result$p.value[3]))
  lines <- readLines("test_loc_data.txt.DIS")
  row <- grep("LOC2[[:space:]]+LOC3", lines, value = TRUE)
  fields <- strsplit(trimws(row), "[[:space:]]+")[[1]]
  values <- as.numeric(tail(fields, 3))
  expect_equal(unname(unlist(result[3, 3:5])), values)
})
