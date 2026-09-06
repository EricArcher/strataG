context("GENEPOP file consistency")

genepop_fixture <- function(pop = c("W", "W", "E"), ids = c("a", "b", "c")) {
  df2gtypes(data.frame(id = ids, pop = pop,
                      L1.1 = c(101, 102, 101), L1.2 = c(101, 102, 102)),
             ploidy = 2)
}

test_that("wrapper settings and return values identify the written file", {
  folder <- tempfile("genepop-test-")
  dir.create(folder)
  old <- setwd(folder)
  on.exit({setwd(old); unlink(folder, recursive = TRUE)})
  # Replace only the external command, not the writer or settings generation.
  wrapper <- genepop
  mock <- new.env(parent = environment(wrapper))
  mock$system <- function(...) 0L
  environment(wrapper) <- mock
  for (filename in c("loc_data.txt", "custom.txt")) {
    result <- wrapper(genepop_fixture(), label = "audit", input.fname = filename,
                       output.ext = ".OUT")
    expect_true(file.exists(filename))
    expect_equal(readLines("settings.txt")[1], paste0("InputFile=", filename))
    expect_equal(unname(result$files["input.fname"]), filename)
    expect_equal(unname(result$files["output.fname"]), paste0(filename, ".OUT"))
    expect_equal(result$locus.names, c(LOC1 = "L1"))
  }
  direct <- genepopWrite(genepop_fixture(), label = "direct")
  expect_equal(direct$fname, "direct_loc_data.txt")
  expect_true(file.exists(direct$fname))
})

test_that("ambiguous population or sample labels fail before writing", {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  writeLines("existing content", path)
  expect_error(genepopWrite(genepop_fixture(c("A,B", "A_B", "A_B")),
                            input.fname = path), "Population labels are not unique")
  expect_equal(readLines(path), "existing content")
  expect_error(genepopWrite(genepop_fixture(c("W", "W", "E"),
                                           c("a,b", "a_b", "c")),
                            input.fname = path), "Sample labels are not unique")
  expect_equal(readLines(path), "existing content")
})
