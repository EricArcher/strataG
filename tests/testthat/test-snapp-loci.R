context("SNAPP locus preservation")

test_that("SNAPP preserves every locus and writes matching sequential rows", {
  for (nloc in c(1L, 2L, 3L, 83L)) {
    x <- data.frame(id = c("a", "b", "c"), pop = "P")
    for (i in seq_len(nloc)) {
      x[[paste0("L", i, ".1")]] <- c("A", "A", "G")
      x[[paste0("L", i, ".2")]] <- c("A", "G", "G")
    }
    # Preserve missing data when another locus keeps the sample present.
    if (nloc > 1L) x[2, c("L1.1", "L1.2")] <- NA
    if (nloc > 2L) x[2, c("L3.1", "L3.2")] <- NA
    g <- df2gtypes(x, ploidy = 2)
    expected <- as.data.frame(g, coded.snps = TRUE)
    labels <- paste(expected$stratum, expected$id, sep = "_")
    expected <- as.matrix(expected[, -(1:2), drop = FALSE])
    rownames(expected) <- labels
    path <- tempfile(fileext = ".nex")
    result <- write.nexus.snapp(g, file = path)
    lines <- readLines(path)
    unlink(path)
    expect_equal(result, expected)
    expect_true(any(grepl(paste0("NCHAR=", nloc, ";"), lines, fixed = TRUE)))
    expect_true(any(grepl("INTERLEAVE=NO", lines, fixed = TRUE)))
    start <- grep("^[[:space:]]*MATRIX[[:space:]]*$", lines)
    end <- which(seq_along(lines) > start & grepl("^[[:space:]]*;", lines))[1]
    rows <- trimws(lines[seq.int(start + 1L, end - 1L)])
    rows <- rows[nzchar(rows)]
    fields <- strsplit(rows, "[[:space:]]+")
    expect_length(fields, nrow(expected))
    expect_equal(vapply(fields, `[`, "", 1), labels)
    sequences <- apply(expected, 1, function(z) {
      z[is.na(z)] <- "?"
      paste0(z, collapse = "")
    })
    expect_equal(vapply(fields, `[`, "", 2), unname(sequences))
  }
})
