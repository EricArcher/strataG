context("Population structure metrics")

test_that("statFst equals fstat", {
  data(msats.g)
  fstat_fst <- unname(hierfstat::fstat(gtypes2genind(msats.g))["Total", "pop"])
  strataG_fst <- unname(statFst(msats.g)$result["estimate"])
  expect_equal(round(fstat_fst, 4), round(strataG_fst, 4))
})

test_that("Hstat equals basic.stats", {
  data(msats.g)
  hierfstat_Hstats <- hierfstat::basic.stats(
    hierfstat::genind2hierfstat(gtypes2genind(msats.g))
  )
  hierfstat_Hstats <- as.matrix(hierfstat_Hstats$perloc)
  hierfstat_Hstats <- round(hierfstat_Hstats[, c("Ho", "Hs", "Ht")], 4)
  strataG_Hstats <- round(t(Hstats(msats.g)), 4)
  expect_equal(hierfstat_Hstats, strataG_Hstats)
})

test_that("Fis equals basic.stats", {
  data(msats.g)
  hierfstat_Hstats <- hierfstat::basic.stats(
    hierfstat::genind2hierfstat(gtypes2genind(msats.g))
  )
  hierfstat_Fis <- round(hierfstat_Hstats$perloc[, "Fis"], 4)
  strataG_Fis <- unname(sapply(
    getLociNames(msats.g), 
    function(loc) unname(statFis(msats.g[, loc, ])$result["estimate"])
  ))
  strataG_Fis <- round(strataG_Fis, 4)
  expect_equal(hierfstat_Fis, strataG_Fis)
})