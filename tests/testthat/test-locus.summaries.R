context("locus summaries")

test_that("allele frequencies are computed and formed correctly", {
  data(msats.g)
  af <- alleleFreqs(msats.g)
  af2 <- alleleFreqs(msats.g, by.strata = TRUE, type = "prop")
  
  expect_that(af, is_a("list"))
  expect_equal(names(af), getLociNames(msats.g))
  expect_equal(unique(sapply(af, class)), "table")
  expect_true(nrow(sapply(af2, dim)) == 2)
  expect_true(all(sapply(af2, dim)[2, ] == 3))
})

test_that("the number of alleles are computed and formed correctly", {
  data(msats.g)
  na <- numAlleles(msats.g)
  na2 <- numAlleles(msats.g, TRUE)
  
  expect_that(na, is_a("data.frame"))
  expect_equal(ncol(na), 2)
  expect_equal(ncol(na2), 3)
  expect_equal(nrow(na), getNumLoci(msats.g))
  expect_equal(sum(na$num.alleles), 68)
  expect_equal(nrow(na2), getNumLoci(msats.g) * getNumStrata(msats.g))
})

test_that("the number of individuals genotyped are computed and formed correctly", {  
  data(msats.g)
  ng <- numGenotyped(msats.g)
  ng2 <- numGenotyped(msats.g, TRUE)
  
  expect_that(ng, is_a("data.frame"))
  expect_equal(ncol(ng), 2)
  expect_equal(ncol(ng2), 3)
  expect_equal(nrow(ng), getNumLoci(msats.g))
  expect_true(all(ng$num.genotyped <= getNumInd(msats.g)))
  expect_equal(nrow(ng2), getNumLoci(msats.g) * getNumStrata(msats.g))
})

test_that("the number of individuals missing genotypes are computed and formed correctly", {  
  data(msats.g)
  nm <- numMissing(msats.g)
  nm2 <- numMissing(msats.g, TRUE)
  
  expect_that(nm, is_a("data.frame"))
  expect_equal(ncol(nm), 2)
  expect_equal(ncol(nm2), 3)
  expect_equal(nrow(nm), getNumLoci(msats.g))
  expect_true(all(nm$num.missing <= getNumInd(msats.g)))
  expect_equal(nrow(nm2), getNumLoci(msats.g) * getNumStrata(msats.g))
  expect_equal(sum(nm$num.missing), sum(is.na(msats.g@data$allele)) / 2)
  expect_equal(sum(nm$num.missing), sum(nm2$num.missing))
})

test_that("the number of individuals missing genotypes are computed and formed correctly", {  
  data(msats.g)
  nm <- numMissing(msats.g)
  nm2 <- numMissing(msats.g, TRUE)
  
  expect_that(nm, is_a("data.frame"))
  expect_equal(ncol(nm), 2)
  expect_equal(ncol(nm2), 3)
  expect_equal(nrow(nm), getNumLoci(msats.g))
  expect_true(all(nm$num.missing <= getNumInd(msats.g)))
  expect_equal(nrow(nm2), getNumLoci(msats.g) * getNumStrata(msats.g))
  expect_equal(sum(nm$num.missing), sum(is.na(msats.g@data$allele)) / 2)
  expect_equal(sum(nm$num.missing), sum(nm2$num.missing))
})

test_that("allelic richness is computed and formed correctly", { 
  data(msats.g)
  ar <- allelicRichness(msats.g)
  ar2 <- allelicRichness(msats.g, TRUE)
  
  expect_that(ar, is_a("data.frame"))
  expect_equal(ncol(ar), 2)
  expect_equal(ncol(ar2), 3)
  expect_equal(nrow(ar), getNumLoci(msats.g))
  expect_equal(nrow(ar2), getNumLoci(msats.g) * getNumStrata(msats.g))
})

test_that("heterozygosity is computed and formed correctly", { 
  data(msats.g)
  data(dloop.g)
  
  het.m.exp <- heterozygosity(msats.g)
  expect_true(min(het.m.exp$exptd.het) > 0.74 & max(het.m.exp$exptd.het) < 0.84)
  
  het.m.exp2 <- heterozygosity(msats.g, TRUE)
  expect_true(min(het.m.exp2$exptd.het) > 0.49 & max(het.m.exp2$exptd.het) < 0.95)

  het.m.obs <- heterozygosity(msats.g, type = "obs")
  expect_true(min(het.m.obs$obsvd.het) > 0.65 & max(het.m.obs$obsvd.het) < 0.77)
  
  het.m.obs2 <- heterozygosity(msats.g, TRUE, type = "o")
  expect_true(min(het.m.obs2$obsvd.het) > 0.51 & max(het.m.obs2$obsvd.het) < 0.95)
  
  het.d.exp <- heterozygosity(dloop.g)
  expect_true(het.d.exp$exptd.het > 0.9 & het.d.exp$exptd.het < 0.92)
  
  het.d.exp2 <- heterozygosity(dloop.g, TRUE)
  expect_true(min(het.d.exp2$exptd.het) > 0.74 & max(het.d.exp2$exptd.het) < 0.97)
})