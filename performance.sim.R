snp.1K <- fastsimcoal(
  fscPopInfo(
    pop.size = rep(10000, 3),
    sample.size = rep(20, 3)
  ),
  fscLocusParams(
    locus.type = "snp",
    num.loci = 1000,
    mut.rate = 1e-5
  ),
  mig.rates = matrix(c(0, 0.1, 0.01, 0.1, 0, 0.001, 0.01, 0.001, 0), ncol = 3)
)

msats.25 <- fastsimcoal(
  fscPopInfo(
    pop.size = rep(10000, 3),
    sample.size = rep(100, 3)
  ),
  fscLocusParams(
    locus.type = "msat",
    num.loci = 25,
    mut.rate = 1e-5
  ),
  mig.rates = matrix(c(0, 0.4, 0.001, 0.001, 0, 0.001, 0.001, 0.4, 0), ncol = 3)
)


snp.pop.struct <- pairwiseTest(snp.1K, stats = c("Chi2", "Fst", "Gst"), 
                               nrep = 1000, num.cores = 4)