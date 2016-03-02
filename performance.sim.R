pop.info <- fscPopInfo(pop.size = rep(10000, 3), sample.size = rep(20, 3))
snp.params <- fscLocusParams(locus.type = "snp", num.loci = 1000, mut.rate = 1e-5)
msat.params <-   fscLocusParams(locus.type = "msat", num.loci = 50, mut.rate = 1e-5)
mig.rates <- matrix(c(0, 0.01, 0.001, 0.01, 0, 0.000001, 0.001, 0.000001, 0), ncol = 3)
hist.ev <- fscHistEv(
  num.gen = c(2000, 2001), source.deme = c(2, 1),
  sink.deme = c(1, 0), prop.migrants = 1
)

snp.1K <- fastsimcoal(pop.info, snp.params, mig.rates, hist.ev)
snp.pop.struct <- pairwiseTest(snp.1K, stats = c("Chi2", "Fst", "Gst"), 
                               nrep = 100, num.cores = 4)

microbenchmark(
  snp.pop.struct <- pairwiseTest(snp.1K, stats = c("Chi2", "Fst", "Gst"), 
                                 nrep = 1000, num.cores = 2), 
  times = 5
)


msats.50 <- fastsimcoal(pop.info, msat.params, mig.rates, hist.ev)
msats.pop.struct <- pairwiseTest(msats.50, stats = c("Chi2", "Fst", "Gst"), 
                                 nrep = 100, num.cores = 4)

microbenchmark(
  msats.pop.struct <- pairwiseTest(msats.50, stats = c("Chi2", "Fst", "Gst"), 
                                   nrep = 1000, num.cores = 2),
  times = 20
)
#15.8s
