## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  comment = NA
)
options(scipen = 99)

## ----eval=FALSE---------------------------------------------------------------
#  if(!require(devtools)) install.packages(devtools)
#  devtools::install_github("ericarcher/strataG", dependencies = TRUE)

## -----------------------------------------------------------------------------
rm(list = ls())
library(strataG)

## -----------------------------------------------------------------------------
deme0 <- fscDeme(deme.size = 1000, sample.size = 4)

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(deme0)

## -----------------------------------------------------------------------------
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 1e-3)

## -----------------------------------------------------------------------------
genetics <- fscSettingsGenetics(msats, num.chrom = 5)

## -----------------------------------------------------------------------------
ex1.params <- fscWrite(demes = demes, genetics = genetics, label = "ex1")

## ----echo=FALSE---------------------------------------------------------------
cat(readLines(file.path(ex1.params$folder, ex1.params$files$input)), sep = '\n')

## -----------------------------------------------------------------------------
print(str(ex1.params))

## -----------------------------------------------------------------------------
ex1.params <- fscRun(ex1.params, num.sim = 1)

## -----------------------------------------------------------------------------
dir(ex1.params$label)

## -----------------------------------------------------------------------------
str(ex1.params[c("run.params", "locus.info")])

## -----------------------------------------------------------------------------
arp.file <- fscReadArp(ex1.params)
str(arp.file)

# The first 6 columns
arp.file[, 1:6]

## -----------------------------------------------------------------------------
# create 3 independent chromosomes with the same structure of four markers
complex.chroms <- fscSettingsGenetics(
  fscBlock_microsat(2, 1e-4),
  fscBlock_dna(4, 1e-5),
  fscBlock_dna(6, 1e-3),
  fscBlock_microsat(2, 1e-5),
  num.chrom = 3
)
complex.params <- fscWrite(
  demes = demes, 
  genetics = complex.chroms,
  label = "complex_chroms"
)
complex.params <- fscRun(complex.params, num.sims = 1)
arp <- fscReadArp(complex.params)
str(arp)

## -----------------------------------------------------------------------------
str(fscReadArp(complex.params, marker = "microsat"))

## -----------------------------------------------------------------------------
fscReadArp(complex.params, marker = "microsat", one.col = TRUE)[, 1:6]

## -----------------------------------------------------------------------------
arp <- fscReadArp(complex.params, chrom = c(1, 3), one.col = TRUE)
str(arp)

## -----------------------------------------------------------------------------
arp <- fscReadArp(complex.params, sep.chrom = TRUE)
str(arp)

## -----------------------------------------------------------------------------
rm(list = ls())
library(strataG)

demes <- fscSettingsDemes(fscDeme(deme.size = 1000, sample.size = 10))
genetics <- fscSettingsGenetics(fscBlock_snp(10, 1e-6), num.chrom = 1000)
p <- fscWrite(demes = demes, genetics = genetics, label = "ex2.snps.1k")
p <- fscRun(p, all.sites = F)
snp.df <- fscReadArp(p)

# an example of the data generated
snp.df[1:6, 1:6]

## -----------------------------------------------------------------------------
snpOccurFreq <- function(mat) {
  # Extract the SNP names from the matrix column names
  snp.name <- colnames(mat[, -(1:2)])
  # Extract the chromosome name (starts with "C" and is followed by numbers) 
  #   from the SNP names
  chrom.names <- regmatches(snp.name, regexpr("^C[[:digit:]]+", snp.name))
  # Count number of occurrences of each chromosome
  chrom.freq <- table(chrom.names)
  # Get frequencies of number of occurrences (how many 1s, 2s, 3s...)
  table(chrom.freq, dnn = NULL)
}

# The occurrence frequencies
snp.occ.freq <- snpOccurFreq(snp.df)
snp.occ.freq

# Convert to proportions
snp.occ.prop <- prop.table(snp.occ.freq)
round(snp.occ.prop, 3)

## -----------------------------------------------------------------------------
sampleOnePerLocus <- function(mat) {
  # Extract the SNP names from the matrix column names
  snp.name <- colnames(mat[, -(1:2)])
  # Extract the chromosome name (starts with "C" and is followed by numbers) 
  #   from the SNP names
  chrom.names <- regmatches(snp.name, regexpr("^C[[:digit:]]+", snp.name))
  # Choose one SNP per chromosome
  one.per.loc <- tapply(colnames(mat[, -(1:2)]), chrom.names, sample, size = 1)
  # Return matrix of 
  mat[, c("id", "deme", one.per.loc)]
}

unlinked.snps <- sampleOnePerLocus(snp.df)
# number of unlinked SNPs
ncol(unlinked.snps) - 2

## -----------------------------------------------------------------------------
genetics <- fscSettingsGenetics(fscBlock_snp(10, 1e-6), num.chrom = 10000)
p <- fscWrite(demes = demes, genetics = genetics, label = "ex2.snps.10k")
p <- fscRun(p, all.sites = F)
snp.df <- fscReadArp(p)

# number of SNPs
ncol(snp.df) - 2

# proportion of n SNPs per locus
round(prop.table(snpOccurFreq(snp.df)), 3)

## -----------------------------------------------------------------------------
genetics <- fscSettingsGenetics(fscBlock_snp(10, 1e-5), num.chrom = 1000)
p <- fscWrite(demes = demes, genetics = genetics, label = "ex2.snps.mut")
p <- fscRun(p, all.sites = F, num.sims = 1)
snp.df <- fscReadArp(p)

# number of SNPs
ncol(snp.df) - 2

# proportion of n SNPs per locus
snp.occ.freq <- snpOccurFreq(snp.df)
round(prop.table(snp.occ.freq), 3)

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(fscDeme(deme.size = 5000, sample.size = 5))
genetics <- fscSettingsGenetics(fscBlock_snp(10, 1e-5), num.chrom = 1000)
p <- fscWrite(demes = demes, genetics = genetics, label = "ex2.inf.sites")
p <- fscRun(p, all.sites = F)
snp.df <- fscReadArp(p)

# number of SNPs recovered
ncol(snp.df) - 2

## -----------------------------------------------------------------------------
p <- fscRun(p, all.sites = F, inf.sites = T, num.sims = 1)
snp.df <- fscReadArp(p)

# number of SNPs recovered
ncol(snp.df) - 2

## ----message=FALSE------------------------------------------------------------
m <- 0.00001
mig.mat <- matrix(c(0, m, m, 0), nrow = 2)
mig.mat

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(fscDeme(1000, 10), fscDeme(1000, 10))
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-6), num.chrom = 1000)
p <- fscWrite(
  demes = demes,
  migration = fscSettingsMigration(mig.mat),
  genetics = genetics,
  label = "ex3.mig.ex"
)
p <- fscRun(p, all.sites = F)

## -----------------------------------------------------------------------------
snp.df <- fscReadArp(p, one.col = F)
snp.g <- df2gtypes(snp.df, ploidy = 2)
overallTest(snp.g, stat = "fst")

## ----results = "hide"---------------------------------------------------------
m.vec <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
fst <- sapply(m.vec, function(m) {
  mig.mat <- matrix(c(0, m, m, 0), nrow = 2)
  p <- fscWrite(
    demes = demes,
    migration = fscSettingsMigration(mig.mat),
    genetics = genetics, 
    label = "mig.test"
  )
  p <- fscRun(p, all.sites = F, inf.sites = T)
  snp.df <- fscReadArp(p, one.col = F)
  snp.g <- df2gtypes(snp.df, ploidy = 2)
  overallTest(snp.g, stat = c("fst"), quietly = T)$result[1, ]
})

## -----------------------------------------------------------------------------
cbind(
  m = m.vec, 
  Nm = 1000 * m.vec, 
  expFst = 1 / ((4 * 1000 * m.vec) + 1), 
  t(fst)
)

## -----------------------------------------------------------------------------
num.demes <- 5
m <- 0.0005

mig.rate <- m / (num.demes - 1)
island.mat <- matrix(rep(mig.rate, num.demes ^ 2), nrow = num.demes)
diag(island.mat) <- 0
island.mat

qgraph::qgraph(island.mat)

## -----------------------------------------------------------------------------
mig.rate <- m / 2
step.mat <- matrix(0, nrow = num.demes, ncol = num.demes)
# set rate for neighbors
for (k in 1:(num.demes - 1)) {
  step.mat[k, k + 1] <- step.mat[k + 1, k] <- mig.rate
}
# demes at ends
step.mat[1, num.demes] <- step.mat[num.demes, 1] <- mig.rate
diag(step.mat) <- 0
step.mat

qgraph::qgraph(step.mat)

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(
  fscDeme(1000, 10), fscDeme(1000, 10), fscDeme(1000, 10), 
  fscDeme(1000, 10), fscDeme(1000, 10)
)
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-5), num.chrom = 1000)
p.island <- fscWrite(
  demes = demes,
  migration = fscSettingsMigration(island.mat),
  genetics = genetics,
  label = "ex3.island"
)
p.island <- fscRun(p.island, all.sites = F)

## -----------------------------------------------------------------------------
p.step <- fscWrite(
  demes = demes,
  migration = fscSettingsMigration(step.mat),
  genetics = genetics,
  label = "ex3.stepping.stone"
)
p.step <- fscRun(p.step, all.sites = F)

## -----------------------------------------------------------------------------
# expected Fst
1 / ((4 * 1000 * m) + 1)

island.g <- df2gtypes(fscReadArp(p.island, one.col = F), ploidy = 2)
statFst(island.g)$result["estimate"]

step.g <- df2gtypes(fscReadArp(p.step, one.col = F), ploidy = 2)
statFst(step.g)$result["estimate"]

## -----------------------------------------------------------------------------
# choose 5 random points in two dimensions
set.seed(50)
deme.pos <- data.frame(
  x = runif(num.demes, 0, 0.5),
  y = runif(num.demes, 0, 0.5)
)
rownames(deme.pos) <- 1:num.demes

plot(deme.pos, type = "n")
text(deme.pos, labels = rownames(deme.pos))

## -----------------------------------------------------------------------------
euc.dist <- dist(deme.pos[, -1], diag = FALSE, upper = TRUE)
scaled.dist <- as.matrix(euc.dist / min(euc.dist))
scaled.dist

## -----------------------------------------------------------------------------
ibd.mat <- 0.05 / scaled.dist
diag(ibd.mat) <- 0
ibd.mat
qgraph::qgraph(ibd.mat)

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(
  fscDeme(1000, 10), fscDeme(1000, 10), fscDeme(1000, 10), 
  fscDeme(1000, 10), fscDeme(1000, 10)
)
p.ibd <- fscWrite(
  demes = demes,
  migration = fscSettingsMigration(ibd.mat),
  genetics = genetics,
  label = "ex3.ibd"
)
p.ibd <- fscRun(p.ibd, all.sites = F)

## ----results = "hide"---------------------------------------------------------
ibd.g <- df2gtypes(fscReadArp(p.ibd, one.col = F), ploidy = 2)
pws <- pairwiseTest(ibd.g, stat = "fst", nrep = 100)$pair.mat$Fst

## -----------------------------------------------------------------------------
pws

## ----message=FALSE------------------------------------------------------------
events <- fscSettingsEvents(
  fscEvent(
    event.time = 2000, 
    source = 1, 
    sink = 2, 
    prop.migrants = 0.05, 
    new.size = 1,
    new.growth = 0,
    migr.mat = 0
  ),
  fscEvent(
    event.time = 2980, 
    source = 1, 
    sink = 1, 
    prop.migrants = 0, 
    new.size = 0.04
  ),
  fscEvent(3000, 1, 0, 1, 1),
  fscEvent(15000, 0, 2, 1, 3)
)

## -----------------------------------------------------------------------------
param.df <- data.frame(
  N = 10 ^ runif(5, 2, 4),
  MIG = 10 ^ runif(5, -8, -5)
)

param.p <- lapply(1:nrow(param.df), function(i) {
  N <- param.df$N[i]
  MIG <- param.df$MIG[i]
  
  p <- fscWrite(
    demes = fscSettingsDemes(fscDeme(N, 5), fscDeme(N, 5)),
    migration = fscSettingsMigration(matrix(c(0, MIG, MIG, 0), nrow = 2)),
    genetics = fscSettingsGenetics(fscBlock_snp(100, 1e-6), num.chrom = 1000),
    label = paste0("param.sim.", i)
  )
  fscRun(p)
})

dir(p$folder, pattern = "param.sim")

## ----warning=FALSE------------------------------------------------------------
p <- fscWrite(
  demes = fscSettingsDemes(fscDeme("N", 5), fscDeme("N", 5)),
  migration = fscSettingsMigration(matrix(c(0, "MIG", "MIG", 0), nrow = 2)),
  genetics = fscSettingsGenetics(fscBlock_snp(100, 1e-6), num.chrom = 1000),
  def = fscSettingsDef(param.df),
  label = "param.sim"
)
p <- fscRun(p)

dir(p$folder, pattern = "param.sim.def")

## -----------------------------------------------------------------------------
obs.p <- fscWrite(
  demes = fscSettingsDemes(fscDeme(7300, 20)),
  events = fscSettingsEvents(
    fscEvent(9800, 0, 0, 0, 3.5),
    fscEvent(9900, 0, 0, 0, 1)
  ),
  genetics = fscSettingsGenetics(fscBlock_snp(10, 2.5e-6), num.chrom = 200000),
  label = "known.1PopBot20Mb"
)
obs.p <- fscRun(obs.p, dna.to.snp = TRUE, no.arl.output = TRUE, num.cores = 3)
obs.sfs <- fscReadSFS(obs.p)
str(obs.sfs)

## -----------------------------------------------------------------------------
demes <- fscSettingsDemes(fscDeme("NCUR", 20))
events <- fscSettingsEvents(
  fscEvent("TBOT", 0, 0, 0, "RESBOT"),
  fscEvent("TENDBOT", 0, 0, 0, "RESENDBOT")
)

## -----------------------------------------------------------------------------
est <- fscSettingsEst(
  fscEstParam("NCUR", is.int = TRUE, distr = "unif", min = 10, max = 100000),
   # default for is.int = TRUE and distr = "unif" 
  fscEstParam("NANC", min = 10, max = 100000),
  fscEstParam("NBOT", min = 10, max = 100000),
  fscEstParam("TBOT", min = 10, max = 10000),
  # these are "complex parameters" (only name and value are given)
  fscEstParam("RESBOT", is.int = FALSE, value = "NBOT/NCUR", output = FALSE),
  fscEstParam("RESENDBOT", is.int = FALSE, value = "NANC/NBOT", output = FALSE),
  fscEstParam("TENDBOT", value = "TBOT+100", output = FALSE),
  obs.sfs = obs.sfs$sfs$marginal[[1]]
)
est

## ----warning=FALSE------------------------------------------------------------
est.p <- fscWrite(
  demes = demes,
  events = events,
  genetics = fscSettingsGenetics(fscBlock_freq(2.5e-6)),
  est = est,
  label = "est.1PopBot20Mb"
)
est.p <- fscRun(est.p, num.sims = 10000, num.cores = 3)

## -----------------------------------------------------------------------------
dir(est.p$folder, pattern = est.p$label)

## -----------------------------------------------------------------------------
param.est <- fscReadParamEst(est.p)
str(param.est)

