## ----echo = FALSE, message = FALSE--------------------------------------------
options(digits = 2)
library(strataG)

## -----------------------------------------------------------------------------
data(msats.g)
msats <- stratify(msats.g, "broad")
msats <- msats[, getLociNames(msats)[1:4], ]

## -----------------------------------------------------------------------------
numAlleles(msats)

## -----------------------------------------------------------------------------
numMissing(msats)

## -----------------------------------------------------------------------------
numMissing(msats, prop = TRUE)

## -----------------------------------------------------------------------------
allelicRichness(msats)

## -----------------------------------------------------------------------------
# observed
heterozygosity(msats, type = "observed")

# expected
heterozygosity(msats, type = "expected")

## -----------------------------------------------------------------------------
propUniqueAlleles(msats)

## -----------------------------------------------------------------------------
theta(msats)

## -----------------------------------------------------------------------------
summarizeLoci(msats)
summarizeLoci(msats, by.strata = TRUE)

## -----------------------------------------------------------------------------
alleleFreqs(msats)
alleleFreqs(msats, by.strata = TRUE)

## -----------------------------------------------------------------------------
# Find samples that share alleles at 2/3rds of the loci
dupGenotypes(msats, num.shared = 0.66)

## -----------------------------------------------------------------------------
library(ape)
data(dolph.seqs)
seq.smry <- summarizeSeqs(as.DNAbin(dolph.seqs))
head(seq.smry)

## -----------------------------------------------------------------------------
bf <- baseFreqs(as.DNAbin(dolph.seqs))

# nucleotide frequencies by site
bf$site.freq[, 1:15]

# overall nucleotide frequencies
bf$base.freqs

## -----------------------------------------------------------------------------
lowFreqSubs(as.DNAbin(dolph.seqs), min.freq = 2)

## -----------------------------------------------------------------------------
data(dolph.haps)
sequenceLikelihoods(as.DNAbin(dolph.haps))

