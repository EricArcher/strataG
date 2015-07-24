## ------------------------------------------------------------------------
library(strataG.devel)
data(msats.g)
smry <- summarizeLoci(msats.g)
head(smry)

## ------------------------------------------------------------------------
# Find samples that share alleles at 2/3rds of the loci
dupGenotypes(msats.g, num.shared = 0.66, num.cores = 2)

## ------------------------------------------------------------------------
data(dolph.seqs)
seq.smry <- summarizeSeqs(as.DNAbin(dolph.seqs))
head(seq.smry)

## ------------------------------------------------------------------------
bf <- baseFreqs(as.DNAbin(dolph.seqs))

# nucleotide frequencies by site
bf$site.freq[, 1:15]

# overall nucleotide frequencies
bf$base.freqs

## ------------------------------------------------------------------------
lowFreqSubs(as.DNAbin(dolph.seqs), min.freq = 2)

## ------------------------------------------------------------------------
data(dolph.haps)
haplotypeLikelihoods(as.DNAbin(dolph.haps))

