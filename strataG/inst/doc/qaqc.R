## ------------------------------------------------------------------------
library(strataGdevel)
data(dolph.strata)
data(dolph.msats)
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
msats <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)
smry <- summarizeLoci(msats)
head(smry)

## ------------------------------------------------------------------------
# Find samples that share alleles at 2/3rds of the loci
dupGenotypes(msats, num.shared = 0.66, num.cores = 2)

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

