## ----echo = FALSE, message = FALSE--------------------------------------------
options(digits = 2)
library(strataG)

## -----------------------------------------------------------------------------
library(ape)
data(dolph.seqs)

i <- sample(1:10, 1)
j <- sample(1:10, 1)
x <- c(rep("n", i), dolph.seqs[[1]], rep("n", j))
x
x.trimmed <- trimNs(as.DNAbin(x))
as.character(as.list(x.trimmed))

## -----------------------------------------------------------------------------
bf <- baseFreqs(dolph.seqs)
bf$site.freqs[, 1:8]
bf$base.freqs

## -----------------------------------------------------------------------------
fs <- fixedSites(dolph.seqs)
fs[1:20]

vs <- variableSites(dolph.seqs)
vs

## -----------------------------------------------------------------------------
fs <- fixedSites(dolph.seqs, bases = c("c", "t"))
fs[1:20]

vs <- variableSites(dolph.seqs, bases = c("c", "t"))
vs

## -----------------------------------------------------------------------------
iupacCode(c("c", "t", "t", "c", "c"))
iupacCode(c("c", "t", "a", "c", "c"))
iupacCode(c("g", "t", "a", "c", "c"))

## -----------------------------------------------------------------------------
validIupacCodes(c("c", "t", "t", "c", "c"))
validIupacCodes(c("c", "t", "a", "c", "c"))
validIupacCodes(c("g", "t", "a", "c", "c"))

## -----------------------------------------------------------------------------
createConsensus(dolph.seqs)

## -----------------------------------------------------------------------------
nd <- nucleotideDiversity(dolph.seqs)
head(nd)

## -----------------------------------------------------------------------------
# create gtypes
data(dolph.seqs)
data(dolph.strata)
rownames(dolph.strata) <- dolph.strata$id
dloop <- df2gtypes(dolph.strata[, c("id", "fine", "id")], ploidy = 1,
             schemes = dolph.strata[, c("fine", "broad")], sequences = dolph.seqs)
dloop <- labelHaplotypes(dloop, "Hap.")

# calculate divergence
nucleotideDivergence(dloop)

## -----------------------------------------------------------------------------
fixedDifferences(dloop)

## -----------------------------------------------------------------------------
x <- as.DNAbin(dolph.seqs)
mostDistantSequences(x, num.seqs = 5)

## -----------------------------------------------------------------------------
mostRepresentativeSequences(x, num.seqs = 5)

