## ------------------------------------------------------------------------
options(digits = 2)
library(strataG.devel)
data(dolph.strata)
data(dolph.msats)
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
msats <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)
rownames(dolph.strata) <- dolph.strata$ids
schemes(msats) <- dolph.strata[, c("fine", "broad")]
msats <- stratify(msats, "broad")
msats <- subset(msats, loci = locNames(msats)[1:4])

## ------------------------------------------------------------------------
numAlleles(msats)

## ------------------------------------------------------------------------
numMissing(msats)

## ------------------------------------------------------------------------
numMissing(msats, prop = TRUE)

## ------------------------------------------------------------------------
allelicRichness(msats)

## ------------------------------------------------------------------------
obsvdHet(msats)
exptdHet(msats)

## ------------------------------------------------------------------------
pctUniqueAlleles(msats)

## ------------------------------------------------------------------------
theta(msats)

## ------------------------------------------------------------------------
summarizeLoci(msats)
summarizeLoci(msats, by.strata = TRUE)

## ------------------------------------------------------------------------
alleleFreqs(msats)
alleleFreqs(msats, by.strata = TRUE)

