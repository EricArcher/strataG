## ------------------------------------------------------------------------
options(digits = 2)
library(strataGdevel)
data(dolph.strata)
data(dolph.msats)
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
msats <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)
rownames(dolph.strata) <- dolph.strata$ids
schemes(msats) <- dolph.strata[, c("fine", "broad")]
msats <- stratify(msats, "fine")
msats <- msats[, locNames(msats)[1:4], ]

## ------------------------------------------------------------------------
statFst(msats)
statGst(msats)

## ------------------------------------------------------------------------
ran.strata <- sample(1:3, size = nInd(msats), replace = TRUE)
statFst(msats, ran.strata)
statGst(msats, ran.strata)

## ------------------------------------------------------------------------
ovl <- overallTest(msats, nrep = 100)

## ------------------------------------------------------------------------
ovl <- overallTest(msats, stats = c(statFst, statChi2), nrep = 100)

## ------------------------------------------------------------------------
pws <- pairwiseTest(msats, stats = c(statFstPrime, statGst), nrep = 100)

## ------------------------------------------------------------------------
pws

## ------------------------------------------------------------------------
popStruct <- popStructTest(msats, stats = c(statFst, statFstPrime), nrep = 100, quietly = TRUE)
popStruct

