## ------------------------------------------------------------------------
# load msat data and create gtypes
library(strataG.devel)
data(dolph.strata)
data(dolph.msats)
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
msats <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)
rownames(dolph.strata) <- dolph.strata$ids
schemes(msats) <- dolph.strata[, c("fine", "broad")]
msats <- stratify(msats, "fine")

# run 10 replicates of STRUCTURE for k = 2
msat.struct <- structureRun(msats, k.range = 2, num.k.rep = 10)

## ------------------------------------------------------------------------
names(msat.struct[[1]])

## ------------------------------------------------------------------------
msat.struct[[1]]$summary

## ------------------------------------------------------------------------
head(msat.struct[[1]]$q.mat)

## ------------------------------------------------------------------------
msat.struct <- structureRun(msats, k.range = 1:5, num.k.rep = 5)
evanno(msat.struct)

## ------------------------------------------------------------------------
msat.clmp <- clumpp(msat.struct, k = 3)
head(msat.clmp)

## ------------------------------------------------------------------------
structurePlot(msat.clmp, horiz = F)

