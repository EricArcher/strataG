rm(list = ls())
library(strataG.devel)
data(dolph.strata)
data(dolph.msats)

msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
msats <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)
rownames(dolph.strata) <- dolph.strata$ids
schemes(msats) <- dolph.strata[, c("fine", "broad")]
msats.g <- stratify(msats, "fine")
save(msats.g, file = "msats.g.rda")

data(dolph.seqs)
dloop.haps <- cbind(dLoop = dolph.strata$id)
rownames(dloop.haps) <- dolph.strata$id
strata.schemes <- dolph.strata[, c("broad", "fine")]
rownames(strata.schemes) <- dolph.strata$id
dloop <- new("gtypes", gen.data = dloop.haps, ploidy = 1,
             schemes = strata.schemes, sequences = dolph.seqs,
             strata = "fine")
dloop.g <- labelHaplotypes(dloop, "Hap.")$gtypes
save(dloop.g, file = "dloop.g.rda")
