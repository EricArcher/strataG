## ----echo = FALSE, message = FALSE--------------------------------------------
options(digits = 2)
library(strataG)

## -----------------------------------------------------------------------------
data(msats.g)
msats <- stratify(msats.g, "fine")
msats <- msats[, getLociNames(msats)[1:4], ]

## -----------------------------------------------------------------------------
statFst(msats)

statGst(msats, nrep = 10, keep.null = TRUE)

## -----------------------------------------------------------------------------
ovl <- overallTest(msats, stats = c("fst", "chi2"), nrep = 1000)

## -----------------------------------------------------------------------------
pws <- pairwiseTest(msats, stats = c("fst.prime", "gst"), nrep = 1000)

## -----------------------------------------------------------------------------
pws

## -----------------------------------------------------------------------------
popStruct <- popStructTest(msats, stats = c("fst", "fst.prime"), nrep = 1000, quietly = TRUE)
popStruct

