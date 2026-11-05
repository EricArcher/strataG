## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = NA,
  fig.height = 8,
  fig.width = 8
)
options(
  scipen = 4, 
  digits = 2,
  rmarkdown.html_vignette.check_title = FALSE
)
library(strataG)
library(tidyverse)

## -----------------------------------------------------------------------------
ind.smry <- summarizeInds(msats.g)

head(ind.smry)

tail(ind.smry)

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
summarizeLoci(msats)
summarizeLoci(msats, by.strata = TRUE)

## -----------------------------------------------------------------------------
thetaH(msats, model = 'IAM')

thetaH(msats, model = 'SMM')

thetaH(msats, model = 'ZourosSMM')

## -----------------------------------------------------------------------------
alleleFreqs(msats)

af <- alleleFreqs(msats, by.strata = TRUE)
str(af)

## -----------------------------------------------------------------------------
# Find samples that share alleles at 2/3rds of the loci
dupes <- dupGenotypes(msats, num.shared = 0.66)
head(dupes)

