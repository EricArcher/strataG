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
dloop.g

## -----------------------------------------------------------------------------
nucleotideDivergence(dloop.g)

## -----------------------------------------------------------------------------
ovl <- overallTest(dloop.g)
ovl

## -----------------------------------------------------------------------------
pws <- pairwiseTest(dloop.g)

str(pws)

## -----------------------------------------------------------------------------
pairwiseSummary(pws)

## -----------------------------------------------------------------------------
pairwiseMatrix(pws, stat = 'Fst')

## -----------------------------------------------------------------------------
pws.msats <- pairwiseTest(msats.g)
pairwiseSummary(pws.msats)

