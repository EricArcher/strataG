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
iupacCodes(c("c", "t", "t", "c", "c"))
iupacCodes(c("c", "t", "a", "c", "c"))
iupacCodes(c("g", "t", "a", "c", "c"))

## -----------------------------------------------------------------------------
iupacCodes(c("c", "t", "t", "c", "c"), all = TRUE)
iupacCodes(c("c", "t", "a", "c", "c"), all = TRUE)
iupacCodes(c("g", "t", "a", "c", "c"), all = TRUE)

## -----------------------------------------------------------------------------
createConsensus(dolph.seqs)

## -----------------------------------------------------------------------------
nd <- nucleotideDiversity(dolph.seqs)
head(nd)

## -----------------------------------------------------------------------------
library(ape)
data(dolph.seqs)
seqs.db <- as.DNAbin(dolph.seqs)

seq.smry <- summarizeSeqs(seqs.db)
head(seq.smry)

## -----------------------------------------------------------------------------
data(dolph.haps)
hap.db <- as.DNAbin(dolph.haps)

hap.lik <- sequenceLikelihoods(hap.db)

hap.lik

## -----------------------------------------------------------------------------
low.freq <- lowFreqSubs(seqs.db, min.freq = 2)

low.freq

## -----------------------------------------------------------------------------
mostDistantSequences(seqs.db, num.seqs = 5)

## -----------------------------------------------------------------------------
mostRepresentativeSequences(seqs.db, num.seqs = 5)

