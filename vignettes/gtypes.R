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
gen.data <- readGenData("msats.csv")
str(gen.data)

## -----------------------------------------------------------------------------
fname <- system.file("extdata/dolph.seqs.fasta", package = "strataG")
dolph.seqs <- read.fasta(fname) 
dolph.seqs

## -----------------------------------------------------------------------------
schemes <- readGenData('schemes.csv')
str(schemes)
head(schemes)
tail(schemes)

## -----------------------------------------------------------------------------
schemes <- readGenData('schemes.csv')
gen.data <- readGenData("msats.csv")

gen.df <- schemes |> 
  select(id, region) |> 
  left_join(gen.data, by = 'id')

str(gen.df)

## -----------------------------------------------------------------------------
msats.g <- df2gtypes(
  gen.df,
  ploidy = 2
)

msats.g

## -----------------------------------------------------------------------------
msats.g <- df2gtypes(
  gen.df,
  ploidy = 2,
  schemes = schemes,
  description = "Example microsatellite data"
)

msats.g

## -----------------------------------------------------------------------------
# grab first allele of first locus to simulate haploid data

#' note that we can use pipes here to read the data directly and 
#' create a gtypes object in the same pipeline

hap.g <- readGenData('schemes.csv') |> 
  select(id, ecotype) |> 
  left_join(
    select(readGenData("msats.csv"), id, D11t.1), 
    by = 'id'
  ) |> 
  df2gtypes(
    ploidy = 1, 
    description = 'haploid example'
  )

hap.g

## -----------------------------------------------------------------------------
seq.df <- schemes |> 
  # we just need the id and strata assignments
  select(id, region) |>  
  # add a new column for the sequence ids
  mutate(dLoop = id) |>
  # remove ids without sequences
  filter(dLoop %in% names(dolph.seqs)) 

dolph.seqs <- read.fasta(
  system.file("extdata/dolph.seqs.fasta", package = "strataG")
) 

seqs.g <- df2gtypes(
  seq.df,
  ploidy = 1,
  sequences = dolph.seqs,
  schemes = schemes,
  description = 'Dolphin mtDNA sequences'
)

seqs.g

## -----------------------------------------------------------------------------
haps.g <- labelHaplotypes(seqs.g)
haps.g

## -----------------------------------------------------------------------------
seqs.g <- sequence2gtypes(dolph.seqs)
seqs.g

## -----------------------------------------------------------------------------
# extract the stratifications, remove unstratified individuals
eco.strata <- schemes |> 
  select(id, ecotype) |> 
  filter(!is.na(ecotype)) |> 
  deframe()

tail(eco.strata)

## -----------------------------------------------------------------------------
# create the gtypes object
seqs.g <- sequence2gtypes(
  dolph.seqs, 
  strata = eco.strata, 
  seq.names = "dLoop",
  description = "dLoop: ecotype stratification"
)

seqs.g

## -----------------------------------------------------------------------------
library(adegenet)
# from example(df2genind)
df <- data.frame(locusA=c("11","11","12","32"),
                 locusB=c(NA,"34","55","15"),
                 locusC=c("22","22","21","22"))
row.names(df) <- .genlab("genotype",4)
obj <- df2genind(df, ploidy=2, ncode=1)
obj

# convert to gtypes
gi.g <- genind2gtypes(obj)
gi.g

## -----------------------------------------------------------------------------
# get 10 random ids
ran.ids <- sample(getNumInd(msats.g), 5)
head(ran.ids)

# index to extract
sub.msats <- msats.g[ran.ids, , ]
sub.msats

## -----------------------------------------------------------------------------
sub.msats <- sub.msats[, c("D11t", "EV37", "TV7"), ]
sub.msats

## -----------------------------------------------------------------------------
sub.msats <- msats.g[, c("Ttr11", "D11t"), "Coastal"]
sub.msats

## -----------------------------------------------------------------------------
# randomly stratify individuals to two populations
msats <- msats.g
new.strata <- sample(c("Pop1", "Pop2"), getNumInd(msats), rep = TRUE)
names(new.strata) <- getIndNames(msats)
setStrata(msats) <- new.strata
msats

## -----------------------------------------------------------------------------
# choose "broad" stratification scheme
msats <- stratify(msats, "ecotype")
msats

## -----------------------------------------------------------------------------
new.schemes <- getSchemes(msats)
new.schemes$inner.outer <- sample(
  c("Inner", "Outer", NA), 
  nrow(new.schemes), 
  rep = TRUE
)
setSchemes(msats) <- new.schemes

## -----------------------------------------------------------------------------
msats.inner.outer <- stratify(msats, "inner.outer")
msats.inner.outer

## -----------------------------------------------------------------------------
msats <- stratify(msats, "offshore")

# original
msats

# permuted
ran.msats <- permuteStrata(msats)
ran.msats

## -----------------------------------------------------------------------------
gen.mat <- as.matrix(msats)
head(gen.mat)

## -----------------------------------------------------------------------------
gen.mat <- as.matrix(msats, one.col = TRUE)
head(gen.mat)

