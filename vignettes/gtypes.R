## ----echo = FALSE, message = FALSE--------------------------------------------
library(strataG)

## -----------------------------------------------------------------------------
gen.data <- readGenData("msats.csv")
str(gen.data)

## -----------------------------------------------------------------------------
fname <- system.file("extdata/dolph.seqs.fasta", package = "strataG")
x <- read.fasta(fname) 
x

## -----------------------------------------------------------------------------
# create a single data.frame with the msat data and stratification
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
str(msats.merge)

# create the gtypes object
msats.fine <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)

## -----------------------------------------------------------------------------
data(dolph.seqs)

seq.df <- dolph.strata[ c("id", "broad", "id")]
colnames(seq.df)[3] <- "D-loop"
dl.g <- df2gtypes(seq.df, ploidy = 1, sequences = dolph.seqs)
dl.g

## -----------------------------------------------------------------------------
dl.haps <- labelHaplotypes(dl.g)
dl.haps

## -----------------------------------------------------------------------------
data(dolph.haps)

haps.g <- sequence2gtypes(dolph.haps)
haps.g

## -----------------------------------------------------------------------------
# extract and name the stratification scheme
strata <- dolph.strata$fine
names(strata) <- dolph.strata$ids

# create the gtypes object
dloop.fine <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop",
  description = "dLoop: fine-scale stratification")
dloop.fine

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
sub.msats <- msats.fine[sample(getNumInd(msats.fine), 10), , ]
sub.msats

## -----------------------------------------------------------------------------
sub.msats <- sub.msats[, c("D11t", "EV37", "TV7"), ]
sub.msats

## -----------------------------------------------------------------------------
sub.msats <- msats.fine[, c("Ttr11", "D11t"), "Coastal"]
sub.msats

## -----------------------------------------------------------------------------
summarizeLoci(msats.fine)
summarizeInds(msats.fine)

## -----------------------------------------------------------------------------
# randomly stratify individuals to two populations
msats <- msats.g
new.strata <- sample(c("Pop1", "Pop2"), getNumInd(msats), rep = TRUE)
names(new.strata) <- getIndNames(msats)
setStrata(msats) <- new.strata
msats

## -----------------------------------------------------------------------------
# choose "broad" stratification scheme
msats <- stratify(msats, "broad")
msats

## -----------------------------------------------------------------------------
new.schemes <- getSchemes(msats)
new.schemes$ran.pop <- sample(c("Pop5", "Pop6"), getNumInd(msats), rep = TRUE)
setSchemes(msats) <- new.schemes

## -----------------------------------------------------------------------------
stratify(msats, "ran.pop")

## -----------------------------------------------------------------------------
# unstratify a random 10 samples
x <- getStrata(msats)
x[sample(getIndNames(msats), 10)] <- NA
msats

## -----------------------------------------------------------------------------
msats <- stratify(msats, "fine")

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

