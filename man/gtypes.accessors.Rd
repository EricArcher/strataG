% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtypes.accessors.R
\name{gtypes.accessors}
\alias{gtypes.accessors}
\alias{accessors}
\alias{getNumInd,gtypes-method}
\alias{getNumInd}
\alias{getNumLoci,gtypes-method}
\alias{getNumLoci}
\alias{getNumStrata}
\alias{getNumStrata,gtypes-method}
\alias{getIndNames}
\alias{getIndNames,gtypes-method}
\alias{getLociNames}
\alias{getLociNames,gtypes-method}
\alias{getAlleleNames}
\alias{getAlleleNames,gtypes-method}
\alias{getStrataNames}
\alias{getStrataNames,gtypes-method}
\alias{getPloidy}
\alias{getPloidy,gtypes-method}
\alias{getStrata}
\alias{getStrata,gtypes-method}
\alias{setStrata<-}
\alias{setStrata<-,gtypes-method}
\alias{setStrata}
\alias{getSchemes}
\alias{getSchemes,gtypes-method}
\alias{setSchemes<-}
\alias{setSchemes<-,gtypes-method}
\alias{setSchemes}
\alias{getSequences}
\alias{getSequences,gtypes-method}
\alias{getDescription}
\alias{getDescription,gtypes-method}
\alias{setDescription<-}
\alias{setDescription<-,gtypes-method}
\alias{setDescription}
\alias{getOther}
\alias{getOther,gtypes-method}
\alias{setOther<-}
\alias{setOther<-,gtypes-method}
\alias{setOther}
\alias{[,gtypes,ANY,ANY,ANY-method}
\alias{index}
\alias{subset}
\title{\code{gtypes} Accessors}
\usage{
\S4method{getNumInd}{gtypes}(x, by.strata = FALSE, ...)

\S4method{getNumLoci}{gtypes}(x, ...)

getNumStrata(x, ...)

\S4method{getNumStrata}{gtypes}(x, ...)

getIndNames(x, ...)

\S4method{getIndNames}{gtypes}(x, by.strata = FALSE, ...)

getLociNames(x, ...)

\S4method{getLociNames}{gtypes}(x, ...)

getAlleleNames(x, ...)

\S4method{getAlleleNames}{gtypes}(x, ...)

getStrataNames(x, ...)

\S4method{getStrataNames}{gtypes}(x, ...)

getPloidy(x, ...)

\S4method{getPloidy}{gtypes}(x, ...)

getStrata(x, ...)

\S4method{getStrata}{gtypes}(x)

setStrata(x) <- value

\S4method{setStrata}{gtypes}(x) <- value

getSchemes(x, ...)

\S4method{getSchemes}{gtypes}(x, ...)

setSchemes(x) <- value

\S4method{setSchemes}{gtypes}(x) <- value

getSequences(x, ...)

\S4method{getSequences}{gtypes}(
  x,
  as.haplotypes = FALSE,
  seqName = NULL,
  as.multidna = FALSE,
  simplify = TRUE,
  ...
)

getDescription(x, ...)

\S4method{getDescription}{gtypes}(x, ...)

setDescription(x) <- value

\S4method{setDescription}{gtypes}(x) <- value

getOther(x, ...)

\S4method{getOther}{gtypes}(x, value = NULL, ...)

setOther(x, name) <- value

\S4method{setOther}{gtypes}(x, name) <- value

\S4method{[}{gtypes,ANY,ANY,ANY}(x, i, j, k, ..., quiet = TRUE, drop = FALSE)
}
\arguments{
\item{x}{a \linkS4class{gtypes} object.}

\item{by.strata}{logical - return results by strata?}

\item{...}{other arguments passed from generics (ignored).}

\item{value}{value being assigned by accessor.}

\item{as.haplotypes}{return sequences as haplotypes? If \code{TRUE}, contents
of \code{@sequences} slot are returned. If \code{FALSE}, one sequence per
individual is returned.}

\item{seqName}{the name (or number) of a set of sequences from the
\code{@sequences} slot to return.}

\item{as.multidna}{return sequences as a \code{multidna} object? If
\code{FALSE}, sequences are returned as a list.}

\item{simplify}{if `getSequences()` would return a single locus, return it 
as a `DNAbin` object (`TRUE`), or a single element named list (`FALSE`).}

\item{name}{name of the value going into the \code{other} list.}

\item{i, j, k}{subsetting slots for individuals (\code{i}), loci (\code{j}), or
strata (\code{k}). See Details for more information.}

\item{quiet}{suppress warnings about unmatched requested individuals, loci,
or strata?}

\item{drop}{if \code{TRUE} the return object will have unused sequences
removed.}
}
\value{
\describe{ \item{nInd}{number of individuals} \item{nLoc}{number of
loci} \item{nStrata}{number of strata} \item{indNames}{vector of
individual/sample names} \item{locNames}{vector of locus names}
\item{strataNames}{vector of strata names for current scheme}
\item{ploidy}{number of alleles at each locus} \item{other}{contents of
\code{@other} slot} \item{strata}{return or modify the current
stratification} \item{schemes}{return or modify the current stratification
schemes} \item{alleleNames}{return a list of alleles at each locus}
\item{sequences}{return the \code{multidna} object in the
\code{@sequences} slot. See \code{\link[apex]{getSequences}} to extract
individual genes or sequences from this object} \item{description}{return the
object's description} }
}
\description{
Accessors for slots in \linkS4class{gtypes} objects.
}
\details{
Indexing a \code{gtypes} object with integers, characters, or
logicals with the \code{[} operator follows the same rules as normal indexing
in R. The order that individuals, loci, and strata are chosen is the order
returned by \code{getIndNames}, \code{getLocNames}, and \code{getStrataNames}
respectively. If unstratified samples are present, they can be selected as a
group either by including \code{NA} in the character or numeric vector of the
\code{k} slot, or by providing a logical vector based on
\code{is.na(strata(g))} to the \code{i} slot.
When setting stratification schemes with \code{setSchemes() <- value}, the
data frame being assigned must have the individual ids in a column named 
'\code{id}' and that column must be the first column in the data frame.
}
\examples{
#--- create a diploid (microsatellite) gtypes object
data(msats.g)
msats.g <- stratify(msats.g, "fine")

getNumStrata(msats.g)
getStrataNames(msats.g)
getNumLoci(msats.g)
getLociNames(msats.g)

# reassign all samples to two randomly chosen strata
new.strata <- sample(c("A", "B"), getNumInd(msats.g), rep = TRUE)
names(new.strata) <- getIndNames(msats.g)
setStrata(msats.g) <- new.strata
msats.g


#--- a sequence example
library(ape)
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
wood.g <- sequence2gtypes(x)
new.strata <- sample(c("A", "B"), getNumInd(wood.g), rep = TRUE)
names(new.strata) <- getIndNames(wood.g)
setStrata(wood.g) <- new.strata
wood.g

# get the multidna sequence object
multi.seqs <- getSequences(wood.g, as.multidna = TRUE)
class(multi.seqs) # "multidna"

# get a list of DNAbin objects
dnabin.list <- getSequences(wood.g)
class(dnabin.list) # "list"

# get a DNAbin object of the first locus
dnabin.1 <- getSequences(wood.g)[[1]]
class(dnabin.1) # "DNAbin"

# getting and setting values in the `other` slot:
getOther(dloop.g)

setOther(dloop.g, "timestamp") <- timestamp()
setOther(dloop.g, "Author") <- "Hoban Washburne"

getOther(dloop.g)
getOther(dloop.g, "timestamp")

setOther(dloop.g, "Author") <- NULL
getOther(dloop.g)


}
\author{
Eric Archer \email{eric.archer@noaa.gov}
}
