% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/privateAlleles.R
\name{privateAlleles}
\alias{privateAlleles}
\title{Private Alleles}
\usage{
privateAlleles(g)
}
\arguments{
\item{g}{a \linkS4class{gtypes} object.}
}
\value{
matrix with the number of private alleles in each strata at each 
  locus. This is the number of alleles only present in one stratum.
}
\description{
The number of private alleles in each strata and locus.
}
\examples{
data(msats.g)

privateAlleles(msats.g)

}
\seealso{
\link{propUniqueAlleles}
}
\author{
Eric Archer \email{eric.archer@noaa.gov}
}
