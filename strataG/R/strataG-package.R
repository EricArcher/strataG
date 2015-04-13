#' strataG
#' 
#' @docType package
#' @name strataG-package
#' @aliases strataG
#' @title Summaries and population structure analyses of haplotypic and genotypic data
#' @keywords package
#' @useDynLib strataG
NULL

#' @docType data
#' @name dolph.strata
#' @title Dolphin Genetic Stratification and Haplotypes
#' @description Assignment of samples to one of two stratifications and mtDNA haplotype designations
#' @usage data(dolph.strata)
#' @format A data.frame of 126 samples and 4 variables
#' @references Lowther et al.
#' @keywords datasets
NULL

#' @docType data
#' @name dolph.msats
#' @title Dolphin Microsatellite Genotypes
#' @description Genotypes for 15 microsatellite loci
#' @usage data(dolph.msats)
#' @format A data.frame of 126 samples and 31 variables
#' @references Lowther et al.
#' @keywords datasets
NULL

#' @docType data
#' @name dolph.seqs
#' @title Dolphin mtDNA Sequences
#' @usage data(dolph.seqs)
#' @format A list of 126 aligned control region sequences
#' @references Lowther et al.
#' @keywords datasets
NULL

#' @docType data
#' @name dolph.haps
#' @title Dolphin mtDNA Haplotype Sequences
#' @usage data(dolph.haps)
#' @format A list of 33 aligned sequences
#' @references Lowther et al.
#' @keywords datasets
NULL

#' @docType data
#' @name bowhead.snps
#' @title Bowhead Whale SNP Genotypes
#' @usage data(bowhead.snps)
#' @format A data.frame of 42 SNPs with sample ids and stratification
#' @references Morin et al.
#' @keywords datasets
NULL

#' @docType data
#' @name bowhead.snp.position
#' @title Bowhead Whale SNP Genotype Groups
#' @usage data(bowhead.snp.position)
#' @format A data.frame of position information for SNPs to be phased
#' @references Morin et al.
#' @keywords datasets
NULL

#' @docType data
#' @name gtype.struct.func
#' @aliases gtype.struct.stat
#' @title Population Structure Statistic Functions
#' @description Functions that calculate population structure statistics. 
#' The function takes a \code{\link{gtypes}} object and perhaps other parameters and returns
#' an object of class \code{gtype.struct.stat} which is a list containing:\cr
#' \tabular{ll}{
#'   \code{stat.name} \tab the name of the statistic.\cr
#'   \code{estimate} \tab the value of the statistic.\cr
#'   \code{strata.freq} \tab a \code{table} of the strata frequencies used in calculating the statistic.
#'This may be different from the frequencies in the original \code{gtypes} object due to exclusion of 
#'samples due to missing data.\cr
#'}
#' @keywords functions
NULL