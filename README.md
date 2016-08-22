![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/strataG)
![](http://cranlogs.r-pkg.org/badges/last-day/strataG?color=red)
![](http://cranlogs.r-pkg.org/badges/last-week/strataG?color=red)
![](http://cranlogs.r-pkg.org/badges/strataG?color=red)
![](http://cranlogs.r-pkg.org/badges/grand-total/strataG?color=red)  
[![DOI](https://zenodo.org/badge/23926/EricArcher/strataG.svg)](https://zenodo.org/badge/latestdoi/23926/EricArcher/strataG)  
[![Travis-CI Build Status](https://travis-ci.org/EricArcher/strataG.svg?branch=master)](https://travis-ci.org/EricArcher/strataG)
[![Build status](https://ci.appveyor.com/api/projects/status/w6tlgm1yl53js2rq?svg=true)](https://ci.appveyor.com/project/EricArcher/stratag-4j93o)  
# strataG

## Description

*strataG* is a toolkit for haploid sequence and multilocus genetic data summaries, and analyses of population structure.

## Installation

To install the stable version from CRAN:

```r
install.packages('strataG')
```

To install the latest version from GitHub:

```r
# make sure you have Rtools installed
if (!require('devtools')) install.packages('devtools')
# install from GitHub
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
```

## Contact

* submit suggestions and bug-reports: <https://github.com/ericarcher/strataG/issues>
* send a pull request: <https://github.com/ericarcher/strataG/>
* e-mail: <eric.archer@noaa.gov>

## version 1.0.6

* Added `read.arlequin` back. Fixed missing function error with `write.arlequin`.
* Added `summarizeSamples`
* Changed `evanno` from base graphics to ggplot2
* Updated logic in `labelHaplotypes` to assign haplotypes if possible alternative site combinations match a present haplotype
* Added Zenodo DOI
* Added shiny app (`strataGUI`) for creating gtypes objects, QA/QC, and population structure analyses
* Added `type` argument to `structurePlot` to select between area and bar charts
* Changed `haplotypeLikelihoods` to `sequenceLikelihoods`
* `neiDa` now creates haplotypes before calculating metric

## version 1.0.5

* Fixed error in dupGenotypes, propSharedLoci, and propSharedIDs where missing genotypes were not being properly counted.
* Added as.data.frame.gtypes.
* Removed gtypes2df.
* Added arguments to as.matrix.gtypes to include id and strata columns in output.
* Removed the jmodeltest function as this functionality is available in the modeltest function in the phangorn package.
* Added conversion functions gtypes2phyDat and phyDat2gtypes to facilitate interoperability with the phangorn package.
* Removed read.arlequin.
* Added alleleNames accessor for gtypes object, which returns list of allele names for each locus.

## version 1.0

* New version with different gtypes format from previous versions. See vignettes for instructions and examples.