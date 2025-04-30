[![CRAN version](http://www.r-pkg.org/badges/version/strataG?color=red)](https://cran.r-project.org/package=strataG)
[![CRAN last day downloads](http://cranlogs.r-pkg.org/badges/last-day/strataG?color=red)](https://cran.r-project.org/package=strataG)
[![CRAN last week downloads](http://cranlogs.r-pkg.org/badges/last-week/strataG?color=red)](https://cran.r-project.org/package=strataG)
[![CRAN last month downloads](http://cranlogs.r-pkg.org/badges/strataG?color=red)](https://cran.r-project.org/package=strataG)
[![CRAN total downloads](http://cranlogs.r-pkg.org/badges/grand-total/strataG?color=red)](https://cran.r-project.org/package=strataG)  
[![Zenodo DOI](https://zenodo.org/badge/23926/EricArcher/strataG.svg)](https://zenodo.org/badge/latestdoi/23926/EricArcher/strataG)  
[![Travis-CI Build Status](https://app.travis-ci.com/EricArcher/strataG.svg?branch=master)](https://app.travis-ci.com/EricArcher/strataG)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/EricArcher/strataG?branch=master&svg=true)](https://ci.appveyor.com/project/EricArcher/strataG)  

# strataG

## Description

*strataG* is a toolkit for haploid sequence and multilocus genetic data summaries, and analyses of population structure. One can select select specific individuals, loci, or strata using standard R '[' indexing methods. . The package contains functions for summarizing haploid and diploid loci (e.g., allelic richness, heterozygosity, haplotypic diversity, etc.), and haploid sequences by locus and by strata as well as functions for computing by-site base frequencies and identifying variable and fixed sites among strata. There are both overall and pairwise standard tests of population structure like PHIst, Fst, Gst, and Jost's D. If individuals are stratified according to multiple schemes, these stratifications can be changed with the `stratify()` function and summaries or tests can be re-run on the new object. The package also includes wrappers for several external programs like *fastsimcoal2*, *STRUCTURE*, and *mafft*. There are also multiple conversion functions for data objects for other population packages such as *adegenet*, *pegas*, and *phangorn*.

## Installation

To install the latest version from GitHub:

```r
# make sure you have devtools installed
if (!require('devtools')) install.packages('devtools')
# install strataG latest version
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
```

## Compilation install error

When installing on a Mac or LINUX system, you may get errors when 
installing during the compilation phase that look like:

```
─  installing *source* package ‘strataG’ ...
   ** using staged installation
   ** libs

...

   ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0' not found
   ld: warning: search path '/opt/gfortran/lib' not found
   ld: library 'gfortran' not found
   clang++: error: linker command failed with exit code 1 (use -v to see invocation)
   make: *** [strataG.so] Error 1
   ERROR: compilation failed for package ‘strataG’
─  removing ‘/private/var/folders/rx/z5h877kx2rx_85q95ct3fdfr0000gn/T/RtmpgBBQQU/Rinstfc0912a21e41/strataG’
         -----------------------------------
   ERROR: package installation failed
Error: Failed to install 'strataG' from GitHub:
  ! System command 'R' failed
```

If you see this, follow these steps to update the compilation-related paths in your `Makeconf` file:

1. Confirm that you have `gfortran` and `gcc` installed. Open Terminal and type

```
<me> ~ % which gfortran
/usr/local/bin/gfortran
<me> ~ % which gcc
/usr/bin/gcc
```

If you see "`gfortran not found`", then you need to install it.
First, if you are on a Mac, make sure you have Xcode Command Line Tools installed:

```
xcode-select --install
```

If `which gfortran` still returns `gfortran not found` and you have `homebrew` then try:

```
brew install gcc
```

2. Locate your `Makeconf` file. Check `/Library/Frameworks/R.framework/Resources/etc`:

```
<me> ~ % ls /Library/Frameworks/R.framework/Resources/etc
Makeconf	Renviron	javaconf	ldpaths		repositories
```

3. Open `Makeconf` in a text editor and find the lines labeled `#Fortran`.

4. Update the `FC` and `F77` paths to match the location given in `which gfortran`:

```
FC = /usr/local/gfortran/bin/gfortran -arch arm64
F77 = /usr/local/gfortran/bin/gfortran
```

5. Search in the same `gfortran` directory for `gcc`, perhaps in `lib`. 
Update `FLIBS` to match the paths for gcc as below:

```
FLIBS =  -L/usr/local/gfortran/lib/gcc/aarch64-apple-darwin23/14.1.0 -L/usr/local/gfortran/lib -lgfortran -lemutls_w -lquadmath
```

For more help, also check out [this](https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r) thread on Stack Overflow.




## Vignettes

Vignettes are available on several topics:

* Creating and manipulating gtypes ("gtypes")
* Genotype and sequence summaries ("summaries")
* Working with sequences ("sequences")
* Tests of population structure ("population.structure")
* Installing external programs ("external.programs")

To see the list of all available vignettes:
```r
browseVignettes("strataG")
```

To open a specific vignette:
```r
vignette("gtypes", "strataG")
```

There is also a tutorial detailing running _fastsimcoal2_ through _strataG_ 
available through the function `fscTutorial()`.

## Citation

The paper can be obtained 
[here](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12559/abstract), 
and is cited as (preferred):   

Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016), *strataG*: An *R* 
package for manipulating, summarizing and analysing population genetic data. 
Mol Ecol Resour. doi:10.1111/1755-0998.12559

If desired, the current release version of the package can be cited as:  

Archer, F. 2025. *strataG*: An *R* package for manipulating, summarizing and 
analysing population genetic data. R package version 1.0.6. 
Zenodo. http://doi.org/10.5281/zenodo.60416  

## Contact

* submit suggestions and bug-reports: <https://github.com/ericarcher/strataG/issues>
* send a pull request: <https://github.com/ericarcher/strataG/>
* e-mail: <eric.ivan.archer@gmail.com>

***  
## Versions

## 2.5.5 (devel)
* added `diagnosability()`
* added `microhaplot2rubias()`
* added `qual2prob()`
* made `rfPermute` and `rmetasim` packages Suggested
* removed melt from structurePlot
* fixed ldNe error when one individual is present
* fixed mafft error and now have mafft .fasta files written to temporary file rather than working directory
* fixed error with `readGenData()` not recognizing `NA`s.
* fixed error with `fs2gtypes()` not formatting multi-block DNA sequence data as gtypes properly

## 2.4.905
* Deleted functions: `alleleFreqFormat`, `as.array.gtypes`
* Changed structure of `gtypes` object, making it __no longer compatible with previous versions__
* Fixed and enhanced `arlequinRead()` so that it will read and parse all .arp files. Added `arp2gtypes()` to create `gtypes` object from parsed .arp files.
* Improved performance of several standard summary functions, most notably `dupGenotypes()`.
* Full rework of _fastsimcoal2_ wrapper. 
* Removed `strataGUI()`.

## 2.1   
* fixed error in ldNe when missing data are present
* added STANDARD marker type to fastsimcoal
* added `na.rm = TRUE` to calculation of mean locus summaries by strata in `summary.gtypes`. This avoids `NaN`s when there is a locus with genotypes missing for all samples.
* explicitly convert `x` to a `data.frame` in `df2gtypes` in case it is a `data.table` or `tibble`.

## 2.0.2

* NOTE: In order to speed up indexing the data in large data sets, this version changes the underlying structure of the `gtypes` object by replacing the `@loci` data.frame slot with a `@data` data.table slot. The data.table has a `id` character column, a `strata` character column, and every column afterwards represents one locus. The `@strata` slot has been removed.
* The `loci` accessor has been removed. 
* Added `as.array` which returns a 3-dimensional array with dimensions of [id, locus, allele].
* The print (show) function for `gtypes` objects no longer shows a by-locus summary. The display was getting too slow for data sets with a large number of loci.
* The `summary` function now includes by-sample results.
* Fixed computational errors in population structure metrics due to incorrect sorting of stratification.
* Added `maf` to return minimum allele frequency for each locus.
* Added `ldNe` to calculate Ne.
* Added `expandHaplotypes` to expand the haplotypes in a `gtypes` object to one sequence per individual.

## 1.0.6 

* Added `read.arlequin` back. Fixed missing function error with `write.arlequin`.
* Added `summarizeSamples`
* Changed `evanno` from base graphics to ggplot2
* Updated logic in `labelHaplotypes` to assign haplotypes if possible alternative site combinations match a present haplotype
* Added Zenodo DOI
* Added shiny app (`strataGUI`) for creating gtypes objects, QA/QC, and population structure analyses
* Added `type` argument to `structurePlot` to select between area and bar charts
* Changed `haplotypeLikelihoods` to `sequenceLikelihoods`
* `neiDa` now creates haplotypes before calculating metric
* Fixed error in `writePhase` that was creating improper input files for PHASE

## 1.0.5

* Fixed error in dupGenotypes, propSharedLoci, and propSharedIDs where missing genotypes were not being properly counted.
* Added as.data.frame.gtypes.
* Removed gtypes2df.
* Added arguments to as.matrix.gtypes to include id and strata columns in output.
* Removed the jmodeltest function as this functionality is available in the modeltest function in the phangorn package.
* Added conversion functions gtypes2phyDat and phyDat2gtypes to facilitate interoperability with the phangorn package.
* Removed read.arlequin.
* Added alleleNames accessor for gtypes object, which returns list of allele names for each locus.

## 1.0

* New version with different gtypes format from previous versions. See vignettes for instructions and examples.
