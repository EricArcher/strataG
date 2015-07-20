#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Hstats_C(NumericMatrix loci, NumericVector strata, int ploidy) {
  // function declarations
  int getMaxInt(NumericVector);
  NumericMatrix table2D(NumericVector, NumericVector);
  NumericVector rowSumC(NumericMatrix);
  NumericVector colMeanC(NumericMatrix);
  
  int nLoc(loci.ncol()), nRow(loci.nrow());
  int nInd(nRow / ploidy), r, c, numAlleles;
  double harmN, meanHet, harmNs;
  NumericVector alleles, genotype, strataN, meanHomFreq;
  NumericMatrix hstats(3, nLoc), locus(nInd, ploidy), homFreq, alleleFreq;
  for(int i = 0; i < nLoc; i++) {
    // Estimate Ho (frequency of all heterozygotes): Equation 5, page 254
    //   translate locus column into matrix
    for(int j = 0; j < nRow; j++) {
      c = floor(j / nInd);
      r = j - (nInd * c);
      locus(r, c) = loci(j, i);
    }
    //   record which homozygote or heterozygote (last column) for each individual
    numAlleles = getMaxInt(loci(_, i));
    genotype = NumericVector(nInd);
    for(int j = 0; j < nInd; j++) {
      alleles = locus(j, _);
      if(any(is_na(alleles))) {
        genotype[j] = NA_REAL;
        continue;
      }
      alleles = unique(alleles);
      if(alleles.size() == 1) {
        genotype[j] = alleles[0];
      } else {
        genotype[j] = numAlleles;
      }
    }
    //   compute homozygote frequency in each population
    homFreq = table2D(strata, genotype);
    strataN = rowSumC(homFreq);
    for(r = 0; r < homFreq.nrow(); r++) {
      homFreq(r, _) = homFreq(r, _) / strataN(r);
    }
    meanHomFreq = NumericVector(homFreq.ncol() - 1);
    for(c = 0; c < homFreq.ncol() - 1; c++) {
      meanHomFreq[c] += mean(homFreq(_, c));
    }
    hstats(0, i) = 1 - sum(meanHomFreq);
    
    // Estimate Hs (expected heterozygosity within strata): Equation 9, page 255
    //   compute homoz frequencies within strata
    alleleFreq = table2D(rep(strata, ploidy), loci(_, i));
    for(r = 0; r < alleleFreq.nrow(); r++) {
      alleleFreq(r, _) = alleleFreq(r, _) / (strataN(r) * ploidy);
    }
    homFreq = NumericMatrix(alleleFreq.nrow(), alleleFreq.ncol());
    for(r = 0; r < homFreq.nrow(); r++) homFreq(r, _) = pow(alleleFreq(r, _), 2);
    meanHet = mean(1 - rowSumC(homFreq));
    harmN = strataN.size() / sum(1 / strataN);
    hstats(1, i) = (harmN / (harmN - 1)) * (meanHet - (hstats(0, i) / 2 / harmN));
    
    // Estimate Ht (expected heterozygosity overall): Equation 11, page 256 
    meanHet = 1 - sum(pow(colMeanC(alleleFreq), 2));
    harmNs = harmN * (sum(strataN * ploidy));
    hstats(2, i) = meanHet + (hstats(1, i) / harmNs) - (hstats(0, i) / 2 / harmNs);
  }
  return hstats;
}