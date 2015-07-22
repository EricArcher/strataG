#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector indGenotype(int nInd, int numAlleles, IntegerMatrix locus) {
  //   record which homozygote or heterozygote (last column) for each individual
  IntegerVector genotype(nInd);
  for(int j = 0; j < nInd; j++) {
    IntegerVector alleles = locus(j, _);
    if(any(is_na(alleles))) {
      genotype[j] = NA_INTEGER;
      continue;
    }
    alleles = unique(alleles);
    if(alleles.size() == 1) {
      genotype[j] = alleles[0];
    } else {
      genotype[j] = numAlleles;
    }
  }
  return genotype;
}

// [[Rcpp::export]]
double HoCalc(int nInd, IntegerVector loci, int ploidy, IntegerVector strata,
              IntegerVector strataN) {
  IntegerMatrix intVecToMat(IntegerVector, int);
  int getMaxInt(IntegerVector);
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector rowSumC(NumericMatrix);
  
  // Estimate Ho (frequency of all heterozygotes): Equation 5, page 254
  //   translate locus column into matrix
  IntegerMatrix locus = intVecToMat(loci, ploidy);
  IntegerVector genotype = indGenotype(nInd, getMaxInt(loci), locus);
  //   compute homozygote frequency in each population
  NumericMatrix homHetFreq = wrap(table2D(strata, genotype));
  for(int r = 0; r < homHetFreq.nrow(); r++) {
    homHetFreq(r, _) = homHetFreq(r, _) / strataN[r];
  }
  NumericVector meanHomFreq(homHetFreq.ncol() - 1);
  for(int c = 0; c < homHetFreq.ncol() - 1; c++) {
    meanHomFreq[c] += mean(homHetFreq(_, c));
  }
  return 1 - sum(meanHomFreq);
}

// [[Rcpp::export]]
double HsCalc(NumericMatrix alleleFreq, int ploidy, IntegerVector strataN, 
              double harmN, double Ho) {    
  NumericVector rowSumC(NumericMatrix);
  
  // Estimate Hs (expected heterozygosity within strata): Equation 9, page 255
  //   compute homoz frequencies within strata
  for(int r = 0; r < alleleFreq.nrow(); r++) {
    alleleFreq(r, _) = alleleFreq(r, _) / (strataN[r] * ploidy);
  }
  NumericMatrix homFreq(alleleFreq.nrow(), alleleFreq.ncol());
  for(int r = 0; r < homFreq.nrow(); r++) homFreq(r, _) = pow(alleleFreq(r, _), 2);
  double meanHet = mean(1 - rowSumC(homFreq));
  return (harmN / (harmN - 1)) * (meanHet - (Ho / 2 / harmN));
}

double HtCalc(NumericMatrix alleleFreq, int ploidy, IntegerVector strataN,
              double harmN, double Ho, double Hs) {
  NumericVector colMeanC(NumericMatrix);
  // Estimate Ht (expected heterozygosity overall): Equation 11, page 256 
  double meanHet = 1 - sum(pow(colMeanC(wrap(alleleFreq)), 2));
  double harmNs = harmN * (sum(strataN * ploidy));
  return meanHet + (Hs / harmNs) - (Ho / 2 / harmNs);
}

// [[Rcpp::export]]
NumericMatrix Hstats_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  IntegerVector calcStrataN(IntegerVector, IntegerVector);
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  
  int nInd(loci.nrow() / ploidy);
  NumericMatrix hstats(3, loci.ncol());
  for(int i = 0; i < loci.ncol(); i++) {
    IntegerVector strataN = calcStrataN(loci(_, i), strata);
    hstats(0, i) = HoCalc(nInd, loci(_, i), ploidy, strata, strataN);
    NumericMatrix alleleFreq = wrap(table2D(rep(strata, ploidy), loci(_, i)));
    NumericVector invN(strataN.size());
    for(int j = 0; j < invN.size(); j++) invN[j] = 1 / double(strataN[j]);
    double harmN = strataN.size() / sum(invN);
    hstats(1, i) = HsCalc(alleleFreq, ploidy, strataN, harmN, hstats(0, i));
    hstats(2, i) = HtCalc(alleleFreq, ploidy, strataN, harmN, hstats(0, i), hstats(1, i));
  }
  return hstats;
}