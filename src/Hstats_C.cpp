#include <Rcpp.h>
using namespace Rcpp;

// Nei, M. and R.K. Chesser. 1983. Estimation of fixation indices and gene diversities. Ann. Hum. Genet. 47:253-259.

// [[Rcpp::export]]
double HoCalc(int nInd, IntegerVector locus, int ploidy, 
              IntegerVector strata, IntegerVector strataN) {
  // function declarations
  IntegerVector idGenotype(IntegerVector, int, int);
  int getMaxInt(IntegerVector);
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  
  // Estimate Ho (frequency of all heterozygotes): Equation 5, page 254
  IntegerVector homGenotype(strata.size());
  int numAlleles(getMaxInt(locus));
  for(int i = 0; i < strata.size(); i++) {
    IntegerVector gt = idGenotype(locus, i, ploidy);
    bool NAs = false;
    for(int j = 0; j < gt.size(); j++) {
      if(IntegerVector::is_na(gt[j])) {
        NAs = true;
        continue;
      }
    }
    if(NAs) { // set to NA if genotype is missing
      homGenotype[i] = NA_INTEGER;
    } else {
      gt = unique(gt);
      if(gt.size() == 1) { // identify allele individual is homozygous for
        homGenotype[i] = gt[0];
      } else { // set to maximum allele value + 1
        homGenotype[i] = numAlleles;
      }
    }
  }
  
  // compute homozygote frequency in each population
  NumericMatrix homHetFreq = wrap(table2D(strata, homGenotype));
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
  // function declarations
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
  // function declarations
  NumericVector colMeanC(NumericMatrix);
  
  // Estimate Ht (expected heterozygosity overall): Equation 11, page 256 
  double meanHet = 1 - sum(pow(colMeanC(wrap(alleleFreq)), 2));
  double harmNs = harmN * (sum(strataN * ploidy));
  return meanHet + (Hs / harmNs) - (Ho / 2 / harmNs);
}

// [[Rcpp::export]]
NumericMatrix Hstats_C(IntegerMatrix loci, IntegerVector strata) {
  // function declarations
  IntegerVector calcStrataN(IntegerVector, IntegerVector, int);
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  double harmonicMean_C(NumericVector);
  
  int ploidy(loci.nrow() / strata.size());
  int nInd(loci.nrow() / ploidy);
  NumericMatrix hstats(3, loci.ncol());
  for(int i = 0; i < loci.ncol(); i++) {
    IntegerVector strataN = calcStrataN(loci(_, i), strata, ploidy);
    hstats(0, i) = HoCalc(nInd, loci(_, i), ploidy, strata, strataN);
    NumericMatrix alleleFreq = wrap(table2D(rep_each(strata, ploidy), loci(_, i)));

    double harmN;
    if(strataN.size() > 1) {
      harmN = harmonicMean_C(wrap(strataN));
    } else {
      harmN = strataN[0];
    }

    hstats(1, i) = HsCalc(alleleFreq, ploidy, strataN, harmN, hstats(0, i));
    hstats(2, i) = HtCalc(alleleFreq, ploidy, strataN, harmN, hstats(0, i), hstats(1, i));
  }
  return hstats;
}