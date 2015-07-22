#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double statChi2_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector rowSumC(NumericMatrix);
  NumericVector colSumC(NumericMatrix);
  NumericMatrix outerC(NumericVector, NumericVector);
  
  strata = rep(strata, ploidy);
  double chi2(0);
  for(int i = 0; i < loci.ncol(); i++) {
    IntegerMatrix obsFreq(table2D(loci(_, i), strata));
    int n = 0;
    for(int r = 0; r < obsFreq.nrow(); r++)
      for(int c = 0; c < obsFreq.ncol(); c++) n += obsFreq(r, c);
    if(obsFreq.nrow() == 0 || obsFreq.ncol() < 2) continue;
    NumericMatrix expFreq = outerC(rowSumC(wrap(obsFreq)), colSumC(wrap(obsFreq)));
    for(int r = 0; r < obsFreq.nrow(); r++) {
      for(int c = 0; c < obsFreq.ncol(); c++) {
        double exp_freq = expFreq(r, c) / n;
        if(exp_freq > 0) chi2 += pow(obsFreq(r, c) - exp_freq, 2) / exp_freq;
      }
    }
  }
  return chi2;
}