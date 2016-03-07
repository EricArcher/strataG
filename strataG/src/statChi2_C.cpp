#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector statChi2_C(IntegerMatrix loci, IntegerMatrix strataMat, int ploidy) {
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector rowSumC(NumericMatrix);
  NumericVector colSumC(NumericMatrix);
  NumericMatrix numOuterC(NumericVector, NumericVector);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    IntegerVector st(rep(strataMat(_, idx), ploidy));
    double chi2(0);
    for(int i = 0; i < loci.ncol(); i++) {
      IntegerMatrix obsFreq(table2D(loci(_, i), st));
      int n = 0;
      for(int r = 0; r < obsFreq.nrow(); r++)
        for(int c = 0; c < obsFreq.ncol(); c++) n += obsFreq(r, c);
      if(obsFreq.nrow() == 0 || obsFreq.ncol() < 2) continue;
      NumericMatrix expFreq = numOuterC(rowSumC(wrap(obsFreq)), colSumC(wrap(obsFreq)));
      for(int r = 0; r < obsFreq.nrow(); r++) {
        for(int c = 0; c < obsFreq.ncol(); c++) {
          double exp_freq = expFreq(r, c) / n;
          if(exp_freq > 0) chi2 += pow(obsFreq(r, c) - exp_freq, 2) / exp_freq;
        }
      }
    }
    
    if(std::isnan(chi2)) chi2 = NA_REAL;
    estVec[idx] = chi2;
  }
  
  return estVec;
}