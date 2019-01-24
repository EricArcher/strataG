#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector statGst_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx)));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _)));
    double est(1 - Hs / Ht);
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}
  
// [[Rcpp::export]]
NumericVector statGstPrime_C(IntegerMatrix loci, IntegerMatrix strataMat, int primeType) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector);
  
  NumericVector estVec(strataMat.ncol());
  int k(unique(strataMat(_, 0)).size());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx)));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), gstMax, est;
    if(primeType == 0) {  // primeType: 0 = Nei's, 1 = Hedrick's
      est = (k * (Ht - Hs)) / ((k * Ht) - Hs);
    } else {
      gstMax = ((k - 1) * (1 - Hs)) / (k - 1 + Hs);
      est = (1 - (Hs / Ht)) / gstMax;
    }
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}

// [[Rcpp::export]]
NumericVector statGstDblPrime_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector);
  
  NumericVector estVec(strataMat.ncol());
  int k(unique(strataMat(_, 0)).size());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx)));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), est;
    est = (k * (Ht - Hs)) / ((k * Ht) - Hs) / (1 - Hs);
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}