#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector statGst_C(IntegerMatrix loci, IntegerMatrix strataMat, int ploidy) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx), ploidy));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _)));
    double est(1 - Hs / Ht);
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}
  
// [[Rcpp::export]]
NumericVector statGstPrime_C(IntegerMatrix loci, IntegerMatrix strataMat, int ploidy, int primeType) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx), ploidy));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), gstMax, est;
    int k(unique(strataMat(_, idx)).size());
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
NumericVector statGstDblPrime_C(IntegerMatrix loci, IntegerMatrix strataMat, int ploidy) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx), ploidy));
    double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), est;
    int k(unique(strataMat(_, idx)).size());
    est = (k * (Ht - Hs)) / ((k * Ht) - Hs) / (1 - Hs);
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}