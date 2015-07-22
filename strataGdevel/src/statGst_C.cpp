#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double statGst_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericMatrix hets(Hstats_C(loci, strata, ploidy));
  double Hs(mean(hets(1, _))), Ht(mean(hets(2, _)));
  double est(1 - Hs / Ht);
  if(isnan(est)) est = NA_REAL;
  return est;
}
  
// [[Rcpp::export]]
double statGstPrime_C(IntegerMatrix loci, IntegerVector strata, int ploidy, int primeType) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericMatrix hets(Hstats_C(loci, strata, ploidy));
  double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), gstMax, est;
  int k(unique(strata).size());
  if(primeType == 0) {  // primeType: 0 = Nei's, 1 = Hedrick's
    est = (k * (Ht - Hs)) / ((k * Ht) - Hs);
  } else {
    gstMax = ((k - 1) * (1 - Hs)) / (k - 1 + Hs);
    est = (1 - (Hs / Ht)) / gstMax;
  }
  if(isnan(est)) est = NA_REAL;
  return est;
}

// [[Rcpp::export]]
double statGstDblPrime_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericMatrix hets(Hstats_C(loci, strata, ploidy));
  double Hs(mean(hets(1, _))), Ht(mean(hets(2, _))), est;
  int k(unique(strata).size());
  est = (k * (Ht - Hs)) / ((k * Ht) - Hs) / (1 - Hs);
  if(isnan(est)) est = NA_REAL;
  return est;
}