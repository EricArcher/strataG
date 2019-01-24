#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector statFis_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector);
  
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    NumericMatrix hets(Hstats_C(loci, strataMat(_, idx)));
    double Ho(mean(hets(0, _))), Hs(mean(hets(1, _)));
    double est((Hs  - Ho) / Hs);
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}