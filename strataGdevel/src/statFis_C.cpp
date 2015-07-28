#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double statFis_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  NumericMatrix Hstats_C(IntegerMatrix, IntegerVector, int);
  
  NumericMatrix hets(Hstats_C(loci, strata, ploidy));
  double Ho(mean(hets(0, _))), Hs(mean(hets(1, _)));
  double est((Hs  - Ho) / Hs);
  if(std::isnan(est)) est = NA_REAL;
  return est;
}