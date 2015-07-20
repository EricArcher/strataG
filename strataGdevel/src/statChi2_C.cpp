#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double statChi2_C(NumericMatrix loci, NumericVector strata) {
  // function declarations
  NumericMatrix table2D(NumericVector, NumericVector);
  NumericVector rowSumC(NumericMatrix);
  NumericVector colSumC(NumericMatrix);
  NumericMatrix outerC(NumericVector, NumericVector);
  
  int ncol(loci.ncol()), n, num_allele, num_strata;
  double chi2(0), exp_freq;
  NumericMatrix obsFreq, expFreq;
  for(int i = 0; i < ncol; i++) {
    obsFreq = table2D(loci(_, i), strata);
    num_allele = obsFreq.nrow();
    num_strata = obsFreq.ncol();
    n = 0;
    for(int r = 0; r < num_allele; r++)
      for(int c = 0; c < num_strata; c++) n += obsFreq(r, c);
    if(obsFreq.nrow() == 0 || obsFreq.ncol() < 2) continue;
    expFreq = outerC(rowSumC(obsFreq), colSumC(obsFreq));
    for(int r = 0; r < num_allele; r++) {
      for(int c = 0; c < num_strata; c++) {
        exp_freq = expFreq(r, c) / n;
        if(exp_freq > 0) chi2 += pow(obsFreq(r, c) - exp_freq, 2) / exp_freq;
      }
    }
  }
  return chi2;
}