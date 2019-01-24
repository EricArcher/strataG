#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector statJostD_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  IntegerVector calcStrataN(IntegerVector, IntegerVector, int);
  double harmonicMean_C(NumericVector);
  
  NumericVector terms(loci.ncol());
  NumericVector estVec(strataMat.ncol());
  int ploidy(loci.nrow() / strataMat.nrow());
  int numStrata(unique(strataMat(_, 0)).size());
  for(int idx = 0; idx < estVec.size(); idx++) {
    for(int i = 0; i < loci.ncol(); i++) {
      IntegerMatrix alleleStrataFreq = table2D(loci(_, i), rep_each(strataMat(_, idx), ploidy));
      NumericMatrix iTerms(2, alleleStrataFreq.nrow());
      for(int a = 0; a < iTerms.ncol(); a++) {
        NumericMatrix jTerms(3, alleleStrataFreq.ncol());
        for(int s = 0; s < jTerms.ncol(); s++) {
          double Nj = sum(alleleStrataFreq(_, s));
          double Nij = alleleStrataFreq(a, s);
          jTerms(0, s) = Nij / Nj;
          jTerms(1, s) = pow(jTerms(0, s), 2);
          jTerms(2, s) = Nij * (Nij - 1) / (Nj * (Nj - 1));
        }
        double aTerm1 = pow(sum(jTerms(0, _)), 2);
        double aTerm2 = sum(jTerms(1, _));
        iTerms(0, a) = (aTerm1 - aTerm2) / (numStrata - 1);
        iTerms(1, a) = sum(jTerms(2, _));
      }
      terms[i] = 1 - sum(iTerms(0, _)) / sum(iTerms(1, _));
    }
    
    for(int i = 0; i < terms.size(); i++) {
      if(terms[i] < 0) terms[i] = 0;
    }
    
    double est;
    if(terms.size() > 1) {
      est = harmonicMean_C(terms);
    } else {
      est = terms[0];
    }
    if(std::isnan(est)) est = NA_REAL;
    estVec[idx] = est;
  }
  return estVec;
}