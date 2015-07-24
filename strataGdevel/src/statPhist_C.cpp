#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ssWPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq,
                        NumericMatrix hapDist) {
  // function declarations
  IntegerMatrix intOuterC(IntegerVector, IntegerVector);
  
  // Calculate sums of squares within strata (Eqn 8a)
  NumericVector ssWPvec(strataFreq.size());
  for(int i = 0; i < ssWPvec.size(); i++) {
    IntegerVector hapFreq = strataHapFreq(_, i);
    IntegerMatrix freqProd = intOuterC(hapFreq, hapFreq);
    for(int r = 0; r < freqProd.nrow(); r++) {
      for(int c = 0; c < freqProd.ncol(); c++) {
        ssWPvec[i] += freqProd(r, c) * hapDist(r, c);
      }
    }
    ssWPvec[i] /= (2 * strataFreq[i]);
  }
 return sum(ssWPvec); 
}

// [[Rcpp::export]]
double ssAPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq, 
                NumericMatrix hapDist) {
  // function declarations
  IntegerMatrix intOuterC(IntegerVector, IntegerVector);
  
  // Calculate sums of squares among strata (Eqn 8b)
  double ssAP(0);
  for(int i = 0; i < strataHapFreq.ncol(); i++) {
    for(int j = 0; j < strataHapFreq.ncol(); j++) {
      IntegerMatrix freqProd = intOuterC(strataHapFreq(_, i), strataHapFreq(_, j));
      for(int r = 0; r < freqProd.nrow(); r++) {
        for(int c = 0; c < freqProd.ncol(); c++) {
          ssAP += freqProd(r, c) * hapDist(r, c);
        }
      }
    }
  }
  return ssAP / (2 * sum(strataFreq));
}

// [[Rcpp::export]]
double statPhist_C(IntegerVector haps, IntegerVector strata, NumericMatrix hapDist) {
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector colSumC(NumericMatrix);
  
  LogicalVector hapsGood = !is_na(haps);
  LogicalVector strataGood = !is_na(strata);
  LogicalVector toUse = hapsGood & strataGood;
  haps = haps[toUse];
  strata = strata[toUse];
  
  // Extract summary values
  IntegerMatrix strataHapFreq = table2D(haps, strata);
  IntegerVector strataFreq = wrap(colSumC(wrap(strataHapFreq)));

  double ssWP = ssWPCalc(strataFreq, strataHapFreq, hapDist);
  double ssAP = ssAPCalc(strataFreq, strataHapFreq, hapDist);
  ssAP = ssAP - ssWP;
  
  // Calculate average sample size correction for among strata variance 
  //   Eqn 9a in paper, but modified as in Table 8.2.1.1 from Arlequin v3.5.1 manual
  //   (denominator is sum{I} - 1)
  int numSamples = sum(strataFreq);
  int numStrata = strataFreq.size();
  NumericVector n2(numStrata);
  for(int i = 0; i < n2.size(); i++) n2[i] = pow(strataFreq[i], 2) / numSamples;
  double n = (numSamples - sum(n2)) / (numStrata - 1);
  
  // Calculate variance components (Table 1)
  //   Set MSD (SSD / df) equal to expected MSD
  double Vc = ssWP / (numSamples - numStrata);
  double Vb = ((ssAP / (numStrata - 1)) - Vc) / n;
  double est(Vb / (Vb + Vc));
  
  if(isnan(est)) est = NA_REAL;
  return est;
}