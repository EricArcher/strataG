#include <Rcpp.h>
using namespace Rcpp;

double ssWPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq,
                        NumericMatrix hapDist) {
  // function declarations
  IntegerMatrix intOuterC(IntegerVector, IntegerVector);
  
  // Calculate sums of squares within strata (Eqn 8a)
  NumericVector ssWPvec(strataFreq.size());
  double d;
  for(int i = 0; i < ssWPvec.size(); i++) {
    IntegerVector hapFreq = strataHapFreq(_, i);
    IntegerMatrix freqProd = intOuterC(hapFreq, hapFreq);
    for(int r = 0; r < freqProd.nrow(); r++) {
      for(int c = 0; c < freqProd.ncol(); c++) {
        d = hapDist(r, c);
        if(std::isnan(d)) return NA_REAL;
        ssWPvec[i] += freqProd(r, c) * d;
      }
    }
    ssWPvec[i] /= (2 * strataFreq[i]);
  }
 return sum(ssWPvec); 
}

double ssAPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq, 
                NumericMatrix hapDist) {
  // function declarations
  IntegerMatrix intOuterC(IntegerVector, IntegerVector);
  
  // Calculate sums of squares among strata (Eqn 8b)
  double ssAP(0);
  double d;
  for(int i = 0; i < strataHapFreq.ncol(); i++) {
    for(int j = 0; j < strataHapFreq.ncol(); j++) {
      IntegerMatrix freqProd = intOuterC(strataHapFreq(_, i), strataHapFreq(_, j));
      for(int r = 0; r < freqProd.nrow(); r++) {
        for(int c = 0; c < freqProd.ncol(); c++) {
          d = hapDist(r, c);
          if(std::isnan(d)) return NA_REAL;
          ssAP += freqProd(r, c) * d;
        }
      }
    }
  }
  return ssAP / (2 * sum(strataFreq));
}

double phistCalc(IntegerVector haps, IntegerVector strata, NumericMatrix hapDist) {  
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector colSumC(NumericMatrix);
  
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
  for(int i = 0; i < n2.size(); i++) n2[i] = std::pow(strataFreq[i], 2.0) / numSamples;
  double n = (numSamples - sum(n2)) / (numStrata - 1);
  
  // Calculate variance components (Table 1)
  //   Set MSD (SSD / df) equal to expected MSD
  double Vc = ssWP / (numSamples - numStrata);
  double Vb = ((ssAP / (numStrata - 1)) - Vc) / n;
  double est(Vb / (Vb + Vc));
  
  if(std::isnan(est)) est = NA_REAL;
  return est;
}

// [[Rcpp::export]]
NumericVector statPhist_C(IntegerMatrix hapMat, IntegerMatrix strataMat, List hapDist) {
  // function declarations
  double harmonicMean_C(NumericVector);
  
  LogicalMatrix hapsGood(hapMat.nrow(), hapMat.ncol());
  for(int gene = 0; gene < hapMat.ncol(); gene++) hapsGood(_, gene) = !is_na(hapMat(_, gene));
  
  NumericVector estVec(strataMat.ncol());
  NumericVector geneVec(hapMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    IntegerVector strata(strataMat(_, idx));
    for(int gene = 0; gene < hapMat.ncol(); gene++) {
      LogicalVector toUse = hapsGood(_, gene);
      IntegerVector haps = hapMat(_, gene);
      geneVec[gene] = phistCalc(haps[toUse], strata[toUse], hapDist[gene]);
    }
    if(geneVec.size() > 1) {
      estVec[idx] = harmonicMean_C(geneVec);
    } else {
      estVec[idx] = geneVec[0];
    }
  }
  return estVec;
}
