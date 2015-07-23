#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp:export]]
double harmonicMean_C(NumericVector x) {
  NumericVector invN(x.size());
  for(int j = 0; j < invN.size(); j++) invN[j] = 1 / x[j];
  return x.size() / sum(invN);
}

// [[Rcpp::export]]
int getMaxInt(IntegerVector x) {
  x = na_omit(x);
  int m(x[0]);
  for(int i = 1; i < x.size(); i++) if(x[i] > m) m = x[i];
  return m;
}

// [[Rcpp::export]]
IntegerMatrix table2D(IntegerVector x, IntegerVector y) {
  int x_max(getMaxInt(x) + 1), y_max(getMaxInt(y) + 1);
  IntegerMatrix tbl(x_max, y_max);
  for(int i = 0; i < x.size(); i++) {
    if(IntegerVector::is_na(x[i]) || IntegerVector::is_na(y[i])) continue;
    tbl(x[i], y[i])++;
  }
  return tbl;
}

// [[Rcpp::export]]
NumericVector rowSumC(NumericMatrix x) {
  double nrow(x.nrow());
  NumericVector rowsum(nrow);
  for(int i = 0; i < nrow; i++) rowsum[i] = sum(x(i, _));
  return rowsum;
}

// [[Rcpp::export]]
NumericVector colSumC(NumericMatrix x) {
  double ncol(x.ncol());
  NumericVector colsum(ncol);
  for(int i = 0; i < ncol; i++) colsum[i] = sum(x(_, i));
  return colsum;
}

// [[Rcpp::export]]
NumericVector colMeanC(NumericMatrix x) {
  double ncol(x.ncol());
  NumericVector colmean(ncol);
  for(int i = 0; i < ncol; i++) colmean[i] = mean(x(_, i));
  return colmean;
}

// [[Rcpp::export]]
IntegerMatrix intOuterC(IntegerVector x, IntegerVector y) {
  int x_len(x.size());
  int y_len(y.size());
  IntegerMatrix tbl(x_len, y_len);
  for(int i = 0; i < x_len; i++) {
    for(int j = 0; j < y_len; j++) {
      if(IntegerVector::is_na(x[i]) || IntegerVector::is_na(y[j])) {
        tbl(i, j) = NA_INTEGER;
      } else {
        tbl(i, j) = x[i] * y[j]; 
      }
    }
  }
  return tbl;
}

// [[Rcpp::export]]
NumericMatrix numOuterC(NumericVector x, NumericVector y) {
  int x_len(x.size());
  int y_len(y.size());
  NumericMatrix tbl(x_len, y_len);
  for(int i = 0; i < x_len; i++) {
    for(int j = 0; j < y_len; j++) {
      if(NumericVector::is_na(x[i]) || NumericVector::is_na(y[j])) {
        tbl(i, j) = NA_REAL;
      } else {
        tbl(i, j) = x[i] * y[j]; 
      }
    }
  }
  return tbl;
}

// [[Rcpp::export]]
IntegerMatrix intVecToMat(IntegerVector x, int ncol) {
  int len(x.size()), r, c;
  int n(len / ncol);
  IntegerMatrix mat(n, ncol);
  for(int i = 0; i < len; i++) {
    c = floor(i / n);
    r = i - (n * c);
    mat(r, c) = x[i];
  }
  return mat;
}

// [[Rcpp::export]]
NumericMatrix numVecToMat(NumericVector x, int ncol) {
  int len(x.size()), r, c;
  int n(len / ncol);
  NumericMatrix mat(n, ncol);
  for(int i = 0; i < len; i++) {
    c = floor(i / n);
    r = i - (n * c);
    mat(r, c) = x[i];
  }
  return mat;
}

// [[Rcpp::export]]
IntegerVector calcStrataN(IntegerVector locus, IntegerVector strata) {
  IntegerVector n(getMaxInt(strata) + 1);
  for(int i = 0; i < strata.size(); i++) {
    bool allelesNA(IntegerVector::is_na(locus[i]));
    bool strataNA(IntegerVector::is_na(strata[i]));
    if(allelesNA || strataNA) continue;
    n[strata[i]]++;
  }
  return n;
}




// [[Rcpp::export]]
double ssWPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq,
                        NumericMatrix hapDist) {
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
  IntegerMatrix intOuterC(IntegerVector, IntegerVector);
  
  // Calculate sums of squares among strata (Eqn 8b)
  double ssAP(0);
  for(int i = 0; i < strataHapFreq.ncol(); i++) {
    for(int j = 0; j < strataHapFreq.ncol(); j++) {
      IntegerMatrix freqProd = intOuterC(strataHapFreq(_, i), strataHapFreq(_, j));
      for(int r = 0; r < freqProd.nrow(); r++) {
        for(int c = 0; c < freqProd.ncol(); c++) {
          ssAP += freqProd(r, c) * hapDist(r, c);
          cout << freqProd(r, c) << " ";
        }
        cout << endl;
      }
      cout << endl;
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
  return Vb / (Vb + Vc);
}




/*** R
library(strataGdevel)
library(ape)
example(gtypes)
dloop <- stratify(dloop, "fine")
hap.dist <- dist.dna(getSequences(sequences(dloop))[levels(loci(dloop)[, 1])], model = "K80", 
                     gamma = FALSE, pairwise.deletion = TRUE, as.matrix = TRUE)


statPhist(dloop)

l <- as.numeric(loci(dloop)[, 1]) - 1
s <- as.numeric(strata(dloop)) - 1
statPhist_C(l, s, hap.dist)

*/