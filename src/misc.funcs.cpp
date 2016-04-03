#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp:export]]
double harmonicMean_C(NumericVector x) {
  x = na_omit(x);
  double result;
  if(all(x > 0).is_true()) {
    double n(x.size());
    NumericVector invX(n);
    for(int i = 0; i < n; i++) invX[i] = 1 / x[i];
    result = n / sum(invX);
  } else {
    double invMeanX = 1 / mean(x);
    double varX = var(x);
    result = 1 / (invMeanX + varX * pow(invMeanX, 3));
  }
  if(std::isnan(result)) result = NA_REAL;
  return result;
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
  double n(len / ncol);
  IntegerMatrix mat(n, ncol);
  for(int i = 0; i < len; i++) {
    c = std::floor(double(i) / n);
    r = i - (int(n) * c);
    mat(r, c) = x[i];
  }
  return mat;
}

// [[Rcpp::export]]
NumericMatrix numVecToMat(NumericVector x, int ncol) {
  int len(x.size()), r, c;
  double n(len / ncol);
  NumericMatrix mat(n, ncol);
  for(int i = 0; i < len; i++) {
    c = std::floor(double(i) / n);
    r = i - (int(n) * c);
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