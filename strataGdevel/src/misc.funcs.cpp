#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int getMaxInt(NumericVector x) {
  int m(x[0]);
  int len(x.size());
  for(int i = 1; i < len; i++) if(x[i] > m) m = x[i];
  return m;
}

// [[Rcpp::export]]
NumericMatrix table2D(NumericVector x, NumericVector y) {
  int len(x.size());
  int x_max(getMaxInt(x));
  int y_max(getMaxInt(y));
  NumericMatrix tbl(x_max, y_max);
  for(int i = 0; i < len; i++) {
    if(NumericVector::is_na(x[i]) || NumericVector::is_na(y[i])) continue;
    tbl(x[i] - 1, y[i] - 1)++;
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
NumericMatrix outerC(NumericVector x, NumericVector y) {
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