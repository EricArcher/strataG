#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix table2D_int(const IntegerVector &x, const IntegerVector &y) {
  IntegerVector unique_x = sort_unique(x);
  IntegerVector unique_y = sort_unique(y);
  IntegerVector x_match = match(x, unique_x) - 1;
  IntegerVector y_match = match(y, unique_y) - 1;
  IntegerMatrix freq(unique_x.size(), unique_y.size());
  rownames(freq) = as<CharacterVector>(unique_x);
  colnames(freq) = as<CharacterVector>(unique_y);
  for(int i = 0; i < x_match.size(); i++) freq(x_match[i], y_match[i])++;
  return freq;
}

// [[Rcpp::export]]
IntegerVector makeLookupVec(const IntegerVector &x, const IntegerVector &y) {
  LogicalVector x_dup = duplicated(x);
  IntegerVector x_vec = x[!x_dup], y_vec = y[!x_dup];
  IntegerVector idx = seq_along(x_vec) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return x_vec[i] < x_vec[j];});
  // And return x in that order
  return y_vec[idx];
}

