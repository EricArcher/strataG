#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


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


// [[Rcpp::export]]
double statChi2_C(NumericMatrix loci, NumericVector strata) {
  int ncol(loci.ncol()), n, num_allele, num_strata;
  double chi2(0), exp_freq;
  NumericMatrix obsFreq;
  NumericMatrix expFreq;
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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x1 <- sample(1:10, 100, rep = T)
x2 <- sample(1:5, 100, rep = T)
y <- table2vec(x1, x2)
z <- outerC(rowSumC(y), colSumC(y)) / length(x1) 


xfreq <- as.vector(table(x1))
yfreq <- as.vector(table(x2))
outerC(xfreq, yfreq)
outer(xfreq, yfreq)


library(strataG.devel)
example(gtypes)
msats <- stratify(msats, "fine")
x <- sapply(loci(msats), as.numeric)
s <- rep(as.numeric(strata(msats)), ploidy(msats))

statChi2_C(x, s)
statChi2(msats)

# lapply(1:ncol(x), function(i) {
#   obs <- table2vec(x[, i], s)
#   rs <- rowSums(obs)
#   cs <- colSums(obs)
#   exptd <- outer(rs, cs)
#   list(obs = obs, exptd = exptd)
# })


dloop <- stratify(dloop, "fine")
x <- sapply(loci(dloop), as.numeric)
s <- rep(as.numeric(strata(dloop)), ploidy(dloop))

statChi2_C(x, s)
statChi2(dloop)


*/
