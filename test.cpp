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

double statJostD_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  IntegerVector calcStrataN(IntegerVector, IntegerVector);
  double harmonicMean_C(NumericVector);
  
  NumericVector est(loci.ncol());
  for(int i = 0; i < loci.ncol(); i++) {
    IntegerMatrix alleleStrataFreq = table2D(loci(_, i), rep(strata, ploidy));
    int numStrata = calcStrataN(loci, strata).size();
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
    est[i] = 1 - sum(iTerms(0, _)) / sum(iTerms(1, _));
  }
  return harmonicMean_C(est);
}

// terms <- sapply(1:ncol(g@loci), function(i) {
//   allele.strata.freq <- table(g@loci[, i], strata, useNA = "no")
//   num.strata <- ncol(allele.strata.freq)
//   i.terms <- sapply(rownames(allele.strata.freq), function(allele) {
//     j.terms <- sapply(colnames(allele.strata.freq), function(strata) {
//       Nj <- sum(allele.strata.freq[, strata])
//       Nij <- allele.strata.freq[allele, strata]
//       a.term1 <- Nij / Nj
//       a.term2 <- a.term1 ^ 2
//       b.term <- Nij * (Nij - 1) / (Nj * (Nj - 1))
//       c(a.term1 = a.term1, a.term2 = a.term2, b.term = b.term)
//     })
//     a.term1 <- sum(j.terms["a.term1", ]) ^ 2
//     a.term2 <- sum(j.terms["a.term2", ])
//     c(a = (a.term1 - a.term2) / (num.strata - 1), 
//       b = sum(j.terms["b.term", ])
//     )
//   })
//   
//   c(a = sum(i.terms["a", ]), b = sum(i.terms["b", ]))
// })
//   est <- 1 - terms["a", ] / terms["b", ]
// 
// c(D = harmonic.mean(est))


/*** R
library(strataGdevel)
example(gtypes)
msats <- stratify(msats, "broad")

statJostD(msats)
statJostD_C(
  sapply(loci(msats), function(x) as.numeric(x) - 1), 
  as.numeric(strata(msats)) - 1,
  ploidy(msats)
)


*/