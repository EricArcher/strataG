#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int getMaxInt(IntegerVector x) {
  x = na_omit(x);
  int m(x[0]);
  for(int i = 1; i < x.size(); i++) if(x[i] > m) m = x[i];
  return m;
}

// [[Rcpp::export]]
NumericMatrix table2D(IntegerVector x, IntegerVector y) {
  int x_max(getMaxInt(x) + 1), y_max(getMaxInt(y) + 1);
  NumericMatrix tbl(x_max, y_max);
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

// // [[Rcpp::export]]
// NumericMatrix numSubMat(NumericMatrix X, LogicalVector T) {
//   NumericVector L(T.size());
//   for(int i = 0; i < T.size(); i++) L[i] = int(T[i]);
//   mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
//   colvec tIdx(L.begin(), T.size(), false); 
//   mat y = Xmat.rows(find(tIdx == 1));
//   return wrap(y);
// }
// 
// 
// // [[Rcpp::export]]
// IntegerMatrix intSubMat(IntegerMatrix X, LogicalVector T) {
//   NumericMatrix Y(X.nrow(), X.ncol());
//   for(int r = 0; r < X.nrow(); r++) {
//     for(int c = 0; c < X.ncol(); c++) Y(r, c) = X(r, c);
//   }
//   return wrap(numSubMat(Y, T));
// }

// [[Rcpp::export]]
NumericMatrix alleleFreqCalc(IntegerVector locVec, IntegerVector strata,
                             int ploidy) {
  // get allele frequencies in each population
  NumericMatrix alleleFreq = table2D(locVec, rep(strata, ploidy));
  NumericVector alleleFreqColSums = colSumC(alleleFreq);
  for(int r = 0; r < alleleFreq.nrow(); r++) {
    for(int c = 0; c < alleleFreq.ncol(); c++) {
      alleleFreq(r, c) = alleleFreq(r, c) / alleleFreqColSums[c];
    }
  }
  return alleleFreq;
}

// [[Rcpp::export]]
NumericMatrix prHetCalc(IntegerVector alleles, IntegerVector nvec,
                        IntegerMatrix locusMat, IntegerVector strata,
                        int ploidy) {
  // get proportion heterozygosity for each allele (rows) and population (cols)
  NumericMatrix prHet(alleles.size(), nvec.size());
  IntegerVector n(nvec.size());
  for(int r = 0; r < locusMat.nrow(); r++) {
    bool allelesNA(any(is_na(locusMat(r, _))));
    bool strataNA(IntegerVector::is_na(strata[r]));
    if(allelesNA || strataNA) continue;
    n[strata[r]] += 1;
    if(unique(locusMat(r, _)).size() == 1) continue;
    for(int c = 0; c < locusMat.ncol(); c++) {
      prHet(locusMat(r, c), strata[r]) += 1;
    }
  }
  for(int c = 0; c < prHet.ncol(); c++) prHet(_, c) = prHet(_, c) / n[c];
  return prHet;
}

// [[Rcpp::export]]
NumericMatrix varCompCalc(IntegerVector nvec, NumericMatrix alleleFreq,
                          NumericMatrix prHet, int r, double nbar, 
                          double rnbar, double nc) {
  // create matrix of Va, Vb, Vc for all alleles
  // rows are variance components (1=Va, 2=Vb, 3=Vc), columns are alleles
  NumericMatrix varcompMat(3, alleleFreq.nrow());
  for(int i = 0; i < alleleFreq.nrow(); i++) {
    double pbar(0), hbar(0);
    NumericVector s(nvec.size());
    for(int k = 0; k < nvec.size(); k++) {
      pbar += nvec[k] * alleleFreq(i, k);
      hbar += nvec[k] * prHet(i, k);
    }
    pbar = pbar / rnbar;
    hbar = hbar / rnbar;
    for(int k = 0; k < nvec.size(); k++) {
      // s =  numerator inside summation for s2 calculation
      s[k] = nvec[k] * pow((alleleFreq(i, k) - pbar), 2);
    }
    double s2 = sum(s) / (r - 1) / nbar;
    
    // used in equations for Va and Vb
    double innerTerm = (pbar * (1 - pbar)) - ((r - 1) * s2 / r);
    double nbarM1 = nbar - 1;
    
    // Va (between strata) - Eqn. 2
    double Va1 = 0.25 * hbar;
    double Va2 = s2 - (innerTerm - Va1) / nbarM1;
    varcompMat(0, i) = nbar * Va2 / nc;
    
    // Vb (between individuals within strata) - Eqn. 3
    double Vb1 = ((2 * nbar - 1) * hbar) / 4 / nbar;
    double Vb2 = innerTerm - Vb1;
    varcompMat(1, i) = nbar * Vb2 / nbarM1;
    
    // Vc (between gametes within individuals) - Eqn. 4
    varcompMat(2, i) = 0.5 * hbar;
  }
  cout << endl;
  return varcompMat;
}


// [[Rcpp::export]]
double statFst_C(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // returns numerator and denominator sums for
  //   theta-w calculation for alleles at each locus (page 1363)
  
  int nInd(loci.nrow() / ploidy);
  LogicalVector toUse(nInd);
  NumericMatrix locusSums(2, loci.ncol());
  for(int i = 0; i < strata.size(); i++) {
    if(!IntegerVector::is_na(strata[i])) strata[i]--;
  }
  
  for(int loc = 0; loc < loci.ncol(); loc++) {
    IntegerVector locVec = loci(_, loc);
    for(int i = 0; i < locVec.size(); i++) {
      if(!IntegerVector::is_na(locVec[i])) locVec[i]--; 
    }
    
    // identify unique alleles and return zeroes if locus is fixed
    IntegerVector alleles = unique(na_omit(locVec));
    if(alleles.size() < 2) continue;
    
    // form matrix of alleles for each individual
    IntegerMatrix locusMat = intVecToMat(locVec, ploidy);

    // calculate variables constant for locus
    IntegerVector nvec(unique(strata).size());
    for(int r = 0; r < locusMat.nrow(); r++) {
      bool allelesNA(any(is_na(locusMat(r, _))));
      bool strataNA(IntegerVector::is_na(strata[r]));
      if(allelesNA || strataNA) continue;
      nvec[strata[r]]++;
    }
    int r = nvec.size();
    if(r < 2) continue;
    double nbar = mean(nvec);
    double rnbar = r * nbar;
    double nc = (rnbar - (sum(pow(nvec, 2)) / rnbar)) / (r - 1);
    
    NumericMatrix alleleFreq = alleleFreqCalc(locVec, strata, ploidy);
    
    NumericMatrix prHet = prHetCalc(alleles, nvec, locusMat, strata, ploidy);
    NumericMatrix varcompMat = varCompCalc(nvec, alleleFreq, prHet, r, nbar, rnbar, nc);
    locusSums(0, loc) = sum(varcompMat(0, _));
    for(int i = 0; i < varcompMat.nrow(); i++) {
      locusSums(1, loc) += sum(varcompMat(i, _));
    }
  }
  
  return sum(locusSums(0, _)) / sum(locusSums(1, _));
}


// #' @rdname popStructStat
// #' @export
// #' 
// statFstPrime <- function(g, strata = NULL, ...) {  
//   if(ploidy(g) == 1) return(c('F\'st' = NA))
//   
//   strata.vec <- if(!is.null(strata)) {
//     rep(strata, length.out = nInd(g)) 
//   } else strata(g)
//   strata.vec <- rep(strata.vec, ploidy(g))
//   fst.max.g <- g
//   for(i in 1:ncol(fst.max.g@loci)) {
//     new.locus <- paste(strata.vec, fst.max.g@loci[, i], sep = ".")
//     new.locus[is.na(g@loci[, i]) | is.na(strata.vec)] <- NA
//     fst.max.g@loci[, i] <- factor(new.locus)
//   }
//   
//   est <- statFst(g, strata, ...) / statFst(fst.max.g, strata, ...)
//   names(est) <- "F'st"
//   est
// }

/*** R
library(microbenchmark)
library(strataGdevel)
example(gtypes)
msats <- stratify(msats, "fine")
microbenchmark(
statFst(msats),
statFst_C(sapply(loci(msats), as.numeric), as.numeric(strata(msats)), ploidy(msats))
)
*/