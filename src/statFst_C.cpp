#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix alleleFreqCalc(IntegerVector locVec, IntegerVector strataRep) {
  IntegerMatrix table2D(IntegerVector, IntegerVector);
  NumericVector colSumC(NumericMatrix);
  
  // get allele frequencies in each population
  NumericMatrix alleleFreq = wrap(table2D(locVec, strataRep));
  NumericVector alleleFreqColSums = colSumC(alleleFreq);
  for(int r = 0; r < alleleFreq.nrow(); r++) {
    for(int c = 0; c < alleleFreq.ncol(); c++) {
      alleleFreq(r, c) = alleleFreq(r, c) / alleleFreqColSums[c];
    }
  }
  
  return alleleFreq;
}

NumericMatrix prHetCalc(IntegerVector locus, int nalleles,
                        IntegerVector strata, IntegerVector nvec, int ploidy) {
  IntegerVector idGenotype(IntegerVector, int, int);
  
  // get proportion heterozygosity for each allele (rows) and population (cols)
  NumericMatrix prHet(nalleles, nvec.size());
  for(int i = 0; i < strata.size(); i++) {
    IntegerVector gt = idGenotype(locus, i, ploidy);
    bool NAs = false;
    for(int j = 0; j < gt.size(); j++) {
      if(IntegerVector::is_na(gt[j])) {
        NAs = true;
        continue;
      }
    }
    if(NAs) continue;
    if(unique(gt).size() == 1) continue;
    for(int a = 0; a < gt.size(); a++) prHet(gt[a], strata[i]) += 1;
  }
  for(int c = 0; c < prHet.ncol(); c++) prHet(_, c) = prHet(_, c) / nvec[c];
  
  return prHet;
}

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
  
  return varcompMat;
}

double fstCalc(IntegerMatrix loci, IntegerVector strata, int ploidy) {
  // function declarations
  IntegerVector calcStrataN(IntegerVector, IntegerVector, int);
  
  // returns numerator and denominator sums for
  //   theta-w calculation for alleles at each locus (page 1363)
  IntegerVector strataRep(rep_each(strata, loci.nrow() / strata.size()));
  NumericMatrix locusSums(2, loci.ncol());
  for(int loc = 0; loc < loci.ncol(); loc++) {
    IntegerVector locVec = loci(_, loc);
    
    // identify unique alleles and return zeroes if locus is fixed
    int nalleles = unique(na_omit(locVec)).size();
    if(nalleles < 2) continue;
    
    IntegerVector nvec = calcStrataN(locVec, strata, ploidy);
    int r = nvec.size();
    if(r < 2) continue;
    double nbar = mean(nvec);
    if(nbar == 1) continue;
    double rnbar = r * nbar;
    double nc = (rnbar - (sum(Rcpp::pow(nvec, 2.0)) / rnbar)) / (r - 1);
    NumericMatrix alleleFreq = alleleFreqCalc(locVec, strataRep);
    NumericMatrix prHet = prHetCalc(locVec, nalleles, strata, nvec, ploidy);
    NumericMatrix varcompMat = varCompCalc(nvec, alleleFreq, prHet, r, nbar, rnbar, nc);
    locusSums(0, loc) = sum(varcompMat(0, _));
    for(int i = 0; i < varcompMat.nrow(); i++) {
      locusSums(1, loc) += sum(varcompMat(i, _));
    }
  }
  double est(sum(locusSums(0, _)) / sum(locusSums(1, _)));
  if(std::isnan(est)) est = NA_REAL;
  return est;
}

IntegerMatrix maxFstLoci(IntegerMatrix loci, IntegerVector strata, int ploidy, IntegerVector maxAllele) {
  IntegerVector st(rep_each(strata, ploidy));
  IntegerMatrix maxLoci(loci.nrow(), loci.ncol());
  
  for(int c = 0; c < maxLoci.ncol(); c++) {
    for(int r = 0; r < maxLoci.nrow(); r++) {
      if(IntegerVector::is_na(maxLoci(r, c))) {
        maxLoci(r, c) = NA_INTEGER;
      } else {
        maxLoci(r, c) = loci(r, c) + (maxAllele[c] * st[r]);
      }
    }
    IntegerVector alleles(unique(maxLoci(_, c)).sort());
    for(int r = 0; r < maxLoci.nrow(); r++) {
      if(IntegerVector::is_na(maxLoci(r, c))) continue;
      for(int a = 0; a < alleles.size(); a++) {
        if(maxLoci(r, c) == alleles[a]) {
          maxLoci(r, c) = a;
          break;
        }
      }
    }
  }
  
  return maxLoci;
}

// [[Rcpp::export]]
NumericVector statFst_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  NumericVector estVec(strataMat.ncol());
  int ploidy(loci.nrow() / strataMat.nrow());
  for(int idx = 0; idx < estVec.size(); idx++) {
    IntegerVector strata(strataMat(_, idx));
    estVec[idx] = fstCalc(loci, strata, ploidy);
  }
  return estVec;
}

// [[Rcpp::export]]
NumericVector statFstPrime_C(IntegerMatrix loci, IntegerMatrix strataMat) {
  IntegerVector maxAllele(loci.ncol());
  for(int c = 0; c < maxAllele.size(); c++) {
    maxAllele[c] = max(na_omit(loci(_, c))) + 1;
  }
  
  int ploidy(loci.nrow() / strataMat.nrow());
  IntegerMatrix maxLoci;
  NumericVector estVec(strataMat.ncol());
  for(int idx = 0; idx < estVec.size(); idx++) {
    IntegerVector strata(strataMat(_, idx));
    maxLoci = maxFstLoci(loci, strata, ploidy, maxAllele);
    estVec[idx] = fstCalc(loci, strata, ploidy) / fstCalc(maxLoci, strata, ploidy);
  }
  return estVec;
}