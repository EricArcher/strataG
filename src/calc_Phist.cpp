#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix outerProd(const NumericVector &x, const NumericVector &y) {
  int nrow(x.size()), ncol(y.size());
  NumericMatrix m(nrow, ncol);
  for(int r = 0; r < nrow; r++) {
    for(int c = 0; c < ncol; c++) m(r, c) = x[r] * y[c];
  } 
  return m;
}

double calc_Phist_base(
  const int &n_st, 
  const IntegerVector &strata, const IntegerMatrix &ind_allele_freq,
  const NumericMatrix &hap_dist, const NumericMatrix &hap_freq
) {
  if(hap_dist == R_NilValue) return(NA_REAL);
  
  int s1, s2;
  double denom, ssd_wp = 0, ssd_ap = 0, phist;
  NumericVector hap_freq1, hap_freq2, st_freq(n_st);
  NumericMatrix freq_prod;
  NumericMatrix outerProd(const NumericVector&, const NumericVector&);
  
  // Calculate sums of squares among strata (Eqn 8b)
  for(s1 = 0; s1 < n_st; s1++) {
    hap_freq1 = hap_freq(_, s1);
    st_freq[s1] = sum(hap_freq1);
    freq_prod = outerProd(hap_freq1, hap_freq1);
    denom = 2 * st_freq[s1];
    ssd_wp += sum(freq_prod * hap_dist) / denom;
    for(s2 = 0; s2 < n_st; s2++) {
      hap_freq2 = hap_freq(_, s2);
      freq_prod = outerProd(hap_freq1, hap_freq2);
      ssd_ap += sum(freq_prod * hap_dist);
    }
  }
  ssd_ap = (ssd_ap / sum(2 * st_freq)) - ssd_wp;
  
  // Calculate average sample size correction for among strata variance
  //  Eqn 9a in paper, but modified as in Table 8.2.1.1 from Arlequin v3.5.1 manual
  //  (denominator is sum{I} - 1)
  double n = sum(st_freq);
  double n_avg = (n - sum(pow(st_freq, 2) / n)) / (n_st - 1);
  
  // Calculate variance components (Table 1)
  //   Set MSD (SSD / df) equal to expected MSD
  double Vc = ssd_wp / (n - n_st);
  double Vb = ((ssd_ap / (n_st - 1)) - Vc) / n_avg;
  
  phist = Vb / (Vb + Vc);
  if(std::isnan(phist)) phist = NA_REAL;
  return phist;
}

// [[Rcpp::export]]
NumericVector calc_Phist(    
  const int &n_ind, const int &n_allele, const int &n_st, 
  const IntegerVector &strata, const IntegerMatrix &ind_allele_freq,
  const NumericMatrix &unit_dist, const NumericMatrix &hap_dist
) {
  int i, h;
  NumericVector result(2);
  NumericMatrix hap_freq(n_allele, n_st);
  
  for(i = 0; i < n_ind; i++) {
    for(h = 0; h < n_allele; h++) {
      hap_freq(h, strata[i]) += ind_allele_freq(i, h);
    }
  }
  
  result[0] = calc_Phist_base(
    n_st, strata, ind_allele_freq, unit_dist, hap_freq
  );
  result[1] = calc_Phist_base(
    n_st, strata, ind_allele_freq, hap_dist, hap_freq
  );
  result.names() = CharacterVector({"Fst", "PHIst"});

  return result;
}