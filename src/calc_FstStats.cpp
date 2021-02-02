#include <Rcpp.h>
using namespace Rcpp;

NumericVector calc_stats(
    const double &Ho, const double &Hs, const double &Ht, 
    const double &num_st, const double &a_term, const double &b_term,
    const double &a, const double &b, const double &c
) {
  
  double Dst = Ht - Hs;
  double Dst_prime = num_st / (num_st - 1) * Dst;
  double Ht_prime = Dst_prime + Hs;
  double Fst = Dst / Ht;
  double Fst_prime = Dst_prime / Ht_prime;
  double Fis = 1 - Ho / Hs;
  double Gst_prime = Fst / ((num_st - 1) * (1 - Hs)) / (num_st - 1 + Hs);
  double Gst_dbl_prime = Fst_prime / (1 - Hs);
  double Dest = Dst_prime / (1 - Hs);
  double Dest_Chao = 1 - (a_term / b_term);
  double wcFit = 1 - (c / (a + b + c));
  double wcFst = a / (a + b + c);
  double wcFis = 1 - (c / (b + c));
  
  return NumericVector::create(
    Ho, Hs, Ht, Ht_prime, Dst, Dst_prime, 
    Fst, Fst_prime, Fis, Gst_prime, Gst_dbl_prime, 
    Dest, Dest_Chao,
    wcFit, wcFst, wcFis
  );
}

double harmonic_mean(const NumericVector &x) {
  if(is_true(any(x == 0))) {
    double inv_mean = 1 / mean(x);
    return 1 / (inv_mean + var(x) * pow(inv_mean, 3));
  } else {
    return x.size() / sum(1 / x);
  }
}

// [[Rcpp::export]]
NumericMatrix calc_FstStats(
    const int &n_ind, const int &n_allele, const int &n_loc, const int &n_st, 
    const double &ploidy, const IntegerVector &strata, 
    const IntegerVector &locus, const IntegerMatrix &ind_allele_freq
) {
  int i, a, l, s;
  double n, lsHo, mean_sum_p2;
  double num_st, mean_n, Ho2mean_n, Ho_l, Hs_l, Ht_l;
  double Ho_ovl = 0, Hs_ovl = 0, Ht_ovl = 0, Nj, Nij;
  double pbar, s2, hbar, inner_term, term_1, term_2;
  double sum_a = 0, sum_b = 0, sum_c = 0;
  NumericVector n_ind_loc, n_gtypes_l, sum_p2(n_loc);
  NumericVector r(n_loc), nbar(n_loc), rnbar(n_loc), nc(n_loc);
  NumericVector var_a(n_loc), var_b(n_loc), var_c(n_loc);
  NumericVector a_term_s(n_st), a_term(n_loc), b_term(n_loc), d(n_loc);
  NumericMatrix n_gtypes(n_loc, n_st), n_hom(n_loc, n_st), sum_p2_ls(n_loc, n_st);
  NumericMatrix allele_count(n_allele, n_st), p_allele(n_allele, n_st);
  NumericMatrix n_het(n_allele, n_st);
  NumericMatrix Ho_ls(n_loc, n_st), Hs_ls(n_loc, n_st), loc_stats(n_loc + 1, 16);
  
  double harmonic_mean(const NumericVector&);
  NumericVector calc_stats(
    const double&, const double&, const double&, 
    const double&, const double&, const double&,
    const double&, const double&, const double&
  );
  
  for(i = 0; i < n_ind; i++) {
    for(a = 0; a < n_allele; a++) {
      allele_count(a, strata[i]) += ind_allele_freq(i, a);
      n_gtypes(locus[a], strata[i]) += ind_allele_freq(i, a);
      if(ind_allele_freq(i, a) == ploidy) {
        n_hom(locus[a], strata[i])++;
      }
      // for wcFst
      if(ind_allele_freq(i, a) >= 1 && ind_allele_freq(i, a) != ploidy) {
        n_het(a, strata[i])++; 
      }
    }
  }
  
  for(s = 0; s < n_st; s++) {
    for(a = 0; a < n_allele; a++) {
      p_allele(a, s) = allele_count(a, s) / n_gtypes(locus[a], s);
      sum_p2_ls(locus[a], s) += pow(p_allele(a, s), 2); 
    }
    for(l = 0; l < n_loc; l++) {
      n = n_gtypes(l, s) / ploidy;
      lsHo = 1 - (n_hom(l, s) / n);
      Ho_ls(l, s) = lsHo;
      Hs_ls(l, s) = (n / (n - 1)) * (1 - sum_p2_ls(l, s) - (lsHo / (2 * n)));
    }
  }
  
  // for wcFst
  for(l = 0; l < n_loc; l++) {
    n_gtypes_l = n_gtypes(l, _) / ploidy;
    r[l] = sum(n_gtypes_l > 0);
    nbar[l] = sum(n_gtypes_l / r[l]);
    rnbar[l] = r[l] * nbar[l];
    nc[l] = (rnbar[locus[l]] - sum(pow(n_gtypes_l, 2) / rnbar[l])) / (r[l] - 1);
  }
  
  for(a = 0; a < n_allele; a++) {
    sum_p2[locus[a]] += pow(mean(p_allele(a, _)), 2);      
    // for Jost's D (D_est_Chao) Jost 2008 eqn 13
    n = sum(n_gtypes(locus[a], _) > 0);
    for(s = 0; s < n_st; s++) {
      Nij = allele_count(a, s);
      Nj = n_gtypes(locus[a], s);
      a_term_s[s] = Nij / Nj;
      b_term[locus[a]] += Nij * (Nij - 1) / (Nj * (Nj - 1));
    }
    a_term[locus[a]] += (pow(sum(a_term_s), 2) - sum(pow(a_term_s, 2))) / (n - 1);
    
    // for wcFst
    pbar = sum(allele_count(a, _) / (rnbar[locus[a]] * ploidy));
    s2 = sum((n_gtypes(locus[a], _) / ploidy) * pow((p_allele(a, _) - pbar), 2) / (rnbar[locus[a]] - nbar[locus[a]]));
    hbar = sum(n_het(a, _) / rnbar[locus[a]]);
    inner_term = (pbar * (1 - pbar)) - (((r[locus[a]] - 1) / r[locus[a]]) * s2);
    term_1 = 1 / (nbar[locus[a]] - 1);
    term_2 = inner_term - (hbar / 4);
    var_a[locus[a]] += (nbar[locus[a]] / nc[locus[a]]) * (s2 - term_1 * term_2);
    term_1 = (((2 * nbar[locus[a]]) - 1) / (4 * nbar[locus[a]])) * hbar;
    var_b[locus[a]] += (nbar[locus[a]] / (nbar[locus[a]] - 1)) * (inner_term - term_1);
    var_c[locus[a]] += hbar / 2;
  }
  
  for(l = 0; l < n_loc; l++) {
    mean_sum_p2 = mean(sum_p2_ls(l, _));
    n_ind_loc = n_gtypes(l, _) / ploidy;
    num_st = sum(n_ind_loc > 0);
    mean_n = num_st / sum(1 / n_ind_loc);
    Ho_l = mean(Ho_ls(l, _));
    Ho2mean_n = Ho_l / 2 / mean_n;
    Hs_l = mean_n / (mean_n - 1) * (1 - mean_sum_p2 - Ho2mean_n);
    Ht_l = 1 - sum_p2[l] + (Hs_l / mean_n / num_st) - (Ho2mean_n / num_st); 
    loc_stats(l + 1, _) = calc_stats(
      Ho_l, Hs_l, Ht_l, num_st, 
      a_term[l], b_term[l],
      var_a[l], var_b[l], var_c[l]
    );
    Ho_ovl += Ho_l;
    Hs_ovl += Hs_l;
    Ht_ovl += Ht_l;
    if(!std::isnan(var_a[l])) sum_a += var_a[l];
    if(!std::isnan(var_b[l])) sum_b += var_b[l];
    if(!std::isnan(var_c[l])) sum_c += var_c[l];
  }
  
  loc_stats(0, _) = calc_stats(
    Ho_ovl / n_loc, Hs_ovl / n_loc, Ht_ovl / n_loc, n_st, 
    harmonic_mean(a_term), harmonic_mean(b_term),
    sum_a, sum_b, sum_c
  );
  
  colnames(loc_stats) = CharacterVector::create(
    "Ho", "Hs", "Ht", "Ht_prime", "Dst", "Dst_prime",
    "Fst", "Fst_prime", "Fis", "Gst_prime", "Gst_dbl_prime", 
    "Dest", "Dest_Chao",
    "wcFit", "wcFst", "wcFis"
  );
  return loc_stats;
}