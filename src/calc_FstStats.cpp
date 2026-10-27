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
  double Gst_prime = Fst * (num_st - 1 + Hs) / ((num_st - 1) * (1 - Hs)); //Fst / ((num_st - 1) * (1 - Hs)) / (num_st - 1 + Hs);
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
  double n, mean_sum_p2;
  double num_st, mean_n, Ho2mean_n, Ho_l, Hs_l, Ht_l;
  double Ho_ovl = 0, Hs_ovl = 0, Ht_ovl = 0, Nj, Nij;
  double pbar;
  double sum_a = 0, sum_b = 0, sum_c = 0;
  NumericVector n_gtypes_l, sum_p2(n_loc);
  NumericVector r(n_loc), nbar(n_loc), rnbar(n_loc), nc(n_loc);
  NumericVector var_a(n_loc), var_b(n_loc), var_c(n_loc);
  NumericVector a_term_s(n_st), a_term(n_loc), b_term(n_loc);
  NumericMatrix n_gtypes(n_loc, n_st), sum_p2_ls(n_loc, n_st);
  NumericMatrix allele_count(n_allele, n_st), p_allele(n_allele, n_st);
  NumericMatrix hom_sum(n_loc, n_st);      // sum of c(c-1) per locus x stratum
  NumericMatrix c2_sum(n_allele, n_st);    // sum of c^2 per allele x stratum
  NumericMatrix Ho_ls(n_loc, n_st), loc_stats(n_loc + 1, 16);
  
  for(i = 0; i < n_ind; i++) {
    for(a = 0; a < n_allele; a++) {
      double c = ind_allele_freq(i, a);
      allele_count(a, strata[i]) += c;
      n_gtypes(locus[a], strata[i]) += c;
      hom_sum(locus[a], strata[i]) += c * (c - 1);
      c2_sum(a, strata[i]) += c * c;
    }
  }
  
  for(s = 0; s < n_st; s++) {
    for(a = 0; a < n_allele; a++) {
      double nt = n_gtypes(locus[a], s);
      if(nt > 0) {
        p_allele(a, s) = allele_count(a, s) / nt;
        sum_p2_ls(locus[a], s) += pow(p_allele(a, s), 2);
      }  // else p_allele stays 0
    }
    for(l = 0; l < n_loc; l++) {
      n = n_gtypes(l, s) / ploidy;
      if(n > 0) Ho_ls(l, s) = 1 - hom_sum(l, s) / (n * ploidy * (ploidy - 1));
    }
  }
  
  // for wcFst
  for(l = 0; l < n_loc; l++) {
    n_gtypes_l = n_gtypes(l, _) / ploidy;
    r[l] = sum(n_gtypes_l > 0);
    nbar[l] = sum(n_gtypes_l / r[l]);
    rnbar[l] = r[l] * nbar[l];
    nc[l] = (rnbar[l] - sum(pow(n_gtypes_l, 2) / rnbar[l])) / (r[l] - 1);
  }
    
  for(a = 0; a < n_allele; a++) {
    int la = locus[a];
    double p_sum = 0;
    int n_pop = 0;
    for(s = 0; s < n_st; s++) {
      Nj = n_gtypes(la, s);
      Nij = allele_count(a, s);
      if(Nj > 0) {
        p_sum += p_allele(a, s);
        n_pop++;
      }
      a_term_s[s] = (Nj > 0) ? Nij / Nj : 0;
      if(Nj > 1) b_term[la] += Nij * (Nij - 1) / (Nj * (Nj - 1));
    }
    sum_p2[la] += pow(p_sum / n_pop, 2);
    a_term[la] += (pow(sum(a_term_s), 2) - sum(pow(a_term_s, 2))) / (n_pop - 1);
    
    // for wcFst
    double N  = rnbar[la];
    double rr = r[la];
    double ssp = 0, ssi = 0, ssg = 0;
    pbar = sum(allele_count(a, _)) / (ploidy * N);
    for(s = 0; s < n_st; s++) {
      double ni = n_gtypes(la, s) / ploidy;
      if(ni > 0) {
        double pi = p_allele(a, s);
        ssp += ploidy * ni * pow(pi - pbar, 2);
        ssi += c2_sum(a, s) / ploidy - ploidy * ni * pi * pi;
        ssg += allele_count(a, s) - c2_sum(a, s) / ploidy;
      }
    }
    double msp = ssp / (rr - 1);
    double msi = ssi / (N - rr);
    double msg = ssg / (N * (ploidy - 1));
    var_c[la] += msg;
    var_b[la] += (msi - msg) / ploidy;
    var_a[la] += (msp - msi) / (ploidy * nc[la]);
  }

  int n_used = 0;
  LogicalVector loc_ok(n_loc);   // FALSE by default
  LogicalVector st_ok(n_st);
  
  for(l = 0; l < n_loc; l++) {
    double sp2_sum = 0, ho_sum = 0, inv_n_sum = 0;
    num_st = 0;
    for(s = 0; s < n_st; s++) {
      n = n_gtypes(l, s) / ploidy;
      if(n > 0) {
        sp2_sum   += sum_p2_ls(l, s);
        ho_sum    += Ho_ls(l, s);
        inv_n_sum += 1 / n;
        num_st++;
      }
    }
    
    if(num_st < 2) {
      loc_stats(l + 1, _) = NumericVector(16, NA_REAL);
      continue;
    }
    loc_ok[l] = true;
    n_used++;
    for(s = 0; s < n_st; s++) if(n_gtypes(l, s) > 0) st_ok[s] = true;
    
    mean_sum_p2 = sp2_sum / num_st;
    Ho_l = ho_sum / num_st;
    mean_n = num_st / inv_n_sum;
    Ho2mean_n = Ho_l * (ploidy - 1) / ploidy / mean_n;
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
  
  if(n_used == 0) stop("No locus has data in at least two strata");
  
  int n_st_used = 0;
  for(s = 0; s < n_st; s++) if(st_ok[s]) n_st_used++;
  
  NumericVector a_ok = a_term[loc_ok];
  NumericVector b_ok = b_term[loc_ok];
  loc_stats(0, _) = calc_stats(
    Ho_ovl / n_used, Hs_ovl / n_used, Ht_ovl / n_used, (double) n_st_used,
    harmonic_mean(a_ok), harmonic_mean(b_ok),
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