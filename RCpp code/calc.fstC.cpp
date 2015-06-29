#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

/*** R
calc.fstC <- function(g) {
  # if haploid data, run 'calc.phist' with no sequences (all haplotypes assumed to be equi-distant)
  if(is.haploid(g)) {
    return(calc.phistC(g))
  }
 est <- calcfst(g$locus.data, g$strata)
  
  if(is.nan(est)) est <- NA
  result <- list(stat.name = "Fst", estimate = est)
  class(result) <- c(class(result), "gtype.struct.stat")
  result
 }
*/

// [[Rcpp::export]]
double calcfst(NumericMatrix locusdata, NumericVector strata){

// function declaration
NumericVector rowSumsC(NumericMatrix, int, int);

// calculations
  int nrow = locusdata.nrow();
  int ncol = locusdata.ncol();
  int numlocus = ncol/2;
  
  NumericMatrix locusSums(2, numlocus);
  NumericVector uniqueS = unique(strata);
  int sizeuniqueS = uniqueS.size();
  uniqueS.sort();

  int locusSumsColumn = 0;
// int column = 0;
// loop through each locus to fill locusSums
  for(int column = 0; column < ncol; column +=2){
    NumericVector a1(nrow), a2(nrow);
    NumericVector uniqueLD(nrow*2);

// create a1 and a2 and uniqueLD
    for(int i = 0; i < nrow; i++){
      a1(i) = locusdata(i, column);
      a2(i) = locusdata(i, column+1);
      uniqueLD(i*2)       = locusdata(i, column); 
      uniqueLD((i*2) + 1) = locusdata(i, column+1);
    }
    
    uniqueLD = unique(uniqueLD);
    uniqueLD.sort();
    

    
// Remove NA from unique alleles
    for (int i = 0; i < uniqueLD.size(); i++){
      if(uniqueLD[i]== -1){
        uniqueLD.erase(i);
        break;
      }
    }
    int sizeuniqueLD = uniqueLD.size();

// if locus fixed for allele break out of loop. consider locusSums = 0 for each locus
    if(sizeuniqueLD < 2){
      break;    // a break out of for loop for locus
    }
    
// counting strata
  NumericVector stratasum(sizeuniqueS);
    for (int i = 0; i < nrow; i++){
      if(a1(i) == -1){
          continue;    
      }
      stratasum(strata(i))+=1;      
    }
    
    if(sizeuniqueS < 2){
      break;    //  a break out of for loop for locus
    }
    
    double nbar = mean(stratasum);
    double rnbar = sizeuniqueS * nbar;
    double nc = 0;
    
    for(int i = 0; i < sizeuniqueS; i++){
      nc += (pow(stratasum(i), 2));
    }
    nc = ((rnbar - (nc/rnbar)) / (sizeuniqueS - 1));

    NumericMatrix varcompmatrix(3, sizeuniqueLD);
// loop through alleles of the locus    
    for(int currentA = 0; currentA < sizeuniqueLD; currentA++){
//    int currentA = 6;
      NumericMatrix allelestats(2, sizeuniqueS);
      
      for(int i = 0; i < nrow; i++){
        if( a1(i) == uniqueLD(currentA)) {
            allelestats(0, strata(i)) += 1; 
            if (a1(i) != a2(i)){
              allelestats(1, strata(i)) += 1;
            }
        }
        if( a2(i) == uniqueLD(currentA)){
            allelestats(0, strata(i)) += 1;
            if (a1(i) != a2(i)){
              allelestats(1, strata(i)) += 1;
            }
        }
      }
      
      double pbar = 0, hbar = 0; 
// doing math on allelestats and math for pbar and hbar 
      for(int currentS = 0; currentS < sizeuniqueS; currentS++){
        allelestats(0, currentS) = (allelestats(0, currentS) / stratasum(currentS) / 2);
        allelestats(1, currentS) = (allelestats(1, currentS) / stratasum(currentS));
        pbar += (allelestats(0, currentS) * stratasum(currentS));
        hbar += (allelestats(1, currentS) * stratasum(currentS));
      }     
      pbar = pbar/rnbar; 
      hbar = hbar/rnbar;

//    calculate s2
      double s2 = 0;
      for(int currentS = 0; currentS < sizeuniqueS; currentS++){
        s2 += (stratasum(currentS) * (pow((allelestats(0, currentS) - pbar),2)));
      }
      s2 = s2 / (sizeuniqueS - 1) / nbar ;  

// calculate innerterm
      double innerterm = (pbar * (1-pbar)) - ((sizeuniqueS - 1) * s2 / sizeuniqueS);
      double nbarm1 = nbar - 1; 

// Va (between strata) - Eqn. 2
      double Va1 = 0.25 * hbar; 
      double Va2 = s2 - (innerterm - Va1) / nbarm1;
      double Va = nbar * Va2 / nc;
// Vb (between individulas within strata) - Eqn. 3
      double Vb1 = (( 2 * nbar - 1) * hbar) / 4 / nbar; 
      double Vb2 = innerterm - Vb1;
      double Vb = nbar * Vb2 / nbarm1;
// Vc (between gametes within individuals) - Eqn. 4
      double Vc = 0.5 * hbar;

      varcompmatrix(0, currentA) = Va;
      varcompmatrix(1, currentA) = Vb; 
      varcompmatrix(2, currentA) = Vc; // rows are variance components (1=Va, 2=Vb, 3=Vc), columns are alleles
    } // for loop through each allele close
    
    NumericVector varcompRowSums = rowSumsC(varcompmatrix, sizeuniqueS, sizeuniqueLD);
    
    locusSums(0, locusSumsColumn) = varcompRowSums(0);
    locusSums(1, locusSumsColumn) = sum(varcompRowSums);

    locusSumsColumn++;
  }  // for loop through each locus close 

  NumericVector locusSumsRowSums = rowSumsC(locusSums, 2, numlocus);
  double est = locusSumsRowSums(0) / locusSumsRowSums(1);

  return(est);  
}

// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix data, int nrow, int ncol){
  NumericVector rowSums(nrow);
    for (int i = 0; i < nrow; i++) {
        double total = 0;
        for (int j = 0; j < ncol; j++) {
            total += data(i, j);
         }
        rowSums[i] = total;
    }
  return(rowSums);
}

