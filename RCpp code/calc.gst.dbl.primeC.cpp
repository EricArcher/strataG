#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*** R
calc.gst.dbl.primeC <- function(g) {
 hets <- calc.est.H.statsC(g)
 est <- calcgstdblprime(g$strata, hets)
 
  if(is.nan(est)) est <- NA
  result <- list(stat.name = "G''st", estimate = est)
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(calc.gst) <- c(class(calc.gst), "gtype.struct.func")
 
*/

// [[Rcpp::export]]
double calcgstdblprime(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowmean(NumericMatrix, int, int);
  double calcgstprime(NumericVector, NumericMatrix);
    
// Calculate primenei
  double primenei = calcgstprime(strata, hets);
      
// double prime calculation
  int ncol = hets.ncol();
  double Hs = rowmean(hets, 1, ncol);
  double est = (primenei / (1 - Hs));
    
  return(est);
}


// [[Rcpp::export]]
double rowmean(NumericMatrix data, int rowToAverage, int ncol){
  double sum = 0;
  for(int i = 0; i < ncol; i++){
    sum += data(rowToAverage, i);
  }
  double mean = sum/ncol;
  return(mean);
}

// [[Rcpp::export]]
double calcgstprime(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowmean(NumericMatrix, int, int);
    
// unique strata calculations
  NumericVector uniqueS = unique(strata);
  uniqueS.sort();
  int sizeuniqueS = uniqueS.size();  //numstrata
    
// count strata frequency
  int numindividuals = strata.size();
  NumericVector StrataFreq(sizeuniqueS);
  for(int i = 0; i < numindividuals; i++){
    StrataFreq(strata[i]) += 1;
  }
    
// gst calculation
  int ncol = hets.ncol();
  double Hs = rowmean(hets, 1, ncol);
  double Ht = rowmean(hets, 2, ncol);
    
  double est = ((sizeuniqueS * (Ht - Hs)) / ((sizeuniqueS * Ht) - Hs));
    
  return(est);
}  
