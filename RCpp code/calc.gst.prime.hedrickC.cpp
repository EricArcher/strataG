#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*** R
calc.gst.prime.hedrickC <- function(g) {
  hets <- calc.est.H.statsC(g)
  gst <- calc.gstC(g)
  est <- calcgstprimeH(g$strata, hets, gst$estimate)
  
  if(is.nan(est)) est <- NA
  result <- list(stat.name = "G'st (Hedrick 2005)", estimate = est)
  class(result) <- c(class(result), "gtype.struct.stat")
  result
  }
class(calc.gst.prime.hedrick) <- c(class(calc.gst.prime.hedrick), "gtype.struct.func")
*/

// [[Rcpp::export]]
double calcgstprimeH(NumericVector strata, NumericMatrix hets, double gstestimate){
//function declarations
  double rowmean(NumericMatrix, int, int);
  
// unique strata calculations 
  NumericVector uniqueS = unique(strata);
  int sizeuniqueS = uniqueS.size();  //numstrata
  
// gst calculation
  int ncol = hets.ncol();
  double Hs = rowmean(hets, 1, ncol);
  double gstmax = (((sizeuniqueS - 1) * (1 - Hs)) / (sizeuniqueS - 1 + Hs));
  double est = gstestimate/gstmax;

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