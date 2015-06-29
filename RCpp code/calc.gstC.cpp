#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*** R
calc.gstC <- function(g) {
  hets <- calc.est.H.statsC(g)
  est <- calcgst(g$strata, hets)
  
  if(is.nan(est)) est <- NA
  result <- list(stat.name = "Gst", estimate = est)
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(calc.gst) <- c(class(calc.gst), "gtype.struct.func")

*/

// [[Rcpp::export]]
double calcgst(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowmean(NumericMatrix, int, int);

// gst calculation
  int ncol = hets.ncol();
  double Hs = rowmean(hets, 1, ncol);
  double Ht = rowmean(hets, 2, ncol);
  double est = 1 - (Hs /Ht);
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