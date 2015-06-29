#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*** R

calc.chi2C <- function(g) {
  stopifnot.gtypes(g)
  if((ncol(g$locus.data)%%2) == 0){
    est <- chi2DN(g$locus.data, g$strata)   
  }
  else {
    est <- chi2N(g$locus.data, g$strata)    
  }
  
  result <- list(stat.name = "Chi2", estimate = est)
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
*/


// [[Rcpp::export]] 
double chi2N(NumericVector locusdata, NumericVector strata) {

// function declarations
  double calcchi2(NumericMatrix, int, int, NumericVector, NumericVector, int);
  NumericVector colSumsC(NumericMatrix, int, int);
  NumericVector rowSumsC(NumericMatrix, int, int);

// create new vectors containing the unique options
  NumericVector uniqueLD = unique(locusdata); 
  NumericVector uniqueS = unique(strata);
  uniqueLD.sort();
  uniqueS.sort();
  
// find size of each vector
  int sizeuniqueLD = uniqueLD.size();
  int sizeuniqueS = uniqueS.size();
  int sizedata = locusdata.size(); // number of individuals

// create blank table "obsfreq" for chi-square test
  NumericMatrix obsfreq(sizeuniqueLD, sizeuniqueS);
    
// counting obsfreq
      for (int i = 0; i < sizedata; i++){
          obsfreq(locusdata[i],strata[i]) += 1; 
      }    
       
// rowSums&colSums of obsfreq 
    int nrow = obsfreq.nrow(), ncol = obsfreq.ncol();
    NumericVector rowSums = rowSumsC(obsfreq, nrow, ncol);
    NumericVector colSums = colSumsC(obsfreq, nrow, ncol);

// Chi2 function call
    double chi2 = calcchi2(obsfreq, nrow, ncol, rowSums, colSums, sizedata);
 
   return(chi2); 
}



// [[Rcpp::export]]
double chi2DN(NumericMatrix locusdata, NumericVector stratadata){

// function declarations
  double calcchi2(NumericMatrix, int, int, NumericVector, NumericVector, int);
  NumericVector colSumsC(NumericMatrix, int, int);
  NumericVector rowSumsC(NumericMatrix, int, int);
 
// Strata calculations
    int n = stratadata.size();
    NumericVector strata(n*2); // new vector used in frequency calculations
    for (int i = 0; i < n; i++){
      strata(i*2)     = stratadata(i);
      strata((i*2)+1) = stratadata(i);
    }



// Calculation for Number of Unique Strata
    NumericVector uniqueS = unique(stratadata);
    uniqueS.sort();
       
// Declared variables for loop  
  double chi2 = 0;
  int nrow = locusdata.nrow(), ncol = locusdata.ncol();

// Loop through columns (Every two columns of locusdata matrix is a different gene) 
  for (int column = 0; column < ncol; column+=2) {
    NumericVector templocusdata(nrow*2);
      for (int i = 0; i < nrow; i++){
        templocusdata(i*2)       = locusdata(i, column); 
        templocusdata((i*2) + 1) = locusdata(i, column+1);
      }

// Calculation for Number of Unique Alleles 
    NumericVector uniqueLD = unique(templocusdata);
    uniqueLD.sort();
    
// Remove NA from unique alleles
    for (int i = 0; i < uniqueLD.size(); i++){
      if(uniqueLD[i]== -1){
        uniqueLD.erase(i);
        break;
      }
    }

// size of each vector
      int sizeuniqueLD = uniqueLD.size();
      int sizeuniqueS = uniqueS.size();
      int sizedata = templocusdata.size();
    
// create blank table "obsfreq" for chi-square test
      NumericMatrix obsfreq(sizeuniqueLD, sizeuniqueS);
 
// counting obsfreq
      for (int i = 0; i < sizedata; i++){
         if(templocusdata[i]== -1){
            continue;
          }
          obsfreq(templocusdata[i], strata[i]) += 1;               
      } 
      
// rowSums&ColSums of obsfreq 
      int obnrow = obsfreq.nrow(), obncol = obsfreq.ncol();
      NumericVector rowSums = rowSumsC(obsfreq, obnrow, obncol);
      NumericVector colSums = colSumsC(obsfreq, obnrow, obncol);

      int sum = 0;
      n = colSums.size();
      for (int i = 0; i< n; i++){
      sum += colSums[i];
      }
      
// Chi2 Calculation
       chi2 += calcchi2(obsfreq, obnrow, obncol, rowSums, colSums, sum);
 
 }
  return(chi2);
}

// [[Rcpp::export]]
double calcchi2(NumericMatrix obsfreq, int nrow, int ncol, NumericVector rowSums, NumericVector colSums, int sizedata){
  double chi2 = 0;
   for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      double expfreq = rowSums[i]*colSums[j] / sizedata;
      double temp = obsfreq(i,j) - expfreq;
      chi2 += temp*temp / expfreq;
    }
  }
  return(chi2);
  
}

// [[Rcpp::export]]
NumericVector colSumsC(NumericMatrix obsfreq, int nrow, int ncol) {
  NumericVector colSums(ncol);
    for (int i = 0; i < ncol; i++) {
        double total = 0;
        for (int j = 0; j < nrow; j++) {
            total += obsfreq(j, i);
         }
        colSums[i] = total;
    }
    return(colSums);
}

// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix obsfreq, int nrow, int ncol){
  NumericVector rowSums(nrow);
    for (int i = 0; i < nrow; i++) {
        double total = 0;
        for (int j = 0; j < ncol; j++) {
            total += obsfreq(i, j);
         }
        rowSums[i] = total;
    }
  return(rowSums);
}

