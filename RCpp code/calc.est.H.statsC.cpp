#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*** R
calc.est.H.statsC <- function(g) { 
  g <- check.gtypes(g, "diploid")
  estHstats(g$locus.data, g$strata)  
  }
*/

// [[Rcpp::export]]
NumericMatrix estHstats(NumericMatrix locusdata, NumericVector strata){
//function declaration
  double colmean(NumericMatrix, int, int); 
  NumericVector rowSumsC(NumericMatrix, int, int);
  
//create estH matrix to store Ho Hs Ht
  int ncol = locusdata.ncol(), numindividuals = locusdata.nrow();  
  NumericMatrix estH(3, ncol/2);
  int estHcol = 0;
  
// strata calculations (unique strata and a double list)  - only needs to be done once
    NumericVector uniqueS = unique(strata);
    uniqueS.sort();
    int sizeuniqueS = uniqueS.size();
    
    NumericVector strataDouble(numindividuals*2); // new vector used in frequency calculations
    for (int i = 0; i < numindividuals; i++){
      strataDouble(i*2)     = strata(i);
      strataDouble((i*2)+1) = strata(i);
    }    

// column to loop through genes ncol = number of genes //counter of homozygotes
  for (int column = 0; column < ncol; column+=2) { 
    
  int numhomozygotes = 0;
// creating a matrix of genotypes (ex: 1/1 = 1.01) with homozygote/heterozygote details in column 2 
    NumericMatrix tempgenotype(numindividuals, 2);
      for (int i = 0; i < numindividuals; i++){   
         tempgenotype(i,0) = (locusdata(i,column)) + (locusdata(i,column+1)/100); 
         if (locusdata(i,column) == locusdata(i,column+1)) {
           tempgenotype(i,1) = 1;
           numhomozygotes += 1;
         }
      }
    
// create unique vector of homozygous genotype options, and remove any NA's 
    NumericVector uniqueHG(numhomozygotes); 
    int counter = 0; //counter for uniqueHG;
    for (int i = 0; i < numindividuals; i++){
      if(tempgenotype(i,1) == 1){
        uniqueHG(counter) = tempgenotype(i,0);
        counter++;
      }
    }
    uniqueHG = unique(uniqueHG);   //remove loop

    for (int i = 0; i < uniqueHG.size(); i++){
      if(uniqueHG[i]== -1.01){
        uniqueHG.erase(i);
        break; }}
    uniqueHG.sort();
    int sizeuniqueHG = uniqueHG.size();
       
  int sizedata = strata.size(); // size of data without NA, will be recalculated for each gene         
// create blank table "genotype.strata.p" for HO test
  NumericMatrix genotypestratap(sizeuniqueS, sizeuniqueHG);
 
// counting genotype.strata.p
      for (int i = 0; i < numindividuals; i++){
        if(tempgenotype(i,0) == -1.01){
            sizedata -= 1;                                  //removing NA's from sizedata
            continue;    
        }
        for(int j = 0; j < sizeuniqueHG; j++){
          if(tempgenotype(i,0) == uniqueHG[j]){
            genotypestratap(strata[i],j) += 1;
            break;
          }
        } 
      }
// counting # from each strata
    NumericVector StrataCount(sizeuniqueS);
    for(int i = 0; i < numindividuals; i++){
      if(tempgenotype(i,0) == -1.01){
        continue;
      }
      else{
        StrataCount(strata[i]) += 1;
      }
    }

//proportion
    for (int j = 0; j < sizeuniqueHG; j++){
      for (int i = 0; i < sizeuniqueS; i++){
        genotypestratap(i,j) = genotypestratap(i,j)/StrataCount(i);
      }
    }



    double sumHomean = 0; 
// Ho calculation
    for(int i = 0; i < uniqueHG.size(); i++){
          sumHomean += colmean(genotypestratap, i, sizeuniqueS); 
    }
    double Ho = 1 - sumHomean;
    estH(0,estHcol) = Ho;
    
// Hs calculations
//create locus data
    NumericVector templocusdata(numindividuals*2);
      for (int i = 0; i < numindividuals; i++){
        templocusdata(i*2)       = locusdata(i, column); 
        templocusdata((i*2) + 1) = locusdata(i, column+1);
      }
// unique locusdata options
    NumericVector uniqueLD = unique(templocusdata);
    uniqueLD.sort();
// Remove NA from unique alleles
    for (int i = 0; i < uniqueLD.size(); i++){
      if(uniqueLD[i]== -1){
        uniqueLD.erase(i);
        break;
      }
    }
    int sizeuniqueLD = uniqueLD.size();
    
// create blank table "allele.strata.p" for Hs test
  NumericMatrix allelestratap(sizeuniqueS, sizeuniqueLD);
 
// counting allele.strata.p
      for (int i = 0; i < numindividuals*2; i++){
        if(templocusdata(i) == -1){
            continue;    
        }
        for(int j = 0; j < sizeuniqueLD; j++){
          if(templocusdata(i) == uniqueLD[j]){
            allelestratap(strataDouble[i],j) += 1;
            break;
          }
        } 
      }  
  
    NumericVector RowSumsASP = rowSumsC(allelestratap, sizeuniqueS, sizeuniqueLD);
//proportion
    for (int j = 0; j < sizeuniqueLD; j++){
      for (int i = 0; i < sizeuniqueS; i++){
        allelestratap(i,j) = allelestratap(i,j)/RowSumsASP(i);
      }
    }  

// calculate for each row
  NumericVector rowsummary(sizeuniqueS);
  for(int i =0; i < sizeuniqueS; i++){
    for(int j = 0; j < sizeuniqueLD; j++){
      rowsummary(i) += (pow(allelestratap(i,j),2));
    }
    rowsummary(i) = 1-rowsummary(i);
  }
  
  double meanrowsumsHS = mean(rowsummary);
// counting strata 
  NumericVector stratasum(sizeuniqueS);
    for (int i = 0; i < numindividuals*2; i++){
      if(templocusdata(i) == -1){
          continue;    
      }
      stratasum(strataDouble(i))+=1;      
    }
    
    double harmonicmeanHS = 0;
//calculate Harmonic mean
  for (int i = 0; i < stratasum.size(); i++){
    harmonicmeanHS += (1/(stratasum(i)/2));
  }
  harmonicmeanHS = stratasum.size()/harmonicmeanHS;

//Calculate Hs
  double Hs = (harmonicmeanHS/(harmonicmeanHS-1)) * (meanrowsumsHS - ((Ho/2)/harmonicmeanHS));
//assign Hs
  estH(1,estHcol) = Hs;
    
//Ht calculation begins 
  double meanalleleHT = 0;
  
  for(int i = 0; i < sizeuniqueLD; i++){
    meanalleleHT +=  (pow(colmean(allelestratap, uniqueLD(i), sizeuniqueS),2));
  }
  meanalleleHT = 1 - meanalleleHT;

  double harmonicNS = harmonicmeanHS * (sizedata*2);
  
  double Ht = meanalleleHT + (Hs / harmonicNS) - ((Ho / 2 ) / harmonicNS);
    estH(2,estHcol) = Ht;
        
    estHcol++;    
  }
  return(estH);
}

// [[Rcpp::export]]
  double colmean(NumericMatrix data, int columnToAverage, int nrow){
    double sum = 0;
    for(int i = 0; i < nrow; i++){
      sum += data(i, columnToAverage);
    }
    double mean = sum/nrow;
    return(mean);
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

