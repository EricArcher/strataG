// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// indGenotype
IntegerVector indGenotype(int nInd, int numAlleles, IntegerMatrix locus);
RcppExport SEXP strataG_indGenotype(SEXP nIndSEXP, SEXP numAllelesSEXP, SEXP locusSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nInd(nIndSEXP);
    Rcpp::traits::input_parameter< int >::type numAlleles(numAllelesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type locus(locusSEXP);
    __result = Rcpp::wrap(indGenotype(nInd, numAlleles, locus));
    return __result;
END_RCPP
}
// HoCalc
double HoCalc(int nInd, IntegerVector loci, int ploidy, IntegerVector strata, IntegerVector strataN);
RcppExport SEXP strataG_HoCalc(SEXP nIndSEXP, SEXP lociSEXP, SEXP ploidySEXP, SEXP strataSEXP, SEXP strataNSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nInd(nIndSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strataN(strataNSEXP);
    __result = Rcpp::wrap(HoCalc(nInd, loci, ploidy, strata, strataN));
    return __result;
END_RCPP
}
// HsCalc
double HsCalc(NumericMatrix alleleFreq, int ploidy, IntegerVector strataN, double harmN, double Ho);
RcppExport SEXP strataG_HsCalc(SEXP alleleFreqSEXP, SEXP ploidySEXP, SEXP strataNSEXP, SEXP harmNSEXP, SEXP HoSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type alleleFreq(alleleFreqSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strataN(strataNSEXP);
    Rcpp::traits::input_parameter< double >::type harmN(harmNSEXP);
    Rcpp::traits::input_parameter< double >::type Ho(HoSEXP);
    __result = Rcpp::wrap(HsCalc(alleleFreq, ploidy, strataN, harmN, Ho));
    return __result;
END_RCPP
}
// Hstats_C
NumericMatrix Hstats_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_Hstats_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(Hstats_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// getMaxInt
int getMaxInt(IntegerVector x);
RcppExport SEXP strataG_getMaxInt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    __result = Rcpp::wrap(getMaxInt(x));
    return __result;
END_RCPP
}
// table2D
IntegerMatrix table2D(IntegerVector x, IntegerVector y);
RcppExport SEXP strataG_table2D(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    __result = Rcpp::wrap(table2D(x, y));
    return __result;
END_RCPP
}
// rowSumC
NumericVector rowSumC(NumericMatrix x);
RcppExport SEXP strataG_rowSumC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(rowSumC(x));
    return __result;
END_RCPP
}
// colSumC
NumericVector colSumC(NumericMatrix x);
RcppExport SEXP strataG_colSumC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(colSumC(x));
    return __result;
END_RCPP
}
// colMeanC
NumericVector colMeanC(NumericMatrix x);
RcppExport SEXP strataG_colMeanC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(colMeanC(x));
    return __result;
END_RCPP
}
// intOuterC
IntegerMatrix intOuterC(IntegerVector x, IntegerVector y);
RcppExport SEXP strataG_intOuterC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    __result = Rcpp::wrap(intOuterC(x, y));
    return __result;
END_RCPP
}
// numOuterC
NumericMatrix numOuterC(NumericVector x, NumericVector y);
RcppExport SEXP strataG_numOuterC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(numOuterC(x, y));
    return __result;
END_RCPP
}
// intVecToMat
IntegerMatrix intVecToMat(IntegerVector x, int ncol);
RcppExport SEXP strataG_intVecToMat(SEXP xSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    __result = Rcpp::wrap(intVecToMat(x, ncol));
    return __result;
END_RCPP
}
// numVecToMat
NumericMatrix numVecToMat(NumericVector x, int ncol);
RcppExport SEXP strataG_numVecToMat(SEXP xSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    __result = Rcpp::wrap(numVecToMat(x, ncol));
    return __result;
END_RCPP
}
// calcStrataN
IntegerVector calcStrataN(IntegerVector locus, IntegerVector strata);
RcppExport SEXP strataG_calcStrataN(SEXP locusSEXP, SEXP strataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type locus(locusSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    __result = Rcpp::wrap(calcStrataN(locus, strata));
    return __result;
END_RCPP
}
// statChi2_C
double statChi2_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statChi2_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statChi2_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// statFis_C
double statFis_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statFis_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statFis_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// alleleFreqCalc
NumericMatrix alleleFreqCalc(IntegerVector locVec, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_alleleFreqCalc(SEXP locVecSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type locVec(locVecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(alleleFreqCalc(locVec, strata, ploidy));
    return __result;
END_RCPP
}
// prHetCalc
NumericMatrix prHetCalc(IntegerVector alleles, IntegerVector nvec, IntegerMatrix locusMat, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_prHetCalc(SEXP allelesSEXP, SEXP nvecSEXP, SEXP locusMatSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type alleles(allelesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvec(nvecSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type locusMat(locusMatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(prHetCalc(alleles, nvec, locusMat, strata, ploidy));
    return __result;
END_RCPP
}
// varCompCalc
NumericMatrix varCompCalc(IntegerVector nvec, NumericMatrix alleleFreq, NumericMatrix prHet, int r, double nbar, double rnbar, double nc);
RcppExport SEXP strataG_varCompCalc(SEXP nvecSEXP, SEXP alleleFreqSEXP, SEXP prHetSEXP, SEXP rSEXP, SEXP nbarSEXP, SEXP rnbarSEXP, SEXP ncSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type nvec(nvecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alleleFreq(alleleFreqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type prHet(prHetSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type nbar(nbarSEXP);
    Rcpp::traits::input_parameter< double >::type rnbar(rnbarSEXP);
    Rcpp::traits::input_parameter< double >::type nc(ncSEXP);
    __result = Rcpp::wrap(varCompCalc(nvec, alleleFreq, prHet, r, nbar, rnbar, nc));
    return __result;
END_RCPP
}
// statFst_C
double statFst_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statFst_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statFst_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// statGst_C
double statGst_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statGst_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statGst_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// statGstPrime_C
double statGstPrime_C(IntegerMatrix loci, IntegerVector strata, int ploidy, int primeType);
RcppExport SEXP strataG_statGstPrime_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP, SEXP primeTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type primeType(primeTypeSEXP);
    __result = Rcpp::wrap(statGstPrime_C(loci, strata, ploidy, primeType));
    return __result;
END_RCPP
}
// statGstDblPrime_C
double statGstDblPrime_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statGstDblPrime_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statGstDblPrime_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// statJostD_C_test
double statJostD_C_test(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statJostD_C_test(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statJostD_C_test(loci, strata, ploidy));
    return __result;
END_RCPP
}
// statJostD_C
double statJostD_C(IntegerMatrix loci, IntegerVector strata, int ploidy);
RcppExport SEXP strataG_statJostD_C(SEXP lociSEXP, SEXP strataSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    __result = Rcpp::wrap(statJostD_C(loci, strata, ploidy));
    return __result;
END_RCPP
}
// ssWPCalc
double ssWPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq, NumericMatrix hapDist);
RcppExport SEXP strataG_ssWPCalc(SEXP strataFreqSEXP, SEXP strataHapFreqSEXP, SEXP hapDistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type strataFreq(strataFreqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type strataHapFreq(strataHapFreqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hapDist(hapDistSEXP);
    __result = Rcpp::wrap(ssWPCalc(strataFreq, strataHapFreq, hapDist));
    return __result;
END_RCPP
}
// ssAPCalc
double ssAPCalc(IntegerVector strataFreq, IntegerMatrix strataHapFreq, NumericMatrix hapDist);
RcppExport SEXP strataG_ssAPCalc(SEXP strataFreqSEXP, SEXP strataHapFreqSEXP, SEXP hapDistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type strataFreq(strataFreqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type strataHapFreq(strataHapFreqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hapDist(hapDistSEXP);
    __result = Rcpp::wrap(ssAPCalc(strataFreq, strataHapFreq, hapDist));
    return __result;
END_RCPP
}
// statPhist_C
double statPhist_C(IntegerVector haps, IntegerVector strata, NumericMatrix hapDist);
RcppExport SEXP strataG_statPhist_C(SEXP hapsSEXP, SEXP strataSEXP, SEXP hapDistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type haps(hapsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hapDist(hapDistSEXP);
    __result = Rcpp::wrap(statPhist_C(haps, strata, hapDist));
    return __result;
END_RCPP
}
