// RegisteringDynamic Symbols

//#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _strataG_calcStrataN(SEXP, SEXP, SEXP);
extern SEXP _strataG_colMeanC(SEXP);
extern SEXP _strataG_colSumC(SEXP);
extern SEXP _strataG_getMaxInt(SEXP);
extern SEXP _strataG_HoCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _strataG_HsCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _strataG_HtCalc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _strataG_Hstats_C(SEXP, SEXP);
extern SEXP _strataG_idGenotype(SEXP, SEXP, SEXP);
extern SEXP _strataG_idStart(SEXP, SEXP);
extern SEXP _strataG_intOuterC(SEXP, SEXP);
extern SEXP _strataG_numOuterC(SEXP, SEXP);
extern SEXP _strataG_rowSumC(SEXP);
extern SEXP _strataG_statChi2_C(SEXP, SEXP);
extern SEXP _strataG_statFis_C(SEXP, SEXP);
extern SEXP _strataG_statFst_C(SEXP, SEXP);
extern SEXP _strataG_statFstPrime_C(SEXP, SEXP);
extern SEXP _strataG_statGst_C(SEXP, SEXP);
extern SEXP _strataG_statGstDblPrime_C(SEXP, SEXP);
extern SEXP _strataG_statGstPrime_C(SEXP, SEXP, SEXP);
extern SEXP _strataG_statJostD_C(SEXP, SEXP);
extern SEXP _strataG_statPhist_C(SEXP, SEXP, SEXP);
extern SEXP _strataG_table2D(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_strataG_calcStrataN",       (DL_FUNC) &_strataG_calcStrataN,       3},
  {"_strataG_colMeanC",          (DL_FUNC) &_strataG_colMeanC,          1},
  {"_strataG_colSumC",           (DL_FUNC) &_strataG_colSumC,           1},
  {"_strataG_getMaxInt",         (DL_FUNC) &_strataG_getMaxInt,         1},
  {"_strataG_HoCalc",            (DL_FUNC) &_strataG_HoCalc,            5},
  {"_strataG_HsCalc",            (DL_FUNC) &_strataG_HsCalc,            5},
  {"_strataG_HtCalc",            (DL_FUNC) &_strataG_HtCalc,            6},
  {"_strataG_Hstats_C",          (DL_FUNC) &_strataG_Hstats_C,          2},
  {"_strataG_idGenotype",        (DL_FUNC) &_strataG_idGenotype,        3},
  {"_strataG_idStart",           (DL_FUNC) &_strataG_idStart,           2},
  {"_strataG_intOuterC",         (DL_FUNC) &_strataG_intOuterC,         2},
  {"_strataG_numOuterC",         (DL_FUNC) &_strataG_numOuterC,         2},
  {"_strataG_rowSumC",           (DL_FUNC) &_strataG_rowSumC,           1},
  {"_strataG_statChi2_C",        (DL_FUNC) &_strataG_statChi2_C,        2},
  {"_strataG_statFis_C",         (DL_FUNC) &_strataG_statFis_C,         2},
  {"_strataG_statFst_C",         (DL_FUNC) &_strataG_statFst_C,         2},
  {"_strataG_statFstPrime_C",    (DL_FUNC) &_strataG_statFstPrime_C,    2},
  {"_strataG_statGst_C",         (DL_FUNC) &_strataG_statGst_C,         2},
  {"_strataG_statGstDblPrime_C", (DL_FUNC) &_strataG_statGstDblPrime_C, 2},
  {"_strataG_statGstPrime_C",    (DL_FUNC) &_strataG_statGstPrime_C,    3},
  {"_strataG_statJostD_C",       (DL_FUNC) &_strataG_statJostD_C,       2},
  {"_strataG_statPhist_C",       (DL_FUNC) &_strataG_statPhist_C,       3},
  {"_strataG_table2D",           (DL_FUNC) &_strataG_table2D,           2},
  {NULL, NULL, 0}
};

void R_init_strataG(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
