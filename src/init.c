// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP strataG_calcStrataN(SEXP, SEXP, SEXP);
extern SEXP strataG_colMeanC(SEXP);
extern SEXP strataG_colSumC(SEXP);
extern SEXP strataG_getMaxInt(SEXP);
extern SEXP strataG_HoCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP strataG_HsCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP strataG_Hstats_C(SEXP, SEXP);
extern SEXP strataG_idGenotype(SEXP, SEXP, SEXP);
extern SEXP strataG_idStart(SEXP, SEXP);
extern SEXP strataG_intOuterC(SEXP, SEXP);
extern SEXP strataG_numOuterC(SEXP, SEXP);
extern SEXP strataG_rowSumC(SEXP);
extern SEXP strataG_statChi2_C(SEXP, SEXP);
extern SEXP strataG_statFis_C(SEXP, SEXP);
extern SEXP strataG_statFst_C(SEXP, SEXP);
extern SEXP strataG_statFstPrime_C(SEXP, SEXP);
extern SEXP strataG_statGst_C(SEXP, SEXP);
extern SEXP strataG_statGstDblPrime_C(SEXP, SEXP);
extern SEXP strataG_statGstPrime_C(SEXP, SEXP, SEXP);
extern SEXP strataG_statJostD_C(SEXP, SEXP);
extern SEXP strataG_statPhist_C(SEXP, SEXP, SEXP);
extern SEXP strataG_table2D(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"strataG_calcStrataN",       (DL_FUNC) &strataG_calcStrataN,       3},
  {"strataG_colMeanC",          (DL_FUNC) &strataG_colMeanC,          1},
  {"strataG_colSumC",           (DL_FUNC) &strataG_colSumC,           1},
  {"strataG_getMaxInt",         (DL_FUNC) &strataG_getMaxInt,         1},
  {"strataG_HoCalc",            (DL_FUNC) &strataG_HoCalc,            5},
  {"strataG_HsCalc",            (DL_FUNC) &strataG_HsCalc,            5},
  {"strataG_Hstats_C",          (DL_FUNC) &strataG_Hstats_C,          2},
  {"strataG_idGenotype",        (DL_FUNC) &strataG_idGenotype,        3},
  {"strataG_idStart",           (DL_FUNC) &strataG_idStart,           2},
  {"strataG_intOuterC",         (DL_FUNC) &strataG_intOuterC,         2},
  {"strataG_numOuterC",         (DL_FUNC) &strataG_numOuterC,         2},
  {"strataG_rowSumC",           (DL_FUNC) &strataG_rowSumC,           1},
  {"strataG_statChi2_C",        (DL_FUNC) &strataG_statChi2_C,        2},
  {"strataG_statFis_C",         (DL_FUNC) &strataG_statFis_C,         2},
  {"strataG_statFst_C",         (DL_FUNC) &strataG_statFst_C,         2},
  {"strataG_statFstPrime_C",    (DL_FUNC) &strataG_statFstPrime_C,    2},
  {"strataG_statGst_C",         (DL_FUNC) &strataG_statGst_C,         2},
  {"strataG_statGstDblPrime_C", (DL_FUNC) &strataG_statGstDblPrime_C, 2},
  {"strataG_statGstPrime_C",    (DL_FUNC) &strataG_statGstPrime_C,    3},
  {"strataG_statJostD_C",       (DL_FUNC) &strataG_statJostD_C,       2},
  {"strataG_statPhist_C",       (DL_FUNC) &strataG_statPhist_C,       3},
  {"strataG_table2D",           (DL_FUNC) &strataG_table2D,           2},
  {NULL, NULL, 0}
};

void R_init_strataG(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
