#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _circglmbayes_atanLF(SEXP, SEXP);
extern SEXP _circglmbayes_circGLMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _circglmbayes_circQuantile(SEXP, SEXP);
extern SEXP _circglmbayes_computeHDI(SEXP, SEXP);
extern SEXP _circglmbayes_computeHDICirc(SEXP, SEXP);
extern SEXP _circglmbayes_computeMeanDirection(SEXP);
extern SEXP _circglmbayes_computeResultantLength(SEXP);
extern SEXP _circglmbayes_estimateMode(SEXP, SEXP);
extern SEXP _circglmbayes_estimateModeCirc(SEXP, SEXP);
extern SEXP _circglmbayes_invAtanLF(SEXP, SEXP);
extern SEXP _circglmbayes_logProbNormal(SEXP, SEXP, SEXP);
extern SEXP _circglmbayes_rvmc(SEXP, SEXP, SEXP);
extern SEXP _circglmbayes_sampleKappa(SEXP, SEXP);
extern SEXP _circglmbayes_truncCauchyPdf(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_circglmbayes_atanLF",                 (DL_FUNC) &_circglmbayes_atanLF,                  2},
    {"_circglmbayes_circGLMC",               (DL_FUNC) &_circglmbayes_circGLMC,               17},
    {"_circglmbayes_circQuantile",           (DL_FUNC) &_circglmbayes_circQuantile,            2},
    {"_circglmbayes_computeHDI",             (DL_FUNC) &_circglmbayes_computeHDI,              2},
    {"_circglmbayes_computeHDICirc",         (DL_FUNC) &_circglmbayes_computeHDICirc,          2},
    {"_circglmbayes_computeMeanDirection",   (DL_FUNC) &_circglmbayes_computeMeanDirection,    1},
    {"_circglmbayes_computeResultantLength", (DL_FUNC) &_circglmbayes_computeResultantLength,  1},
    {"_circglmbayes_estimateMode",           (DL_FUNC) &_circglmbayes_estimateMode,            2},
    {"_circglmbayes_estimateModeCirc",       (DL_FUNC) &_circglmbayes_estimateModeCirc,        2},
    {"_circglmbayes_invAtanLF",              (DL_FUNC) &_circglmbayes_invAtanLF,               2},
    {"_circglmbayes_logProbNormal",          (DL_FUNC) &_circglmbayes_logProbNormal,           3},
    {"_circglmbayes_rvmc",                   (DL_FUNC) &_circglmbayes_rvmc,                    3},
    {"_circglmbayes_sampleKappa",            (DL_FUNC) &_circglmbayes_sampleKappa,             2},
    {"_circglmbayes_truncCauchyPdf",         (DL_FUNC) &_circglmbayes_truncCauchyPdf,          3},
    {NULL, NULL, 0}
};

void R_init_circglmbayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
