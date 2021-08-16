#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP IHSEP_mloglik1c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP IHSEP_mloglik1d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP IHSEP_mloglik1e(SEXP, SEXP, SEXP, SEXP, SEXP);

/* .External calls */
extern SEXP mloglik1a(SEXP);
extern SEXP mloglik1b(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"IHSEP_mloglik1c", (DL_FUNC) &IHSEP_mloglik1c, 5},
    {"IHSEP_mloglik1d", (DL_FUNC) &IHSEP_mloglik1d, 5},
    {"IHSEP_mloglik1e", (DL_FUNC) &IHSEP_mloglik1e, 5},
    {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalEntries[] = {
    {"mloglik1a", (DL_FUNC) &mloglik1a, 10},
    {"mloglik1b", (DL_FUNC) &mloglik1b,  8},
    {NULL, NULL, 0}
};

void R_init_IHSEP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
}
