#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _JointGWAS_matrix_multiply_toDense(SEXP, SEXP);
extern SEXP _JointGWAS_partial_matrix_multiply_toDense(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_JointGWAS_matrix_multiply_toDense",         (DL_FUNC) &_JointGWAS_matrix_multiply_toDense,         2},
    {"_JointGWAS_partial_matrix_multiply_toDense", (DL_FUNC) &_JointGWAS_partial_matrix_multiply_toDense, 3},
    {NULL, NULL, 0}
};

void R_init_JointGWAS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
