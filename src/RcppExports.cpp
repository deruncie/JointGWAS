// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "JointGWAS_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matrix_multiply_toDense
MatrixXd matrix_multiply_toDense(SEXP X_, SEXP Y_);
RcppExport SEXP _JointGWAS_matrix_multiply_toDense(SEXP X_SEXP, SEXP Y_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y_(Y_SEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_multiply_toDense(X_, Y_));
    return rcpp_result_gen;
END_RCPP
}
// partial_matrix_multiply_toDense
MatrixXd partial_matrix_multiply_toDense(SEXP X_, SEXP Y_, VectorXi j);
RcppExport SEXP _JointGWAS_partial_matrix_multiply_toDense(SEXP X_SEXP, SEXP Y_SEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(partial_matrix_multiply_toDense(X_, Y_, j));
    return rcpp_result_gen;
END_RCPP
}
