// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// tabulate
List tabulate(Reference solResponse, Reference solveRequest);
RcppExport SEXP _iceR_tabulate(SEXP solResponseSEXP, SEXP solveRequestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Reference >::type solResponse(solResponseSEXP);
    Rcpp::traits::input_parameter< Reference >::type solveRequest(solveRequestSEXP);
    rcpp_result_gen = Rcpp::wrap(tabulate(solResponse, solveRequest));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iceR_tabulate", (DL_FUNC) &_iceR_tabulate, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_iceR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
