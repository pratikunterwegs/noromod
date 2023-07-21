// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// norovirus_model_cpp
Rcpp::List norovirus_model_cpp(const double& t, Eigen::VectorXd& state, Rcpp::List parameters);
RcppExport SEXP _noromod_norovirus_model_cpp(SEXP tSEXP, SEXP stateSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(norovirus_model_cpp(t, state, parameters));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_noromod_norovirus_model_cpp", (DL_FUNC) &_noromod_norovirus_model_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_noromod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
