// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// intraLaplacian
NumericMatrix intraLaplacian(NumericVector mod);
RcppExport SEXP nettools_intraLaplacian(SEXP modSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type mod(modSEXP );
        NumericMatrix __result = intraLaplacian(mod);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// directProd
NumericMatrix directProd(NumericMatrix Li, int n);
RcppExport SEXP nettools_directProd(SEXP LiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type Li(LiSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        NumericMatrix __result = directProd(Li, n);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// IntraAdj
NumericMatrix IntraAdj(NumericVector mod);
RcppExport SEXP nettools_IntraAdj(SEXP modSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type mod(modSEXP );
        NumericMatrix __result = IntraAdj(mod);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}