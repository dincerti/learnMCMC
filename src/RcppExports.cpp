// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// blmC
List blmC(int sims, arma::vec y, arma::mat X, arma::vec b0, arma::mat B0inv, double sigma2, double c0, double d0);
RcppExport SEXP learnMCMC_blmC(SEXP simsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP b0SEXP, SEXP B0invSEXP, SEXP sigma2SEXP, SEXP c0SEXP, SEXP d0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type sims(simsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B0inv(B0invSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< double >::type d0(d0SEXP);
    __result = Rcpp::wrap(blmC(sims, y, X, b0, B0inv, sigma2, c0, d0));
    return __result;
END_RCPP
}
// rtnormInv
double rtnormInv(const double lower, const double upper, const double mean, const double sd);
RcppExport SEXP learnMCMC_rtnormInv(SEXP lowerSEXP, SEXP upperSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double >::type sd(sdSEXP);
    __result = Rcpp::wrap(rtnormInv(lower, upper, mean, sd));
    return __result;
END_RCPP
}