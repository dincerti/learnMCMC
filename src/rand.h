# ifndef RTNORM_H
# define RTNORM_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double rtnormInv(const double lower,
            const double upper,
            const double mean,
            const double sd) ;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

# endif
