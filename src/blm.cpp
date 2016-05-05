// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rand.h"
using namespace Rcpp;

// [[Rcpp::export]]
List blmC(int sims, arma::vec y, arma::mat X, arma::vec b0, arma::mat B0inv,
          double sigma2, double c0 = .1, double d0 = .1){
  int n = y.n_elem;
  int k = X.n_cols;
  arma::mat beta_store(sims, k);
  arma::vec sigma2_store(sims);

  // Crossproducts
  arma::mat XtX = arma::trans(X) * X;
  arma::vec Xty = arma::trans(X) * y;

  // Gibbs Sampler
  for(int i = 0; i < sims; ++i){
    // Draw beta
    arma::mat B =  arma::inv(B0inv + XtX/sigma2);
    arma::mat b = B * (B0inv * b0  + Xty/sigma2);
    arma::mat beta = mvrnormArma(1, b, B);

    // Draw sigma2
    arma::vec epsilon = y - X * arma::trans(beta);
    double crossprod = arma::dot(epsilon, epsilon);
    sigma2 = 1/R::rgamma((n + c0)/2, 1/((d0 + crossprod)/2));

    // Store
    beta_store.row(i) = beta;
    sigma2_store(i) = sigma2;
  }

  // Return
  return List::create(Named("beta") = beta_store,
                      Named("sigma2") = sigma2_store);
}

