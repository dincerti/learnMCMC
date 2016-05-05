// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rand.h"
using namespace Rcpp;

// // Simulate from multivariate normal distribution
// // [[Rcpp::export]]
// arma::mat mvrnormArma2(int n, arma::vec mu, arma::mat sigma) {
//   int ncols = sigma.n_cols;
//   arma::mat Y = arma::randn(n, ncols);
//   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
// }

// [[Rcpp::export]]
List bprobitC(int sims, arma::vec y, arma::mat X, arma::vec b0, arma::mat B0inv,
          arma::vec beta_start){
  // initialize
  int n = y.n_elem;
  int k = X.n_cols;
  arma::vec ystar(n);
  arma::mat beta_store(sims, k);
  arma::mat b_store(sims, k);
  arma::mat ystar_store(n, sims);

  // starting values
  arma::mat beta(k, 1);
  beta.col(0) = beta_start;

  // Crossproducts
  arma::mat B = arma::inv(B0inv + arma::trans(X) * X);

  // Gibbs Sampler
  for (int i = 0; i < sims; ++i){
    // Draw ystar
    for (int j = 0; j < n; j++){
      if(y(j) == 0){
        ystar(j) = rtnormInv(R_NegInf, 0, arma::dot(X.row(j), beta), 1);
      }
      else if (y(j) == 1){
        ystar(j) = rtnormInv(0, R_PosInf, arma::dot(X.row(j), beta), 1);
      }
    }

    // Draw beta
    arma::vec b = B * (B0inv * b0  + arma::trans(X) * ystar);
    beta = mvrnormArma(1, b, B);

    // Store
    beta_store.row(i) = beta;
    ystar_store.col(i) = ystar;
    b_store.row(i) = arma::trans(b);
  }

  // Return
  return List::create(Named("beta") = beta_store,
                      Named("ystar") = ystar_store,
                      Named("b") = b_store);
}

// [[Rcpp::export]]
arma::mat simulate(int sims, arma::vec y, arma::vec xb){
  int n = y.n_elem;
  arma::vec ystar(n);
  arma::mat ystar_store(n, sims);
  for (int i = 0; i < sims; ++i){
    for (int j = 0; j < n; j++){
      if(y(j) == 0){
        ystar(j) = rtnormInv(R_NegInf, 0, xb(j), 1);
      }
      else if (y(j) == 1){
        ystar(j) = rtnormInv(0, R_PosInf, xb(j), 1);
      }
    }
  ystar_store.col(i) = ystar;
  }
  return ystar_store;
}
