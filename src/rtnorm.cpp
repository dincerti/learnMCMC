#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rtnormInv(const double lower = -5, const double upper = 5,
                  const double mean= 0, const double sd = 1) {
  double pl, pu, u, rand;
  if (lower == R_NegInf){
      double pl = 0;
  }
  else {
    pl = R::pnorm(lower, mean, sd, 1, 0);
  }
  if (upper == R_PosInf){
      pu = 1;
  }
  else {
    pu = R::pnorm(upper, mean, sd, 1, 0);
  }
  u =  R::runif(pl, pu);
  rand = R::qnorm(u, mean, sd, 1, 0);
  return rand;
}




