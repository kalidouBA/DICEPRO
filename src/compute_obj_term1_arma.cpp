#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
double compute_obj_term1_arma(const arma::mat& W,
                              const arma::mat& H,
                              const arma::mat& B,
                              double sigma) {
  // Compute tcrossprod equivalent: W * H.t()  (dimensions: nrow(W) x nrow(H))
  arma::mat R = W * H.t() - B;
  // Compute elementwise asinh, square, and sum the result.
  // Note: arma::asinh() is available in Armadillo.
  double sum_val = arma::accu(arma::asinh(R) % arma::asinh(R));
  return sum_val / (2.0 * sigma * sigma);
}
