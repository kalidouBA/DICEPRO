#include <RcppEigen.h>
#include <cmath>
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
double compute_obj_term1_eigen_fast(const Eigen::MatrixXd& W,
                                    const Eigen::MatrixXd& H,
                                    const Eigen::MatrixXd& B,
                                    double sigma) {
  // Compute R = W * H.transpose() - B
  Eigen::MatrixXd R = W * H.transpose() - B;

  double sum_val = 0.0;
  const int n = R.size();
  // Get a pointer to the data of R for fast iteration.
  const double* r_ptr = R.data();

  // Loop through all elements of R without creating an extra matrix
  for (int i = 0; i < n; i++) {
    double asinh_val = std::asinh(r_ptr[i]);
    sum_val += asinh_val * asinh_val;
  }

  return sum_val / (2.0 * sigma * sigma);
}
