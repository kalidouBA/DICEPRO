#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd compute_grad_eigen_fast(const Eigen::MatrixXd& W,
                                            const Eigen::MatrixXd& H,
                                            const Eigen::MatrixXd& B,
                                            double sigma,
                                            const Eigen::VectorXd& lambda,
                                            double gamma) {
  // Dimensions:
  int n_gene = W.rows();    // rows of W
  int n_cells = W.cols();   // columns of W (and H)
  int n_sample = H.rows();  // rows of H (samples)

  double sigma2 = sigma * sigma;
  double sigma3 = sigma2 * sigma;

  // Compute the residual matrix: W * H^T - B
  //   Dimensions: n_gene x n_sample
  Eigen::MatrixXd res = W * H.transpose() - B;

  // Compute element-wise asinh(res)
  Eigen::MatrixXd asr = res.unaryExpr([](double x){ return std::asinh(x); });

  // Compute element-wise denominator: sqrt(1 + res^2)
  Eigen::MatrixXd denom = (1.0 + res.array().square()).sqrt().matrix();

  // Compute temp = asr / (denom * sigma^2)
  Eigen::MatrixXd temp = asr.array() / (denom.array() * sigma2);

  // Compute grad_sigma = (-1/sigma^3)*sum(asr^2) + (n_gene*n_sample)/sigma
  double grad_sigma = (-1.0 / sigma3) * asr.array().square().sum() + (n_gene * n_sample) / sigma;

  // Compute grad_W_last = temp * H.col(last)
  //   Dimensions: n_gene x 1 (since temp is n_gene x n_sample and H.col(last) is n_sample x 1)
  Eigen::VectorXd grad_W_last = temp * H.col(n_cells - 1);

  // Compute grad_H = temp^T * W
  //   Dimensions: n_sample x n_cells
  Eigen::MatrixXd grad_H = temp.transpose() * W;

  // Compute per-sample adjustments:
  // For each sample j: adjustment = lambda(j) + gamma*(sum(H.row(j)) - 1)
  Eigen::VectorXd H_row_sum = H.rowwise().sum();
  Eigen::VectorXd adjust = lambda + gamma * (H_row_sum.array() - 1.0).matrix();

  // Add the adjustment to each row of grad_H (broadcasting the adjustment vector)
  grad_H = grad_H.colwise() + adjust;

  // Build final gradient vector by concatenating:
  // 1. grad_W_last (length n_gene)
  // 2. grad_H flattened in column-major order (length n_sample*n_cells)
  // 3. grad_sigma (scalar)
  int total_len = n_gene + n_sample * n_cells + 1;
  Eigen::VectorXd grad(total_len);

  // Place grad_W_last at the beginning
  grad.segment(0, n_gene) = grad_W_last;

  // Flatten grad_H in column-major order and insert
  Eigen::Map<Eigen::VectorXd> grad_H_vec(grad_H.data(), grad_H.size());
  grad.segment(n_gene, n_sample * n_cells) = grad_H_vec;

  // Set the final scalar grad_sigma
  grad(total_len - 1) = grad_sigma;

  return grad;
}
