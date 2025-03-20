#' NMF with L-BFGS-B Optimization
#'
#' Perform non-negative matrix factorization (NMF) using L-BFGS-B optimization. This method uses an iterative approach
#' to solve a regularized least-squares problem for the factorization of a matrix into non-negative matrices W and H.
#' It incorporates a regularization term to ensure that the solution meets certain constraints, such as cell type composition.
#'
#' @param r_dataset A list containing the dataset for the factorization, including:
#'   - `B`: A matrix of size (N_gene x N_sample), the input matrix to factorize (gene expression data).
#'   - `P_cb`: A matrix of size (N_sample x N_cellsType), cell type proportions for the samples.
#'   - `W_cb`: A matrix of size (N_gene x N_cellsType), initial matrix W (gene-cell type matrix).
#' @param W_prime Initial value for W (optional, defaults to 0). It is a scalar.
#' @param p_prime Initial value for p_cb (optional, defaults to 0). It is a scalar.
#' @param lambda_ Regularization parameter (optional, defaults to 10). It is a scalar.
#' @param gamma_par Penalty factor for the sum of cell type proportions (optional, defaults to 100). It is a scalar.
#' @param N_unknownCT Number of unknown cell types to model (optional, defaults to 1). It is a scalar.
#' @param con A list of control parameters for the optimization process passed to `optim` (optional).
#'   - `fnscale`: Scaling factor for the objective function (optional, defaults to 1).
#'   - `maxit`: Maximum number of iterations (optional, defaults to 5000).
#'   - `tmax`: Maximum number of function evaluations (optional, defaults to 50).
#'
#' @return A list containing:
#'   - `w`: The optimized matrix W, a vector of length N_gene.
#'   - `H`: The optimized matrix H (cell type proportions), a matrix of size (N_sample x N_cellsType).
#'   - `loss`: The final value of the objective function (scalar).
#'   - `constraint`: The final value of the constraint term (penalty for sum of cell type proportions), a scalar.
#'
#' @details This function applies the L-BFGS-B optimization method to solve the NMF problem with a regularized objective
#'   function. The regularization includes a penalty for the sum of cell type proportions and a scaling factor.
#'   The optimization minimizes the objective function subject to constraints on the cell type proportions.
#'
#' @importFrom stats optim
#'
#' @export

nmf_lbfgsb <- function(r_dataset, W_prime = 0, p_prime = 0, lambda_ = 10, gamma_par = 100,
                       N_unknownCT = 1, con = list(fnscale = 1, maxit = 5e3, tmax = 50)) {

  B <- r_dataset$B
  p_cb <- r_dataset$P_cb
  W_cb <- r_dataset$W_cb
  N_gene <- nrow(B); N_sample <- ncol(B); N_cellsType <- ncol(W_cb) + N_unknownCT

  p_cb_C <- matrix(p_prime, nrow = N_sample, ncol = N_unknownCT)
  W_cb_C <- matrix(W_prime, nrow = N_gene, ncol = N_unknownCT)
  colnames(p_cb_C) <- colnames(W_cb_C) <- paste0("Uknown_", 1:N_unknownCT)

  H <- cbind(p_cb, p_cb_C)
  H <- H / rowSums(H)

  W <- cbind(W_cb, W_cb_C)

  residuals_vec <- as.vector(tcrossprod(W,H) - B)
  transformed_residuals_vec <- asinh(residuals_vec)
  sigma_par <- sd(transformed_residuals_vec)
  theta <- c(W_cb_C, H, sigma_par)

  lambda_par <- rep(lambda_, N_sample)
  IterateIndex <- 1
  listH <- list()
  all_loss <- list()

  obj_fun <- function(theta) {
    W[, N_cellsType] <- theta[1:N_gene]
    H <- matrix(theta[(N_gene + 1):(length(theta) - 1)], nrow = N_sample, ncol = N_cellsType)
    listH[[IterateIndex]] <<- theta[(N_gene + 1):(length(theta) - 1)]
    sigma_par <- theta[length(theta)]
    obj_term1 <- compute_obj_term1_arma(W,H,as.matrix(B),sigma_par)
    obj_term2 <- (N_sample * N_gene / 2) * log(2 * pi * sigma_par^2)

    h_H <- rowSums(H) - 1
    obj_term3 <- lambda_par %*% h_H
    obj_term4 <- (gamma_par / 2) * sum(h_H^2)
    obj_value <- obj_term1 + obj_term2
    penalty <- obj_term3 + obj_term4
    f_val <- obj_value + penalty
    all_loss[[IterateIndex]] <<- f_val
    IterateIndex <<- IterateIndex + 1
    return(f_val)
  }

  grad_obj_fun <- function(theta) {
    W[, N_cellsType] <- theta[1:N_gene]
    H <- matrix(theta[(N_gene + 1):(length(theta) - 1)], nrow = N_sample, ncol = N_cellsType)
    sigma_par <- theta[length(theta)]

    residual <- tcrossprod(W,H) - B

    asr <- asinh(residual)
    temp <- asr / sqrt(1 + residual^2) / sigma_par^2
    grad_W <- temp %*% H
    grad_H <- crossprod(temp, W)

    h_H <- rowSums(H) - 1
    grad_H <- grad_H + lambda_par + gamma_par * h_H

    grad_sigma <- (-1 / sigma_par^3) * sum(asr^2) + (N_sample * N_gene) / sigma_par

    grad <- c(as.vector(grad_W[, N_cellsType]), grad_H, grad_sigma)

    return(grad)
  }

  result <- optim(par = theta, fn = obj_fun, gr = grad_obj_fun,
                  lower = c(rep(0, length(theta) - 1), 1e-6),
                  upper = c(rep(Inf, N_gene * N_cellsType), rep(1, N_sample * N_cellsType), Inf),
                  control = con, method = "L-BFGS-B")

  if (is.null(result$par) || !is.finite(result$value)) {
    return(NULL)
  }

  theta <- result$par
  W_opt <- W
  W_opt[, N_cellsType] <- theta[1:N_gene]
  H_opt <- matrix(theta[(N_gene + 1):(length(theta) - 1)], nrow = N_sample, ncol = N_cellsType)

  dimnames(H_opt) <- dimnames(H)
  h_H <- abs(rowSums(H_opt) - 1)
  constraints = sum(h_H)

  H <- as.data.frame(H_opt[, 1:(N_cellsType - N_unknownCT)])
  p_prime <- H_opt[, N_unknownCT]
  results <- list(w = as.vector(W_opt[, N_cellsType]), H = H, p_prime = p_prime, loss = result$value, constraint = 1 - constraints)
  return(results)
}
