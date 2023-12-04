#' NMF with Projected Conjugate Gradient Optimization
#'
#' Perform non-negative matrix factorization (NMF) using projected conjugate
#' gradient optimization. Define the NMF function with projected conjugate gradient optimization
#'
#' @param V Input matrix for factorization.
#' @param k Number of components.
#' @param W Initial matrix W (optional, defaults to random non-negative values).
#' @param H Initial matrix H (optional, defaults to random non-negative values).
#' @param max_iter Maximum number of iterations (default is 100).
#' @return A list containing factorized matrices W and H.
#'
#' @importFrom stats optim
#'
#' @examples
#' set.seed(123)
#' data_matrix <- matrix(runif(100), nrow = 10, ncol = 10)
#' k <- 2
#' result <- nmf_conjugate_gradient(data_matrix, k)
#' print(result$W)
#' print(result$H)
#'
#' @export

nmf_conjugate_gradient <- function(V, k, W = NULL, H = NULL, max_iter = 100) {

  # Initialize W and H with random non-negative values
  if(is.null(W) || is.null(H)){
    W <- matrix(runif(nrow(V) * k), nrow = nrow(V), ncol = k)
    H <- matrix(runif(k * ncol(V)), nrow = k, ncol = ncol(V))
  }
  # Flatten matrices for use in optimization
  theta <- c(W, H)

  # Define the objective function and its gradient
  obj_fun <- function(theta) {
    W <- matrix(theta[1:(nrow(V) * k)], nrow = nrow(V), ncol = k)
    H <- matrix(theta[(nrow(V) * k + 1):length(theta)], nrow = k, ncol = ncol(V))

    # Calculate Frobenius norm squared
    obj_value <- sum((W %*% H - V)^2)

    return(obj_value)
  }

  grad_obj_fun <- function(theta) {
    W <- matrix(theta[1:(nrow(V) * k)], nrow = nrow(V), ncol = k)
    H <- matrix(theta[(nrow(V) * k + 1):length(theta)], nrow = k, ncol = ncol(V))

    # Calculate gradients
    grad_W <- 2 * as.matrix(W %*% H - V) %*% t(H)
    grad_H <- 2 * t(W) %*% as.matrix(W %*% H - V)

    return(c(grad_W, grad_H))
  }

  # Define the projection function to enforce non-negativity
  project_to_non_negative <- function(matrix) {
    matrix[matrix < 0] <- 0
    return(matrix)
  }

  # Perform projected conjugate gradient optimization
  result <- optim(par = theta, fn = obj_fun, gr = grad_obj_fun, method = "L-BFGS-B", control = list(maxit = max_iter))

  # Extract factorized matrices W and H and enforce non-negativity
  # W <- project_to_non_negative(matrix(result$par[1:(nrow(V) * k)], nrow = nrow(V), ncol = k))
  # H <- project_to_non_negative(matrix(result$par[(nrow(V) * k + 1):length(result$par)], nrow = k, ncol = ncol(V)))

  W <- matrix(result$par[1:(nrow(V) * k)], nrow = nrow(V), ncol = k)
  H <- matrix(result$par[(nrow(V) * k + 1):length(result$par)], nrow = k, ncol = ncol(V))

  n <- length(V)
  squared_diff <- (V - W%*%H)^2
  mse <- sum(squared_diff) / n
  rmse_nmf <- sqrt(mse)

  # Return the factorized matrices W and H and rmse_nmf
  list(W = W, H = H, Error = rmse_nmf)
}
