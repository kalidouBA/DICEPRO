#' NMF with Projected Conjugate Gradient Optimization
#'
#' Perform non-negative matrix factorization (NMF) using projected conjugate
#' gradient optimization. Define the NMF function with projected conjugate gradient optimization
#'
#' @param V Input matrix for factorization.
#' @param k Number of components. Default is 1
#' @param W Initial matrix W (optional, defaults to random non-negative values).
#' @param H Initial matrix H (optional, defaults to random non-negative values).
#' @return A list containing factorized matrices W and H.
#'
#' @importFrom stats optim
#' @importFrom parallel detectCores
#'
#' @export

nmf_conjugate_gradient <- function(V, k = 1, W = NULL, H = NULL) {

  # Initialize W and H with random non-negative values
  if(is.null(W) || is.null(H)){
    W <- matrix(runif(nrow(V) * k), nrow = nrow(V), ncol = k)
    H <- matrix(runif(k * ncol(V)), nrow = k, ncol = ncol(V))
  }

  par_nmf <- function(i){
    W <- as.matrix(cbind.data.frame(W, unk = 0))
    H <- as.matrix(cbind.data.frame(H, unk = i))
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

    # Perform projected conjugate gradient optimization
    result <- optim(par = theta, fn = obj_fun, gr = grad_obj_fun, method = "L-BFGS-B",
                    lower = rep(0, length(theta)), control = list(maxit = 1000, trace = FALSE, factr = 1e-8))

    W <- matrix(result$par[1:(nrow(V) * k)], nrow = nrow(V), ncol = k)
    H <- matrix(result$par[(nrow(V) * k + 1):length(result$par)], nrow = k, ncol = ncol(V))

    n <- length(V)
    squared_diff <- (V - W%*%H)^2
    mse <- sum(squared_diff) / n
    rmse_nmf <- sqrt(mse)

    list(W = W, H = t(H), Error = rmse_nmf)
  }
 cl <-  min(5, detectCores())
  result_list <- mclapply(seq(0,1,0.25), function(i) {
    par_nmf(i)
  }, mc.cores = cl)

  ind <- which.min(lapply(result_list, `[[`, 3))
  result_list[[ind]]
}
