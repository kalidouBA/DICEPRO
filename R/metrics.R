#' Compute Performance Metrics for Comparing Two Data Matrices of abundances
#'
#' This function computes performance metrics to compare two data matrices.
#' The function calculates various performance metrics, including adjusted R-squared (R2_adj),
#' Root Mean Square Error Relative RMSE (RRMSE).
#'
#' @param x A matrix of cell type proportions from the i iteration
#' @param y A matrix of cell type proportions from the i+1 iteration.
#' @param metric A metric select between adjusted R-squared (R2_adj), Root Mean Square Error Relative RMSE (RRMSE).
#' @export
#'
#' @importFrom Metrics rmse
#' @importFrom reshape2 melt
#'
#'
#' @return A data frame containing performance metrics for each variable in the comparison.
#' @details The function calculates and returns the following performance metrics for each variable:
#' \itemize{
#'   \item \code{R2_adj}: Adjusted R-squared.
#'   \item \code{RRMSE}: Relative Root Mean Square Error (RRMSE).
#'}
#'
#' @examples
#' if(interactive()){
#'  x <- matrix(rnorm(50), ncol = 5)
#'  y <- matrix(rnorm(50), ncol = 5)
#'  colnames(x) <- colnames(y) <- paste0("C_", 1:5)
#'  result <- computPerf(x, y, "RRMSE")
#'}

computPerf <- function(x, y, metric) {

  # Calculate Root Mean Squared Error Relative RMSE (RRMSE)
  if(metric == "RRMSE"){
    num = sum((x - y)^2)
    den = sum(y^2)
    squared_error = num/den
    rrmse_loss = sqrt(squared_error)
    rrmse_loss <- round(rrmse_loss,2)
    perf <- rrmse_loss
  }

  # Fit a linear model to calculate R-squared and Adjusted R-squared
  else{
    model <- lm(y ~ x)
    R2 <- summary(model)$r.squared
    R2_adj <- summary(model)$adj.r.squared
    perf <- R2_adj
  }

  return(perf)
}


#' Compute precision of estimated
#'
#' This function computes the precision of estimated by evaluating
#' the Frobenius norm of the difference between the gene expression matrix (Y)
#' and the product of the reference matrix and the estimated abundance matrix (Y_tilde),
#' normalized by the Frobenius norm of the the gene expression matrix (Y).
#'
#' @param Y The gene expression matrix.
#' @param Y_tilde The product of the reference matrix and the estimated abundance matrix.
#' @export
#'
#' @return The precision of sequential simulation.
#'
#' @examples
#' Y <- matrix(c(3, 5, 7, 2, 6, 4, 0, 2, 8), nrow=3, ncol=3, byrow=TRUE)
#' Y_tilde <- Y
#' precision <- compute_precision(Y, Y_tilde)
#' cat("Precision of sequential simulation:", precision, "\n")

compute_precision <- function(Y, Y_tilde) {
  frobenius_norm <- function(matrix) {
    sqrt(sum(matrix^2))
  }

  diff_matrix <- Y - Y_tilde
  norm_diff <- frobenius_norm(diff_matrix)
  norm_Y <- frobenius_norm(Y)

  precision <- norm_diff / norm_Y

  return(precision)
}
