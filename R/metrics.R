#' Compute Performance Metrics for Comparing Two Data Matrices of abundances
#'
#' This function computes performance metrics to compare two data matrices.
#' The function calculates various performance metrics, including adjusted R-squared (R2_adj),
#' Root Mean Square Error Relative RMSE (RRMSE).
#'
#' @param outDec_1 A matrix of cell type proportions from the i iteration
#' @param outDec_2 A matrix of cell type proportions from the i+1 iteration.
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
#'  outDec_1 <- matrix(rnorm(50), ncol = 5)
#'  outDec_2 <- matrix(rnorm(50), ncol = 5)
#'  colnames(outDec_1) <- colnames(outDec_2) <- paste0("C_", 1:5)
#'  result <- computPerf(outDec_1, outDec_2, "RRMSE")
#'}

computPerf <- function(outDec_1, outDec_2, metric) {

  # Normalize outDec_1 and outDec_2 matrices to sum to 1 along each row
  outDec_1 <- sweep(outDec_1, 1, rowSums(outDec_1), FUN = "/")
  outDec_2 <- sweep(outDec_2, 1, rowSums(outDec_2), FUN = "/")
  x <- melt(outDec_1, id.vars=NULL)$value
  y <- melt(outDec_2, id.vars=NULL)$value

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

  return(data.frame(metric = perf))
}

