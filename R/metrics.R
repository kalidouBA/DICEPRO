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
#'  result <- computPerf(outDec_1, outDec_2)
#'}

computPerf <- function(outDec_1, outDec_2, metric) {

  # Normalize outDec_1 and outDec_2 matrices to sum to 1 along each row
  outDec_1 <- sweep(outDec_1, 1, rowSums(outDec_1), FUN = "/")
  outDec_2 <- sweep(outDec_2, 1, rowSums(outDec_2), FUN = "/")

  CellType <- colnames(outDec_1)
  # Initialize a data frame to store performance metrics
  perfAll <- NULL

  # Iterate over each cell type
  for (ct in CellType) {
    perf <- c()
    x <- as.vector(outDec_1[, ct])
    y <- as.vector(outDec_2[, ct])


    # Calculate Root Mean Squared Error Relative RMSE (RRMSE)
    if(metric == "RRMSE"){
      num = sum((x - y)^2)
      den = sum(y^2)
      squared_error = num/den
      rrmse_loss = sqrt(squared_error)
      rrmse_loss <- round(rrmse_loss,2)
      perf <- append(perf, rrmse_loss)
    }

    # Fit a linear model to calculate R-squared and Adjusted R-squared
    else{
      model <- lm(y ~ x)
      R2 <- summary(model)$r.squared
      R2_adj <- summary(model)$adj.r.squared
      perf <- append(perf, R2_adj)
    }
    perfAll <- rbind(perfAll,perf)
  }

  perfAll <- cbind.data.frame(metric = perfAll, CellType = as.factor(CellType))

  return(perfAll)
}

