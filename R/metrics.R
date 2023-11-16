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
  colnameIntersect <- intersect(colnames(outDec_1), colnames(outDec_2))

  # Normalize outDec_1 and outDec_2 matrices to sum to 1 along each row
  outDec_1 <- sweep(outDec_1[, colnameIntersect], 1, rowSums(outDec_1[, colnameIntersect]), FUN = "/")
  outDec_2 <- sweep(outDec_2[, colnameIntersect], 1, rowSums(outDec_2[, colnameIntersect]), FUN = "/")

  # Initialize a data frame to store performance metrics
  perfAll <- NULL

  # Iterate over each cell type
  for (v_ in colnameIntersect) {
    perf <- c()
    x <- as.vector(outDec_1[, v_])
    y <- as.vector(outDec_2[, v_])


    # Calculate Root Mean Squared Error Relative RMSE (RRMSE)
    if(metric %in% c("RRMSE", "ALL")){
      num = sum((x - y)^2)
      den = sum(y^2)
      squared_error = num/den
      rrmse_loss = sqrt(squared_error)
      rrmse_loss <- round(rrmse_loss,2)
      perf <- append(perf, rrmse_loss)
    }

    # Fit a linear model to calculate R-squared and Adjusted R-squared
    if(metric %in% c("R2_adj", "ALL")){
      model <- lm(y ~ x)
      R2 <- summary(model)$r.squared
      R2_adj <- summary(model)$adj.r.squared
      perf <- append(perf, R2_adj)
    }

    perfAll <- rbind(perfAll, colMeans(perf))
  }

  perfAll <- cbind(perfAll, colnameIntersect)

  return(perfAll)
}

