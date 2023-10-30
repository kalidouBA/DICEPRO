#' Compute Distances Between Columns of a Matrix
#'
#' This function calculates various distance metrics between columns of a matrix.
#' The supported distance metrics include Euclidean distance, Manhattan distance,
#' Concordance Correlation Coefficient (CCC), Intra-class Correlation (ICC), and
#' cross-product.
#'
#' @param matrix_input The matrix for which distances are to be calculated.
#'
#' @export
#'
#' @importFrom DescTools CCC
#' @importFrom psych ICC
#'
#' @return A list of distance matrices and a vector of fold comparisons.
#' @details The function calculates and returns the following distance metrics:
#' \itemize{
#'   \item \code{Euclidean}: Euclidean distance matrix.
#'   \item \code{Manhattan}: Manhattan distance matrix.
#'   \item \code{ccc}: Concordance Correlation Coefficient (CCC) matrix.
#'   \item \code{icc}: Intra-class Correlation (ICC) matrix.
#'   \item \code{CrossProd}: Cross-product matrix.
#'}
#'
#' @examples
#' matrixdata <- matrix(rnorm(100), ncol = 5)
#' result <- compute_distances(matrixdata)

compute_distances <- function(matrix_input) {
  num_cols <- ncol(matrix_input)
  foldName <- colnames(matrix_input)

  euclidean_dist <- manhattan_dist <- ccc_dist <- icc_dist <- cp_dist <- matrix(0, nrow = num_cols, ncol = num_cols)
  alComparisionFold <- c()
  for (i in 1:(num_cols-1)) {
    for (j in (i+1):num_cols) {
      col1 <- matrix_input[, i]
      col2 <- matrix_input[, j]

      # Calculate Euclidean Distance
      euclidean_dist[i, j] <- sqrt(sum((col1 - col2)^2))

      # Calculate Manhattan Distance
      manhattan_dist[i, j] <- sum(abs(col1 - col2))

      # Calculate CCC
      ccc_dist[i, j] <- CCC(col1, col2, ci = "z-transform",conf.level = 0.95, na.rm = FALSE)[['rho.c']][['est']]
      # Calculate ICC
      dat <- data.frame(col1, col2)
      temp <- ICC(dat, missing=FALSE, alpha=.05, lmer=FALSE)
      res_ICC <- temp$results$ICC[3]
      icc_dist[i, j] <- res_ICC

      # calculate crossProd
      resCrossProd <- crossprod(col1, col1)
      cp_dist[i, j] <- resCrossProd

      alComparisionFold <- append(alComparisionFold, paste0(foldName[i], " vs ", foldName[j]))
    }
  }

  result <- list(
    Euclidean = euclidean_dist[upper.tri(euclidean_dist, diag=FALSE)],
    Manhattan = manhattan_dist[upper.tri(manhattan_dist, diag=FALSE)],
    ccc = ccc_dist[upper.tri(ccc_dist, diag=FALSE)],
    icc = icc_dist[upper.tri(icc_dist, diag=FALSE)],
    CrossProd = cp_dist[upper.tri(cp_dist, diag=FALSE)]
  )
  return(list(result, alComparisionFold))
}


#' Compute Performance Metrics for Comparing Two Data Matrices of abundances
#'
#' This function computes performance metrics to compare two data matrices.
#' The function calculates various performance metrics, including R-squared, adjusted R-squared,
#' Concordance Correlation Coefficient (CCC), Intra-class Correlation (ICC), Root Mean Square Error (RMSE),
#' Relative RMSE (RRMSE), and the name of the variable being compared.
#'
#' @param truth The true data matrix.
#' @param estimated The estimated data matrix to compare against the true data.
#' @param it_ An identifier or label for the comparison.
#'
#' @export
#'
#' @importFrom DescTools CCC
#' @importFrom psych ICC
#' @importFrom Metrics rmse
#'
#'
#' @return A data frame containing performance metrics for each variable in the comparison.
#' @details The function calculates and returns the following performance metrics for each variable:
#' \itemize{
#'   \item \code{R2}: R-squared.
#'   \item \code{R2_adj}: Adjusted R-squared.
#'   \item \code{ccc}: Concordance Correlation Coefficient (CCC).
#'   \item \code{icc}: Intra-class Correlation (ICC).
#'   \item \code{RRMSE}: Relative Root Mean Square Error (RRMSE).
#'   \item \code{RMSE}: Root Mean Square Error (RMSE).
#'   \item \code{Variable}: Name of the variable.
#'   \item \code{Comparison}: Identifier for the comparison.
#'}
#'
#' @examples
#'  truth_data <- matrix(rnorm(50), ncol = 5)
#'  estimated_data <- matrix(rnorm(50), ncol = 5)
#'  colnames(truth_data) <- colnames(estimated_data) <- paste0("C_", 1:5)
#'  result <- computPerf(truth_data, estimated_data, "Comparison_1")

computPerf <- function(truth, estimated, it_){
  colnameIntersect <- intersect(colnames(truth), colnames(estimated))
  num_cols <- length(colnameIntersect)
  truth <- sweep(truth[,colnameIntersect], 1, rowSums(truth[,colnameIntersect]), FUN = "/")
  estimated <- sweep(estimated[,colnameIntersect], 1, rowSums(estimated[,colnameIntersect]), FUN = "/")

  perfAll <- NULL
  for (v_ in colnameIntersect) {
    x <- truth[,v_]
    y <- as.vector(estimated[,v_])

    s2 <- (x - y)^2 / 2
    ms2 <- mean(s2)
    mx <- mean(x)
    my <- mean(y)
    ccc <- CCC(y, x,ci = "z-transform",conf.level = 0.95, na.rm = FALSE)[['rho.c']][['est']]
    dat <- data.frame(x, y)
    temp <- ICC(dat, missing=FALSE, alpha=.05, lmer=FALSE)
    res_ICC <- temp$results$ICC[3]
    RMSE <- rmse(actual = x, predicted = y)
    RRMSE <- RMSE/sd(x)
    model <- lm(y ~ x)
    R2 <- summary(model)$r.squared
    R2_adj <- summary(model)$adj.r.squared
    perfAll <- rbind(perfAll, c(R2, R2_adj, ccc, res_ICC, RRMSE, RMSE, v_, it_))
  }
  return(perfAll)
}



#' Compute Normalized Frobenius Norm of Absolute Error Matrix
#'
#' This function calculates the normalized Frobenius norm of the absolute error matrix
#' with respect to a reference bulk data matrix.
#'
#' @param absolute_error_matrix The absolute error matrix to be normalized.
#' @param bulkData The reference bulk data matrix.
#'
#' @export
#'
#' @return The normalized Frobenius norm of the absolute error matrix.
#'
#' @details The function calculates the Frobenius norm of the absolute error matrix and
#' normalizes it by dividing it by the Frobenius norm of the reference bulk data matrix.
#' The result is a measure of the error relative to the magnitude of the reference data.
#'
#' @examples
#' error_matrix <- matrix(rnorm(100), ncol = 5)
#' reference_data <- matrix(rnorm(100), ncol = 5)
#' norm_result <- normFrob(error_matrix, reference_data)

normFrob <- function(absolute_error_matrix, bulkData){
  frobenius_norm <- norm(as.matrix(absolute_error_matrix), "F")
  frobenius_norm_bulk <- norm(as.matrix(bulkData), "F")
  frobenius_norm_standard <- frobenius_norm / frobenius_norm_bulk
  return(frobenius_norm_standard)
}
