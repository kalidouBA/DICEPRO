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
#' @importFrom stats lm
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


#' Compare Two Data Frames Column by Column
#'
#' This function compares two numeric data frames of equal dimensions by computing
#' statistical metrics (NRMSE, adjusted R^2, ICC3, and CCC) for each column,
#' after Z-score normalization. It returns per-column results and a combined summary.
#'
#' @param x A numeric data frame.
#' @param y A numeric data frame with the same dimensions as `x`.
#' @param method Aggregation method for combined metrics: "weighted" (default, weights by SD) or "unweighted".
#'
#' @return A list with:
#' \describe{
#'   \item{by_column}{Named list of per-column metrics: NRMSE, R^2 adjusted, ICC3, and CCC.}
#'   \item{combined}{Aggregated metrics across all columns (weighted or unweighted).}
#'   \item{method}{Description of the aggregation method used.}
#' }
#'
#' @details
#' - Z-score normalization is applied to each column.
#' - If standard deviation is zero, the column is set to zero.
#' - ICC3 is computed if the \pkg{irr} package is available.
#' - CCC (Concordance Correlation Coefficient) is computed using the formula from Lin (1989).
#'
#' @references
#' Lin, L. I. (1989). A concordance correlation coefficient to evaluate reproducibility. *Biometrics*, 45(1), 255-268.
#'
#' @examples
#' \dontrun{
#' df1 <- data.frame(a = rnorm(100), b = runif(100))
#' df2 <- data.frame(a = df1$a + rnorm(100, sd = 0.1), b = df1$b + rnorm(100, sd = 0.1))
#' compare_abundances(df1, df2)
#' }
#'
#' @importFrom stats cov var weighted.mean setNames
#' @importFrom utils read.table
#' @export
compare_abundances <- function(x, y, method = "weighted") {
  tryCatch({
    # 1. Verifications initiales
    if (!identical(dim(x), dim(y))) stop("Dimensions incoherentes")
    if (any(is.na(x)) || any(is.na(y))) stop("NA detectes")
    if (any(is.infinite(as.matrix(x)))) stop("Inf dans x")
    if (any(is.infinite(as.matrix(y)))) stop("Inf dans y")

    # 2. Normalisation Z-score
    normalize <- function(v) {
      if (sd(v) == 0) return(v * 0)
      (v - mean(v)) / sd(v)
    }

    x_norm <- as.data.frame(lapply(x, normalize))
    y_norm <- as.data.frame(lapply(y, normalize))

    # 3. Calcul des metriques par colonne
    metrics <- lapply(colnames(x), function(col) {
      vx <- x_norm[[col]]
      vy <- y_norm[[col]]

      # NRMSE
      nrmse <- sqrt(mean((vx - vy)^2)) / sd(vy)

      # R^2 ajuste
      r2_adj <- summary(lm(vy ~ vx))$adj.r.squared

      # ICC3
      icc_val <- if (requireNamespace("irr", quietly = TRUE)) {
        irr::icc(data.frame(vx, vy), model = "twoway", type = "agreement")$value
      } else NA

      # Concordance Correlation Coefficient (CCC)
      mean_x <- mean(vx)
      mean_y <- mean(vy)
      var_x <- var(vx)
      var_y <- var(vy)
      cov_xy <- cov(vx, vy)
      ccc_val <- (2 * cov_xy) / (var_x + var_y + (mean_x - mean_y)^2)

      c(NRMSE = nrmse, R2_adj = r2_adj, ICC3 = icc_val, CCC = ccc_val)
    })

    # 4. Ponderation
    weights <- if (method == "weighted") sapply(x, sd) else rep(1, ncol(x))
    combined <- list(
      NRMSE = weighted.mean(sapply(metrics, `[`, "NRMSE"), weights),
      R2_adj = weighted.mean(sapply(metrics, `[`, "R2_adj"), weights),
      ICC3 = weighted.mean(sapply(metrics, `[`, "ICC3"), weights, na.rm = TRUE),
      CCC = weighted.mean(sapply(metrics, `[`, "CCC"), weights)
    )

    # 5. Resultat final
    list(
      by_column = setNames(metrics, colnames(x)),
      combined = combined,
      method = paste("Methode:", method, if (method == "weighted") "(ponderee par SD)" else "")
    )

  }, error = function(e) {
    warning("Erreur: ", e$message)
    return(NULL)
  })
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
