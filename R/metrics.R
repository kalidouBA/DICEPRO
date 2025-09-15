#' Compute Precision of Estimation
#'
#' This function computes the precision of estimation by evaluating
#' the Frobenius norm of the difference between the gene expression matrix (Y)
#' and the product of the reference matrix and the estimated abundance matrix (Y_tilde),
#' normalized by the Frobenius norm of the gene expression matrix (Y).
#'
#' @param Y The gene expression matrix
#' @param Y_tilde The product of the reference matrix and the estimated abundance matrix
#'
#' @return The precision value as a numeric scalar
#'
#' @examples
#' Y <- matrix(c(3, 5, 7, 2, 6, 4, 0, 2, 8), nrow = 3, ncol = 3, byrow = TRUE)
#' Y_tilde <- Y
#' precision <- compute_precision(Y, Y_tilde)
#' cat("Precision of estimation:", precision, "\n")
#'
#' @export
compute_precision <- function(Y, Y_tilde) {
  frobenius_norm <- function(matrix) {
    sqrt(sum(matrix^2, na.rm = TRUE))
  }

  diff_matrix <- Y - Y_tilde
  norm_diff <- frobenius_norm(diff_matrix)
  norm_Y <- frobenius_norm(Y)

  precision <- norm_diff / norm_Y

  return(precision)
}


#' Compute Performance Metrics for Comparing Two Data Matrices
#'
#' This function computes performance metrics to compare two data matrices.
#' It calculates either Normalized Root Mean Square Error (NRMSE) or
#' R-squared (R2) based on the specified metric.
#'
#' @param x A matrix of values (e.g., cell type proportions from iteration i)
#' @param y A matrix of values (e.g., cell type proportions from iteration i+1)
#' @param metric A character string specifying the metric to compute:
#'               "NRMSE" for Normalized Root Mean Square Error or
#'               "R2" for R-squared
#' @param method Character string specifying the normalization method for NRMSE:
#'               "sd" (standard deviation), "mean", "maxmin" (range), or "iq" (interquartile range)
#' @param transformation Character string specifying data transformation for NRMSE:
#'               "none", "sqrt", "4thrt", "log", "log10", "log2", "log1p", "arcsine", or "other"
#' @param trans_function Character string specifying custom transformation function
#'               (required when transformation = "other")
#'
#' @return A numeric value representing the computed performance metric
#'
#' @details
#' For NRMSE: Uses the nrmse function with specified parameters
#' For R2: Uses the r2_1dim function
#'
#' @examples
#' if(interactive()){
#'   x <- matrix(rnorm(50), ncol = 5)
#'   y <- matrix(rnorm(50), ncol = 5)
#'   colnames(x) <- colnames(y) <- paste0("C_", 1:5)
#'   result <- computPerf(x, y, "NRMSE")
#' }
#'
#' @export
computPerf <- function(x, y, metric, method = "sd", transformation = "none", trans_function = "none") {

  # Calculate Normalized Root Mean Squared Error (NRMSE)
  if(metric == "NRMSE"){
    # Convert matrices to vectors for calculation
    x_vec <- as.vector(x)
    y_vec <- as.vector(y)

    # Remove NA values
    excl <- is.na(x_vec) | is.na(y_vec)
    x_vec <- x_vec[!excl]
    y_vec <- y_vec[!excl]

    if (length(x_vec) < 2) {
      stop("Not enough observations to calculate NRMSE")
    }

    # Use the nrmse function with specified parameters
    perf <- nrmse(x_vec, y_vec, method = method, transformation = transformation, trans_function = trans_function)
    perf <- round(perf, 2)
  }

  # Calculate R-squared
  else if(metric == "R2"){
    # Convert matrices to vectors for calculation
    x_vec <- as.vector(x)
    y_vec <- as.vector(y)

    # Remove NA values
    excl <- is.na(x_vec) | is.na(y_vec)
    x_vec <- x_vec[!excl]
    y_vec <- y_vec[!excl]

    if (length(x_vec) < 2) {
      stop("Not enough observations to calculate R2")
    }

    # Use the r2_1dim function
    perf <- r2_1dim(x_vec, y_vec)
    perf <- round(perf, 2)
  }
  else {
    stop("Invalid metric specified. Choose 'NRMSE' or 'R2'")
  }

  return(perf)
}

#' Compare Two Data Frames Column by Column
#'
#' This function compares two numeric data frames of equal dimensions by computing
#' statistical metrics (NRMSE, R^2, ICC, and CCC) for each column,
#' after Z-score normalization. It returns per-column results and a combined summary.
#'
#' @param x A numeric data frame
#' @param y A numeric data frame with the same dimensions as `x`
#' @param method Aggregation method for combined metrics: "weighted" (default, weights by SD) or "unweighted"
#' @param nrmse_method Character string specifying the normalization method for NRMSE:
#'               "sd" (standard deviation), "mean", "maxmin" (range), or "iq" (interquartile range)
#' @param transformation Character string specifying data transformation for NRMSE:
#'               "none", "sqrt", "4thrt", "log", "log10", "log2", "log1p", "arcsine", or "other"
#' @param trans_function Character string specifying custom transformation function
#'               (required when transformation = "other")
#'
#' @return A list with:
#' \describe{
#'   \item{by_column}{Named list of per-column metrics: NRMSE, R^2, ICC, and CCC}
#'   \item{combined}{Aggregated metrics across all columns (weighted or unweighted)}
#'   \item{method}{Description of the aggregation method used}
#' }
#'
#' @details
#' - Z-score normalization is applied to each column
#' - If standard deviation is zero, the column is set to zero
#' - ICC is computed using full_icc_gaussian function
#' - CCC (Concordance Correlation Coefficient) is computed using the formula from Lin (1989)
#'
#' @references
#' Lin, L. I. (1989). A concordance correlation coefficient to evaluate reproducibility.
#' *Biometrics*, 45(1), 255-268.
#'
#' @examples
#' \dontrun{
#' df1 <- data.frame(a = rnorm(100), b = runif(100))
#' df2 <- data.frame(a = df1$a + rnorm(100, sd = 0.1), b = df1$b + rnorm(100, sd = 0.1))
#' compare_abundances(df1, df2)
#' }
#'
#' @importFrom stats cov var weighted.mean setNames sd
#' @importFrom utils read.table
#' @export
compare_abundances <- function(x, y, method = "weighted",
                               nrmse_method = "sd", transformation = "none",
                               trans_function = "none") {
  tryCatch({
    # 1. Initial verifications
    if (!identical(dim(x), dim(y))) stop("Inconsistent dimensions")
    if (any(is.na(x)) || any(is.na(y))) stop("NA values detected")
    if (any(is.infinite(as.matrix(x)))) stop("Infinite values in x")
    if (any(is.infinite(as.matrix(y)))) stop("Infinite values in y")

    # 2. Z-score normalization
    normalize <- function(v) {
      if (sd(v, na.rm = TRUE) == 0) return(v * 0)
      (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    }

    x_norm <- as.data.frame(lapply(x, normalize))
    y_norm <- as.data.frame(lapply(y, normalize))

    # 3. Calculate metrics by column
    metrics <- lapply(colnames(x), function(col) {
      vx <- x_norm[[col]]
      vy <- y_norm[[col]]

      # NRMSE with specified parameters
      nrmse_val <- nrmse(vx, vy, method = nrmse_method,
                         transformation = transformation,
                         trans_function = trans_function)

      # R^2
      r2_val <- r2_1dim(vx, vy)

      # ICC using full_icc_gaussian
      icc_val <- tryCatch({
        full_icc_gaussian(matrix(vx, ncol = 1), matrix(vy, ncol = 1))
      }, error = function(e) NA)

      # Concordance Correlation Coefficient (CCC)
      mean_x <- mean(vx, na.rm = TRUE)
      mean_y <- mean(vy, na.rm = TRUE)
      var_x <- var(vx, na.rm = TRUE)
      var_y <- var(vy, na.rm = TRUE)
      cov_xy <- cov(vx, vy, use = "complete.obs")
      ccc_val <- (2 * cov_xy) / (var_x + var_y + (mean_x - mean_y)^2)

      c(NRMSE = nrmse_val, R2 = r2_val, ICC = icc_val, CCC = ccc_val)
    })

    # 4. Weighting
    weights <- if (method == "weighted") sapply(x, sd, na.rm = TRUE) else rep(1, ncol(x))
    combined <- list(
      NRMSE = weighted.mean(sapply(metrics, `[`, "NRMSE"), weights, na.rm = TRUE),
      R2 = weighted.mean(sapply(metrics, `[`, "R2"), weights, na.rm = TRUE),
      ICC = weighted.mean(sapply(metrics, `[`, "ICC"), weights, na.rm = TRUE),
      CCC = weighted.mean(sapply(metrics, `[`, "CCC"), weights, na.rm = TRUE)
    )

    # 5. Final result
    list(
      by_column = setNames(metrics, colnames(x)),
      combined = combined,
      method = paste("Method:", method, if (method == "weighted") "(weighted by SD)" else ""),
      nrmse_params = paste("NRMSE method:", nrmse_method, "| Transformation:", transformation)
    )

  }, error = function(e) {
    warning("Error: ", e$message)
    return(NULL)
  })
}

#' Calculate Normalized Root Mean Square Error (NRMSE)
#'
#' This function calculates the Normalized Root Mean Square Error (NRMSE)
#' between predicted and observed values with various normalization options.
#'
#' @param pred Numeric vector of predicted values
#' @param obs Numeric vector of observed values
#' @param method Character string specifying the normalization method:
#'               "sd" (standard deviation), "mean", "maxmin" (range), or "iq" (interquartile range)
#' @param transformation Character string specifying data transformation:
#'               "none", "sqrt", "4thrt", "log", "log10", "log2", "log1p", "arcsine", or "other"
#' @param trans_function Character string specifying custom transformation function
#'               (required when transformation = "other")
#'
#' @return The NRMSE value as a numeric scalar
#'
#' @details
#' The function supports various data transformations and normalization methods.
#' Back-transformation is applied before calculating NRMSE when using transformations.
#'
#' @examples
#' pred <- rnorm(100)
#' obs <- pred + rnorm(100, sd = 0.1)
#' nrmse_value <- nrmse(pred, obs, method = "sd", transformation = "none")
#'
#' @importFrom stats sd quantile
#' @export
nrmse <- function(pred, obs, method = "sd", transformation = "none", trans_function = "none") {

  # Data input validation
  if (length(pred) != length(obs))
    stop(paste("The observation and prediction vectors do not have the same length.",
               "Obs:", length(obs), "Pred:", length(pred)))
  if (!method %in% c("mean", "sd", "maxmin", "iq"))
    stop("Invalid method specified. Choose from: 'mean', 'sd', 'maxmin', 'iq'")
  if (!transformation %in% c("none", "sqrt", "4thrt", "log", "log10",
                             "log2", "log1p", "arcsine", "other"))
    stop("Invalid transformation specified")
  if (transformation == "other" & trans_function == "none")
    stop("trans_function must be specified when transformation = 'other'")

  # Backtransform if needed
  if (transformation == "sqrt") {
    obs <- obs^2
    pred <- pred^2
  }
  if (transformation == "4thrt") {
    obs <- obs^4
    pred <- pred^4
  }
  if (transformation == "log") {
    obs <- exp(obs)
    pred <- exp(pred)
  }
  if (transformation == "log10") {
    obs <- 10^obs
    pred <- 10^pred
  }
  if (transformation == "log2") {
    obs <- 2^obs
    pred <- 2^pred
  }
  if (transformation == "log1p") {
    obs <- expm1(obs)
    pred <- expm1(pred)
  }
  if (transformation == "arcsine") {
    obs <- sin(obs)^2
    pred <- sin(pred)^2
  }
  if (transformation == "other") {
    x <- obs
    obs <- eval(parse(text = trans_function))
    x <- pred
    pred <- eval(parse(text = trans_function))
  }

  excl <- is.na(obs) | is.na(pred)
  obs_in <- obs[!excl]
  pred_in <- pred[!excl]

  if (length(obs_in) < 2) {
    out <- NA
    message("Not enough observations/predictions to calculate NRMSE, returning NA")
  } else {
    # Calculation of NRMSE
    sq_sums <- sum((obs_in - pred_in)^2, na.rm = TRUE)
    mse <- sq_sums / length(obs_in)
    rmse <- sqrt(mse)

    # Normalization
    if (method == "sd") {
      out <- rmse / sd(obs_in, na.rm = TRUE)
    } else if (method == "mean") {
      out <- rmse / mean(obs_in, na.rm = TRUE)
    } else if (method == "maxmin") {
      out <- rmse / (max(obs_in, na.rm = TRUE) - min(obs_in, na.rm = TRUE))
    } else if (method == "iq") {
      out <- rmse / (quantile(obs_in, 0.75, na.rm = TRUE) - quantile(obs_in, 0.25, na.rm = TRUE))
      names(out) <- NULL
    }

    out <- abs(out)
  }

  return(out)
}

#' Calculate R-squared for One-Dimensional Data
#'
#' This function calculates the R-squared value between two numeric vectors,
#' with an option for adjusted R-squared.
#'
#' @param x Numeric vector of predictor values
#' @param y Numeric vector of response values
#' @param adjusted Logical indicating whether to calculate adjusted R-squared (default: FALSE)
#'
#' @return The R-squared value as a numeric scalar
#'
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100, sd = 0.1)
#' r2 <- r2_1dim(x, y)
#' r2_adj <- r2_1dim(x, y, adjusted = TRUE)
#'
#' @export
r2_1dim <- function(x, y, adjusted = FALSE) {

  res <- cor(x, y, use = "complete.obs")^2

  if (adjusted) {
    res <- 1 - (1 - res) * (length(x) - 1) / (length(x) - 2)
  }

  return(res)
}

#' Calculate Intraclass Correlation Coefficient (ICC) for Gaussian Data
#'
#' This function calculates the ICC(3,1) (two-way mixed effects, absolute agreement)
#' for multivariate Gaussian data using linear mixed effects models.
#'
#' @param x Matrix of observed values (subjects x variables)
#' @param y Matrix of predicted values (same dimensions as x)
#'
#' @return The ICC value as a numeric scalar
#'
#' @details
#' Requires the 'lme4' and 'insight' packages. The function uses a linear mixed
#' effects model with random effects for population and subject within population.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100), ncol = 5)
#' y <- x + matrix(rnorm(100, sd = 0.1), ncol = 5)
#' icc_val <- full_icc_gaussian(x, y)
#' }
#'
#' @importFrom lme4 lmer VarCorr
#' @export
full_icc_gaussian <- function(x, y) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Please install 'lme4'")
  if (!requireNamespace("insight", quietly = TRUE)) stop("Please install 'insight'")

  n <- nrow(x)
  p <- ncol(x)

  df <- data.frame(
    value = c(as.vector(x), as.vector(y)),
    subject = factor(rep(1:n, times = 2 * p)),
    population = factor(rep(rep(colnames(x), each = n), times = 2)),
    method = factor(rep(c("obs", "pred"), each = n * p))
  )

  model <- lme4::lmer(
    value ~ 1 + (1 | population/subject),
    data = df
  )

  var <- lme4::VarCorr(model)

  var_r1 <- var$population[1, 1]
  var_r2 <- var$`subject:population`[1, 1]
  var_e <- attributes(var)$sc^2

  icc3 <- (var_r1 + var_r2) / ((var_r1 + var_r2) + var_e)
  return(as.numeric(icc3))
}
