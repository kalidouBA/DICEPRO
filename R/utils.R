#' Convert a Matrix to a Data Frame
#'
#' This function converts a given matrix into a data frame while preserving its row names and column names.
#' If the input is not a matrix, it is returned as is.
#'
#' @param mat A matrix to be converted into a data frame.
#'
#' @return A data frame with the same content, row names, and column names as the input matrix.
#'         If the input is not a matrix, it is returned unchanged.

convert_matrix_to_df <- function(mat) {
  if (is.matrix(mat)) {
    df <- as.data.frame(mat)
    rownames(df) <- rownames(mat)
    colnames(df) <- colnames(mat)
    return(df)
  }
  return(mat)
}


#' @importFrom stats cov var weighted.mean setNames
#' @importFrom utils read.table

utils::globalVariables(c(
  "lambda_", "gamma", "frobNorm", "constraint", "abs_constraint",
  "p_prime", "duration", "penalty", "objectiveValue"
))


#' Generate experiment paths
#' @param output_dir Output directory
#' @param bulkName Bulk name
#' @param refName Reference name
#' @param hspaceTechniqueChoose Hyperparameter space technique
#' @param algo_select Algorithm selection
#' @return List of paths
#' @export
generate_experiment_paths <- function(output_dir, bulkName, refName,
                                      hspaceTechniqueChoose, algo_select) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  data_dir <- file.path(output_dir, paste0(bulkName, "_", refName, "_",
                                           hspaceTechniqueChoose, "_", algo_select))
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  optim_dir <- file.path(data_dir, "optim")
  dir.create(optim_dir, showWarnings = FALSE, recursive = TRUE)

  report_dir <- file.path(optim_dir, "report")
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

  exp_name <- "optim"
  config_path <- file.path(optim_dir, paste0(exp_name, ".config.json"))
  trials_json_path <- file.path(optim_dir, paste0(exp_name, ".trials.json"))
  report_path <- file.path(report_dir, paste0(exp_name, ".report.png"))

  list(
    data_dir = data_dir,
    optim_dir = optim_dir,
    report_dir = report_dir,
    config_path = config_path,
    trials_json_path = trials_json_path,
    report_path = report_path
  )
}

#' Safe "or default" operator
#'
#' Returns `b` if `a` is NULL, otherwise returns `a`.
#'
#' @param a Object to check
#' @param b Default value if `a` is NULL
#' @return `a` if not NULL, otherwise `b`
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


#' Hyperparameter Optimization Tools
#'
#' A complete R implementation of hyperparameter optimization tools for
#' non-negative matrix factorization and other machine learning tasks.
#'
#' @details
#' This package provides functions for hyperparameter optimization including:
#' \itemize{
#'   \item Configuration parsing and validation
#'   \item Search space definition and sampling
#'   \item Optimization algorithms (random search, TPE, etc.)
#'   \item Result visualization and reporting
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/yourusername/hyperopt-r}
#'   \item Report bugs at \url{https://github.com/yourusername/hyperopt-r/issues}
#' }
#'
#' @name hyperopt-package
"_PACKAGE"
