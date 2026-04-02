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
#' @return List of paths
#' @export
.generate_experiment_paths <- function(output_dir, bulkName, refName) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  data_dir <- file.path(output_dir, paste0("dicepro_",bulkName, "_", refName))
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)


  config_path <- file.path(data_dir, "optim.config.json")

  list(
    data_dir = data_dir,
    config_path = config_path
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
#' @name or-default-operator
#' @rdname or-default-operator
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

