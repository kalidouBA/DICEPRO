#' Convert a Matrix to a Data Frame
#'
#' This function converts a given matrix into a data frame while preserving its row names and column names.
#' If the input is not a matrix, it is returned as is.
#'
#' @param mat A matrix to be converted into a data frame.
#'
#' @return A data frame with the same content, row names, and column names as the input matrix.
#'         If the input is not a matrix, it is returned unchanged.
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3, dimnames = list(c("A", "B", "C"), c("X", "Y", "Z")))
#' df <- convert_matrix_to_df(mat)
#' print(df)
#'

convert_matrix_to_df <- function(mat) {
  if (is.matrix(mat)) {
    df <- as.data.frame(mat)
    rownames(df) <- rownames(mat)
    colnames(df) <- colnames(mat)
    return(df)
  }
  return(mat)
}
