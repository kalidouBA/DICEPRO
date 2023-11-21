#' Plot Cell Abundance Heatmap and Error Plot
#'
#' This function combines the results of a cell abundance heatmap and an error plot for visualization.
# It takes a list `x` containing the results of both plots and generates the combined plot.
# You can pass additional arguments to customize the appearance of the plots.
#'
#' @param x A list containing the results of both the cell abundance heatmap and error plot between folds.
#' @param ... Additional arguments to customize the appearance of the plots.
#'
#' @return A combined plot showing the cell abundance heatmap and error plot.
#'
#'@method plot DICEPRO
#'@import patchwork
#'
#' @seealso Functions for generating individual plots: \code{\link{heatmap_abundances}}, \code{\link{metric_plot}}.
#'
#' @export

plot.DICEPRO <- function(x, ...){

  if(length(unique(x$Matrix_prediction$Iterate)) > 1){
    # Combine the heatmap and error between folds
    heatmap_abundances(x$Matrix_prediction) /
      metric_plot(x$performs2plot)
  } else {
    warning("Only one unique iteration found. No plot generated.")
  }

}

