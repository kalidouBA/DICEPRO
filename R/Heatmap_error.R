#' Create a Heatmap of Cell Abundances
#'
#' This function generates a heatmap to visualize cell abundances across different iterations and cell types.
#'
#' @param res2plot A data frame containing the results to be plotted, including Iteration, Cell_Type, and Abundances.
#'
#' @return A heatmap of cell abundances.
#'
#' @details This function uses the ggplot2 package to create a heatmap. It plots cell abundances across iterations and cell types, with colors representing the abundance levels.
#'
#' @seealso Other functions for data visualization: \code{\link{ggplot}}, \code{\link{geom_raster}}.
#'
#' @export
#'
#' @importFrom reshape2 melt
#' @import ggplot2

heatmap_abundances <- function(res2plot){
  value <- L1 <- Iterate <- sem <- CI_lower <- CI_upper <- Iteration <- Cell_Type <- Abundances <- NULL
  data2plot <- reshape2::melt(res2plot, id.vars = "Iterate")
  colnames(data2plot) <- c("Iteration",  "Cell_Type", "Abundances")

  # Create a heatmap using ggplot2
  p <- ggplot(data2plot, aes(Iteration,Cell_Type, fill=Abundances)) +
    geom_raster() +
    scale_fill_viridis_c(direction=-1) +
    theme_bw(base_size = 9, base_family = "Helvetica") %+replace%
    theme(
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9,angle = 90),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9,hjust=1),
      axis.ticks =  element_line(colour = "black"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9,angle=90),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      plot.margin = unit(c(0.5,  1, 1, 1), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size=9),
      legend.key.size = unit(0.5, 'cm')
    )
  return(p)
}

#' Create an Error/performance Plot
#'
#' This function generates an metric plot to visualize metric values between folds across different iterations.
#'
#' @param perf2plot A data frame containing the metric values to be plotted, including Iterate, and metric
#'
#' @return An metric plot showing the mean metric values with confidence intervals.
#'
#' @details This function calculates the mean and confidence intervals of metric values between folds and creates a line plot with shaded confidence intervals.
#' It is useful for visualizing the variation in metric across iterations.
#'
#' @seealso Other functions for data visualization: \code{\link{ggplot}}, \code{\link{geom_line}}, \code{\link{geom_ribbon}}, \code{\link{facet_wrap}}.
#'
#' @export
#'
#' @import ggplot2

metric_plot <- function(perf2plot){
  metric <- It <- CI_lower <- CI_upper <- sem <- NULL
  # Calculate mean and confidence intervals of metric values
  perf2plot_mean <- perf2plot %>%
    filter(!is.na(metric)) %>%
    group_by(It) %>%
    summarise(n = n(),
              mean = mean(metric),
              median = median(metric),
              sd = sd(metric)) %>%
    mutate(sem = sd / sqrt(n - 1),
           CI_lower = mean + qt((1 - 0.95) / 2, n - 1) * sem,
           CI_upper = mean - qt((1 - 0.95) / 2, n - 1) * sem)

  # Create the metric plot using ggplot2
  p <- ggplot(perf2plot_mean, aes(x=It, y=mean)) +
    geom_line() +
    geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper),color="grey70",alpha=0.4)+
    labs(x = "Iteration", y = "Error between folds") +
    theme_bw(base_size = 9, base_family = "Helvetica") %+replace%
    theme(
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9,angle = 90),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9,hjust=1),
      axis.ticks =  element_line(colour = "black"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9,angle=90),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      plot.margin = unit(c(0.5,  1, 1, 1), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size=9),
      legend.key.size = unit(0.5, 'cm')
    )
  return(p)
}


