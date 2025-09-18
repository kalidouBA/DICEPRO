#' Plot Cell Abundance Heatmap and Error Plot
#'
#' This function combines the results of a cell abundance heatmap and an error plot for visualization.
#' It takes a list `x` containing the results of both plots and generates the combined plot.
#' You can pass additional arguments to customize the appearance of the plots.
#'
#' @param x A list containing the results of both the cell abundance heatmap and error plot between folds.
#' @param ... Additional arguments to customize the appearance of the plots.
#'
#' @return A combined plot showing the cell abundance heatmap and error plot.
#'
#' @method plot DICEPRO
#' @import patchwork
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

#' Plot Hyperparameter Optimization Report
#'
#' This function visualizes hyperparameter optimization results as a scatter matrix
#' and violin plots for the best parameters.
#'
#' @param exp Path to the experiment results folder
#' @param params Vector of parameters to visualize
#' @param metric Metric for point size (default: "loss")
#' @param loss_metric Loss metric (default: "loss")
#' @param loss_behaviour Loss behaviour ("min" or "max")
#' @param not_log Parameters not to log-scale
#' @param categorical Categorical parameters
#' @param max_deviation Maximum deviation for outlier detection
#' @param title Plot title
#' @return Combined ggplot object with scatter matrix and violin plots
#' @export
#' @import ggplot2
#' @importFrom purrr map map_dbl compact
#' @importFrom jsonlite read_json
#' @importFrom gridExtra grid.arrange
plot_hyperopt_report_v2 <- function(exp, params, metric = "loss", loss_metric = "loss",
                                    loss_behaviour = "min", not_log = NULL,
                                    categorical = NULL, max_deviation = NULL,
                                    title = NULL) {

  # ---- Internal functions ----
  get_results <- function(exp) {
    report_path <- file.path(exp, "results")
    if (!dir.exists(report_path)) stop(paste("The folder", report_path, "does not exist"))

    files <- list.files(report_path, full.names = TRUE)
    results <- list()

    for (file in files) {
      if (file.info(file)$isdir) next
      res <- tryCatch(jsonlite::read_json(file), error = function(e) {
        warning(paste("Error reading file", file, ":", e$message))
        NULL
      })
      if (!is.null(res)) results[[length(results) + 1]] <- res
    }

    if (length(results) == 0) stop("No results found")
    return(results)
  }

  scale_01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }

  filter_outliers <- function(values, max_dev) {
    mean_val <- mean(values, na.rm = TRUE)
    abs(values - mean_val) < (mean_val + max_dev)
  }

  # ---- Load results ----
  results <- get_results(exp)

  loss <- sapply(results, function(x) x$returned_dict[[loss_metric]] %||% NA_real_)
  scores <- sapply(results, function(x) x$returned_dict[[metric]] %||% NA_real_)
  keep <- !is.na(scores) & !is.na(loss)
  scores <- scores[keep]
  loss <- loss[keep]

  values <- list()
  for (p in params) {
    values[[p]] <- sapply(results, function(x) {
      if (!is.null(x$current_params[[p]])) x$current_params[[p]] else NA_real_
    })
  }

  # ---- Outlier filtering ----
  if (!is.null(max_deviation)) {
    keep <- filter_outliers(loss, max_deviation)
    loss <- loss[keep]
    scores <- scores[keep]
    for (p in params) values[[p]] <- values[[p]][keep]
  }

  categorical <- categorical %||% character()
  if (length(categorical) > 0) {
    for (p in categorical) values[[p]] <- as.character(values[[p]])
  }

  # ---- Sorting ----
  all_numerical <- lapply(values[setdiff(params, categorical)], unlist)
  all_categorical <- lapply(values[categorical], unlist)
  sorted_idx <- do.call(order, c(all_numerical, all_categorical, list(loss, scores)))

  loss <- loss[sorted_idx]
  scores <- scores[sorted_idx]
  for (p in params) values[[p]] <- values[[p]][sorted_idx]

  scores_scaled <- scale_01(scores)

  # ---- Best points ----
  if (loss_behaviour == "min") {
    lmaxs <- loss > min(loss, na.rm = TRUE)
  } else {
    lmaxs <- loss < max(loss, na.rm = TRUE)
  }

  top_pct <- ceiling(length(scores) * 0.05)
  smaxs <- order(scores, decreasing = TRUE)[1:top_pct]
  cmaxs <- rep(NA, length(scores))
  cmaxs[smaxs] <- scale_01(scores[smaxs])

  # ---- Prepare data frame ----
  df <- as.data.frame(values)
  df$loss <- loss
  df$scores <- scores_scaled
  df$smaxs <- seq_along(scores) %in% smaxs
  df$lmaxs <- lmaxs
  df$cmaxs <- cmaxs

  # ---- Plot function ----
  make_plot <- function(x, y, diag = FALSE) {
    if (diag) {
      p <- ggplot(df, aes(x = .data[[x]], y = loss)) +
        geom_point(aes(size = scores, color = lmaxs), alpha = 0.7) +
        geom_point(data = df[df$smaxs, ], aes(size = scores, fill = cmaxs),
                   shape = 21, color = "black") +
        scale_size_continuous(range = c(1, 10)) +
        scale_color_manual(values = c("TRUE" = "orange", "FALSE" = "red")) +
        scale_fill_viridis_c() +
        theme_minimal() + theme(legend.position = "none")
      if (!(x %in% not_log)) p <- p + scale_x_log10()
    } else {
      p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
        geom_point(aes(size = scores, color = loss), alpha = 0.7) +
        geom_point(data = df[df$smaxs, ], aes(size = scores, fill = cmaxs),
                   shape = 21, color = "black") +
        scale_size_continuous(range = c(1, 10)) +
        scale_color_viridis_c() +
        scale_fill_viridis_c() +
        theme_minimal() + theme(legend.position = "none")
      if (!(x %in% not_log)) p <- p + scale_x_log10()
      if (!(y %in% not_log)) p <- p + scale_y_log10()
    }
    return(p)
  }

  # ---- Scatter matrix ----
  plots <- list()
  for (i in seq_along(params)) {
    row_plots <- list()
    for (j in seq_along(params)) {
      p1 <- params[i]
      p2 <- params[j]
      row_plots[[j]] <- make_plot(p2, p1, diag = p1 == p2)
    }
    plots[[i]] <- row_plots
  }

  plots_flat <- unlist(plots, recursive = FALSE)
  plot_matrix <- gridExtra::grid.arrange(
    grobs = plots_flat,
    nrow = length(params),
    ncol = length(params)
  )

  # ---- Violin / bar plots for top points ----
  violin_plots <- list()
  for (i in seq_along(params)) {
    p <- params[i]
    df_best <- df[df$smaxs, ]
    if (p %in% categorical) {
      violin_plots[[i]] <- ggplot(df_best, aes(x = .data[[p]])) +
        geom_bar(fill = "forestgreen", alpha = 0.3) +
        theme_minimal() + labs(x = p, y = "Count")
    } else {
      p_plot <- ggplot(df_best, aes(x = 1, y = .data[[p]])) +
        geom_violin(fill = "forestgreen", alpha = 0.3) +
        geom_boxplot(width = 0.1, fill = "orange", alpha = 0.7) +
        geom_point(aes(color = cmaxs), position = position_jitter(width = 0.1)) +
        scale_color_viridis_c() +
        theme_minimal() + labs(x = p, y = "") +
        theme(axis.text.x = element_blank())
      if (!(p %in% not_log)) p_plot <- p_plot + scale_y_log10()
      violin_plots[[i]] <- p_plot
    }
  }

  violin_grid <- gridExtra::grid.arrange(
    grobs = violin_plots,
    nrow = 1,
    ncol = length(params)
  )

  # ---- Combine plots ----
  final_plot <- gridExtra::grid.arrange(
    plot_matrix, violin_grid,
    nrow = 2, heights = c(3, 1),
    top = title
  )

  return(final_plot)
}

