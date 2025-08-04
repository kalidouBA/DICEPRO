#' Generate a Kraljic Pareto plot without ground-truth-based metrics
#'
#' This function generates an interactive Plotly figure showing the Pareto frontier
#' based on Frobenius norm and constraint deviation, excluding performance indicators
#' calculated against ground truth (e.g., R2, NRMSE).
#'
#' The function reads input data from a folder derived from dataset and algorithm names,
#' processes the data, computes the Pareto frontier using the \code{KraljicMatrix} package,
#' and saves an interactive HTML file. It also returns the frontier solution with the lowest
#' Frobenius norm and constraint deviation closest to zero (i.e., constraint closest to 1).
#'
#' @param outputDir Name of the directory saving output.
#'
#' @return A named list corresponding to the best solution on the Pareto frontier.
#' @export
#'
#' @importFrom dplyr mutate row_number
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom htmlwidgets saveWidget
#' @importFrom KraljicMatrix get_frontier
#' @importFrom purrr map compact
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' plot_kraljic(outputDir)
#' }
#'
#'
plot_kraljic <- function(outputDir) {
  json_dir <- paste0(outputDir, "/optim/results")

  data2plot <- load_all_json(json_dir)
  rownames(data2plot) <- gsub("^.*\\.(\\d+)$", "\\1", rownames(data2plot))

  data2plot <- data2plot %>%
    dplyr::mutate(
      log_lambda = log10(lambda_ + 1),
      log_gamma = log10(gamma + 1),
      log_frob = log10(frobNorm),
      scaled_constraint = ifelse(abs(constraint - 1) <= 0.05, 1, 0)
    )
  data2plot$abs_constraint <- abs(1 - data2plot$constraint)

  # Compute Pareto frontier
  frontier_result <- KraljicMatrix::get_frontier(
    data = data2plot,
    x = frobNorm,
    y = abs_constraint,
    quadrant = "bottom.left",
    decreasing = FALSE
  )
  data2plotfrontier_result <- data2plot[rownames(frontier_result), ]

  # Plot
  fig_kraljic <- plotly::plot_ly() %>%
    plotly::add_trace(
      data = data2plot,
      x = ~frobNorm,
      y = ~abs_constraint,
      type = 'scatter',
      mode = 'markers',
      marker = list(color = 'gray', size = 6),
      name = "Valid Solutions",
      text = ~paste(
        "<b>\u03BB:</b>", formatC(lambda_, format = "e", digits = 2), "<br>",
        "<b>\u03B3:</b>", formatC(gamma, format = "e", digits = 2), "<br>",
        "<b>Frobenius Norm:</b>", formatC(frobNorm, format = "e", digits = 2), "<br>",
        "<b>Constraint:</b>", round(constraint, 4), "<br>",
        "<b>p_prime:</b>", signif(p_prime, 3), "<br>",
        "<b>Duration:</b>", signif(duration, 3), "<br>",
        "<b>Penalty:</b>", signif(penalty, 3), "<br>",
        "<b>Objective:</b>", signif(objectiveValue, 3)
      ),
      hoverinfo = "text"
    ) %>%
    plotly::add_trace(
      data = data2plotfrontier_result[order(data2plotfrontier_result$frobNorm), ],
      x = ~frobNorm,
      y = ~abs_constraint,
      type = 'scatter',
      mode = 'markers+lines',
      marker = list(color = 'red', size = 8, symbol = 'diamond'),
      line = list(color = 'red', width = 2),
      name = "Pareto Frontier (Kraljic)",
      text = ~paste(
        "<b>\u03BB:</b>", formatC(lambda_, format = "e", digits = 2), "<br>",
        "<b>\u03B3:</b>", formatC(gamma, format = "e", digits = 2), "<br>",
        "<b>Frobenius Norm:</b>", formatC(frobNorm, format = "e", digits = 2), "<br>",
        "<b>Constraint:</b>", round(constraint, 4), "<br>",
        "<b>p_prime:</b>", signif(p_prime, 3), "<br>",
        "<b>Duration:</b>", signif(duration, 3), "<br>",
        "<b>Penalty:</b>", signif(penalty, 3), "<br>",
        "<b>Objective:</b>", signif(objectiveValue, 3)
      ),
      hoverinfo = "text"
    ) %>%
    plotly::layout(
      title = "<b>Frobenius Norm vs Constraint (KraljicMatrix)</b>",
      xaxis = list(title = "<b>Frobenius Norm</b> <i>(-log scale)</i>", type = "log", autorange = "reversed"),
      yaxis = list(title = "<b>abs(1-Constraint)</b> <i>(-log scale)</i>", type = "log", autorange = "reversed"),
      legend = list(x = 0.02, y = 0.98)
    )

  configured_fig <- fig_kraljic %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "pareto_front_log_scale",
        width = 500,
        height = 300,
        scale = 1.5
      ),
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("lasso2d", "select2d")
    )

  htmlwidgets::saveWidget(
    widget = configured_fig,
    file = file.path(sprintf("%s/optim/report/front_pareto_kraljic.html", outputDir)),
    libdir = "lib",
    selfcontained = FALSE
  )

  # ---- Return best solution as named list ----
  best_idx <- which.min(
    data2plotfrontier_result$frobNorm + 10 * abs(1 - data2plotfrontier_result$constraint)
  )
  bestHP <- list('frontPreto' = data2plotfrontier_result)
  bestHP$best_idx <- best_idx
  return(bestHP)
}




#' Process a single JSON file into a data frame
#'
#' This function reads a JSON file containing optimization results,
#' extracts key fields such as loss, constraint, status, duration,
#' and flattens the current parameters into separate columns.
#' It returns a data frame of these values.
#'
#' If the file cannot be parsed, it returns NULL with a warning.
#'
#' @param file_path Character string specifying the path to the JSON file.
#'
#' @return A data frame with columns for loss, constraint, status, duration,
#'   and one column for each parameter in \code{current_params},
#'   or \code{NULL} if an error occurs.
#'
#' @importFrom jsonlite fromJSON
#'
#' @examples
#' \dontrun{
#' df <- process_single_json("path/to/file.json")
#' }
#'
process_single_json <- function(file_path) {
  tryCatch({
    data <- jsonlite::fromJSON(file_path)

    # Extract main elements
    main_data <- data.frame(
      loss = data$returned_dict$loss,
      constraint = data$returned_dict$constraint,
      status = data$returned_dict$status,
      duration = data$returned_dict$duration,
      stringsAsFactors = FALSE
    )

    # Extract current_params and add as columns
    params <- t(data.frame(unlist(data$returned_dict$current_params)))

    # Combine main data and parameters
    cbind(main_data, params)
  }, error = function(e) {
    warning(paste("Error reading file", file_path, ":", e$message))
    NULL
  })
}

#' Load and combine all JSON files from a directory into a data frame
#'
#' This function searches the specified directory for all JSON files,
#' applies \code{process_single_json} to each,
#' combines the results into a single data frame,
#' and adds a unique integer \code{id} column.
#'
#' Files that cause errors during parsing are silently skipped with a warning.
#'
#' @param directory Character string specifying the folder containing JSON files.
#'
#' @return A data frame combining all processed JSON files with an added \code{id} column.
#'
#' @importFrom dplyr bind_rows mutate row_number
#' @importFrom purrr map compact
#'
#' @examples
#' \dontrun{
#' combined_df <- load_all_json("path/to/json_folder")
#' }
#'
load_all_json <- function(directory) {
  # List all JSON files in the directory (case-insensitive)
  json_files <- list.files(
    path = directory,
    pattern = "\\.json$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Read and process each JSON file using process_single_json (assumed defined elsewhere)
  df_list <- purrr::map(json_files, process_single_json)

  # Remove NULL elements and combine into a single data frame
  dplyr::bind_rows(purrr::compact(df_list)) %>%
    # Add unique identifier column "id" at the start
    dplyr::mutate(id = dplyr::row_number(), .before = 1)
}
