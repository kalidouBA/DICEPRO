#' Generate a Kraljic Pareto plot without ground-truth-based metrics
#'
#' This function generates an interactive Plotly figure showing the Pareto frontier
#' based on Frobenius norm and constraint deviation, excluding performance indicators
#' calculated against ground truth. It also identifies optimal points based on
#' different constraint criteria and includes corresponding TSV data.
#'
#' The function reads input data from a folder derived from dataset and algorithm names,
#' processes the data, computes the Pareto frontier using the \code{KraljicMatrix} package,
#' and saves an interactive HTML file. It returns the frontier points and optimal solutions.
#'
#' @param outputDir Character string. Name of the directory saving output.
#' @param mixName Character string. Name of the mixture dataset (for labeling).
#' @param refName Character string. Name of the reference dataset (for labeling).
#'
#' @return A list containing:
#'   \item{frontier_points}{Dataframe containing all points on the Pareto frontier}
#'   \item{optimal_points}{List with two optimal points:
#'     \itemize{
#'       \item{DiceproOptCstrt: Point with constraint closest to 1}
#'       \item{DiceproOptCstrt_0.1: Point minimizing Frobenius norm with constraint within ±10%}
#'     }
#'   }
#'   \item{optimal_tsv_data}{List with TSV data corresponding to optimal points}
#'   \item{plot}{The generated plotly object}
#'
#' @details
#' The function performs the following steps:
#' 1. Reads and preprocesses JSON performance data
#' 2. Identifies the Pareto frontier using KraljicMatrix::get_frontier()
#' 3. Finds optimal points based on constraint criteria
#' 4. Loads corresponding TSV files for optimal points
#' 5. Generates an interactive plotly visualization
#' 6. Saves HTML plots to the specified directory
#'
#' @examples
#' \dontrun{
#' result <- plot_kraljic(
#'   outputDir = "results/",
#'   mixName = "mix1",
#'   refName = "ref1"
#' )
#'
#' # Access optimal points
#' optimal_DiceproOptCstrt <- result$optimal_points$DiceproOptCstrt
#' optimal_DiceproOptCstrt_0.1 <- result$optimal_points$DiceproOptCstrt_0.1
#'
#' # Access corresponding TSV data
#' tsv_data_DiceproOptCstrt <- result$optimal_tsv_data$DiceproOptCstrt
#' tsv_data_DiceproOptCstrt_0.1 <- result$optimal_tsv_data$DiceproOptCstrt_0.1
#' }
#'
#' @importFrom dplyr mutate filter slice row_number
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom htmlwidgets saveWidget
#' @importFrom KraljicMatrix get_frontier
#' @importFrom purrr map compact
#' @importFrom dplyr bind_rows
#' @importFrom readr read_tsv
#' @export
plot_kraljic <- function(outputDir, mixName, refName) {
  json_dir <- paste0(outputDir, "/optim/results")
  tsv_dir <- paste0(outputDir, "/optim")  # Directory containing TSV files

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

  # Finding the optimal points
  optimal_points <- list()
  optimal_tsv_data <- list()  # To store corresponding TSV data

  # Point with constraint closest to 1
  optimal_points$DiceproOptCstrt <- data2plotfrontier_result %>%
    filter(abs_constraint == min(abs_constraint)) %>%
    slice(1) %>%
    mutate(
      type = "DiceproOptCstrt",
      mix = mixName,
      ref = refName
    )

  # Point that minimizes the Frobenius norm with a constraint of ±10%
  optimal_points$DiceproOptCstrt_0.1 <- data2plotfrontier_result %>%
    filter(scaled_constraint == 1) %>%
    filter(constraint >= 0.9 & constraint <= 1.1) %>%
    filter(frobNorm == min(frobNorm)) %>%
    slice(1) %>%
    mutate(
      type = "DiceproOptCstrt_0.1",
      mix = mixName,
      ref = refName
    )

  # Load corresponding TSV files for optimal points
  load_tsv_data <- function(optimal_point, point_name) {
    # Construct the expected TSV filename pattern
    lambda_val <- optimal_point$lambda_
    gamma_val <- optimal_point$gamma

    # Format the values to match the filename pattern
    lambda_formatted <- formatC(lambda_val, format = "e", digits = 2)
    gamma_formatted <- formatC(gamma_val, format = "e", digits = 2)

    # Create the filename pattern to search for
    pattern <- sprintf("H_lambda_%s_gamma_%s",
                       gsub("\\.", "_", lambda_formatted),
                       gsub("\\.", "_", gamma_formatted))

    # Find matching TSV files
    tsv_files <- list.files(
      path = tsv_dir,
      pattern = paste0("^", pattern, ".*\\.tsv$"),
      full.names = TRUE
    )

    if (length(tsv_files) > 0) {
      # Read the first matching TSV file
      tsv_data <- readr::read_tsv(tsv_files[1], show_col_types = FALSE)
      return(tsv_data)
    } else {
      warning(paste("No TSV file found for", point_name,
                    "with lambda =", lambda_val, "gamma =", gamma_val))
      return(NULL)
    }
  }

  # Load TSV data for both optimal points
  optimal_tsv_data$DiceproOptCstrt <- load_tsv_data(optimal_points$DiceproOptCstrt, "DiceproOptCstrt")
  optimal_tsv_data$DiceproOptCstrt_0.1 <- load_tsv_data(optimal_points$DiceproOptCstrt_0.1, "DiceproOptCstrt_0.1")

  # Plot (unchanged from your original code)
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

  # Return results including TSV data
  return(list(
    frontier_points = data2plotfrontier_result,
    optimal_points = optimal_points,
    optimal_tsv_data = optimal_tsv_data,
    plot = configured_fig
  ))
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
