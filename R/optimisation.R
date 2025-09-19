#' Redirect output to a log file
#'
#' Redirects standard output and messages to a specified log file.
#'
#' @param log_file Character. Path to the log file.
#' @return None. Output is redirected to the log file.
#' @export
redirect_output_to_file <- function(log_file) {
  sink(log_file, append = TRUE, type = "output")
  sink(log_file, append = TRUE, type = "message")
}

#' Reset output to console
#'
#' Resets standard output and messages to the R console.
#'
#' @return None. Output is reset.
#' @export
reset_output <- function() {
  sink(type = "output")
  sink(type = "message")
}

#' Check if a value contains NaN or Inf
#'
#' Checks whether a numeric value, vector, matrix, or data frame contains NaN or Inf.
#'
#' @param value Numeric, vector, matrix, or data frame to check.
#' @return Logical. TRUE if NaN or Inf is present, FALSE otherwise.
#' @export
contains_nan_or_inf <- function(value) {
  if (is.numeric(value)) {
    return(any(is.nan(value) | is.infinite(value)))
  } else if (is.matrix(value) || is.data.frame(value)) {
    return(any(is.nan(as.matrix(value)) | is.infinite(as.matrix(value))))
  }
  return(FALSE)
}

#' Run hyperparameter optimization experiment
#'
#' Performs hyperparameter optimization, saves trials and results, and generates reports.
#'
#' @param dataset List. Dataset containing matrices B, W, and P.
#' @param W_prime Numeric. Additional parameter for objective function.
#' @param bulkName Character. Name identifier for bulk data.
#' @param refName Character. Name identifier for reference data.
#' @param hp_max_evals Integer. Maximum number of hyperparameter evaluations.
#' @param algo_select Character. Algorithm selection, e.g., "random".
#' @param output_base_dir Character. Base directory for saving outputs.
#' @param hspaceTechniqueChoose Character. Hyperparameter space technique ("all" or "restrictionEspace").
#' @return A list with results, plots, and TSV data.
#' @export
run_experiment <- function(dataset, W_prime = 0, bulkName = "", refName = "",
                           hp_max_evals = 100, algo_select = "random",
                           output_base_dir = ".", hspaceTechniqueChoose = "restrictionEspace") {

  paths <- generate_experiment_paths(output_base_dir, bulkName, refName,
                                     hspaceTechniqueChoose, algo_select)
  list2env(paths, envir = environment())

  # ---- Hyperparameter configuration ----
  if (hspaceTechniqueChoose == "all") {
    hyperopt_config <- list(
      exp = paths$optim_dir,
      hp_max_evals = hp_max_evals,
      hp_method = algo_select,
      seed = 4,
      hp_space = list(
        lambda_ = c("loguniform", 1, 1e+8),
        gamma = c("loguniform", 1, 1e+8),
        p_prime = c("loguniform", 1e-6, 1)
      )
    )
    hp_params <- c("gamma", "lambda_", "p_prime")
  } else if (hspaceTechniqueChoose == "restrictionEspace") {
    hyperopt_config <- list(
      exp = paths$optim_dir,
      hp_max_evals = hp_max_evals,
      hp_method = algo_select,
      seed = 4,
      hp_space = custom_space()
    )
    hp_params <- c("gamma_factor", "lambda_", "p_prime")
  } else {
    stop(sprintf("Technique '%s' not available. Choose 'all' or 'restrictionEspace'.",
                 hspaceTechniqueChoose))
  }

  # ---- Save config JSON ----
  jsonlite::write_json(hyperopt_config, paths$config_path, auto_unbox = TRUE, pretty = TRUE)


  # ---- Run hyperparameter search ----
  research_v2(objective_opt, dataset, config_path = paths$config_path, hp_space = hyperopt_config$hp_space,
              report_path = paths$report_dir, W_prime = W_prime)

  # ---- Plot Kraljic ----
  fp_plot <- plot_kraljic(outputDir = paths$data_dir, mixName = bulkName, refName = refName)

  # Points optimaux
  optimal_DiceproOptCstrt <- fp_plot$optimal_points$DiceproOptCstrt
  optimal_DiceproOptCstrt_0.1 <- fp_plot$optimal_points$DiceproOptCstrt_0.1

  # TSV correspondants
  tsv_data_DiceproOptCstrt <- fp_plot$optimal_tsv_data$DiceproOptCstrt
  tsv_data_DiceproOptCstrt_0.1 <- fp_plot$optimal_tsv_data$DiceproOptCstrt_0.1

  # ---- Generate hyperparameter optimization report ----
  fig_constraint <- plot_hyperopt_report_v2(
    exp = paths$optim_dir,
    params = hp_params,
    metric = "constraint",
    loss_metric = "constraint",
    loss_behaviour = "min",
    not_log = NULL,
    categorical = NULL,
    max_deviation = NULL,
    title = paste0("Hyperparameter Optimization - ", bulkName)
  )

  # Save report plot
  ggsave(
    filename = file.path(paths$report_dir, paste0("hyperopt_report_", bulkName, ".png")),
    plot = fig_constraint,
    width = 14, height = 10
  )

  return(list(
    optimal_points = list(
      DiceproOptCstrt = optimal_DiceproOptCstrt,
      DiceproOptCstrt_0.1 = optimal_DiceproOptCstrt_0.1
    ),
    optimal_tsv_data = list(
      DiceproOptCstrt = tsv_data_DiceproOptCstrt,
      DiceproOptCstrt_0.1 = tsv_data_DiceproOptCstrt_0.1
    ),
    pareto_frontier = fp_plot$frontier_points,
    plot = fp_plot$plot,
    hyperopt_plot = fig_constraint
  ))
}

#' Define custom hyperparameter search space
#'
#' Returns a restricted hyperparameter space compatible with hyperopt/research_v2.
#'
#' @return List of vectors for lambda_, gamma_factor, and p_prime.
#' @export
custom_space <- function() {
  list(
    lambda_      = c("loguniform", 1, 1e5),
    gamma_factor = c("loguniform", 2, 1e2),
    p_prime      = c("loguniform", 1e-1, 1)
  )
}
