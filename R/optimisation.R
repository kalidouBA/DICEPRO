#' Run hyperparameter optimization experiment
#'
#' Performs hyperparameter optimization, saves trials and results, and generates reports.
#'
#' @param dataset List. Dataset containing matrices `B`, `W`, and `P`.
#' @param bulkName Character. Name identifier for bulk data.
#' @param refName Character. Name identifier for reference data.
#' @param hp_max_evals Integer. Maximum number of hyperparameter evaluations.
#' @param algo_select Character. Algorithm selection, e.g., "random".
#' @param output_base_dir Character. Base directory for saving outputs.
#' @param hspaceTechniqueChoose Character. Hyperparameter space technique ("all" or "restrictionEspace").
#' @return None. Saves JSON files and plots to disk.
#' @export
run_experiment <- function(dataset, W_prime = 0, bulkName = "", refName = "", hp_max_evals = 100,
                           algo_select = "random", output_base_dir = ".", hspaceTechniqueChoose = "restrictionEspace") {

  paths <- generate_experiment_paths(output_base_dir, bulkName, refName, hspaceTechniqueChoose, algo_select)
  list2env(paths, envir = environment())

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
  } else if (hspaceTechniqueChoose == "restrictionEspace") {
    hyperopt_config <- list(
      exp = paths$optim_dir,
      hp_max_evals = hp_max_evals,
      hp_method = algo_select,
      seed = 4
    )
    hp_space <- custom_space()
  } else {
    stop(sprintf("Technique '%s' not available. Choose 'all' or 'restrictionEspace'.", hspaceTechniqueChoose))
  }

  jsonlite::write_json(hyperopt_config, paths$config_path, auto_unbox = TRUE, pretty = TRUE)

  research_v2(objective, dataset, config_path = paths$config_path, hp_space = hyperopt_config$hp_space,
              report_path = paths$report_dir, W_prime = W_prime)

  fp_plot <- plot_kraljic(outputDir = paths$data_dir, mixName = bulkName, refName = refName)

  # Access optimal points
  optimal_DiceproOptCstrt <- fp_plot$optimal_points$DiceproOptCstrt
  optimal_DiceproOptCstrt_0.1 <- fp_plot$optimal_points$DiceproOptCstrt_0.1

  # Access corresponding TSV data
  tsv_data_DiceproOptCstrt <- fp_plot$optimal_tsv_data$DiceproOptCstrt
  tsv_data_DiceproOptCstrt_0.1 <- fp_plot$optimal_tsv_data$DiceproOptCstrt_0.1

  if (hspaceTechniqueChoose == "all") {
    fig_constraint <- plot_hyperopt_report_v2(optim_dir, c("gamma", "lambda_", "p_prime"), metric = "constraint")
  } else {
    fig_constraint <- plot_hyperopt_report_v2(optim_dir, c("gamma_factor", "lambda_", "p_prime"), metric = "constraint")
  }

  ggsave(report_path, fig_constraint)

  # Return a comprehensive list of results
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
#' Returns a list of functions generating random hyperparameter values.
#'
#' @return List of functions for `lambda_`, `gamma_factor`, and `p_prime`.
#' @export
custom_space <- function() {
  list(
    lambda_ = function() exp(runif(1, log(1), log(1e5))),
    gamma_factor = function() exp(runif(1, log(2), log(1e2))),
    p_prime = function() exp(runif(1, log(1e-1), log(1)))
  )
}
