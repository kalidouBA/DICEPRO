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

#' NMF LBFGSB Hyperparameter Optimization
#'
#' Calls the R function `nmf_lbfgsb` to perform NMF optimization with given hyperparameters.
#'
#' @param dataset List containing matrices `B`, `W`, and `P`.
#' @param W_prime Optional numeric matrix. Initial W'.
#' @param p_prime Optional numeric matrix. Initial P'.
#' @param lambda_ Numeric. Regularization parameter lambda.
#' @param gamma Numeric. Regularization parameter gamma.
#' @param path2save Character. Path to save results.
#' @return List. Output from `nmf_lbfgsb` including loss, constraints, and matrices.
#' @import DICEPRO
#' @export
nmf_lbfgsb_hyperOpt <- function(dataset, W_prime = NULL, p_prime = NULL, lambda_ = 10, gamma = 100, path2save = "") {
  B <- as.data.frame(dataset$B)
  W_cb <- as.data.frame(dataset$W)
  P_cb <- as.data.frame(dataset$P)

  # Ensure p_prime is a numeric vector or matrix
  if (is.null(p_prime)) {
    N_sample <- ncol(B)
    N_unknownCT <- 1
    p_prime <- matrix(0.1, nrow = N_sample, ncol = N_unknownCT)
  } else if (is.numeric(p_prime) && is.vector(p_prime)) {
    N_sample <- ncol(B)
    N_unknownCT <- 1
    p_prime <- matrix(p_prime, nrow = N_sample, ncol = N_unknownCT)
  }

  if (nrow(P_cb) != ncol(B)) {
    P_cb <- t(P_cb)
  }

  r_dataset <- list(
    B = B,
    W_cb = W_cb,
    P_cb = P_cb
  )

  result <- nmf_lbfgsb(r_dataset = r_dataset, W_prime = W_prime, p_prime = p_prime, lambda_ = lambda_, gamma_par = gamma, path2save = path2save,
                       N_unknownCT = N_unknownCT)

  return(result)
}


#' Objective function for hyperparameter optimization
#'
#' Defines the objective function to minimize for hyperparameter tuning.
#'
#' @param dataset List. Dataset containing matrices `B`, `W`, and `P`.
#' @param config List. Configuration for the experiment, e.g., directories.
#' @param lambda_ Numeric. Regularization parameter lambda.
#' @param gamma_factor Numeric. Optional multiplier for gamma.
#' @param gamma Numeric. Regularization parameter gamma.
#' @param p_prime Numeric. Optional matrix for initialization.
#' @param W_prime Numeric. Optional matrix for initialization.
#' @return List. Contains loss, constraint, status, and current parameter values.
#' @export
objective_opt <- function(dataset, config = list(), lambda_ = NULL, gamma_factor = NULL, gamma = NULL, p_prime = NULL, W_prime = 0) {

  exp_dir <- ifelse(!is.null(config$exp), config$exp, ".")

  if (!is.null(gamma_factor)) {
    gamma <- lambda_ * gamma_factor
  }

  saveHpath <- exp_dir
  result <- nmf_lbfgsb_hyperOpt(dataset = dataset, W_prime = W_prime, p_prime = p_prime, lambda_ = lambda_, gamma = gamma, path2save = saveHpath)

  result_dict <- lapply(result, function(x) {
    if (is.numeric(x)) {
      x
    } else if (is.matrix(x) || is.data.frame(x)) {
      num_cols <- sapply(x, is.numeric)
      as.numeric(as.matrix(x[, num_cols, drop = FALSE]))
    } else {
      x
    }
  })



  if (contains_nan_or_inf(result_dict$loss) || contains_nan_or_inf(result_dict$constraint)) {
    return(list(
      loss = Inf,
      constraint = Inf,
      status = "fail",
      p_prime_estm = NULL,
      current_params = list(
        gamma = gamma,
        lambda_ = lambda_,
        p_prime = p_prime,
        frobNorm = Inf,
        constNorm = Inf,
        c1 = Inf,
        c2 = Inf,
        objectiveValue = Inf,
        penalty = Inf
      )
    ))
  }

  return(list(
    loss = result_dict$objectiveValue,
    constraint = result_dict$constraint,
    status = "OK",
    p_prime_estm = result_dict$p_prime,
    current_params = list(
      gamma = gamma,
      lambda_ = lambda_,
      p_prime = p_prime[1, 1],
      frobNorm = result_dict$frobNorm,
      constNorm = result_dict$constNorm,
      c1 = result_dict$c1,
      c2 = result_dict$c2,
      objectiveValue = result_dict$objectiveValue,
      penalty = result_dict$penalty
    )
  ))
}

#' Objective wrapper for hyperparameter optimization
#'
#' Wraps an objective function with error handling and automatic JSON result saving.
#'
#' @param objective_opt Function. Objective function to optimize. Must accept arguments `(dataset, config, params)` and return a list containing at least a `loss` element.
#' @param dataset List. Dataset used for optimization (e.g., bulk matrix, reference signatures, precomputed proportions).
#' @param config List. Parsed configuration object from JSON.
#' @param report_path Character. Directory path where JSON results are saved.
#' @param params List. Current hyperparameter set to evaluate.
#' @return List containing objective function output, status, start time, duration, and error message (if any).
#' @importFrom utils glob2rx
#' @keywords internal
objective_wrapper <- function(objective_opt, dataset, config, report_path, params, W_prime = NULL) {
  tryCatch({
    start_time <- Sys.time()

    # Pass W_prime to objective
    returned_dict <- objective_opt(dataset = dataset, config = config,
                               lambda_ = params$lambda_,
                               gamma_factor = params$gamma_factor,
                               gamma = params$gamma,
                               p_prime = params$p_prime,
                               W_prime = W_prime)

    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

    returned_dict$start_time <- as.numeric(start_time)
    returned_dict$duration <- duration

    save_file <- sprintf("%.7f_hyperopt_results", returned_dict$loss)
    params$p_prime <- params$p_prime[1,1]
    # Save results
    json_dict <- list(returned_dict = returned_dict, current_params = params)
    save_file_path <- file.path(report_path, save_file)
    existing_files <- list.files(report_path, pattern = glob2rx(paste0(save_file, "*")))
    save_file_path <- paste0(save_file_path, "_", length(existing_files) + 1, "call.json")

    jsonlite::write_json(json_dict, save_file_path, auto_unbox = TRUE, pretty = TRUE)

    return(returned_dict)
  }, error = function(e) {
    start_time <- Sys.time()
    returned_dict <- list(
      status = "FAIL",
      start_time = as.numeric(start_time),
      error = conditionMessage(e)
    )

    save_file <- paste0("ERR", as.numeric(start_time), "_hyperopt_results")
    save_file_path <- file.path(report_path, save_file)
    jsonlite::write_json(list(returned_dict = returned_dict), save_file_path, auto_unbox = TRUE)

    return(returned_dict)
  })
}

#' Hyperparameter optimization for DICEPRO
#'
#' Performs hyperparameter optimization using a specified objective function.
#' Supports user-defined search spaces or configuration-defined search spaces.
#'
#' @param objective_opt Function. Objective function to optimize. Must return a list containing a numeric `loss`.
#' @param dataset List. Dataset containing matrices such as bulk expression, reference signatures, and precomputed proportions.
#' @param config_path Character. Path to JSON configuration file specifying experiment settings (exp path, hp_max_evals, hp_method, hp_space, seed).
#' @param hp_space List (optional). User-defined hyperparameter search space. If `NULL`, parsed from configuration file.
#' @param report_path Character (optional). Path to store intermediate results. If `NULL`, defaults to `config$exp/results`.
#' @param W_prime Description of W_prime
#'
#' @return List with:
#' \describe{
#'   \item{best}{Best hyperparameter set found.}
#'   \item{trials}{List of all trial results, each containing `params` and `result`.}
#' }
#'
#' @details
#' 1. Load configuration from JSON.
#' 2. Create or parse hyperparameter search space.
#' 3. Loop for `hp_max_evals` iterations:
#'    - Sample a parameter set from the search space.
#'    - Evaluate objective function with `objective_wrapper()` (handles errors and saves JSON).
#'    - Track the best parameters based on the `loss`.
#' 4. Return the best parameters and full trial history.
#'
#' @examples
#' \dontrun{
#' config_path <- "path/to/optim_config.json"
#' dataset <- list(B = bulk_matrix, W = reference_matrix, P = precomputed_proportions)
#' result <- research_v2(my_objective_function, dataset, config_path)
#' }
#'
#' @export
research_v2 <- function(objective_opt, dataset, config_path, hp_space = NULL, report_path = NULL, W_prime = NULL) {
  # Load configuration
  config <- .get_conf_from_json(config_path)
  report_path <- .get_report_path(exp_dir = config$exp)

  # Create search space
  if (is.null(hp_space)) {
    search_space <- list()
    for (arg in names(config$hp_space)) {
      search_space[[arg]] <- .parse_hyperopt_searchspace(arg, config$hp_space[[arg]])
    }
  } else {
    search_space <- hp_space
  }

  # Optimization loop
  trials <- list()
  set.seed(config$seed %||% 42)

  # Progress bar
  pb <- progress::progress_bar$new(
    format = "Hyperopt [:bar] :percent eta: :eta",
    total = config$hp_max_evals,
    clear = FALSE,
    width = 60
  )

  for (i in 1:config$hp_max_evals) {
    params <- .sample_from_space(space = search_space)

    # Ensure p_prime is numeric and a matrix if needed
    if (!is.null(params$p_prime) && is.numeric(params$p_prime) && is.vector(params$p_prime)) {
      N_sample <- ncol(dataset$B)
      N_unknownCT <- ifelse(!is.null(W_prime), length(W_prime), 1)
      params$p_prime <- matrix(params$p_prime, nrow = N_sample, ncol = N_unknownCT)
    }

    # Pass W_prime to objective_wrapper
    result <- objective_wrapper(
      objective_opt   = objective_opt,
      dataset     = dataset,
      config      = config,
      report_path = report_path,
      params      = params,
      W_prime     = W_prime
    )

    trials[[i]] <- list(params = params, result = result)

    # Update progress bar
    pb$tick()
  }
}



#' Parse configuration from JSON file
#'
#' Internal function to parse hyperparameter optimization configuration
#'
#' @param config Configuration list
#' @return Parsed configuration list
#' @keywords internal
.parse_config <- function(config) {
  required_args <- c("exp", "hp_max_evals", "hp_method")

  for (arg in required_args) {
    if (is.null(config[[arg]])) {
      stop(paste("No", arg, "argument found in configuration file."))
    }
  }

  valid_methods <- c("tpe", "random", "atpe", "anneal")
  if (!config$hp_method %in% valid_methods) {
    stop(paste("Unknown hyperopt algorithm:", config$hp_method,
               "Available algorithms: 'random', 'tpe', 'atpe', 'anneal'."))
  }

  return(config)
}

#' Get configuration from JSON file
#'
#' Load and parse hyperparameter optimization configuration
#'
#' @param confpath Path to JSON configuration file
#' @return Configuration list
#' @keywords internal
.get_conf_from_json <- function(confpath) {
  if (!file.exists(confpath)) {
    stop(paste("Training conf", confpath, "not found."))
  }

  config <- jsonlite::read_json(confpath)
  return(.parse_config(config))
}

#' Parse hyperopt search space
#'
#' Internal function to parse hyperparameter search space definitions
#'
#' @param arg Parameter name
#' @param specs Specification list
#' @return Search space definition
#' @keywords internal
.parse_hyperopt_searchspace <- function(arg, specs) {
  type <- specs[[1]]

  switch(type,
         "choice" = list(type = "choice", choices = specs[-1]),
         "randint" = list(type = "randint", low = specs[[2]], high = specs[[3]]),
         "uniform" = list(type = "uniform", low = specs[[2]], high = specs[[3]]),
         "quniform" = list(type = "quniform", low = specs[[2]], high = specs[[3]], q = specs[[4]]),
         "loguniform" = list(type = "loguniform", low = specs[[2]], high = specs[[3]]),
         "qloguniform" = list(type = "qloguniform", low = specs[[2]], high = specs[[3]], q = specs[[4]]),
         "normal" = list(type = "normal", mu = specs[[2]], sigma = specs[[3]]),
         "qnormal" = list(type = "qnormal", mu = specs[[2]], sigma = specs[[3]], q = specs[[4]]),
         "lognormal" = list(type = "lognormal", mu = specs[[2]], sigma = specs[[3]]),
         "qlognormal" = list(type = "qlognormal", mu = specs[[2]], sigma = specs[[3]], q = specs[[4]]),
         stop(paste("Unknown search space type:", type))
  )
}

#' Get report path
#'
#' Create and return report directory path
#'
#' @param exp_name Experiment name
#' @param base_path Base path
#' @return Report path
#' @keywords internal
.get_report_path <- function(exp_dir) {
  # exp_dir: directory where experiment results are saved (config$exp)
  report_path <- file.path(exp_dir, "results")

  # Create all parent directories if needed
  if (!dir.exists(report_path)) {
    dir.create(report_path, recursive = TRUE)
  }
  return(report_path)
}
