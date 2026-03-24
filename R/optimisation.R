#' Run a DICEPRO hyperparameter optimisation experiment
#'
#' Builds the hyperparameter search space from \code{hspaceTechniqueChoose},
#' writes the configuration to disk, runs \code{research_hyperOpt()}, and returns
#' the collected trials.
#'
#' @param dataset              List containing at least \code{$B}, \code{$W},
#'   and \code{$P} matrices.
#' @param W_prime              Numeric matrix (or \code{0}). Passed verbatim to
#'   \code{research_hyperOpt()} as the initial \eqn{W} for NMF.
#' @param bulkName             Character scalar. Identifier for the bulk dataset
#'   (used in output path construction).
#' @param refName              Character scalar. Identifier for the reference
#'   dataset (used in output path construction).
#' @param hp_max_evals         Positive integer. Number of hyperparameter trials
#'   to run.
#' @param algo_select          Character scalar. Sampling algorithm passed to
#'   the config (e.g. \code{"random"}).
#' @param output_base_dir      Character scalar. Root directory for all outputs
#'   (default \code{"."}).
#' @param hspaceTechniqueChoose Character scalar. Search-space strategy:
#'   \describe{
#'     \item{\code{"all"}}{Full independent log-uniform grid for
#'       \code{lambda_}, \code{gamma}, \code{p_prime}.}
#'     \item{\code{"restrictionEspace"}}{Restricted space where \code{gamma}
#'       is the base variable and \code{lambda_} is derived as
#'       \code{lambda_ = gamma * lambda_factor} with
#'       \code{lambda_factor} in (2, 100).}
#'   }
#'
#' @return The list returned by \code{\link{research_hyperOpt}}: \code{trials},
#'   \code{W}, and \code{H}.
#'
#' @export
run_experiment <- function(dataset,
                           W_prime = 0,
                           bulkName,
                           refName,
                           hp_max_evals,
                           algo_select,
                           output_base_dir= ".",
                           hspaceTechniqueChoose) {

  # ---- Output paths ----------------------------------------------------------
  paths <- .generate_experiment_paths(
    output_base_dir,
    bulkName,
    refName
  )

  # ---- Search-space configuration --------------------------------------------
  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("all", "restrictionEspace")
  )

  base_config <- list(
    exp          = paths$data_dir,
    hp_max_evals = hp_max_evals,
    hp_method    = algo_select,
    seed         = 4L
  )

  # FIX: both spaces are now defined as named lists compatible with
  # .parse_hyperopt_searchspace() / .sample_from_space() in hypersearch.R.
  #
  # "all"               — lambda_ and gamma are independent (no constraint).
  # "restrictionEspace" — gamma is the base variable; lambda_ is derived
  #                       inside .sample_from_space() as:
  #                           lambda_ = gamma * lambda_factor
  #                       lambda_factor in [2, 100] ensures lambda_ >= 2*gamma.
  raw_space <- switch(
    hspaceTechniqueChoose,
    all = list(
      lambda_ = c("loguniform", 1,    1e8),
      gamma   = c("loguniform", 1,    1e8),
      p_prime = c("loguniform", 1e-6, 1)
    ),
    restrictionEspace = .custom_space()   # list(gamma, lambda_factor, p_prime)
  )

  # Convert the raw atomic-vector specs to the named-list format that
  # .sample_from_space() in hypersearch.R expects.
  parsed_space <- lapply(
    stats::setNames(names(raw_space), names(raw_space)),
    function(arg) .parse_hyperopt_searchspace(arg, raw_space[[arg]])
  )

  hyperopt_config <- c(base_config, list(hp_space = raw_space))

  # ---- Persist configuration -------------------------------------------------
  jsonlite::write_json(hyperopt_config,
                       path       = paths$config_path,
                       auto_unbox = TRUE,
                       pretty     = TRUE)

  # ---- Run search ------------------------------------------------------------
  # Pass parsed_space directly so research_hyperOpt() skips re-parsing from
  # the JSON file (which would see the raw atomic-vector format).
  research_hyperOpt(
    objective_opt = objective_opt,
    dataset       = dataset,
    config_path   = paths$config_path,
    hp_space      = parsed_space,          # FIX: pass already-parsed space
    W_prime       = W_prime
  )
}


#' Define custom (restricted) hyperparameter search space
#'
#' Returns the \code{restrictionEspace} search space where \code{gamma} is
#' the base variable and \code{lambda_} is derived as
#' \code{lambda_ = gamma * lambda_factor}.
#'
#' Compatible with both:
#' \itemize{
#'   \item \code{.sample_from_space()} in \code{lambda_gamma_initialize.R}
#'     (expects atomic vectors).
#'   \item \code{run_experiment()} which converts these to named lists before
#'     passing to \code{research_hyperOpt()}.
#' }
#'
#' @return Named list of atomic vectors \code{c(type, low, high)}.
#' @keywords internal
.custom_space <- function() {
  list(
    gamma         = c("loguniform", 1,    1e5),   # base variable
    lambda_factor = c("loguniform", 2,    1e2),   # min=2: lambda_ >= 2*gamma
    p_prime       = c("loguniform", 1e-1, 1)      # log-scale, median ~0.31
  )
}
