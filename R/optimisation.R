#' Run a DICEPRO hyperparameter optimisation experiment
#'
#' Builds the hyperparameter search space from \code{hspaceTechniqueChoose},
#' and runs \code{research_hyperOpt()}, returning the collected trials.
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
                           output_base_dir = ".",
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

  # MODIF: config construit comme objet R — plus besoin de jsonlite pour
  # l'écrire sur disque puis le relire aussitôt après.
  # exp pointe toujours vers paths$data_dir pour que objective_opt() puisse
  # s'en servir comme répertoire de sauvegarde intermédiaire.
  base_config <- list(
    exp          = paths$data_dir,
    hp_max_evals = hp_max_evals,
    hp_method    = algo_select,
    seed         = 4L
  )

  raw_space <- switch(
    hspaceTechniqueChoose,
    all = list(
      lambda_ = c("loguniform", 1,    1e8),
      gamma   = c("loguniform", 1,    1e8),
      p_prime = c("loguniform", 1e-6, 1)
    ),
    restrictionEspace = .custom_space()
  )

  parsed_space <- lapply(
    stats::setNames(names(raw_space), names(raw_space)),
    function(arg) .parse_hyperopt_searchspace(arg, raw_space[[arg]])
  )

  # MODIF: on assemble la config finale en mémoire (raw_space inclus pour
  # validation dans .parse_config(), mais on ne la sérialise plus sur disque).
  hyperopt_config <- c(base_config, list(hp_space = raw_space))

  # MODIF: on passe `config` directement à research_hyperOpt() au lieu de
  # config_path. Le fichier JSON intermédiaire n'existe plus.
  research_hyperOpt(
    objective_opt = objective_opt,
    dataset       = dataset,
    config        = hyperopt_config,   # <-- objet R, plus de fichier JSON
    hp_space      = parsed_space,
    W_prime       = W_prime
  )
}


#' Define custom (restricted) hyperparameter search space
#'
#' Returns the \code{restrictionEspace} search space where \code{gamma} is
#' the base variable and \code{lambda_} is derived as
#' \code{lambda_ = gamma * lambda_factor}.
#'
#' @return Named list of atomic vectors \code{c(type, low, high)}.
#' @keywords internal
.custom_space <- function() {
  list(
    gamma         = c("loguniform", 1,    1e5),
    lambda_factor = c("loguniform", 2,    1e2),
    p_prime       = c("loguniform", 1e-1, 1)
  )
}
