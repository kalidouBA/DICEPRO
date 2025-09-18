# Déclarer les variables globales pour éviter les notes "no visible binding"
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "lambda_", "gamma", "frobNorm", "constraint", "abs_constraint",
    "p_prime", "duration", "penalty", "objectiveValue",
    "scaled_constraint", "optim_dir", "report_path",
    "extract_final_results", "plot_hyperopt_report_v2"
  ))
}

# Importer explicitement les fonctions de stats qui étaient non reconnues
# (évite les notes sur cor, rpois, etc.)
#' @importFrom stats cor rpois cov var weighted.mean setNames
#' @importFrom utils read.table
NULL
