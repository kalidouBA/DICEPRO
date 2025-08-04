#' DICEPRO: Deconvolution and Inference of Cell Proportions with Optimization
#'
#' This function performs cell type deconvolution of bulk gene expression data using
#' a combination of a deconvolution method and Non-Negative Matrix Factorization
#' (NMF) with L-BFGS-B optimization. It can integrate known cell type signatures
#' and estimate the proportions of unknown cell types.
#'
#' @param reference A numeric matrix of reference gene expression data, where rows are genes
#'                  and columns are cell types (genes x cellTypes).
#' @param bulk A numeric matrix of bulk gene expression data, where rows are genes
#'             and columns are samples (genes x samples).
#' @param methodDeconv A character string specifying the deconvolution method to be used.
#'                     Options include: "CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism".
#' @param cibersortx_email (Optional) Email address for the CIBERSORTx web-based tool.
#' @param cibersortx_token (Optional) Token for the CIBERSORTx web-based tool.
#' @param W_prime A numeric value or matrix for the initial matrix W for unknown cell types.
#'                Default is 0.
#' @param bulkName A string specifying the name of the bulk dataset.
#' @param refName A string specifying the name of the reference dataset.
#' @param hp_max_evals A numeric value specifying the maximum number of evaluations
#'                     for the optimization algorithm. Default is 100.
#' @param N_unknownCT A numeric value indicating the number of unknown cell types to estimate. Default is 1.
#' @param algo_select A character string selecting the optimization algorithm. Options: "random", "tpe", "atpe". Default is "random".
#' @param output_path Optional character string specifying the output path to save results.
#'                    If NULL, results will be saved in the current working directory.
#' @param hspaceTechniqueChoose A character string to control the hyperparameter search space.
#'                              Options: "restrictionEspace", "all". Default is "restrictionEspace".
#' @param out_Decon Optional. A matrix of precomputed cell type proportions. If provided,
#'                  the `running_method()` step will be skipped, and this matrix will be used directly.
#'
#' @return A list of class `"DICEPRO"` containing:
#' \item{Prediction}{Estimated cell type proportions (excluding unknown cell types).}
#' \item{New_signature}{Optimized gene expression signatures including unknown cell types.}
#' \item{W_prime}{Optimized W matrix for unknown cell types.}
#' \item{P_prime}{Estimated proportions for unknown cell types.}
#' \item{constraint}{Final constraint value from the optimization.}
#'
#' @details
#' This function integrates deconvolution and optimization steps. If `out_Decon` is NULL,
#' the selected deconvolution method is run first via `running_method()`. If `out_Decon`
#' is provided, it is used directly as the predicted proportions matrix.
#' Then the optimization is performed using Python (via `reticulate`) through the script
#' `optimisation.py` bundled with the DICEPRO package.
#'
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rpois runif rnorm
#' @importFrom reticulate import_from_path r_to_py
#' @importFrom utils read.table
#'
## usethis namespace: start
#' @import RcppEigen
#' @useDynLib DICEPRO, .registration = TRUE
## usethis namespace: end
#'
#' @examples
#' \dontrun{
#' # Simulated data example
#' set.seed(1234)
#' sim_data <- simulation(loi = "gauss", scenario = "", bias = TRUE, nSample = 20,
#'                        prop = NULL, nGenes = 50, nCellsType = 30)
#'
#' result <- DICEPRO(reference = sim_data$reference,
#'                   bulk = sim_data$bulk,
#'                   methodDeconv = "CSx",
#'                   bulkName = "SimBulk",
#'                   refName = "SimRef")
#'
#' # Using a precomputed deconvolution matrix
#' precomputed_P <- matrix(runif(100), nrow = 10, ncol = 10)
#' rownames(precomputed_P) <- colnames(sim_data$bulk)
#' colnames(precomputed_P) <- colnames(sim_data$reference)
#'
#' result2 <- DICEPRO(reference = sim_data$reference,
#'                    bulk = sim_data$bulk,
#'                    out_Decon = precomputed_P,
#'                    bulkName = "SimBulk",
#'                    refName = "SimRef")
#' }
#'
#' @export


DICEPRO <- function(reference, bulk, methodDeconv = "DCQ", cibersortx_email = NULL, cibersortx_token = NULL,
                    W_prime = 0, bulkName = "Bulk", refName = "Reference", hp_max_evals = 100, N_unknownCT = 1,
                    algo_select = "random", output_path = NULL, hspaceTechniqueChoose = "restrictionEspace",
                    out_Decon = NULL) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  cellTypeName <- colnames(reference)
  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  if(is.matrix(reference)){
    stopifnot(is.numeric(reference))
    reference <- reference[geneIntersect, ]
  } else {
    reference <- apply(reference[geneIntersect, ], 2, as.numeric)
  }

  if(is.matrix(bulk)){
    stopifnot(is.numeric(bulk))
    bulk <- bulk[geneIntersect, ]
  } else {
    bulk <- apply(bulk[geneIntersect, ], 2, as.numeric)
  }

  rownames(reference) <- rownames(bulk) <- geneIntersect

  # === Choose whether to use running_method() or a provided out_Decon ===
  if (is.null(out_Decon)) {
    out_Dec <- t(running_method(bulk, reference, methodDeconv, cibersortx_email, cibersortx_token))
  } else {
    message("Utilisation de l'objet out_Decon fourni, sans appel à running_method().")
    out_Dec <- out_Decon
    if (!is.matrix(out_Dec)) {
      stop("L'objet 'out_Decon' doit être une matrice.")
    }
  }

  bulk_df <- as.data.frame(bulk)
  reference_df <- as.data.frame(reference)
  out_Dec_df <- as.data.frame(out_Dec)

  bulk_py <- reticulate::r_to_py(bulk_df)
  W_cb_py <- reticulate::r_to_py(reference_df)
  out_Dec_cb_py <- reticulate::r_to_py(out_Dec_df)

  dataset <- list('B' = bulk_py, 'W' = W_cb_py, 'P' = out_Dec_cb_py)

  script_path <- system.file("python/optimisation.py", package = "DICEPRO")
  suppressWarnings({
    optimisation <- reticulate::import_from_path("optimisation", path = dirname(script_path))
  })

  dirName <- paste0(bulkName, "_", refName, "_", hspaceTechniqueChoose, "_", algo_select)
  if (is.null(output_path)) {
    output_path <- getwd()
  }

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  # === Run optimization ===
  optimisation$run_experiment(dataset, bulkName, refName, hp_max_evals, algo_select, output_path, hspaceTechniqueChoose)

  # === Call the plotting function ===
  outputDir <- sprintf("%s/%s", output_path, dirName)
  bestHP <- try({
    plot_kraljic(outputDir)
  }, silent = TRUE)

  class(bestHP) <- "DICEPRO"
  return(bestHP)
}

