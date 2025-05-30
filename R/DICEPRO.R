#' DICEPRO: Deconvolution and Inference of Cell Proportions with Optimization
#'
#' This function performs cell type deconvolution of bulk gene expression data using
#' a combination of a deconvolution method and Non-Negative Matrix Factorization
#' (NMF) with L-BFGS-B optimization. It can integrate known cell type signatures
#' and estimate the proportions of unknown cell types.
#'
#' @param reference A matrix of reference gene expression data, where rows are genes
#'                  and columns are cell types (cellType x gene).
#' @param bulk A matrix of bulk gene expression data, where rows are genes
#'             and columns are samples (gene x sample).
#' @param methodDeconv A character string specifying the deconvolution method to be used.
#'                     Options include: "CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism".
#' @param cibersortx_email An optional email address for the CIBERSORTx web-based tool.
#' @param cibersortx_token An optional token for the CIBERSORTx web-based tool.
#' @param W_prime A numeric value for the initial matrix W for unknown cell types,
#'                default is 0.
#' @param bulkName A string indicating the name of the bulk samples dataset.
#' @param refName A string indicating the name of the reference dataset.
#' @param hp_max_evals A numeric value indicating the maximum number of evaluations for optimization.
#'                     Default is 100.
#' @param N_unknownCT A numeric value for the number of unknown cell types, default is 1.
#' @param algo_select (character) Choice of optimization algorithm. Available options: `"random"` or `"tpe"` or `"atpe"`. Par défaut, `"random"`.
#' @param output_path Optional. A file path (character) to save the output results.
#'                    If NULL, results are written to the current working directory (`getwd()`).
#' @param hspaceTechniqueChoose Control the search space for hyperparameters. Available options: `"restrictionEspace"` or `"all"`. Par défaut, `"restrictionEspace"`.
#'
#' @return A list containing the following components:
#' \item{Prediction}{A matrix of the estimated cell type proportions for the bulk samples, excluding
#'                  the unknown cell type.}
#' \item{New_signature}{A matrix of the optimized gene-cell type matrix for the cell types, including
#'                     the unknown cell type.}
#' \item{W_prime}{The optimized matrix W for the unknown cell types.}
#' \item{P_prime}{The optimized proportions for the unknown cell types.}
#' \item{constraint}{The final constraint value from the optimization process.}
#'
#' @details The function performs cell type deconvolution using a specified method (e.g.,
#'          "CSx") followed by optimization using non-negative matrix factorization (NMF).
#'          The NMF part is performed using the `nmf_lbfgsb` function, where the parameters
#'          for the optimization (e.g., `W_prime`, `p_prime`, `lambda_`, etc.) can be customized.
#'          The optimization process is executed via a Python script using `reticulate` for
#'          interfacing with Python.
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
#' # Example usage of DICEPRO function
#' set.seed(2101)
#' data_simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE, nSample = 20, prop = NULL,
#'                               nGenes = 50, nCellsType = 30)
#'
#' result <- DICEPRO(reference = data_simulation$reference,
#'                   bulk = data_simulation$bulk,
#'                   methodDeconv = "CSx",
#'                   W_prime = 0,
#'                   bulkName = "Bulk_Sample",
#'                   refName = "Reference_Data",
#'                   hp_max_evals = 100,
#'                   N_unknownCT = 1)
#' }
#'
#' @export

DICEPRO <- function(reference, bulk, methodDeconv = "DCQ", cibersortx_email = NULL, cibersortx_token = NULL,
                    W_prime = 0, bulkName = "Bulk", refName = "Reference", hp_max_evals = 100, N_unknownCT = 1,
                    algo_select = "random", output_path = NULL, hspaceTechniqueChoose = "restrictionEspace") {

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

  out_Dec <- t(running_method(bulk, reference, methodDeconv, cibersortx_email, cibersortx_token))

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
  outputDir <- sprintf("%s/%s", output_path,dirName)
  bestHP <- try({
    plot_kraljic(outputDir)
  }, silent = TRUE)

  bestHP$bestEstimation <- read.table(
    sprintf("%s/optim/H_lambda_%s_gamma_%s.tsv", outputDir,
            round(bestHP$lambda_, 2), round(bestHP$gamma, 2)), header = T,row.names = 1
  )

  class(bestHP) <- "DICEPRO"
  return(bestHP)
}
