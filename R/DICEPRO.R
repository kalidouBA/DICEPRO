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
#' @importFrom stats optim
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rpois runif rnorm
#' @importFrom reticulate import_from_path
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

DICEPRO <- function(reference, bulk, methodDeconv = "CSx", cibersortx_email = NULL, cibersortx_token = NULL,
                    W_prime = 0, bulkName = "", refName = "", hp_max_evals = 100, N_unknownCT = 1) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  cellTypeName <- colnames(reference)
  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  if(is.matrix(reference)){
    stopifnot(is.numeric(reference))
    reference <- reference[geneIntersect, ]
  }else{
    reference <- apply(reference[geneIntersect, ], 2, as.numeric)
    reference <- as.matrix(reference)
  }

  if(is.matrix(bulk)){
    stopifnot(is.numeric(bulk))
    bulk <- bulk[geneIntersect, ]
  }else{
    bulk <- apply(bulk[geneIntersect, ], 2, as.numeric)
    bulk <- as.matrix(bulk)
  }

  rownames(reference) <- rownames(bulk) <- geneIntersect

  out_Dec <- t(running_method(bulk, reference, methodDeconv, cibersortx_email, cibersortx_token))

  # Python-based optimization part
  script_path <- system.file("python/optimisation.py", package = "DICEPRO")
  optimisation <- import_from_path("optimisation", path = dirname(script_path))

  bulk <- convert_matrix_to_df(bulk)
  print(dimnames(bulk))
  reference <- convert_matrix_to_df(reference)
  print(dimnames(reference))
  out_Dec <- convert_matrix_to_df(out_Dec)
  print(dimnames(out_Dec))

  bulk_py <- r_to_py(bulk)
  W_cb_py <- r_to_py(reference)
  out_Dec_cb_py <- r_to_py(out_Dec)

  optimisation$run_experiment(bulk_py, W_cb_py, out_Dec_cb_py, bulkName, refName, hp_max_evals)
  pathDir <- paste0("dataTestNMFOptHyper/", bulkName, "_", refName)

  best_hyperOpt <- paste0(pathDir, "/best_hyperparameters_", bulkName, "_", refName, ".json")
  r_dataset <- list('B' = bulk, "P_cb" = out_Dec, "W_cb" = reference)
  if(file.exists(best_hyperOpt)){
    bestHP <- fromJSON(best_hyperOpt)
    resultDicepro <- nmf_lbfgsb(r_dataset, W_prime, p_prime = bestHP$p_prime, lambda_ = bestHP$lambda_, gamma_par = bestHP$gamma)
    results <- list("Prediction" = out_Dec_Update, "W_prime" = resultDicepro$w, "P_prime" = resultDicepro$p_prime, "constraint" = resultDicepro$constraint)
    class(results) <- "DICEPRO"
    return(results)
  }
}
