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
#' @param p_prime A numeric value for the initial proportions of unknown cell types,
#'                default is 0.
#' @param lambda_ A numeric value for the regularization parameter, default is 10.
#' @param gamma_par A numeric value for the penalty parameter, default is 100.
#' @param N_unknownCT A numeric value for the number of unknown cell types, default is 1.
#'
#' @return A list containing the following components:
#' \item{Prediction}{The estimated cell type proportions for the bulk samples, excluding
#'                  the unknown cell type.}
#' \item{New_signature}{The optimized gene-cell type matrix for the cell types, including
#'                     the unknown cell type.}
#'
#' @details The function performs cell type deconvolution using a specified method (e.g.,
#'          "CSx") followed by optimization using non-negative matrix factorization (NMF).
#'          The NMF part is performed using the `nmf_lbfgsb` function, where the parameters
#'          for the optimization (e.g., `W_prime`, `p_prime`, `lambda_`, etc.) can be customized.
#'
#' @importFrom stats optim
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rpois runif rnorm
#'
## usethis namespace: start
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
#'                   p_prime = 0,
#'                   lambda_ = 10,
#'                   gamma_par = 100,
#'                   N_unknownCT = 1)
#' }
#'
#' @export

DICEPRO <- function(reference, bulk, methodDeconv = "CSx", cibersortx_email = NULL, cibersortx_token = NULL,
                    W_prime = 0, p_prime = 0, lambda_ = 10, gamma_par = 100, N_unknownCT = 1) {

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
  k_CT <- ncol(reference) + 1

  res <- nmf_lbfgsb(r_dataset = list(B = bulk, P_cb = out_Dec, W_cb = reference),
                    W_prime = W_prime, p_prime = p_prime, lambda_ = lambda_,
                    gamma_par = gamma_par, N_unknownCT = N_unknownCT)

  out_Dec_Update <- res$H
  dimnames(out_Dec_Update) <- dimnames(out_Dec)

  results <- list("Prediction" = out_Dec_Update, "W_prime" = res$w, "P_prime" = res$p_prime, "constraint" = res$constraint)

  class(results) <- "DICEPRO"
  return(results)
}
