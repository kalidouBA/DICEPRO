#' Semi-Supervised Deconvolution of Bulk RNA-Seq Data
#'
#' This function performs semi-supervised deconvolution of bulk RNA-seq data using the CIBERSORTx or
#' Digital Cell Quantification method or others. It iteratively estimates cell type proportions and refines
#' reference data using non-negative matrix factorization (NMF).
#'
#' @param reference The reference data matrix containing cell type (column) gene expression (row) profiles.
#' @param bulk The bulk RNA-seq data matrix. In the columns we have the samples and in the rows the genes.
#' @param methodDeconv A character vector specifying the deconvolution method.
#'        Supported values are "CSx" or "DCQ" or others. The method to use for deconvolution. Options include \code{CSx},
#'        \code{DCQ}, \code{CDSeq}, \code{DeconRNASeq}, \code{FARDEEP} and \code{BayesPrism}.
#'
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @seealso Other functions used within deconvolution: \code{\link{running_method}}, \code{\link{nmf}}.
#'
#' @return A list data simulated.
#'
#' @details The function calculates and returns the simulations:
#' \itemize{
#'   \item \code{Prediction}: A matrix or data frame containing the estimated cell type proportions for each sample in the bulk RNA-Seq data.
#'   \item \code{Matrix_prediction}: A matrix or data frame containing the cell type proportions for each sample in the bulk RNA-Seq data, with iteration information.
#'   \item \code{New_signature}: The refined reference data matrix with unknown components.
#'   \item \code{Optimal_iteration}: The optimal iteration at which the convergence criteria were met.
#'   \item \code{performs2plot}: The variation performance between iteration.
#' }
#'
#' @details This function performs deconvolution of bulk RNA-Seq data using either
#' the CSx or DCQ method. It first prepares the data and runs the selected method
#' for cell type proportion estimation. The function then calculates error norms,
#' performs cross-validation.
#' It estimates unknown components using Non-Negative Matrix Factorization (NMF)
#' and calculates distances between unknown components and known cell types.
#' The function returns the estimated cell type proportions and their associated matrices.
#'
#'
#' The function performs several steps for semi-supervised deconvolution, including:
#' \enumerate{
#'   \item Running CIBERSORTx or Digital Cell Quantification or others for cell type proportion estimation.
#'   \item Estimating unknown components using non-negative matrix factorization (NMF).
#'}
#'
#' @export
#'
#' @importFrom NMF nmf basis
#' @importFrom ComICS dcq
#' @importFrom parallel mclapply
#' @import dplyr
#'
#' @examples
#' if(interactive()){
#' simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE,
#' nSample = 10, prop = NULL, nGenes = 50, nCellsType = 5)
#' cellTypeOut <- sample(1:ncol(simulation$reference), 2)
#' refDataIncomplet <- simulation$reference[,-cellTypeOut]
#' results <- DICEPRO(reference = refDataIncomplet, bulk = simulation$bulk, methodDeconv = "DCQ")
#' print(results)
#' }

DICEPRO <- function(reference, bulk, methodDeconv = "CSx", cibersortx_email = NULL, cibersortx_token = NULL) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  cellTypeName <- colnames(reference)
  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  reference <- apply(reference[geneIntersect, ], 2, as.numeric)
  bulk <- as.data.frame(apply(bulk[geneIntersect, ], 2, as.numeric))
  rownames(reference) <- rownames(bulk) <- geneIntersect

  matrixAbundances <- performs <- normFrobs <- performs2plot <- opt <- NULL

  out_Dec <- running_method(bulk, W, methodDeconv, cibersortx_email, cibersortx_token)
  B_Deconv <- as.matrix(W) %*% t(out_Dec)

  k_CT <- ncol(W)

  res <- nmf_conjugate_gradient(V = bulk, W = W, H = out_Dec, k_CT+1)
  W <- res$W

  colnames(W) <- c(cellTypeName, paste0("Unknown"))

  out_Dec_Update <- res$H[,1:ncol(out_Dec)]
  dimnames(out_Dec_Update) <- dimnames(out_Dec)
  rownames(performs2plot) <- NULL

  results <- list("Prediction" = out_Dec_Update, "New_signature" = W)

  class(results) <- "DICEPRO"
  return(results)
}

