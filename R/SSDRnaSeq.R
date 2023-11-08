#' Semi-Supervised Deconvolution of Bulk RNA-Seq Data
#'
#' This function performs semi-supervised deconvolution of bulk RNA-seq data using the CIBERSORTx or
#' Digital Cell Quantification method or others. It iteratively estimates cell type proportions and refines
#' reference data using non-negative matrix factorization (NMF).
#'
#' @param reference The reference data matrix containing cell type gene expression profiles.
#' @param bulk The bulk RNA-seq data matrix.
#' @param k_folds The number of cross-validation folds.
#' @param nIteration The number of iterations to perform for refining the reference data.
#' @param methodDeconv A character vector specifying the deconvolution method.
#'        Supported values are "CSx" or "DCQ" or others. The method to use for deconvolution. Options include \code{CSx},
#'        \code{DCQ}, \code{CDSeq}, \code{DeconRNASeq}, \code{FARDEEP} and \code{BayesPrism}.
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @seealso Other functions used within deconvolution: \code{\link{running_method}}, \code{\link{createFolds}}, \code{\link{nmf}}, \code{\link{compute_distances}}.
#'
#' @return A list data simulated.
#'
#' @details The function calculates and returns the simulations:
#' \itemize{
#'   \item \code{Prediction}: A matrix or data frame containing the estimated cell type proportions for each sample in the bulk RNA-Seq data.
#'   \item \code{Matrix_prediction}: A matrix or data frame containing the cell type proportions for each sample in the bulk RNA-Seq data, with iteration information.
#'}
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
#'   \item Calculating and storing performance metrics.
#'   \item Estimating unknown components using non-negative matrix factorization (NMF).
#'   \item Computing distances between unknown components and known cell types.
#'   \item Iteratively refining the reference data with unknown components.
#'}
#'
#' @export
#'
#' @importFrom caret createFolds
#' @importFrom NMF nmf basis
#' @importFrom reshape2 melt
#' @importFrom ComICS dcq
#' @import dplyr
#'
#' @examples
#' if(interactive()){
#'
#' simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE, nSample = 10, prop = NULL,
#'                          nGenes = 50, nCellsType = 5)
#' cellTypeOut <- sample(1:ncol(simulation$reference), 1)
#' refDataIncomplet <- simulation$reference[,-cellTypeOut]
#' results <- SSDRnaSeq(reference = refDataIncomplet, bulk = simulation$bulk, k_folds = 2,
#'                      nIteration = 2, methodDeconv = "DCQ")
#' print(results)
#'
#'}

SSDRnaSeq <- function(reference, bulk, k_folds = 5, nIteration = 1, methodDeconv = "CSx",
                      cibersortx_email=NULL, cibersortx_token = NULL) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))
  stopifnot(nIteration > 0)
  stopifnot(ncol(reference) > k_folds)

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  cellTypeName <- colnames(reference)
  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  reference <- apply(reference[geneIntersect, ], 2, as.numeric)
  bulk <- apply(bulk[geneIntersect, ], 2, as.numeric)
  rownames(reference) <- rownames(bulk) <- geneIntersect

  # Initialize variables to store results and metrics
  errorFrob <- list()
  matrixAbundances <- BetweenUnknown <- NULL

  if(nIteration > 0){
    for (iterate_ in 1:nIteration) {
      out_Dec <- running_method(bulk, reference, methodDeconv,  cibersortx_email, cibersortx_token)
      # Calculate error norms
      bulkDeconv <- as.matrix(reference) %*% out_Dec
      diff_bulk <- bulk - bulkDeconv
      errorFrob[[iterate_]] <- normFrob(diff_bulk, bulk)

      # Create folds for cross-validation
      flds <- createFolds(1:ncol(bulk), k = k_folds, list = TRUE, returnTrain = FALSE)

      # Initialize a matrix to store unknown components
      matUnknown <- matrix(ncol = k_folds, nrow = nrow(reference))

      for (indFld in 1:k_folds) {
        resNMF <- nmf(as.matrix(abs(diff_bulk)[, flds[[indFld]]]), rank = 1)
        matUnknown[, indFld] <- as.vector(basis(resNMF))
      }

      colnames(matUnknown) <- names(flds)

      # Compute distances between unknown components and known cell types
      resDist <- compute_distances(matUnknown)

      # Store distance results
      resDist <- melt(resDist[[1]]) %>% cbind.data.frame(
        Folds = as.factor(resDist[[2]]),
        Iterate = iterate_
      )

      BetweenUnknown <- rbind(BetweenUnknown, resDist)
      matrixAbundances <- rbind(matrixAbundances,
                                cbind.data.frame(t(out_Dec)[,cellTypeName],
                                                 "Iterate" = iterate_))
      # Estimate unknown components using NMF
      resNMF <- nmf(x = abs(diff_bulk), rank = 1)
      unknownMat <- as.data.frame(basis(resNMF))
      colnames(unknownMat) <- paste0("Unknown_", iterate_)
      reference <- cbind(reference, unknownMat)
    }
  }
  results <- list("Prediction" = out_Dec,
                  "Matrix_prediction" = matrixAbundances,
                  "Error_folds" = BetweenUnknown)

  class(results) <- "SSDRnaSeq"
  return(results)
}
