#' Semi-Supervised Deconvolution of Bulk RNA-Seq Data
#'
#' This function performs semi-supervised deconvolution of bulk RNA-seq data using the CIBERSORTx method.
# It iteratively estimates cell type proportions and refines reference data using non-negative matrix factorization (NMF).
#
#' @param reference The reference data matrix containing cell type gene expression profiles.
#' @param bulk The bulk RNA-seq data matrix.
#' @param k_folds The number of cross-validation folds.
#' @param nIteration The number of iterations to perform for refining the reference data.
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @return abundances of cell types estimated.
#' @details The function performs several steps for semi-supervised deconvolution, including:
#' \enumerate{
#'   \item Running CIBERSORTx for cell type proportion estimation.
#'   \item Calculating and storing performance metrics.
#'   \item Estimating unknown components using non-negative matrix factorization (NMF).
#'   \item Computing distances between unknown components and known cell types.
#'   \item Iteratively refining the reference data with unknown components.
#'}
#'
#' @examples
#' reference_data <- matrix(rnorm(100), ncol = 5)
#' bulk_data <- matrix(rnorm(100), ncol = 5)
#' cibersortx_email <- "your_email@example.com"
#' cibersortx_token <- "your_token"
#' results <- SSDRnaSeq(reference_data, bulk_data, 5, 50, cibersortx_email, cibersortx_token)
#'
#' @import base
#' @import utils
#' @importFrom readr read.delim2
#' @importFrom caret createFolds
#' @importFrom NMF nmf
#' @export

SSDRnaSeq <- function(reference, bulk, k_folds, nIteration, cibersortx_email, cibersortx_token) {
  # Prepare bulk data for CIBERSORTx
  bulkData <- cbind.data.frame(Gene=rownames(bulk), bulk)
  write.table(bulkData, file = paste0("bulk.txt"), sep = "\t", row.names = F, quote = F)

  # Initialize variables to store results and metrics
  typeOFabundances <- list()
  matAbundances <- list()
  allError <- list()
  errorFrob <- list()
  BetweenUnknown <- list()

  for (iterate_ in 1:nIteration) {
    # Prepare reference data for CIBERSORTx
    sc_ref <- cbind.data.frame(GeneSymbol = rownames(reference), reference)
    write.table(sc_ref, file = "SignCSx.txt", sep = "\t", row.names = F, quote = F)

    # Run CIBERSORTx for cell type proportion estimation
    out_Dec <- run_CSx(cibersortx_email, cibersortx_token)

    # Store iteration-specific results
    typeOFabundances[[iterate_]] <- rep(paste0("it_", iterate_), nrow(truthAbUse))
    matAbundances[[iterate_]] <- t(out_Dec)[, colnames(truthAbUse)]

    # Calculate and store performance metrics
    resPerf <- computPerf(truth = truthAbUse, estimated = t(out_Dec), it_ = iterate_)
    allError[[iterate_]] <- resPerf

    # Calculate error norms
    bulkCSx <- as.matrix(reference) %*% out_Dec
    diff_bulk <- bulk - bulkCSx
    errorFrob[[iterate_]] <- normFrob(diff_bulk, bulk)

    # Create folds for cross-validation
    flds <- caret::createFolds(1:ncol(bulk), k = k_folds, list = TRUE, returnTrain = FALSE)

    # Initialize a matrix to store unknown components
    matUnknown <- matrix(ncol = k_folds, nrow = nrow(reference))

    for (indFld in 1:k_folds) {
      resNMF <- nmf(x = as.matrix(abs(diff_bulk)[, flds[[indFld]]], rank = 1))
      matUnknown[, indFld] <- as.vector(basis(resNMF))
    }

    colnames(matUnknown) <- names(flds)

    # Compute distances between unknown components and known cell types
    resDist <- compute_distances(matrix_ = matUnknown)

    # Store distance results
    resDist <- reshape2::melt(resDist[[1]]) %>% cbind.data.frame(
      Folds = as.factor(resDist[[2]]),
      Iterate = iterate_
    )

    BetweenUnknown <- rbind(BetweenUnknown, resDist)

    # Estimate unknown components using NMF
    resNMF <- nmf(x = abs(diff_bulk), rank = 1)
    unknownMat <- as.data.frame(basis(resNMF))
    colnames(unknownMat) <- paste0("Unknown_", iterate_)
    reference <- cbind(reference, unknownMat)
  }

  # Return the final CIBERSORTx results
  return(out_Dec)
}
