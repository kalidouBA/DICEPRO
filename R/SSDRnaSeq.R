#' Semi-Supervised Deconvolution of Bulk RNA-Seq Data
#'
#' This function performs semi-supervised deconvolution of bulk RNA-seq data using the CIBERSORTx or
#' Digital Cell Quantification method. It iteratively estimates cell type proportions and refines
#' reference data using non-negative matrix factorization (NMF).
#'
#' @param reference The reference data matrix containing cell type gene expression profiles.
#' @param bulk The bulk RNA-seq data matrix.
#' @param k_folds The number of cross-validation folds.
#' @param nIteration The number of iterations to perform for refining the reference data.
#' @param methodDeconv A character vector specifying the deconvolution method.
#'        Supported values are "CSx" or "DCQ".
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @return abundances of cell types estimated.
#'
#' @details This function performs deconvolution of bulk RNA-Seq data using either
#' the CSx or DCQ method. It first prepares the data and runs the selected method
#' for cell type proportion estimation. The function then calculates error norms,
#' performs cross-validation, and estimates unknown components using NMF.
#'
#' The function performs several steps for semi-supervised deconvolution, including:
#' \enumerate{
#'   \item Running CIBERSORTx or Digital Cell Quantification for cell type proportion estimation.
#'   \item Calculating and storing performance metrics.
#'   \item Estimating unknown components using non-negative matrix factorization (NMF).
#'   \item Computing distances between unknown components and known cell types.
#'   \item Iteratively refining the reference data with unknown components.
#'}
#'
#' @export
#' @importFrom caret createFolds
#' @importFrom NMF nmf basis
#' @importFrom reshape2 melt
#' @importFrom ComICS dcq
#' @import dplyr

SSDRnaSeq <- function(reference, bulk, k_folds, nIteration, methodDeconv = c("DCQ", "CSx"),
                      cibersortx_email=NULL, cibersortx_token = NULL) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ"))

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  # Initialize variables to store results and metrics
  errorFrob <- list()
  BetweenUnknown <- NULL

  if(methodDeconv == "CSx"){
    # Prepare bulk data for CIBERSORTx
    bulkData <- cbind.data.frame(Gene=rownames(bulk), bulk)
    write.table(bulkData, file = paste0("bulk.txt"), sep = "\t", row.names = F, quote = F)
  }

  for (iterate_ in 1:nIteration) {
    if(methodDeconv == "CSx"){
      # Prepare reference data for CIBERSORTx
      sc_ref <- cbind.data.frame(GeneSymbol = rownames(reference), reference)
      write.table(sc_ref, file = "SignCSx.txt", sep = "\t", row.names = F, quote = F)

      # Run CIBERSORTx for cell type proportion estimation
      out_Dec <- run_CSx(cibersortx_email, cibersortx_token)
      }
    else{
      out_Dec = t(ComICS::dcq(reference_data = reference, mix_data = bulk,
                              marker_set = as.data.frame(row.names(reference)))$average)
      out_Dec = apply(out_Dec,2,function(x) ifelse(x < 0, 0, x))
      }

    # Calculate error norms
    bulkDeconv <- as.matrix(reference) %*% out_Dec
    diff_bulk <- bulk - bulkCSx
    errorFrob[[iterate_]] <- normFrob(diff_bulk, bulk)

    # Create folds for cross-validation
    flds <- createFolds(1:ncol(bulk), k = k_folds, list = TRUE, returnTrain = FALSE)

    # Initialize a matrix to store unknown components
    matUnknown <- matrix(ncol = k_folds, nrow = nrow(reference))

    for (indFld in 1:k_folds) {
      resNMF <- nmf(x = as.matrix(abs(diff_bulk)[, flds[[indFld]]], rank = 1))
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

    # Estimate unknown components using NMF
    resNMF <- nmf(x = abs(diff_bulk), rank = 1)
    unknownMat <- as.data.frame(basis(resNMF))
    colnames(unknownMat) <- paste0("Unknown_", iterate_)
    reference <- cbind(reference, unknownMat)
  }
  return(out_Dec)
}
