#' Semi-Supervised Deconvolution of Bulk RNA-Seq Data
#'
#' This function performs semi-supervised deconvolution of bulk RNA-seq data using the CIBERSORTx or
#' Digital Cell Quantification method or others. It iteratively estimates cell type proportions and refines
#' reference data using non-negative matrix factorization (NMF).
#'
#' @param reference The reference data matrix containing cell type (column) gene expression (row) profiles.
#' @param bulk The bulk RNA-seq data matrix. In the columns we have the samples and in the rows the genes.
#' @param nIteration The number of iterations to perform for refining the reference data.
#' @param methodDeconv A character vector specifying the deconvolution method.
#'        Supported values are "CSx" or "DCQ" or others. The method to use for deconvolution. Options include \code{CSx},
#'        \code{DCQ}, \code{CDSeq}, \code{DeconRNASeq}, \code{FARDEEP} and \code{BayesPrism}.
#' @param metric A metric used to detect the convergence method.
#'        Supported values are "R2_adj", "ICC", "RRMSE", "simRatio".
#'        Choose the metric to evaluate the convergence method. Options include \code{R2_adj},
#'        which represents Adjusted R-squared and \code{RRMSE} for Relative Root Mean Squared Error.
#'        The convergence method will be determined based on the selected metric.
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
#'   \item Calculating and storing performance metrics.
#'   \item Estimating unknown components using non-negative matrix factorization (NMF).
#'   \item Computing distances between unknown components and known cell types.
#'   \item Iteratively refining the reference data with unknown components.
#'}
#'
#' @export
#'
#' @importFrom NMF nmf basis
#' @importFrom ComICS dcq
#' @import dplyr
#'
#' @examples
#' if(interactive()){
#'
#' simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE,
#' nSample = 10, prop = NULL, nGenes = 50, nCellsType = 5)
#' cellTypeOut <- sample(1:ncol(simulation$reference), 2)
#' refDataIncomplet <- simulation$reference[,-cellTypeOut]
#' results <- SSDRnaSeq(reference = refDataIncomplet, bulk = simulation$bulk,
#' nIteration = 5, methodDeconv = "DCQ")
#' print(results)
#' }

SSDRnaSeq <- function(reference, bulk, nIteration = 50, methodDeconv = "CSx", metric = "RRMSE",
                      cibersortx_email = NULL, cibersortx_token = NULL) {

  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))
  stopifnot(metric %in% c("RRMSE", "R2_adj"))
  stopifnot(nIteration > 0)

  if(length(methodDeconv) > 1)
    methodDeconv <- methodDeconv[1]

  cellTypeName <- colnames(reference)
  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  reference <- apply(reference[geneIntersect, ], 2, as.numeric)
  bulk <- apply(bulk[geneIntersect, ], 2, as.numeric)
  rownames(reference) <- rownames(bulk) <- geneIntersect

  # Initialize variables to store results and metrics
  matrixAbundances <- performs <- performs2plot <- opt <- NULL

  for (iterate_ in 0:nIteration) {
    message("Current iteration ++++++++++++++++++++++++++++++++ ", iterate_)

    out_Dec <- running_method(bulk, reference, methodDeconv,  cibersortx_email, cibersortx_token)

    bulkDeconv <- as.matrix(reference) %*% out_Dec
    diff_bulk <- bulk - bulkDeconv

    matrixAbundances <- rbind(matrixAbundances, cbind.data.frame(t(out_Dec)[,cellTypeName], "Iterate" = iterate_))

    if(iterate_ > 0){
      perform_it <- computPerf(outDec_1 = matrixAbundances[matrixAbundances$Iterate == iterate_-1, cellTypeName],
                               outDec_2 = matrixAbundances[matrixAbundances$Iterate == iterate_, cellTypeName], metric)

      performs <- c(performs, perform_it)
      performs2plot <- rbind.data.frame(performs2plot, data.frame(metric = perform_it, Iterate = iterate_))

      if (length(performs) > 1 &&
          ((metric == "R2_adj" && performs[iterate_] > 0.99) ||
           (metric == "RRMSE" && performs[iterate_-1] - performs[iterate_] < 0 ) ||
           iterate_ == nIteration)) {
        opt <- ifelse(iterate_ == nIteration, iterate_, iterate_ - 1)
        message("Convergence criteria Done with optimal criteria: ", opt, "\nBreaking the loop.")

        break
      }
    }

    # Estimate one unknown component using NMF
    resNMF <- nmf(x = abs(diff_bulk), rank = 1)
    unknownMat <- as.data.frame(basis(resNMF))
    colnames(unknownMat) <- paste0("Unknown_", iterate_)
    reference <- cbind(reference, unknownMat)
  }

  rownames(performs2plot) <- NULL
  results <- list("Prediction" = out_Dec, "Matrix_prediction" = matrixAbundances,
                  "New_signature" = reference, "Optimal_iteration" = opt,
                  "performs2plot" = performs2plot)

  class(results) <- "SSDRnaSeq"
  return(results)
}
