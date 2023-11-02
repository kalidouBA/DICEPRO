#' Run Deconvolution Method
#'
#' This function performs cell deconvolution using various methods.
#'
#' @param bulk A matrix of bulk RNA-seq data.
#' @param reference A matrix of reference cell type gene expression profiles.
#' @param methodDeconv The method to use for deconvolution. Options include "CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", and "BayesPrism".
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @return A matrix containing the estimated cell proportions for each cell type.
#'
#' @details This function allows you to run cell deconvolution using different methods, including CSx, DCQ, CDSeq, DeconRNASeq, FARDEEP, and BayesPrism.
#' The choice of method should be specified using the `methodDeconv` parameter.
#'
#' @seealso Documentation of the methods used for cell deconvolution: \code{\link{run_CSx}}, \code{\link{dcq}}, \code{\link{CDSeq}}, \code{\link{DeconRNASeq}},
#' \code{\link{fardeep}}, \code{\link{new.prism}}, \code{\link{run.prism}}, \code{\link{get.fraction}}.
#'
#' @importFrom ComICS dcq
#' @importFrom CDSeq CDSeq
#' @importFrom DeconRNASeq DeconRNASeq
#' @importFrom FARDEEP fardeep
#' @importFrom parallel detectCores
#' @importFrom BayesPrism new.prism run.prism get.fraction

running_method <- function(bulk, reference, methodDeconv = "CSx", cibersortx_email, cibersortx_token){
  nCellType <- ncol(reference)
  common <- intersect(row.names(bulk), row.names(reference))

  bulkIntersect <- bulk[common,]
  refIntersect <- reference[common,]
  geneLenght <- nrow(bulk)
  row.names(bulkIntersect) <- row.names(refIntersect) <- common

  switch(methodDeconv,
         CSx={
             out_Dec <- run_CSx(bulk, reference, cibersortx_email, cibersortx_token)
             break
           },
         DCQ={
             out_Dec = t(dcq(reference_data = refIntersect, mix_data = bulkIntersect,
                                     marker_set = as.data.frame(row.names(refIntersect)),
                             alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 30)$average)
             out_Dec = apply(out_Dec, 2, function(x) ifelse(x < 0, 0, x))
             break
           },
         DeconRNASeq = {
           out_Dec = DeconRNASeq(datasets = bulk, signatures = reference, proportions = NULL, checksig = FALSE,
                                 known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all
           rownames(out_Dec) = colnames(bulk)
           break
           },
         FARDEEP = {
           out_Dec <- fardeep(X = reference, Y = bulk, nn = TRUE, intercept = TRUE, permn = 100, QN = FALSE)$abs.beta
           break
           }
         )
  return(out_Dec)
}
