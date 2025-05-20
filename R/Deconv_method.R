#' Run Deconvolution Method
#'
#' This function performs cell deconvolution using various methods, including CSx, DCQ, DeconRNASeq, FARDEEP, etc.
#'
#' @param bulk A matrix of bulk RNA-seq data where rows are genes and columns are samples.
#' @param reference A matrix of reference cell type gene expression profiles where rows are genes and columns are cell types.
#' @param methodDeconv The method to use for deconvolution. Options include "CSx", "DCQ", "DeconRNASeq", "FARDEEP", and "BayesPrism".
#' @param cibersortx_email The CIBERSORTx account email (used if `methodDeconv` is "CSx").
#' @param cibersortx_token The CIBERSORTx account token (used if `methodDeconv` is "CSx").
#'
#' @return A matrix containing the estimated cell proportions for each cell type.
#'
#' @details This function allows you to run cell deconvolution using different methods. The choice of method is specified by the `methodDeconv` parameter.
#'          Each method runs its corresponding deconvolution procedure and returns the estimated proportions of the cell types in the bulk data.
#'
#' @seealso Documentation of the methods used for cell deconvolution: \code{\link{run_CSx}}, \code{\link{dcq}}, \code{\link{DeconRNASeq}},
#' \code{\link{fardeep}}.
#'
#' @importFrom ComICS dcq
#' @importFrom DeconRNASeq DeconRNASeq
#' @importFrom FARDEEP fardeep
#' @importFrom parallel detectCores
#'
#' @export
running_method <- function(bulk, reference, methodDeconv = "CSx", cibersortx_email, cibersortx_token){
  dimname_ <- dimnames(reference)
  reference <- apply(reference, 2, as.numeric)
  dimnames(reference) <- dimname_

  dimname_ <- dimnames(bulk)
  bulk <- apply(bulk, 2, as.numeric)
  dimnames(bulk) <- dimname_

  nCellType <- ncol(reference)
  common <- intersect(row.names(bulk), row.names(reference))

  bulkIntersect <- as.data.frame(bulk[common,])
  refIntersect <- reference[common,]
  geneLenght <- nrow(bulk)
  row.names(bulkIntersect) <- row.names(refIntersect) <- common
  colnames(bulkIntersect) <- colnames(bulk)

  if (methodDeconv == "CSx")
    out_Dec <- run_CSx(bulk, reference, cibersortx_email, cibersortx_token)
  else if(methodDeconv == "DCQ"){
    out_Dec = t(dcq(reference_data = refIntersect, mix_data = bulkIntersect,
                    marker_set = as.data.frame(row.names(refIntersect)),
                    alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 30)$average)
    out_Dec = apply(out_Dec, 2, function(x) ifelse(x < 0, 0, x))
  }
  else if(methodDeconv == "DeconRNASeq"){
    out_Dec = DeconRNASeq(datasets = as.data.frame(bulk), signatures = as.data.frame(reference),
                          proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE,
                          fig = FALSE)$out.all
    rownames(out_Dec) = colnames(bulk)
    out_Dec <- t(out_Dec)
  }
  else if(methodDeconv == "FARDEEP")
    out_Dec <- t(fardeep(X = reference, Y = bulk, nn = TRUE, intercept = TRUE, permn = 100, QN = FALSE)$abs.beta)

  return(out_Dec)
}
