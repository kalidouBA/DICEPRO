#' Run CIBERSORTx Deconvolution Method
#'
#' This function runs the CIBERSORTx deconvolution method using Docker with the provided email and token.
#' The method is used to estimate cell type proportions in bulk RNA-seq data.
#'
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @export
#'
#' @import base
#' @import utils
#'
#' @return A matrix containing the estimated cell type proportions.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Runs CIBERSORTx with the provided email and token.
#'   \item Reads the output results from CIBERSORTx.
#'   \item Extracts cell type proportions and returns them as a matrix.
#'}
#'
#' @examples
#'if(interactive()){
#' cibersortx_email <- "your_email@example.com"
#' cibersortx_token <- "your_token"
#' estimated_proportions <- run_methods_deconv(cibersortx_email, cibersortx_token)
#' }


run_CSx <- function(cibersortx_email, cibersortx_token){
  verbose <- FALSE
  in_dir <- out_dir <- getwd()
  command_to_run <- paste0("docker run -v ", in_dir, ":/src/data -v ", out_dir,
                           ":/src/outdir cibersortx/fractions --username ", cibersortx_email,
                           " --token ",  cibersortx_token," --mixture bulk.txt --sigmatrix SignCSx.txt --QN FALSE")
  code <- system(command_to_run)
  if(code == 0){
    out_CSx <- read.delim2(paste0(out_dir, "/CIBERSORTx_Results.txt"),row.names = 1)
    rownames_ <- rownames(out_CSx)
    out_CSx <- apply(out_CSx[,!names(out_CSx) %in% c("P.value", "Correlation", "RMSE")],2, as.numeric)
    rownames(out_CSx) <- rownames_
  }
  return(t(out_CSx))
}
