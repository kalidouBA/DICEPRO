#' Run CIBERSORTx Deconvolution Method
#'
#' This function runs the CIBERSORTx deconvolution method using Docker with the provided email and token.
#' The method is used to estimate cell type proportions in bulk RNA-seq data.
#'
#' @param bulk A matrix of bulk RNA-seq data.
#' @param reference A matrix of reference cell type gene expression profiles.
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#'
#' @importFrom utils read.delim2 write.table
#'
#' @return A matrix containing the estimated cell type proportions.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Runs CIBERSORTx with the provided email and token.
#'   \item Reads the output results from CIBERSORTx.
#'   \item Extracts cell type proportions and returns them as a matrix.
#'}
#' @export

run_CSx <- function(bulk, reference, cibersortx_email, cibersortx_token){
  verbose <- FALSE
  in_dir <- out_dir <- getwd()
  sc_ref <- cbind.data.frame(GeneSymbol = rownames(reference), reference)
  write.table(sc_ref, file = "SignCSx.txt", sep = "\t", row.names = F, quote = F)

  bulk_ <- cbind.data.frame(Gene=rownames(bulk),bulk)
  write.table(bulk_, file = "bulk.txt", sep = "\t", row.names = F, quote = F)

  command_to_run <- paste0("docker run -v ", in_dir, ":/src/data -v ", out_dir,
                           ":/src/outdir cibersortx/fractions --username ", cibersortx_email,
                           " --token ",  cibersortx_token," --mixture bulk.txt --sigmatrix SignCSx.txt --QN FALSE")
  code <- system(command_to_run)
  if(code == 0){
    out_CSx <- read.delim2(paste0(out_dir, "/CIBERSORTx_Results.txt"),row.names = 1)
    rownames_ <- rownames(out_CSx)
    out_CSx <- as.data.frame(apply(out_CSx[,!names(out_CSx) %in% c("P.value", "Correlation", "RMSE")],2, as.numeric))
    rownames(out_CSx) <- rownames_
  }
  return(t(out_CSx))
}
