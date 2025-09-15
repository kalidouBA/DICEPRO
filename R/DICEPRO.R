#' DICEPRO: Deconvolution and Inference of Cell Proportions with Optimization
#'
#' This function performs cell type deconvolution of bulk gene expression data
#' using a combination of a deconvolution method and Non-Negative Matrix
#' Factorization (NMF) with L-BFGS-B optimization. It can integrate known
#' cell type signatures and estimate the proportions of unknown cell types.
#'
#' @param reference Numeric matrix of reference gene expression data (genes x cell types).
#' @param bulk Numeric matrix of bulk gene expression data (genes x samples).
#' @param methodDeconv Deconvolution method: "CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism".
#' @param cibersortx_email Email for CIBERSORTx (optional).
#' @param cibersortx_token Token for CIBERSORTx (optional).
#' @param W_prime Initial W matrix for unknown cell types (default 0).
#' @param bulkName Name of bulk dataset.
#' @param refName Name of reference dataset.
#' @param hp_max_evals Max evaluations for optimization (default 100).
#' @param N_unknownCT Number of unknown cell types to estimate (default 1).
#' @param algo_select Optimization algorithm: "random", "tpe", "atpe" (default "random").
#' @param output_path Path to save results (default current working directory).
#' @param hspaceTechniqueChoose Hyperparameter search space: "restrictionEspace", "all" (default "restrictionEspace").
#' @param out_Decon Precomputed deconvolution matrix (optional).
#'
#' @return List of class "DICEPRO" containing:
#' \item{Prediction}{Estimated cell type proportions (excluding unknowns).}
#' \item{New_signature}{Optimized gene expression signatures including unknown cell types.}
#' \item{W_prime}{Optimized W matrix for unknown cell types.}
#' \item{P_prime}{Estimated proportions for unknown cell types.}
#' \item{constraint}{Final constraint value from optimization.}
#' \item{optimization_results}{Complete optimization results including both optimal points and TSV data.}
#'
#' @export
DICEPRO <- function(reference, bulk,
                    methodDeconv = "DCQ",
                    cibersortx_email = NULL,
                    cibersortx_token = NULL,
                    W_prime = 0,
                    bulkName = "Bulk",
                    refName = "Reference",
                    hp_max_evals = 100,
                    N_unknownCT = 1,
                    algo_select = "random",
                    output_path = NULL,
                    hspaceTechniqueChoose = "restrictionEspace",
                    out_Decon = NULL) {

  # --- Validate deconvolution method ---
  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism"))
  methodDeconv <- methodDeconv[1]
  if (methodDeconv == "CSx" && (is.null(cibersortx_token) || is.null(cibersortx_email))) {
    stop("CIBERSORTx token is required for method 'CSx'. Please provide 'cibersortx_token'.")
  }

  # --- Match genes between reference and bulk ---
  geneIntersect <- intersect(rownames(reference), rownames(bulk))
  if (is.matrix(reference)) {
    stopifnot(is.numeric(reference))
    reference <- reference[geneIntersect, ]
  } else {
    reference <- apply(reference[geneIntersect, ], 2, as.numeric)
  }

  if (is.matrix(bulk)) {
    stopifnot(is.numeric(bulk))
    bulk <- bulk[geneIntersect, ]
  } else {
    bulk <- apply(bulk[geneIntersect, ], 2, as.numeric)
  }

  rownames(reference) <- rownames(bulk) <- geneIntersect

  # --- Run deconvolution if out_Decon not provided ---
  if (is.null(out_Decon)) {
    out_Dec <- t(running_method(bulk, reference, methodDeconv, cibersortx_email, cibersortx_token))
  } else {
    message("Using provided 'out_Decon', skipping running_method().")
    out_Dec <- out_Decon
    if (!is.matrix(out_Dec)) stop("The 'out_Decon' object must be a matrix.")
  }

  # Build dataset list
  dataset <- list(B = bulk, W = reference, P = out_Dec)

  # Determine output directory
  dirName <- paste0(bulkName, "_", refName, "_", hspaceTechniqueChoose, "_", algo_select)
  if (is.null(output_path)) output_path <- getwd()
  output_dir <- file.path(output_path, dirName)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Run hyperparameter optimization (R-native)
  optimization_results <- run_experiment(
    dataset = dataset,
    W_prime = W_prime,
    bulkName = bulkName,
    refName = refName,
    hp_max_evals = hp_max_evals,
    algo_select = algo_select,
    output_base_dir = output_path,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  )

  # Extract both optimal points
  optimal_DiceproOptCstrt <- optimization_results$optimal_points$DiceproOptCstrt
  optimal_DiceproOptCstrt_0.1 <- optimization_results$optimal_points$DiceproOptCstrt_0.1

  # Extract both TSV datasets
  tsv_DiceproOptCstrt <- optimization_results$optimal_tsv_data$DiceproOptCstrt
  tsv_DiceproOptCstrt_0.1 <- optimization_results$optimal_tsv_data$DiceproOptCstrt_0.1

  # Extract final results from the primary optimal point (DiceproOptCstrt)
  final_results_DiceproOptCstrt <- if (!is.null(tsv_DiceproOptCstrt)) {
    extract_final_results(tsv_DiceproOptCstrt)
  } else {
    list()
  }

  # Extract final results from the secondary optimal point (DiceproOptCstrt_0.1)
  final_results_DiceproOptCstrt_0.1 <- if (!is.null(tsv_DiceproOptCstrt_0.1)) {
    extract_final_results(tsv_DiceproOptCstrt_0.1)
  } else {
    list()
  }

  # Build the comprehensive results list
  final_results <- list(
    # Primary results from DiceproOptCstrt
    DiceproOptCstrt = list(Prediction = out_Dec, New_signature = final_results_DiceproOptCstrt$W_final %||% reference,
                           W_prime = final_results_DiceproOptCstrt$W_prime_final %||% W_prime, P_prime = final_results_DiceproOptCstrt$P_prime_final,
                           constraint = optimal_DiceproOptCstrt$constraint,
                           lambda_ = optimal_DiceproOptCstrt$lambda_,
                           gamma = optimal_DiceproOptCstrt$gamma,
                           p_prime = optimal_DiceproOptCstrt$p_prime,
                           frobNorm = optimal_DiceproOptCstrt$frobNorm),

    # Secondary results from DiceproOptCstrt_0.1
    DiceproOptCstrt_0.1 = list(secondary_Prediction = out_Dec,
                              secondary_New_signature = final_results_DiceproOptCstrt_0.1$W_final %||% reference,
                              secondary_W_prime = final_results_DiceproOptCstrt_0.1$W_prime_final %||% W_prime,
                              secondary_P_prime = final_results_DiceproOptCstrt_0.1$P_prime_final,
                              secondary_constraint = optimal_DiceproOptCstrt_0.1$constraint,
                              secondary_lambda_ = optimal_DiceproOptCstrt_0.1$lambda_,
                              secondary_gamma = optimal_DiceproOptCstrt_0.1$gamma,
                              secondary_p_prime = optimal_DiceproOptCstrt_0.1$p_prime,
                              secondary_frobNorm = optimal_DiceproOptCstrt_0.1$frobNorm),

    # Complete optimization results
    optimization_results = optimization_results
  )

  class(final_results) <- "DICEPRO"
  return(final_results)
}
