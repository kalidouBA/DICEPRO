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
  stopifnot(methodDeconv %in% c("CSx", "DCQ", "CDSeq", "DeconRNASeq", "FARDEEP", "BayesPrism", "abis"))
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

  # Use the private helper function to extract results with safe null handling
  final_results <- list(
    DiceproOptCstrt = if (!is.null(optimal_DiceproOptCstrt)) {
      .get_optimal_results(optimal_DiceproOptCstrt, tsv_DiceproOptCstrt, out_Dec, reference, W_prime)
    } else {
      list()
    },
    DiceproOptCstrt_0.1 = if (!is.null(optimal_DiceproOptCstrt_0.1)) {
      .get_optimal_results(optimal_DiceproOptCstrt_0.1, tsv_DiceproOptCstrt_0.1, out_Dec, reference, W_prime, prefix = "secondary_")
    } else {
      list()
    },
    optimization_results = optimization_results
  )


  class(final_results) <- "DICEPRO"
  return(final_results)
}


#' @title Extract final results from a TSV dataset
#' @description Private helper function to extract relevant final results from a TSV dataset.
#'   This function assumes the TSV contains columns named `W_final`, `W_prime_final`, and `P_prime_final`.
#'
#' @param tsv_data A data frame read from the TSV file corresponding to an optimal point.
#'
#' @return A list with the following elements:
#'   - `W_final`: final W matrix (or NULL if not present)
#'   - `W_prime_final`: final W_prime values (or NULL if not present)
#'   - `P_prime_final`: final P_prime values (or NULL if not present)
#'
#' @keywords internal
.extract_final_results <- function(tsv_data) {
  list(
    W_final = tsv_data$W_final %||% NULL,
    W_prime_final = tsv_data$W_prime_final %||% NULL,
    P_prime_final = tsv_data$P_prime_final %||% NULL
  )
}


# NULL coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' @title Extract Results for a Single Optimal Point (Internal)
#' @description
#' Private helper function to extract final results for a given optimal point.
#' It safely handles missing TSV data and optionally prefixes the names for secondary results.
#'
#' @param opt_point A list or data frame row representing the optimal point.
#'   Must contain `constraint`, `lambda_`, `gamma`, `p_prime`, and `frobNorm`.
#' @param tsv_data A data frame read from the corresponding TSV file, or `NULL`.
#' @param out_Dec The deconvolution result matrix.
#' @param reference The reference matrix.
#' @param W_prime Initial W matrix for unknown cell types.
#' @param prefix Optional string to prefix the result names (default `""`).
#'
#' @return A named list containing:
#'   - `Prediction`: Default prediction object (`out_Dec`).
#'   - `New_signature`: New signature matrix (from TSV or fallback `reference`).
#'   - `W_prime`: Optimized W_prime value.
#'   - `P_prime`: Optimized P_prime value.
#'   - `constraint`, `lambda_`, `gamma`, `p_prime`, `frobNorm`: Values from the optimal point.
#'
#' @keywords internal
#' @noRd
.get_optimal_results <- function(opt_point, tsv_data, out_Dec, reference, W_prime, prefix = "") {
  final_res <- if (!is.null(tsv_data)) .extract_final_results(tsv_data) else list()

  res <- list(
    Prediction = out_Dec,
    New_signature = final_res$W_final %||% reference,
    W_prime = final_res$W_prime_final %||% W_prime,
    P_prime = final_res$P_prime_final,
    constraint = opt_point$constraint,
    lambda_ = opt_point$lambda_,
    gamma = opt_point$gamma,
    p_prime = opt_point$p_prime,
    frobNorm = opt_point$frobNorm
  )

  if (prefix != "") names(res) <- paste0(prefix, names(res))
  res
}
