# =============================================================================
# DICEPRO Main Function
#
# Public API  : DICEPRO()
# Private fns : .normalize_zscore_per_gene()
# =============================================================================


# -----------------------------------------------------------------------------
# .normalize_zscore_per_gene  [private]
# -----------------------------------------------------------------------------

#' Z-score normalisation per gene
#'
#' Normalises each row (gene) of a numeric matrix by centering on its mean
#' and scaling by its standard deviation. Genes with SD = 0 or \code{NA}
#' are silently removed.
#'
#' @param mat A numeric matrix or data frame with genes as rows and samples
#'   as columns.
#'
#' @return A numeric matrix with the same columns as \code{mat} but
#'   potentially fewer rows (genes with SD = 0 removed). Each remaining
#'   row has mean = 0 and SD = 1.
#'
#' @importFrom stats sd
#' @keywords internal
#' @noRd
.normalize_zscore_per_gene <- function(mat) {

  mat <- as.matrix(mat)

  gene_sd    <- apply(mat, 1L, sd, na.rm = TRUE)
  keep_genes <- !is.na(gene_sd) & gene_sd > 0
  mat        <- mat[keep_genes, , drop = FALSE]

  t(apply(mat, 1L, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
}


# -----------------------------------------------------------------------------
# DICEPRO  [public]
# -----------------------------------------------------------------------------

#' Semi-supervised bulk RNA-seq deconvolution with NMF hyperparameter
#' optimisation
#'
#' \code{DICEPRO} performs cell-type deconvolution of bulk RNA-seq data by
#' combining a supervised deconvolution step (estimating proportions of
#' \emph{known} cell types) with an unsupervised NMF step (discovering and
#' quantifying \emph{unknown} cell types). Hyperparameters
#' \eqn{(\lambda, \gamma, p')} are selected automatically via Pareto-frontier
#' analysis and knee-point detection.
#'
#' @section Pipeline:
#' \enumerate{
#'   \item \strong{Normalisation} (optional): z-score per gene applied to
#'     both \code{reference} and \code{bulk} via
#'     \code{.normalize_zscore_per_gene()}.
#'   \item \strong{Gene intersection}: only genes present in both matrices
#'     are retained.
#'   \item \strong{Supervised deconvolution}: \code{\link{running_method}}
#'     estimates known cell-type proportions.
#'   \item \strong{Hyperparameter search}: \code{\link{run_experiment}}
#'     evaluates \code{hp_max_evals} configurations using the chosen
#'     \code{algo_select} strategy.
#'   \item \strong{Best configuration}: \code{\link{best_hyperParams}}
#'     selects the knee point on the Pareto frontier.
#'   \item \strong{Report}: plots are saved under
#'     \code{output_path/DICEPRO_<bulkName>_<refName>/report/}.
#' }
#'
#' @param reference    Numeric matrix of reference gene expression profiles,
#'   dimensions \eqn{G \times K} (genes \eqn{\times} known cell types).
#'   Row names must be gene symbols.
#' @param bulk         Numeric matrix of bulk RNA-seq expression,
#'   dimensions \eqn{G \times N} (genes \eqn{\times} samples).
#'   Row names must be gene symbols matching \code{reference}.
#' @param methodDeconv Character scalar. Deconvolution method to use for
#'   the supervised step. One of \code{"CSx"}, \code{"DCQ"},
#'   \code{"CDSeq"}, or \code{"FARDEEP"} (default \code{"CSx"}).
#' @param cibersortx_email Character scalar. Registered email for the
#'   CIBERSORTx web service. Required when \code{methodDeconv = "CSx"};
#'   ignored otherwise.
#' @param cibersortx_token Character scalar. API token for the CIBERSORTx
#'   web service. Required when \code{methodDeconv = "CSx"};
#'   ignored otherwise.
#' @param W_prime      Either \code{0} (default, initialise a single unknown
#'   cell-type column automatically) or a numeric matrix of dimensions
#'   \eqn{G \times U} providing initial signatures for \eqn{U} unknown cell
#'   types.
#' @param bulkName     Character scalar. Label for the bulk dataset, used in
#'   the output directory name and plot titles (default \code{"Bulk"}).
#' @param refName      Character scalar. Label for the reference dataset,
#'   used in the output directory name and plot titles
#'   (default \code{"Reference"}).
#' @param hp_max_evals Positive integer. Number of hyperparameter
#'   configurations to evaluate during the search
#'   (default \code{100}). Increase to \eqn{\geq 200} for production use.
#' @param N_unknownCT  Positive integer. Number of unknown cell types to
#'   estimate when \code{W_prime = 0} (default \code{1}).
#' @param algo_select  Character scalar. Hyperparameter sampling strategy.
#'   \describe{
#'     \item{\code{"random"}}{Pure random search (default).}
#'     \item{\code{"tpe"}}{Tree-structured Parzen Estimator.}
#'     \item{\code{"atpe"}}{Accepted alias for \code{"tpe"}.}
#'     \item{\code{"anneal"}}{Accepted in config; currently falls back to
#'       \code{"random"}.}
#'   }
#' @param output_path  Character scalar. Root directory for all outputs.
#'   Defaults to \code{getwd()} when \code{NULL}.
#' @param hspaceTechniqueChoose Character scalar. Hyperparameter space
#'   parameterisation (\code{"all"} or \code{"restrictionEspace"}).
#' @param out_Decon    Optional numeric matrix. Pre-computed deconvolution
#'   result. When provided, \code{\link{running_method}} is skipped.
#' @param normalize    Logical scalar. When \code{TRUE} (default), both
#'   matrices are z-score normalised per gene before deconvolution.
#'
#' @return An object of class \code{"DICEPRO"} (a named list), or
#'   \code{invisible(NULL)} when no valid hyperparameter configuration is
#'   found. The list contains:
#' \describe{
#'   \item{hyperparameters}{Named list with \code{lambda} and \code{gamma}.}
#'   \item{metrics}{Named list with \code{loss} and \code{constraint}.}
#'   \item{trials}{data.frame of all successful trial results.}
#'   \item{W}{Signature matrix of the selected trial.}
#'   \item{H}{Proportion matrix of the selected trial.}
#'   \item{plot}{ggplot2 figure of the Pareto frontier.}
#'   \item{plot_hyperopt}{ggplot2 scatter-matrix of all hyperparameter
#'     configurations.}
#' }
#'
#' @seealso
#'   \code{\link{running_method}}, \code{\link{run_experiment}},
#'   \code{\link{best_hyperParams}}, \code{\link{plot_hyperopt}}.
#'
#' @export
DICEPRO <- function(reference, bulk,
                    methodDeconv          = "CSx",
                    cibersortx_email      = NULL,
                    cibersortx_token      = NULL,
                    W_prime               = 0,
                    bulkName              = "Bulk",
                    refName               = "Reference",
                    hp_max_evals          = 100,
                    N_unknownCT           = 1,
                    algo_select           = "random",
                    output_path           = NULL,
                    hspaceTechniqueChoose = "all",
                    out_Decon             = NULL,
                    normalize             = TRUE) {

  # ---- Input validation ------------------------------------------------------
  if (!is.logical(normalize) || length(normalize) != 1L)
    stop("'normalize' must be a single logical value (TRUE or FALSE).")

  algo_select <- match.arg(
    tolower(algo_select),
    c("random", "tpe", "atpe", "anneal")
  )

  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("restrictionEspace", "all")
  )

  # ---- Normalisation ---------------------------------------------------------
  if (normalize) {
    reference <- .normalize_zscore_per_gene(reference)
    bulk      <- .normalize_zscore_per_gene(bulk)
  } else {
    reference <- as.matrix(reference)
    bulk      <- as.matrix(bulk)
    message("normalize = FALSE: matrices used as-is.")
  }

  # ---- Gene intersection -----------------------------------------------------
  geneIntersect <- intersect(rownames(reference), rownames(bulk))
  if (length(geneIntersect) == 0L)
    stop("No common genes between 'reference' and 'bulk'.")

  reference <- reference[geneIntersect, , drop = FALSE]
  bulk      <- bulk[geneIntersect, , drop = FALSE]
  rownames(reference) <- rownames(bulk) <- geneIntersect
  message(sprintf("Gene intersection: %d genes retained.", length(geneIntersect)))

  # ---- Supervised deconvolution ----------------------------------------------
  if (is.null(out_Decon)) {
    methodDeconv <- match.arg(methodDeconv, c("CSx", "DCQ", "CDSeq", "FARDEEP"))

    if (methodDeconv == "CSx" &&
        (is.null(cibersortx_token) || is.null(cibersortx_email))) {
      stop("CIBERSORTx credentials are required for methodDeconv = 'CSx'. ",
           "Please provide both 'cibersortx_email' and 'cibersortx_token'.")
    }

    out_Dec <- t(running_method(
      bulk, reference, methodDeconv,
      cibersortx_email, cibersortx_token
    ))

  } else {

    out_Decon <- as.matrix(out_Decon)

    if (ncol(out_Decon) != ncol(reference))
      stop(sprintf(
        "'out_Decon' has %d cell type(s) but 'reference' has %d.",
        ncol(out_Decon), ncol(reference)
      ))

    if (nrow(out_Decon) != ncol(bulk))
      stop(sprintf(
        "'out_Decon' has %d sample(s) but 'bulk' has %d columns.",
        nrow(out_Decon), ncol(bulk)
      ))

    message(sprintf(
      "Using provided 'out_Decon' -- skipping running_method(). Dimensions: %d samples x %d cell types.",
      nrow(out_Decon), ncol(out_Decon)
    ))

    out_Dec <- t(out_Decon)
  }

  # ---- Build dataset list ----------------------------------------------------
  dataset <- list(B = bulk, W = reference, P = out_Dec)

  # ---- Output directory ------------------------------------------------------
  dirName    <- paste0("DICEPRO_", bulkName, "_", refName)
  if (is.null(output_path)) output_path <- getwd()
  output_dir <- file.path(output_path, dirName)

  # ---- Hyperparameter column names for the report plot ----------------------
  hp_params <- switch(
    hspaceTechniqueChoose,
    all               = c("gamma", "lambda_",      "p_prime"),
    restrictionEspace = c("gamma", "lambda_factor", "p_prime")
  )

  # ---- Hyperparameter optimisation -------------------------------------------
  out <- run_experiment(
    dataset               = dataset,
    W_prime               = W_prime,
    bulkName              = bulkName,
    refName               = refName,
    hp_max_evals          = hp_max_evals,
    algo_select           = algo_select,
    output_base_dir       = output_path,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  ) |>
    (\(res) best_hyperParams(
      trials_df = res$trials,
      W         = res$W,
      H         = res$H,
      savePaths = output_dir
    ))()

  # ---- Guard: all trials failed ----------------------------------------------
  if (is.null(out)) {
    warning(
      "DICEPRO: no valid hyperparameter configuration found. ",
      "Try increasing 'hp_max_evals' or checking dataset dimensions."
    )
    return(invisible(NULL))
  }

  # ---- Hyperparameter optimisation plot --------------------------------------
  out$plot_hyperopt <- plot_hyperopt(
    x      = structure(out, class = "DICEPRO"),
    params = hp_params
  )

  # ---- Persist plots ---------------------------------------------------------
  report_dir <- file.path(output_dir, "report")
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

  ggplot2::ggsave(
    filename = file.path(report_dir, "hyperopt_report.pdf"),
    plot     = out$plot_hyperopt,
    width    = 4L * length(hp_params),
    height   = 4L * length(hp_params) + 2L,
    units    = "in",
    device   = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    filename = file.path(report_dir, "pareto_frontier.pdf"),
    plot     = out$plot,
    width    = 8L,
    height   = 6L,
    units    = "in",
    device   = grDevices::cairo_pdf
  )

  message(sprintf("Results saved to: %s", report_dir))

  structure(out, class = "DICEPRO")
}
