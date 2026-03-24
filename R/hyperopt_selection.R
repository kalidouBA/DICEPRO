# =============================================================================
# NMF Hyperparameter Selection -- Pareto Frontier & Knee Point
#
# Public API  : best_hyperParams()
# Private fns : .select_knee_pareto()
#               .plot_pareto_kraljic()
#               .make_preddata_fname()
#
# Pipeline:
#   1. Filter valid hyperparameter trials
#   2. Compute constraint deviation
#   3. Extract Pareto frontier
#   4. Detect the knee point (best compromise)
#   5. Optionally generate an interactive Plotly visualisation
# =============================================================================


# -----------------------------------------------------------------------------
# .select_knee_pareto
# -----------------------------------------------------------------------------

#' Select the knee point on a Pareto frontier
#'
#' Identifies the configuration that maximises the perpendicular distance to
#' the line connecting the *utopia* \eqn{(0,0)} and *nadir* \eqn{(1,1)} points
#' in a bi-objective space. Both objectives are min-max normalised before
#' scoring, making the procedure fully scale-free.
#'
#' @param df             data.frame of Pareto-optimal solutions.
#' @param loss_col       Character scalar. Column name of the loss objective
#'   (minimised).
#' @param constraint_col Character scalar. Column name of the constraint
#'   deviation objective (minimised).
#' @param eps            Numeric scalar. Guard against division by zero when
#'   the range of an objective is effectively zero (default \code{1e-12}).
#' @param return_scores  Logical. If \code{TRUE} (default), returns a list with
#'   diagnostic scores; otherwise returns only the selected row.
#'
#' @return
#' When \code{return_scores = TRUE}, a named list:
#' \describe{
#'   \item{best_row}{One-row data.frame of the selected configuration.}
#'   \item{best_index}{Integer row index within \code{df}.}
#'   \item{knee_score}{Numeric vector of perpendicular-distance scores.}
#'   \item{f1_norm}{Normalised loss values.}
#'   \item{f2_norm}{Normalised constraint-deviation values.}
#' }
#' When \code{return_scores = FALSE}, returns \code{best_row} directly.
#'
#' @keywords internal
#' @noRd
.select_knee_pareto <- function(df,
                                loss_col       = "loss",
                                constraint_col = "abs_constraint",
                                eps            = 1e-12,
                                return_scores  = TRUE) {

  stopifnot(is.data.frame(df))
  if (!loss_col       %in% names(df)) stop("Column not found: ", loss_col)
  if (!constraint_col %in% names(df)) stop("Column not found: ", constraint_col)

  norm01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < eps) return(rep(0, length(x)))
    (x - rng[1L]) / diff(rng)
  }

  f1_norm <- norm01(df[[loss_col]])
  f2_norm <- norm01(df[[constraint_col]])

  # Utopia = (0,0), Nadir = (1,1) after normalisation.
  # Perpendicular distance to the line y = x: d = |f1 - f2| / sqrt(2).
  knee_score <- abs(f1_norm - f2_norm) / sqrt(2)
  best_index <- which.max(knee_score)
  best_row   <- df[best_index, , drop = FALSE]

  if (return_scores) {
    list(
      best_row   = best_row,
      best_index = best_index,
      knee_score = knee_score,
      f1_norm    = f1_norm,
      f2_norm    = f2_norm
    )
  } else {
    best_row
  }
}


#' Generate an interactive Pareto frontier plot
#'
#' Builds a Plotly scatter plot showing all valid configurations, the Pareto
#' frontier, and the selected knee-point solution.
#'
#' @param data_all           data.frame of all valid trials.
#' @param data_frontier      data.frame of Pareto-frontier points.
#' @param bestPareto         One-row data.frame of the selected configuration.
#' @param constraint_threshold Numeric scalar. Constraint threshold shown in
#'   plot title metadata.
#'
#' @return A Plotly figure object.
#'
#' @importFrom plotly plot_ly add_trace layout config
#' @keywords internal
#' @noRd
.plot_pareto_kraljic <- function(data_all,
                                 data_frontier,
                                 bestPareto,
                                 constraint_threshold = 0.1) {

  .hover <- function(df) {
    with(df, paste0(
      "<b>\u03BB:</b> ",     formatC(lambda_,       format = "e", digits = 2), "<br>",
      "<b>\u03B3:</b> ",     formatC(gamma,         format = "e", digits = 2), "<br>",
      "<b>Frobenius:</b> ",  formatC(frobNorm,      format = "e", digits = 2), "<br>",
      "<b>Constraint:</b> ", round(constraint, 4),                             "<br>",
      "<b>p_prime:</b> ",    signif(p_prime,   3),                             "<br>",
      "<b>Duration:</b> ",   signif(duration,  3),                             "<br>",
      "<b>Penalty:</b> ",    signif(penalty,   3),                             "<br>",
      "<b>Objective:</b> ",  signif(objectiveValue, 3)
    ))
  }

  frontier_ordered <- data_frontier[order(data_frontier$frobNorm), ]

  plotly::plot_ly(data = data_all, x = ~frobNorm, y = ~abs_constraint) |>
    plotly::add_trace(
      type      = "scatter",
      mode      = "markers",
      marker    = list(color = "gray", size = 6L, opacity = 0.6),
      name      = "Valid solutions",
      text      = .hover(data_all),
      hoverinfo = "text"
    ) |>
    plotly::add_trace(
      data      = frontier_ordered,
      x         = ~frobNorm,
      y         = ~abs_constraint,
      type      = "scatter",
      mode      = "markers+lines",
      marker    = list(color = "red", size = 8L, symbol = "diamond"),
      line      = list(color = "red", width = 2L),
      name      = "Pareto frontier",
      text      = .hover(frontier_ordered),
      hoverinfo = "text"
    ) |>
    plotly::add_trace(
      data      = bestPareto,
      x         = ~frobNorm,
      y         = ~abs_constraint,
      type      = "scatter",
      mode      = "markers",
      marker    = list(color = "blue", size = 14L, symbol = "star",
                       line = list(color = "white", width = 1L)),
      name      = "Best pareto",
      text      = .hover(bestPareto),
      hoverinfo = "text"
    ) |>
    plotly::layout(
      title  = "<b>Frobenius Norm vs Constraint Deviation \u2014 Pareto Frontier</b>",
      xaxis  = list(title = "<b>Frobenius Norm</b>",        type = "log",
                    autorange = "reversed"),
      yaxis  = list(title = "<b>|1 \u2212 Constraint|</b>", type = "log",
                    autorange = "reversed"),
      legend = list(x = 0.02, y = 0.98)
    ) |>
    plotly::config(
      toImageButtonOptions = list(
        format   = "svg",
        filename = "pareto_front_kraljic",
        width    = 500L,
        height   = 300L,
        scale    = 1.5
      ),
      displaylogo            = FALSE,
      modeBarButtonsToRemove = c("lasso2d", "select2d")
    )
}

# -----------------------------------------------------------------------------
# .make_preddata_fname
# -----------------------------------------------------------------------------

#' Build standard hyperparameter result filenames
#'
#' @param lambda Numeric vector of \eqn{\lambda} values.
#' @param gamma  Numeric vector of \eqn{\gamma} values.
#'
#' @return Character vector of filenames (same length as \code{lambda}).
#'
#' @keywords internal
#' @noRd
.make_preddata_fname <- function(lambda, gamma) {
  paste0(
    "hyperopt_results_lambda_", round(lambda, 2),
    "_gamma_",                  round(gamma,  2),
    ".RData"
  )
}


# -----------------------------------------------------------------------------
# best_hyperParams  (public API)
# -----------------------------------------------------------------------------

#' Select optimal hyperparameters using a Pareto frontier
#'
#' Identifies the best \eqn{(\lambda, \gamma)} pair from a set of
#' optimisation trials by:
#' \enumerate{
#'   \item Filtering trials that violate the constraint threshold.
#'   \item Computing constraint deviation (\code{abs_constraint}).
#'   \item Extracting the Pareto frontier (Frobenius norm vs. loss).
#'   \item Selecting the knee point -- the best loss/constraint trade-off.
#'   \item Generating an interactive Pareto plot saved under
#'         \code{savePaths/report/}.
#' }
#'
#' The \code{W} and \code{H} matrices produced by the winning trial are
#' returned directly, avoiding any file I/O for the result matrices.
#'
#' @param trials_df          data.frame of optimisation trials as returned by
#'   \code{\link{research_hyperOpt}} (\code{$trials}).
#' @param W                  List of \eqn{W} matrices, one per trial, as
#'   returned by \code{research_hyperOpt} (\code{$w_unknown}).
#' @param H                  List of deconvolution outputs (\eqn{H} matrices),
#'   one per trial, as returned by \code{research_hyperOpt}
#'   (\code{$out_deconvolution}).
#' @param savePaths          Character scalar. Root directory for all outputs
#'   (plot saved under \code{savePaths/report/}).
#' @param constraint_threshold Numeric scalar. Maximum allowed value of
#'   \eqn{|1 - \text{constraint}|}; trials exceeding this are discarded
#'   (default \code{0.1}).
#'
#' @return
#' A named list, or \code{invisible(NULL)} when no valid configuration
#' survives filtering:
#' \describe{
#'   \item{hyperparameters}{List with \code{lambda} and \code{gamma} scalars.}
#'   \item{metrics}{List with \code{loss} and \code{constraint} scalars.}
#'   \item{W}{The \eqn{W} matrix of the selected trial.}
#'   \item{H}{The \eqn{H} (deconvolution) matrix of the selected trial.}
#'   \item{plot}{Plotly figure of the Pareto frontier.}
#' }
#'
#' @importFrom dplyr mutate filter
#' @importFrom KraljicMatrix get_frontier
#' @export
best_hyperParams <- function(trials_df,
                             W,
                             H,
                             savePaths,
                             constraint_threshold = 0.1) {
  frob_H <- loss <- abs_constraint <- lambda_ <- gamma <- NULL

  # ---- Guard: empty input --------------------------------------------------
  if (nrow(trials_df) == 0L) {
    warning("trials_df is empty -- no configurations to evaluate.")
    return(invisible(NULL))
  }

  # ---- Attach trial index before any filtering -----------------------------
  # Preserves the link to the corresponding W / H list elements.
  trials_df$trial_index <- seq_len(nrow(trials_df))

  # ---- Derive auxiliary columns and apply constraint filter ----------------
  trials_df <- trials_df |>
    dplyr::mutate(
      log_lambda     = log10(lambda_ + 1),
      log_gamma      = log10(gamma   + 1),
      log_frob       = log10(frobNorm),
      abs_constraint = abs(1 - constraint)
    ) |>
    dplyr::filter(abs_constraint <= constraint_threshold)

  if (nrow(trials_df) == 0L) {
    warning(sprintf(
      "No trials pass the constraint filter (abs_constraint <= %g).",
      constraint_threshold
    ))
    return(invisible(NULL))
  }

  # ---- Pareto frontier -----------------------------------------------------
  frontier_result <- KraljicMatrix::get_frontier(
    data     = trials_df,
    x        = frob_H,
    y        = loss,
    quadrant = "bottom.right"
  )
  frontier_df <- trials_df[rownames(frontier_result), , drop = FALSE]

  # ---- Knee-point selection ------------------------------------------------
  knee_res <- .select_knee_pareto(
    df             = frontier_df,
    loss_col       = "frob_H",
    constraint_col = "var_H"
  )

  best_row   <- knee_res$best_row
  best_trial <- best_row$trial_index   # index into the original W / H lists

  # ---- Interactive plot ----------------------------------------------------
  plot_res <- .plot_pareto_kraljic(
    data_all             = trials_df,
    data_frontier        = frontier_df,
    bestPareto           = best_row,
    constraint_threshold = constraint_threshold
  )

  # ---- Assemble output -----------------------------------------------------
  list(
    hyperparameters = list(
      lambda = best_row$lambda_,
      gamma  = best_row$gamma
    ),
    metrics = list(
      loss       = best_row$loss,
      constraint = best_row$constraint
    ),
    trials = trials_df,
    W      = W[[best_trial]],
    H      = H[[best_trial]][,-1],
    plot = plot_res
  )
}
