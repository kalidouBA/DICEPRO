# FIX: declare NSE symbol used in ggplot2::aes() to silence R CMD CHECK NOTE.
utils::globalVariables("bound_label")

# -----------------------------------------------------------------------------
# .custom_space  [private]
# -----------------------------------------------------------------------------

#' Default hyperparameter search space for gamma-lambda visualisation
#'
#' Returns the \code{restrictionEspace} search space where \code{gamma} is
#' the base variable and \code{lambda_} is derived as
#' \code{lambda_ = gamma * lambda_factor}.
#'
#' @return Named list of atomic vectors \code{c(type, low, high)}.
#' @keywords internal
#' @noRd
.custom_space <- function() {
  list(
    gamma         = c("loguniform", 1,    1e5),
    lambda_factor = c("loguniform", 2,    1e2),
    p_prime       = c("loguniform", 1e-1, 1)
  )
}


# -----------------------------------------------------------------------------
# create_gamma_lambda_plot  [public]
# -----------------------------------------------------------------------------

#' Create a gamma-lambda hyperparameter plot
#'
#' Generates a ggplot2 scatter plot showing the relationship between
#' \eqn{\gamma} and \eqn{\lambda} sampled from the search space defined by
#' \code{hspaceTechniqueChoose}, overlaid with the feasibility bounds when
#' applicable.
#'
#' @param hspaceTechniqueChoose Character scalar. Search-space strategy:
#'   \describe{
#'     \item{\code{"all"}}{Full independent log-uniform grid for
#'       \code{lambda_}, \code{gamma}, \code{p_prime}.
#'       No constraint between \eqn{\lambda} and \eqn{\gamma} — feasibility
#'       bounds are not shown.}
#'     \item{\code{"restrictionEspace"}}{Restricted space where \code{gamma}
#'       is the base variable and \code{lambda_} is derived as
#'       \code{lambda_ = gamma * lambda_factor} with
#'       \code{lambda_factor} in \eqn{[2, 100]}.
#'       The lower (\eqn{\lambda = 2\gamma}) and upper
#'       (\eqn{\lambda = 100\gamma}) feasibility bounds are drawn.}
#'   }
#' @param n_samples Positive integer. Number of configurations to draw.
#'   Default: \code{200L}.
#'
#' @return A ggplot2 figure object (invisibly also saved as
#'   \code{gamma_lambda_plot.pdf} in the working directory).
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_x_log10
#'   scale_y_log10 labs theme_bw ggsave
#' @importFrom stats runif
#' @export
#' @examples
#' \dontrun{
#' # Independent sampling — no constraint lines
#' create_gamma_lambda_plot(hspaceTechniqueChoose = "all")
#'
#' # Restricted space — lambda_ = gamma * lambda_factor, factor in [2, 100]
#' create_gamma_lambda_plot(hspaceTechniqueChoose = "restrictionEspace")
#' }
create_gamma_lambda_plot <- function(hspaceTechniqueChoose = c("all", "restrictionEspace"),
                                     n_samples = 200L) {

  hspaceTechniqueChoose <- match.arg(hspaceTechniqueChoose)

  # ---- Build the raw space (same logic as run_experiment) ------------------
  raw_space <- switch(
    hspaceTechniqueChoose,
    all = list(
      lambda_ = c("loguniform", 1,    1e8),
      gamma   = c("loguniform", 1,    1e8),
      p_prime = c("loguniform", 1e-6, 1)
    ),
    restrictionEspace = .custom_space()
  )

  # ---- Sample n_samples configurations ------------------------------------
  samples <- lapply(seq_len(n_samples), function(i) .sample_from_space(raw_space))

  plot_df <- data.frame(
    gamma   = vapply(samples, function(s) as.numeric(s$gamma),   numeric(1L)),
    lambda_ = vapply(samples, function(s) as.numeric(s$lambda_), numeric(1L))
  )

  # ---- Build the plot -------------------------------------------------------
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = gamma, y = lambda_)) +
    ggplot2::geom_point(colour = "royalblue", size = 1.5, alpha = 0.45) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  if (hspaceTechniqueChoose == "restrictionEspace") {
    gamma_range <- exp(seq(log(min(plot_df$gamma)), log(max(plot_df$gamma)),
                           length.out = 100L))
    bounds_df <- data.frame(
      gamma       = rep(gamma_range, 2L),
      lambda_     = c(2 * gamma_range, 100 * gamma_range),
      bound_label = rep(
        c("\u03bb = 2\u03b3  (lower bound)", "\u03bb = 100\u03b3  (upper bound)"),
        each = 100L
      )
    )

    p <- p +
      ggplot2::geom_line(
        data    = bounds_df,
        mapping = ggplot2::aes(x = gamma, y = lambda_, colour = bound_label),
        linetype = "dashed",
        linewidth = 0.8
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "\u03bb = 2\u03b3  (lower bound)"   = "red",
          "\u03bb = 100\u03b3  (upper bound)" = "green4"
        ),
        name = NULL
      )

    title_str <- paste0(
      "\u03b3 vs \u03bb  \u2014  restrictionEspace  ",
      "(\u03bb = \u03b3 \u00d7 \u03bb_factor,  \u03bb_factor \u2208 [2, 100])"
    )
  } else {
    title_str <- "\u03b3 vs \u03bb  \u2014  all  (independent log-uniform sampling)"
  }

  p <- p +
    ggplot2::labs(
      title = title_str,
      x = expression(gamma ~ "(log scale)"),
      y = expression(lambda ~ "(log scale)")
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave("gamma_lambda_plot.pdf", plot = p, width = 7, height = 5)
  p
}


#' Sample hyperparameters from a defined search space
#'
#' Each hyperparameter is defined as an atomic vector \code{c(type, low, high)}.
#' When \code{lambda_factor} is present, \code{lambda_} is automatically
#' derived as \code{lambda_ = gamma * lambda_factor}.
#'
#' @param space Named list. Each element is \code{c(type, low, high)}.
#' @return Named list of sampled values, including derived \code{lambda_}.
#' @keywords internal
#' @importFrom stats runif
.sample_from_space <- function(space) {
  params <- list()

  for (param_name in names(space)) {
    spec <- space[[param_name]]
    type <- spec[1]
    low  <- as.numeric(spec[2])
    high <- as.numeric(spec[3])

    params[[param_name]] <- switch(type,
                                   "randint"    = sample(seq.int(low, high - 1L), 1L),
                                   "uniform"    = runif(1, low, high),
                                   "loguniform" = exp(runif(1, log(low), log(high))),
                                   stop(paste("Unknown type:", type))
    )
  }

  # gamma is the base variable — lambda_ is derived, not sampled
  if ("lambda_factor" %in% names(params) && "gamma" %in% names(params)) {
    params$lambda_ <- params$gamma * params$lambda_factor
  }

  return(params)
}
