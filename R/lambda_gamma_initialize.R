#' Create gamma-lambda interactive plot
#'
#' Generates an interactive plot showing the relationship between
#' gamma and lambda hyperparameters
#'
#' @return plotly object
#' @export
#' @examples
#' \dontrun{
#' create_gamma_lambda_plot()
#' }
create_gamma_lambda_plot <- function() {
  space <- custom_space()

  # Sample hyperparameters
  samples <- purrr::map(1:200, ~{
    params <- .sample_from_space(space)
    list(
      lambda_      = params$lambda_,
      gamma_factor = params$gamma_factor,
      gamma        = params$lambda_ * params$gamma_factor,
      p_prime      = params$p_prime
    )
  })

  gamma_samples  <- purrr::map_dbl(samples, "gamma")
  lambda_samples <- purrr::map_dbl(samples, "lambda_")

  gamma_factor_min <- min(purrr::map_dbl(samples, "gamma_factor"))
  gamma_factor_max <- max(purrr::map_dbl(samples, "gamma_factor"))

  gamma_range <- seq(min(gamma_samples), max(gamma_samples), length.out = 100)

  lambda_min <- gamma_range / gamma_factor_max
  lambda_max <- gamma_range / gamma_factor_min

  # Create plotly plot
  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = gamma_samples,
      y = lambda_samples,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 3, color = 'royalblue', opacity = 0.4),
      name = 'Samples'
    ) %>%
    plotly::add_trace(
      x = gamma_range,
      y = lambda_min,
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'red'),
      name = '\u03bb = \u03b3 / max(\u03b3_factor)'
    ) %>%
    plotly::add_trace(
      x = gamma_range,
      y = lambda_max,
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'green'),
      name = '\u03bb = \u03b3 / min(\u03b3_factor)'
    ) %>%
    plotly::layout(
      title = 'Echantillons (\u03b3, \u03bb) avec contrainte \u03b3 = \u03bb * \u03b3_factor',
      xaxis = list(title = '\u03b3', type = 'log'),
      yaxis = list(title = '\u03bb', type = 'log'),
      hovermode = 'closest'
    )

  htmlwidgets::saveWidget(p, "gamma_lambda_interactive_plot.html")
  return(p)
}


#' Sample hyperparameters from a defined search space
#'
#' Internal function to sample a set of hyperparameters from a given search space.
#' Each hyperparameter is defined as a vector: c(type, low, high), where `type`
#' can be "uniform", "loguniform", or "randint".
#'
#' If `gamma_factor` exists in the search space, `gamma` will be automatically
#' computed as `gamma = lambda_ * gamma_factor`.
#'
#' @param space List. Each element is a vector of the form c(type, low, high).
#'              For example: list(lambda_ = c("loguniform", 1, 1e5))
#' @return List of sampled hyperparameters. Includes `gamma` if `gamma_factor` exists.
#' @keywords internal
.sample_from_space <- function(space) {
  params <- list()

  for (param_name in names(space)) {
    spec <- space[[param_name]]

    # spec must be a vector: c(type, low, high)
    type <- spec[1]
    low  <- as.numeric(spec[2])
    high <- as.numeric(spec[3])

    params[[param_name]] <- switch(type,
                                   "randint"    = sample(low:(high-1), 1),
                                   "uniform"    = runif(1, low, high),
                                   "loguniform" = exp(runif(1, log(low), log(high))),
                                   stop(paste("Unknown type:", type))
    )
  }

  # Compute gamma if gamma_factor exists
  if ("gamma_factor" %in% names(params)) {
    params$gamma <- params$lambda_ * params$gamma_factor
  }

  return(params)
}

