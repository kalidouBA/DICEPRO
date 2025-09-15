#' Sample from search space
#'
#' Internal function to sample parameters from search space
#'
#' @param space Search space definition (list of hyperparameters with type, low, high)
#' @return List of sampled parameters
#' @keywords internal
.sample_from_space <- function(space) {
  params <- list()
  for (param_name in names(space)) {
    spec <- space[[param_name]]

    # Assumer que spec = c(type, low, high)
    type <- spec[1]
    low  <- as.numeric(spec[2])
    high <- as.numeric(spec[3])

    switch(type,
           "choice" = { stop("Choice type not implemented in vector format") },
           "randint" = { params[[param_name]] <- sample(low:(high-1), 1) },
           "uniform" = { params[[param_name]] <- runif(1, low, high) },
           "loguniform" = { params[[param_name]] <- exp(runif(1, log(low), log(high))) },
           stop(paste("Unknown type:", type))
    )
  }
  return(params)
}


#' Custom search space for gamma, lambda_, and p_prime parameters
#'
#' Defines a search space for NMF hyperparameter optimization.
#'
#' @return List of hyperparameter specifications
#' @export
#' @examples
#' space <- custom_space()
custom_space <- function() {
  list(
    gamma    = list(type = "loguniform", low = 1, high = 1e5),
    lambda_  = list(type = "loguniform", low = 2, high = 1e2),
    p_prime  = list(type = "loguniform", low = 1e-1, high = 1)
  )
}

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
  samples <- purrr::map(1:200, ~{
    params <- .sample_from_space(space)
    list(
      gamma = params$gamma,
      lambda_ = params$gamma * params$lambda_factor,
      p_prime = params$p_prime
    )
  })

  gamma_samples <- purrr::map_dbl(samples, "gamma")
  lambda_samples <- purrr::map_dbl(samples, "lambda_")

  gamma_range <- 10^seq(0, 5, length.out = 100)
  lambda_min <- 2 * gamma_range
  lambda_max <- 1e2 * gamma_range

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
      name = 'λ = 2 × γ'
    ) %>%
    plotly::add_trace(
      x = gamma_range,
      y = lambda_max,
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'green'),
      name = 'λ = 10⁴ × γ'
    ) %>%
    plotly::layout(
      title = "Échantillons (γ, λ) avec contrainte λ = k·γ, k ∈ [10, 10⁴]",
      xaxis = list(title = 'γ', type = 'log'),
      yaxis = list(title = 'λ', type = 'log'),
      hovermode = 'closest'
    )

  htmlwidgets::saveWidget(p, "gamma_lambda_interactive_plot.html")
  return(p)
}
