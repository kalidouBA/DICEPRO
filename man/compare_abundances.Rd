% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{compare_abundances}
\alias{compare_abundances}
\title{Compare Two Data Frames Column by Column}
\usage{
compare_abundances(x, y, method = "weighted")
}
\arguments{
\item{x}{A numeric data frame.}

\item{y}{A numeric data frame with the same dimensions as \code{x}.}

\item{method}{Aggregation method for combined metrics: "weighted" (default, weights by SD) or "unweighted".}
}
\value{
A list with:
\describe{
\item{by_column}{Named list of per-column metrics: NRMSE, R^2 adjusted, ICC3, and CCC.}
\item{combined}{Aggregated metrics across all columns (weighted or unweighted).}
\item{method}{Description of the aggregation method used.}
}
}
\description{
This function compares two numeric data frames of equal dimensions by computing
statistical metrics (NRMSE, adjusted R^2, ICC3, and CCC) for each column,
after Z-score normalization. It returns per-column results and a combined summary.
}
\details{
\itemize{
\item Z-score normalization is applied to each column.
\item If standard deviation is zero, the column is set to zero.
\item ICC3 is computed if the \pkg{irr} package is available.
\item CCC (Concordance Correlation Coefficient) is computed using the formula from Lin (1989).
}
}
\examples{
\dontrun{
df1 <- data.frame(a = rnorm(100), b = runif(100))
df2 <- data.frame(a = df1$a + rnorm(100, sd = 0.1), b = df1$b + rnorm(100, sd = 0.1))
compare_abundances(df1, df2)
}

}
\references{
Lin, L. I. (1989). A concordance correlation coefficient to evaluate reproducibility. \emph{Biometrics}, 45(1), 255-268.
}
