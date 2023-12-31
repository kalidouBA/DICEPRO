% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSx.R
\name{run_CSx}
\alias{run_CSx}
\title{Run CIBERSORTx Deconvolution Method}
\usage{
run_CSx(bulk, reference, cibersortx_email, cibersortx_token)
}
\arguments{
\item{bulk}{A matrix of bulk RNA-seq data.}

\item{reference}{A matrix of reference cell type gene expression profiles.}

\item{cibersortx_email}{The CIBERSORTx account email.}

\item{cibersortx_token}{The CIBERSORTx account token.}
}
\value{
A matrix containing the estimated cell type proportions.
}
\description{
This function runs the CIBERSORTx deconvolution method using Docker with the provided email and token.
The method is used to estimate cell type proportions in bulk RNA-seq data.
}
\details{
The function performs the following steps:
\enumerate{
\item Runs CIBERSORTx with the provided email and token.
\item Reads the output results from CIBERSORTx.
\item Extracts cell type proportions and returns them as a matrix.
}
}
