% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_em_loglikelihood.R
\name{plot_em_loglikelihood}
\alias{plot_em_loglikelihood}
\title{Plot EM log-likelihood over iterations}
\usage{
plot_em_loglikelihood(log_likelihood, outpath, ID = "")
}
\arguments{
\item{log_likelihood}{A matrix or data.frame with columns "iteration" and "logLikelihood".}

\item{outpath}{Output directory.}

\item{ID}{Optional identifier to include in the filename.}
}
\value{
A saved PDF plot.
}
\description{
Plots the log-likelihood values from EM iterations and saves to a PDF.
}
