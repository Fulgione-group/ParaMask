% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_distance_cutoff.R
\name{estimate_distance_cutoff}
\alias{estimate_distance_cutoff}
\title{Estimate distance cutoff from multicopy SNPs}
\usage{
estimate_distance_cutoff(
  het,
  dist_em_rep = 100,
  min_dist = 5,
  max_dist = 1000,
  cdist = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{het}{Data frame including Chromosome, Position, and EM_class.}

\item{dist_em_rep}{Number of EM replicates (e.g. 100).}

\item{min_dist}{Minimum allowed cutoff.}

\item{max_dist}{Maximum allowed cutoff.}

\item{cdist}{Logical. If TRUE, use conservative cutoff.}

\item{verbose}{Logical, print progress.}

\item{em_algorithm_geom}{Function performing EM on geometric mixture.}
}
\value{
A list with \code{dist_cutoff}, \code{p1}, \code{p2}, \code{dist_cutoff_samples}, and \code{distances}.
}
\description{
Runs EM algorithm on distances between multicopy SNPs to identify cutoff separating clusters.
}
