% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_type_assignments.R
\name{average.mat.rows}
\alias{average.mat.rows}
\title{Average Matrix Rows}
\usage{
average.mat.rows(m, ids, f = colMeans)
}
\arguments{
\item{m}{matrix}

\item{ids}{ids of the samples you care about}

\item{f}{aggregation function}
}
\value{
a list of binary scores
}
\description{
function for computing matrix averages
}
