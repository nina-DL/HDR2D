% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All_HDR_measures.R
\name{CDF.e.2d}
\alias{CDF.e.2d}
\title{This function computes M3:e-CDF distance in Definition 5 (concentration measure)}
\usage{
CDF.e.2d(sn, eps = NULL)
}
\arguments{
\item{sn}{data sample of dimension n x d, with d = 2}

\item{eps}{a bidimensional hyperparameter with positive elements, representing the half-length of the rectangular neighborhood}
}
\value{
the M3:e-CDF distance for each point in the sample
}
\description{
This function computes M3:e-CDF distance in Definition 5 (concentration measure)
}
\examples{
# Generate some bivariate data
R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
CDF.e.2d(sn = draws_2d)
}
