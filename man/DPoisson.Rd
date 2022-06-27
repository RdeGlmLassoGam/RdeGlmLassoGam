\name{Double Poisson}

\alias{Double Poisson}
\alias{dDPois}
\alias{pDPois}
\alias{qDPois}
\alias{rDPois}

\title{The Double Poisson Distribution}

\description{Density, distribution function, quantile function and random generation for the double Poisson distribution with parameters \code{mu} and \code{theta}.}

\usage{
dDPois(x, mu, theta, correction = TRUE)
pDPois(q, mu, theta)
qDPois(p, mu, theta)
rDPois(n, mu, theta)
}


\arguments{
\item{x}{vector of (non-negative integer) quantiles.}

\item{q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of random values to return.}

\item{mu}{vector of (non-negative) means.}

\item{theta}{vector of (non-negative) dispersions, such that the variance equals V(mu)/theta.}

\item{correction}{boolean indicating whether the normalizing constant should be used or not. Default value is TRUE.}
}

\value{
\code{dDPois} gives the density,
\code{pDPois} gives the distribution function,
\code{qDPois} gives the quantile function, and
\code{rDPois} generates random deviates.
}

\examples{
library(RdeGlmLassoGam)
dDPois(x=5, mu=4, theta=1, correction = TRUE)
pDPois(q=3, mu=2, theta=2)
qDPois(p=0.6, mu=3, theta=0.5)
rDPois(10, mu=6, theta=1.5)
}


