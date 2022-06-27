\name{Double Binomial}

\alias{Double Binomial}
\alias{dDBinom}
\alias{pDBinom}
\alias{qDBinom}
\alias{rDBinom}

\title{The Double Binomial Distribution}

\description{Density, distribution function, quantile function and random generation for the double binomial distribution with parameters \code{n}, \code{mu} and \code{theta}.}

\usage{
dDBinom(x, n, mu, theta, correction = TRUE)
pDBinom(q, n, mu, theta)
qDBinom(p, n, mu, theta)
rDBinom(m, n, mu, theta)
}

\arguments{
\item{x}{vector of quantiles.}

\item{q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of trials.}

\item{m}{number of random values to return.}

\item{mu}{vector of (non-negative) means.}

\item{theta}{vector of (non-negative) dispersions, such that the variance equals V(mu)/theta.}

\item{correction}{boolean indicating whether the normalizing constant should be used or not. Default value is TRUE.}
}

\value{
\code{dDBinom} gives the density,
\code{pDBinom} gives the distribution function,
\code{qDBinom} gives the quantile function, and
\code{rDBinom} generates random deviates.
}


\examples{
library(RdeGlmLassoGam)
dDBinom(x=0.9, n=10, mu=0.8, theta=2, correction = TRUE)
pDBinom(q=0.2, n=5, mu=0.5, theta=1)
qDBinom(p=0.6, n=4, mu=0.3, theta=0.5)
rDBinom(m=10, n=9, mu=0.6, theta=1.5)
}


