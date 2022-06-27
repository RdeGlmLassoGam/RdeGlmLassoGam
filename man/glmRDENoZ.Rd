\name{glmRDENoZ}
\alias{glmRDENoZ}

\title{Estimating the coefficients of the Robust Double Exponential (RDE) estimator in case there is no dispersion part}

\description{This function estimates the RDE coefficients in case there is no dispersion part.}

\usage{glmRDENoZ(y, X, m = NULL, family, beta.ini, weights.on.xz1 = "none",
  contX = NULL, weightFunction1 = "Huber",
  optionList = list(huberC = 2, tukeyC = 3, tol = 1e-4, maxIt = 100, alpha = 0.75))}

\arguments{
\item{y}{the response vector.}

\item{X}{the model matrix for the mean.}

\item{m}{in case \code{family} equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.}

\item{family}{a character string indicating the family. This can be "poisson" or "binomial".}

\item{beta.ini}{optional initial value for the beta coefficients.}

\item{weights.on.xz1}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. The tuning parameter \code{alpha} can be provided in the \code{optionList}.}

\item{contX}{indices of the columns in \code{X} representing the the continuous variables.  When \code{weights.on.xz1} is equal to "covMcd", this parameter should be specified. Default value is NULL.}

\item{weightFunction1}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".}

\item{optionList}{list that should contain the tuning parameter \code{huberC} in case \code{weightFunction} equals "Huber", or the tuning parameter \code{tukeyC} in case \code{weightFunction} equals "Tukey". Furthermore, the tuning parameter \code{alpha} for the function robustbase::covMcd can be provided, with default value 0.7 . Finally, \code{maxIt} which is the maximum number of iterations and \code{tol}, which is the tolerance can be provided here. Default value is list(huberC = 2, tukeyC = 3, tol = 1e-4, maxIt = 100, alpha = 0.75).}
}

\value{
\code{glmRDENoZ} returns a list that contains the following items:
\item{betaEstimate}{the estimated beta coefficients.}
\item{thetaStart}{the beta coefficients from which the algorithm starts.}
\item{cutOff1}{the tuning parameter that is used in the function to diminish the effect of large residuals in the model for the mean.}
\item{weightF1}{the function that is used to diminish the effect of large residuals in the model for the mean.}
\item{weights.xz1}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean.}
\item{M}{M-matrix for this model.}
\item{Q}{Q-matrix for this model.}
\item{ASVar}{Asymptotic variance matrix.}
}

\examples{
library(RdeGlmLassoGam)
library(tidyr)
d <- as.data.frame(UCBAdmissions)
d <- spread(d, Admit, Freq)
d[order(d$Dept), ]

Data <- d
attach(Data)

response <- Admitted / (Admitted +  Rejected)

Data$total <- (Admitted +  Rejected)

family <- "binomial"

designX = cbind(rep(1,nrow(Data)), (Dept=="B"), (Dept=="C"), (Dept=="D"), (Dept=="E"),
  (Dept=="F"), (Gender=="Female"))
colnames(designX) = c("Intercept", "DeptB", "DeptC", "DeptD", "DeptE", "DeptF",
  "GenderFemale")

optionList <- list(huberC = 1.85,  tukeyC = 4.75, tol = 1e-4,  maxIt = 15)

rob.est <- glmRDENoZ(y = response ,X=designX, weights.on.xz1 = "none",
  weightFunction1 = "Tukey", m = Data$total, family = "binomial",
  optionList = optionList)
}


