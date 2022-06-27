\name{CalcTuningParam}

\alias{CalcTuningParam}

\title{Determine the optimal tuning parameter for the RDE estimator}

\description{This function determines the smallest tuning parameter for the RDE estimator such that the desired AMSE is obtained.}

\usage{CalcTuningParam(cValues, AMSECrit, y, X, Z, m, family, weights.on.xz1 = "none",
  weights.on.xz2=NULL, weightFunction1, weightFunction2 = NULL)}

\arguments{
\item{cValues}{the vector with values for the tuning parameter that should be checked.}

\item{AMSECrit}{the critical AMSE value.}

\item{y}{the response vector.}

\item{X}{the model matrix for the mean.}

\item{Z}{the model matrix for the dispersion.}

\item{m}{in case \code{family} equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.}

\item{family}{a character string indicating the family. This can be "poisson" or "binomial".}

\item{weights.on.xz1}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd..}

\item{weights.on.xz2}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the dispersion. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. Default value is NULL, meaning that the same weights as for the dispersion model are used.}

\item{weightFunction1}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".}

\item{weightFunction2}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the dispersion. This can be "Huber" or "Tukey". Default value is NULL, meaning that the same function as for the mean model is used.}
}

\value{
\code{CalcTuningParam} returns a list that contains the following items:
\item{tunPar}{vector of three elements, which are the lowest tuning parameters such that \code{AMSECrit} is obtained respectively for beta, gamma and alpha.}
\item{AMSEBeta}{vector of AMSE values obtained for each possible tuning paramater in \code{cValues} focussing on beta.}
\item{AMSEGamma}{vector of AMSE values obtained for each possible tuning paramater in \code{cValues} focussing on gamma.}
\item{AMSEAlpha}{vector of AMSE values obtained for each possible tuning paramater in \code{cValues} focussing on alpha.}
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
designZ = as.matrix(rep(1,nrow(Data)))
colnames(designZ) = "Intercept"

CalcTuningParam(cValues = seq(1.5,2.1,0.1), AMSECrit = 0.9, y = response ,X=designX,
  Z=designZ, family= "binomial", m = Data$total, weightFunction1 = "Huber")
}
