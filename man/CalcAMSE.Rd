\name{CalcAMSE}

\alias{CalcAMSE}

\title{Calculate AMSE}

\description{This function calculates the Asymptotic Mean Squared Error (AMSE) of the Robust Double Exponential (RDE) estimator.}

\usage{CalcAMSE(y, X, Z,family, m = NULL, weights.xz1, weights.xz2, alpha, cutOff1,
  weightF1, cutOff2, weightF2)}

\arguments{
\item{y}{the response vector.}

\item{X}{the model matrix for the mean.}

\item{Z}{the model matrix for the dispersion.}

\item{family}{a character string indicating the family. This can be "poisson" or "binomial".}

\item{m}{in case \code{family} equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.}

\item{weights.xz1}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean.}

\item{weights.xz2}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the dispersion.}

\item{alpha}{a numeric vector containing the beta and gamma estimates.}

\item{cutOff1}{the tuning parameter that is used in the function to diminish the effect of large residuals in the model for the mean.}

\item{weightF1}{the function that is used to diminish the effect of large residuals in the model for the mean.}

\item{cutOff2}{the tuning parameter that is used in the function to diminish the effect of large residuals in the model for the dispersion.}

\item{weightF2}{the function that is used to diminish the effect of large residuals in the model for the dispersion.}
}

\value{\code{calcAMSE} returns a vector of three elements, which are the Asymptoctic Mean Squared Errors of the RDE estimator beta, gamma and alpha.}


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

optionList <- list(huberC = 1.85, tukeyC = 4.75, tol = 1e-4, maxIt = 15)

rob.est <- RDE(y = response ,X=designX, Z=designZ, weights.on.xz1 = "none",
  weights.on.xz2 = "none", weightFunction1 = "Huber", weightFunction2 = "Huber",
  muEstim = 'glm', thetaEstim = 'glm',
  mBin = Data$total, family = "binomial", optionList = optionList)

CalcAMSE(y = response ,X=designX, Z=designZ, family= "binomial", m = Data$total,
  weights.xz1=rob.est$weights.xz1, weights.xz2 = rob.est$weights.xz2,
  alpha = c(rob.est$beta, rob.est$gamma), cutOff1 = rob.est$c1,
  weightF1 = rob.est$v_fun1, cutOff2 = rob.est$c2, weightF2 = rob.est$v_fun2)
}
