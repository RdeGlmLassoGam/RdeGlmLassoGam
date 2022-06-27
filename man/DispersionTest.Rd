\name{DispersionTest}

\alias{DispersionTest}

\title{Hypothesis testing for RDE estimator}

\description{This function performs the hypothesis test that \code{nzeroBeta} specific beta coefficients and \code{nzeroGamma} specific gamma coefficients of the RDE estimator are zero.}

\usage{DispersionTest(y, X, Z, m = NULL, alphaHAlt, alphaHNull, nzeroBeta, nzeroGamma,
  family, weightFunction1, weightFunction2 = NULL, weights.on.xz1, weights.on.xz2,
  optionList, M, Q)}

\arguments{
\item{y}{the response vector.}

\item{X}{the model matrix for the mean. The order of the variables should be the same as in \code{alphaHNull}.}

\item{Z}{the model matrix for the dispersion. The order of the variables should be the same as in \code{alphaHNull}.}

\item{m}{in case \code{family} equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.}

\item{alphaHAlt}{the estimated coefficients under the alternative hypothesis. The order of the coefficients should  be the same as specified in \code{alphaHNull}.}

\item{alphaHNull}{the estimated coefficients under the null hypothesis. In this vector, the beta coefficients are followed by the gamma coefficients.  Note that for both parts, we first have the non-zero coefficients followed by the zero coeffcients; e.g. \eqn{(\beta_1, \beta_2, ..., 0,..., 0, \gamma_1, \gamma_2, ..., 0, 0)}.}

\item{nzeroBeta}{number of beta coefficients that are zero under the null hypothesis.}

\item{nzeroGamma}{number of gamma coefficients that are zero under the null hypothesis.}

\item{family}{a character string indicating the family. This can be "poisson" or "binomial".}

\item{weightFunction1}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".}

\item{weightFunction2}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the dispersion. This can be "Huber" or "Tukey". Default value is NULL, meaning that the same function as for the mean model is used.}

\item{weights.on.xz1}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean.}

\item{weights.on.xz2}{a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the dispersion.}

\item{optionList}{list that should contain the tuning parameter \code{huberC} in case \code{weightFunction} equals "Huber", or the tuning parameter \code{tukeyC} in case \code{weightFunction} equals "Tukey".}

\item{M}{M-matrix from the alternative hypothesis.}

\item{Q}{Q-matrix from the alternative hypothesis.}
}

\value{
\code{DispersionTest} returns a list that contains the following items:
\item{DispTestValue}{the test value of the hypothesis test.}
\item{DispTestEigValue}{the positive eigenvalues discussed in the paper.}
\item{DispTestpValue}{the p-value of the hypothesis test.}
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

optionList <- list(huberC = 1.85,  tukeyC = 4.75, tol = 1e-4,  maxIt = 15)

rob.est <- RDE(y = response ,X=designX, Z=designZ, weights.on.xz1 = "none",
  weights.on.xz2 = "none", weightFunction1 = "Tukey", weightFunction2 = "Huber",
  mBin = Data$total, family = "binomial", optionList = optionList)

designX <- designX[,c(1,3,4,5,6,2,7)]

rob.est_Final <- RDE(y = response ,X=designX, Z=designZ,
  weights.on.xz1 = "none", weights.on.xz2 = "none", weightFunction1 = "Tukey",
  weightFunction2 = "Tukey", mBin = Data$total, family = family, optionList = optionList)

rob.est_Final_NULL <- RDE(y = response ,X=designX[,1:5], Z=designZ,
  weights.on.xz1 = "none", weights.on.xz2 = "none", weightFunction1 = "Tukey",
  weightFunction2 = "Tukey", mBin = Data$total, family = family, optionList = optionList)

DispTest <- DispersionTest(y = response ,X = designX, Z = designZ, m = Data$total,
  alphaHAlt = c(rob.est_Final$betaEstimate,rob.est_Final$gammaEstimate),
  alphaHNull = c(rob.est_Final_NULL$betaEstimate,0,0,rob.est_Final_NULL$gammaEstimate),
  nzeroBeta = 2, nzeroGamma = 0, family=family,
  weights.on.xz1 = rob.est_Final$weights.xz1,
  weights.on.xz2 = rob.est_Final$weights.xz2, weightFunction1 = "Tukey",
  optionList = optionList, M = rob.est_Final$M, Q = rob.est_Final$Q)
}
