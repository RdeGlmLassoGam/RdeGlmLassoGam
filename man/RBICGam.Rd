\name{RBICGam}

\alias{RBICGam}

\title{Estimating the RBIC value for the dispersion parameters.}

\description{This function estimates the RBIC value for the mean parameters. This can be used for the tuning of the smoothing parameter in the dispersion model of the RDE estimator in the GAM setting.}

\usage{RBICGam(y, X, Z, mBin = NULL, family,
                muEstim = 'glm', thetaEstim = 'glm',
                p1=3, p2=NULL, K1=30, K2=NULL, sp1=-1, sp2=NULL,
                smooth.basis1="cr", smooth.basis2="cr",
                intercept = TRUE, standardize = TRUE,
                beta.ini = NULL,gamma.ini = NULL,
                lambdaBeta=0, lambdaGamma=0,
                weights.on.xz1 = "none", weights.on.xz2 = NULL,
                rowLev = TRUE, contX = NULL, contZ = NULL,
                weightFunction1 = "Huber", weightFunction2 = NULL,
                optionList = list(huberC = 2,
                                  tukeyC = 3,
                                  tol = 1e-4,
                                  maxIt = 100,
                                  alpha = 0.75
                ))}

\arguments{
\item{y}{the response vector.}

\item{X}{the model matrix for the mean.}

\item{Z}{the model matrix for the dispersion.}

\item{mBin}{in case \code{family} equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.}

\item{family}{a character string indicating the family. This can be "poisson" or "binomial".}

\item{muEstim}{a character string indicating the estimating method for the mean. This can be "glm" in case of a generalised linear model, "pen" in case of a penalised generalised linear model or "gam" in case of a generalised additive model.}

\item{thetaEstim}{a character string indicating the estimating method for the dispersion. This can be "glm" in case of a generalised linear model, "pen" in case of a penalised generalised linear model or "gam" in case of a generalised additive model.}

\item{p1}{order of the basis for the mean model in case \code{muEstim} equals "gam". It depends on the option of smooth.basis1. Default value is 3.}

\item{K1}{number of knots of the basis for the mean model in case \code{muEstim} equals "gam"; dependent on the option of smooth.basis1. Default value is 30.}

\item{sp1}{a vector of smoothing parameter for the mean model in case \code{muEstim} equals "gam". If only one value is specified, it will be used for all smoothing parameters in this model.}

\item{smooth.basis1}{the specification of basis for the mean model in case \code{muEstim} equals "gam". Four choices are available: "tp" = thin plate regression spline, "cr" = cubic regression spline, "ps" = P-splines, "tr" = truncated power spline.}

\item{p2}{order of the basis for the dispersion model in case \code{thetaEstim} equals "gam". It depends on the option of smooth.basis2. Default value is 3.}

\item{K2}{number of knots of the basis for the dispersion model in case \code{thetaEstim} equals "gam"; dependent on the option of smooth.basis2. Default value is 30.}

\item{sp2}{a vector of smoothing parameters for the dispersion model in case \code{thetaEstim} equals "gam". If only one value is specified, it will be used for all smoothing parameters in this model.}

\item{smooth.basis2}{the specification of basis for the dispersion model in case \code{thetaEstim} equals "gam". Four choices are available: "tp" = thin plate regression spline, "cr" = cubic regression spline, "ps" = P-splines, "tr" = truncated power spline. Default value is "cr".}

\item{intercept}{Logical indicating wether an intercept should be present. Default value is TRUE.}

\item{standardize}{Logical indicating whether design matrices X and Z are robustly standardized. If FALSE, the continuous predictors are robustly standardized in the algorithm (minus median, divided by median absolute deviation). This is strongly advised when using a penalized regression. Default value is TRUE.}

\item{beta.ini}{optional initial value for the beta coefficients.}

\item{gamma.ini}{optional initial value for the gamma coefficients.}

\item{lambdaBeta}{(Non-negative) regularization parameter for lasso on mean GLM. Default value is 0, which means no regularization.}

\item{lambdaGamma}{(Non-negative) regularization parameter for lasso on dispersion GLM. Default value is 0, which means no regularization.}

\item{weights.on.xz1}{only of importance when there are more observations than variabels. It is a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. The tuning parameter \code{alpha} can be provided in the \code{optionList}.}

\item{weights.on.xz2}{only of importance when there are more observations than variabels. It is a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the dispersion. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. The tuning parameter \code{alpha} can be provided in the \code{optionList}. Default value is NULL, meaning that the same weights as for the dispersion model are used.}

 \item{rowLev}{logical, only important when there are more variables than observations, which is TRUE (default) when observations get a low weight once the whole observation is detected as a leverage.}

\item{contX}{indices of the columns in \code{X} representing the the continuous variables.  When \code{weights.on.xz1} or \code{weights.on.xz2} is equal to "covMcd", this parameter or \code{contZ} should be specified. Default value is NULL.}

\item{contZ}{indices of the columns in \code{Z} representing the the continuous variables. When \code{weights.on.xz1} or \code{weights.on.xz2} is equal to "covMcd", this parameter or \code{contX} should be specified. Default value is NULL.}

\item{weightFunction1}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".}

\item{weightFunction2}{a character string indicating which weight function is used to diminish the effect of large residuals in the model for the dispersion. This can be "Huber" or "Tukey". Default value is NULL, meaning that the same function as for the mean model is used.}

\item{optionList}{list that should contain the tuning parameter \code{huberC} in case \code{weightFunction} equals "Huber", or the tuning parameter \code{tukeyC} in case \code{weightFunction} equals "Tukey". Furthermore, the tuning parameter \code{alpha} for the function robustbase::covMcd can be provided, with default value 0.75 . Finally, \code{maxIt} which is the maximum number of iterations and \code{tol}, which is the tolerance can be provided here. Default value is list(huberC = 2, tukeyC = 3, tol = 1e-4, maxIt = 100, alpha = 0.75).}
}

\value{
\code{RBICGam} returns a list that contains the following items:
\item{RBICGam}{the estimated RBIC value for the dispersion parameters.}
\item{fit}{the fit of which the RBIC value for the dispersion parameters was determined.}
}


\examples{
library(RdeGlmLassoGam)
n = 250
true.family <- poisson(link="log")
    
x <- cbind(runif(n,-3,3))
true.eta <- 2*sin(x+1)+2
true.mu <- true.family$linkinv(true.eta)

z <- cbind(runif(n,0,1))
true.eta.z <-  cos(2*z[,1])
true.theta <- exp(true.eta.z)
    
y=rep(NA,n)
for(j in 1:n){
  y[j] <- rDPois(n = 1, mu = true.mu[j], theta = true.theta[j]) # generate response vector
}
                
res <- RBICGam(X=x, Z = z, y=y, family=true.family, p1=3,
            muEstim = 'gam', thetaEstim = 'gam', K1=10,
            sp1=5,sp2=5,
            smooth.basis1="cr",smooth.basis2="cr",
            contX = 1, contZ = 1, weights.on.xz1="none", weights.on.xz2="none",
            weightFunction1 = "Tukey", weightFunction2 = NULL,
            optionList = list(tukeyC = 4.5,huberC = 1.8,
                              tol = 1e-4,
                              maxIt = 100,
                              alpha = 0.7
            ))
res$RBICGam
}
