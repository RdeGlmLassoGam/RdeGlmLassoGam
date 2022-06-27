#calculate EBIC for gamma
EBICGam = function(y, X, Z, mBin = NULL, family,
                   muEstim = 'glm', thetaEstim = 'glm',
                   p1=3, p2=NULL, K1=30, K2=NULL, sp1=-1, sp2=NULL,
                   smooth.basis1="cr", smooth.basis2="cr",
                   intercept = TRUE, standardize = TRUE,
                   beta.ini= NULL,gamma.ini = NULL,
                   lambdaBeta=0, lambdaGamma=0,
                   weights.on.xz1 = "none", weights.on.xz2 = NULL,
                   rowLev = TRUE, contX = NULL, contZ = NULL,
                   weightFunction1 = "Huber", weightFunction2 = NULL,
                   optionList = list(huberC = 2,
                                     tukeyC = 3,
                                     tol = 1e-4,
                                     maxIt = 100,
                                     alpha = 0.75
                   )){
  
 fit = RDE(y, X, Z, mBin, family, muEstim, thetaEstim, p1, p2, K1, K2, sp1, sp2,
            smooth.basis1, smooth.basis2, intercept, standardize, beta.ini,gamma.ini,
            lambdaBeta, lambdaGamma, weights.on.xz1, weights.on.xz2, rowLev,
            contX, contZ, weightFunction1, weightFunction2, optionList)
  
  tau = 0.5
  n = length(y)
  p2 = ncol(Z)
  
  if(is.character(family)){
    if(family == 'poisson'){
      family = poisson("log")
    }else if(family == "binomial"){
      family = binomial("logit")
    }
  }
  
  if(family[[1]] == "binomial"){
    if(is.null(mBin)){
      stop('Please provide the variable mBin.')
    }
  }else{
    mBin = rep(1,n)
  }

  absGam = sum(fit$gammaEstimate !=0)

  cutOff = fit$c2
  mu = fit$fitted.mu.values
  theta = fit$fitted.theta.values

  v_fun = fit$v_fun2
  nu_fun = fit$nu_fun2
  weights.xz = fit$weights.xz2
  varFun = family$variance(mu) / (mBin * theta)
  PearsonResid <- (y - mu) / sqrt(varFun)
  ind = 1:length(y)
  if(family[[1]]=="poisson"){
    yNotZero = ind[which(y!=0)]
    UTheta = 1/(2*theta) - mu
    UTheta[yNotZero] = UTheta[yNotZero] + y[yNotZero]*log(mu[yNotZero]*exp(1)/y[yNotZero])
  }else if(family[[1]]=="binomial"){
    yNotOne = ind[which(y< 1 - 1e-12)]
    yNotZero = ind[which(y>1e-12)]

    UTheta = 1/(2*theta)
    UTheta[yNotZero] = UTheta[yNotZero] + mBin[yNotZero]*y[yNotZero]*log(mu[yNotZero]/y[yNotZero])
    UTheta[yNotOne] = UTheta[yNotOne] + mBin[yNotOne]*(1-y[yNotOne])*log((1-mu[yNotOne])/(1-y[yNotOne]))
  }else {
    stop("Requested family not yet implemented.")
  }
  nu = v_fun(PearsonResid, cutOff)*UTheta
  expNu = rep(NA, n)
  expNuUTheta = rep(NA, n)
  if (family[[1]] == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmax(0, floor(mu - cutOff * (1 / theta) * mu))
    maxValue <- ceiling(mu + cutOff * (1/theta) * mu)
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)

    for (i in 1:n) {
      j <- minValue[i]:maxValue[i]
      jZero <- which(j == 0)
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      weight <- v_fun(R, cutOff)

      LLTerm = 1/(2*theta[i]) - mu[i] + j*log(mu[i]*exp(1)/j)
      LLTerm[jZero] = 1/(2*theta[i]) - mu[i]
      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUTheta[i] <- sum(weight * LLTerm^2* Probs)
    }
  } else if (family[[1]] == "binomial") {
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = mBin, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = mBin, prob = mu)
    for (i in 1:n) {
      j <- minValue[i]:maxValue[i] / mBin[i]
      Probs <-  dDBinom(j, n = mBin[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (mBin[i] * theta[i]))

      indJ = 1:length(j)
      jNotOne = indJ[which(j < 1 - 1e-12)]
      jNotZero = indJ[which(j > 1e-12)]

      LLTerm = rep(1/(2*theta[i]), length(j))
      LLTerm[jNotZero] = LLTerm[jNotZero] + mBin[i]*j[jNotZero]*log(mu[i]/j[jNotZero])
      LLTerm[jNotOne] = LLTerm[jNotOne] + mBin[i]*(1-j[jNotOne])*log((1-mu[i])/(1-j[jNotOne]))

      weight <- v_fun(R, cutOff)
      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUTheta[i] <- sum(weight * LLTerm^2* Probs)
    }
  } else {
    stop("Requested family not yet implemented.")
  }
  pseuResp =  log(theta) + (nu - expNu)/(expNuUTheta * theta)
  weightsAM =  diag(c(expNuUTheta * weights.xz * theta^2))

  EBICGam = t(pseuResp- log(theta)) %*% weightsAM %*% (pseuResp- log(theta)) +(log(n) +tau*log(p2))*absGam/n
  return(list(EBICGam=EBICGam, fit=fit))
}