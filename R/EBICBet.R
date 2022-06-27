#calculate EBIC for beta
EBICBet = function(y, X, Z, mBin = NULL, family,
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
  p1 = ncol(X)
  absBet = sum(fit$betaEstimate !=0)
  
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

  mu = fit$fitted.mu.values
  theta = fit$fitted.theta.values

  v_fun = fit$v_fun1
  nu_fun = fit$nu_fun1
  weights.xz = fit$weights.xz1
  
  c =fit$c1

  varFun = family$variance(mu) / (mBin * theta)
  PearsonResid <- (y - mu) / sqrt(varFun)
  UMu <- PearsonResid / sqrt(varFun)
  eta = family$linkfun(mu)

  nu = v_fun(PearsonResid, c)*UMu

  expNu = rep(NA, n)
  expNuUMu = rep(NA, n)

  if (family[[1]] == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmax(0, floor(mu - c * (1 / theta) * mu))
    maxValue <- ceiling(mu + c * (1/theta) * mu)
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)

    for (i in 1:n) {
      j <- minValue[i]:maxValue[i]
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      LLTerm <- R / (sqrt(mu[i]/theta[i]))
      weight <- v_fun(R, c)

      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUMu[i] <- sum(weight * LLTerm^2 * Probs)
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
      LLTerm <- R / sqrt(mu[i] * (1 - mu[i]) / (mBin[i] * theta[i]))
      weight <- v_fun(R, c)

      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUMu[i] <- sum(weight * LLTerm^2 * Probs)
    }
  } else {
    print("stop")
    stop("Requested family not yet implemented.")
  }
  pseuResp = eta +(nu - expNu)/(expNuUMu*family$mu.eta(eta))
  weightsAM = diag(c(expNuUMu * weights.xz * family$mu.eta(eta)^2))

  EBICBet = t(pseuResp- eta) %*% weightsAM %*% (pseuResp- eta) +(log(n) +tau*log(p1))*absBet/n
  return(list(EBICBet=EBICBet, fit=fit))
}


