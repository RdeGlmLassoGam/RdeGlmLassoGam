#calculate RBIC for gamma
RBICGam = function(y, X, Z, mBin = NULL, family,
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
                   )){
  
  fit = RDE(y, X, Z, mBin, family, muEstim, thetaEstim, p1, p2, K1, K2, sp1, sp2,
            smooth.basis1, smooth.basis2, intercept, standardize, beta.ini,gamma.ini,
            lambdaBeta, lambdaGamma, weights.on.xz1, weights.on.xz2, rowLev,
            contX, contZ, weightFunction1, weightFunction2, optionList)
  
  n = length(y)
  p1 = ncol(fit$B1)
  p2 = ncol(fit$B2)
  
  beta = fit$betaEstimate
  gamma = fit$gammaEstimate
  
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
  
  weightF1 = fit$v_fun1
  weightF2 = fit$v_fun2
  
  weights.xz1 = fit$weights.xz1
  weights.xz2 = fit$weights.xz2
  
  c1 =fit$c1
  c2 =fit$c2
  
  nu_fun = fit$nu_fun2
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
  
  nu = weightF2(PearsonResid, c2)*UTheta
  expNu = rep(NA, n)
  expNuUTheta = rep(NA, n)
  
  if (family[[1]] == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmax(0, floor(mu - c2 * (1 / theta) * mu))
    maxValue <- ceiling(mu + c2 * (1/theta) * mu)
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
      weight <- weightF2(R, c2)
      
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
      
      weight <- weightF2(R, c2)
      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUTheta[i] <- sum(weight * LLTerm^2* Probs)
    }
  } else {
    stop("Requested family not yet implemented.")
  }
  pseuResp =  log(theta) + (nu - expNu)/pmax((expNuUTheta * theta), 10^(-12))
  weightsAM =  diag(c(expNuUTheta * weights.xz2 * theta^2))
  
  rS = fit$sD2
  
  if (family[[1]] == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmin(pmax(0, floor(mu - c1 * (1 / theta) * mu)),pmax(0, floor(mu - c2 * (1 / theta) * mu)))
    maxValue <- pmax(ceiling(mu + c1 * (1/theta) * mu),ceiling(mu + c2 * (1/theta) * mu))
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)
    
    BnumMat = 0
    BdenMat1 = 0
    BdenMat2 = 0
    BdenMat3 = 0
    GnumMat = 0
    GdenMat1 = 0
    GdenMat2 = 0
    GdenMat3 = 0
    BGnumMat12 = 0
    M12 =0
    M21 =0
    Q12 =0
    Q21 =0
    a1 = 0
    a2 = 0
    
    for (i in 1:n) {
      j <- minValue[i]:maxValue[i]
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      weight1 <- weightF1(R, c1)
      weight1[which(R == 0)] <- 1
      weight2 <- weightF2(R, c2)
      weight2[which(R == 0)] <- 1
      
      UMu <- (j-mu[i])/(mu[i]/theta[i])
      UTheta <- 1 / (2 * theta[i]) - mu[i] + j * log(exp(1) * mu[i] / j)
      if(j[1]==0){
        UTheta[1] <- 1 / (2 * theta[i]) - mu[i]
      }
      xxt = tcrossprod(fit$B1[i,])
      zzt = tcrossprod(fit$B2[i,])
      xzt = tcrossprod(fit$B1[i,],fit$B2[i,])
      zxt = tcrossprod(fit$B2[i,],fit$B1[i,])
      BnumMat = BnumMat + sum(UMu^2 * Probs)*mu[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*mu[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*mu[i]*t(fit$B1[i,])
      BdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*mu[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * mu[i] * fit$B1[i,]
      
      GnumMat = GnumMat + sum(UTheta^2 * Probs)*theta[i]^2*zzt
      GdenMat1 = GdenMat1 + sum(weight2 * UTheta^2 * Probs)*weights.xz2[i]*theta[i]^2*zzt
      GdenMat2 = GdenMat2 + sum(UTheta * Probs)*theta[i]*t(fit$B2[i,])
      GdenMat3 = GdenMat3 + sum(weight2^2 * UTheta^2 * Probs)*weights.xz2[i]^2*theta[i]^2*zzt
      a2 = a2 + sum(weight2 * UTheta * Probs) * weights.xz2[i] * theta[i] * fit$B2[i,]
      
      BGnumMat12 = BGnumMat12 + sum(UMu * UTheta * Probs)*mu[i]*theta[i]*xzt
      M12 <- M12 + sum(weight1 * UMu * UTheta * Probs)*weights.xz1[i]*mu[i]*theta[i]*xzt
      M21 <- M21 + sum(weight2 * UTheta * UMu * Probs)*weights.xz2[i]*mu[i]*theta[i]*zxt
      Q12 <- Q12 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*mu[i]*theta[i]*xzt
      Q21 <- Q21 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*mu[i]*theta[i]*zxt
    }
    K = matrix(0,p1+p2,p1+p2); Q = matrix(0,p1+p2,p1+p2); M = matrix(0,p1+p2,p1+p2)
    K[1:p1,1:p1]=BnumMat; K[1:p1,(p1+1):(p1+p2)]=BGnumMat12; K[(p1+1):(p1+p2),1:p1]=t(BGnumMat12); K[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GnumMat
    Q[1:p1,1:p1]=BdenMat3-1/n*tcrossprod(a1); Q[(p1+1):(p1+p2),1:p1]=Q21-1/n*tcrossprod(a2,a1); Q[1:p1,(p1+1):(p1+p2)]=Q12-1/n*tcrossprod(a1,a2); Q[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p1,1:p1]=BdenMat1-1/n*a1%*%BdenMat2; M[(p1+1):(p1+p2),1:p1]=M21-1/n*a2%*%BdenMat2; M[1:p1,(p1+1):(p1+p2)]=M12-1/n*a1%*%GdenMat2; M[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GdenMat1-1/n*a2%*%GdenMat2
    Binv = solve(1/n*BdenMat1-1/n^2*a1%*%BdenMat2,diag(rep(1,p1)))
    Bden = sum(diag(Binv%*%(1/n*BdenMat3-1/n^2*tcrossprod(a1))%*%t(Binv)))
    Ginv = solve(1/n*GdenMat1-1/n^2*a2%*%GdenMat2,diag(rep(1,p2)))
    Gden = sum(diag(Ginv%*%(1/n*GdenMat3-1/n^2*tcrossprod(a2))%*%t(Ginv)))
    Pinv = solve((1/n*(M[(p1+1):(p1+p2),(p1+1):(p1+p2)] + t(rS)%*%rS) ),diag(rep(1,p2)))
  }else if(family[[1]] == "binomial"){
    linpred <- as.numeric(fit$B1 %*% beta)
    mu <- as.numeric(1 / (1 + exp(-linpred)))
    theta <- as.numeric(exp(fit$B2 %*% gamma))
    der <- exp(-linpred)/(1+exp(-linpred))^2
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = mBin, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = mBin, prob = mu)
    
    BnumMat = 0
    BdenMat1 = 0
    BdenMat2 = 0
    BdenMat3 = 0
    GnumMat = 0
    GdenMat1 = 0
    GdenMat2 = 0
    GdenMat3 = 0
    BGnumMat12 = 0
    M12 =0
    M21 =0
    Q12 =0
    Q21 =0
    a1 = 0
    a2 = 0
    
    for (i in 1:n) {
      j <- minValue[i]:maxValue[i] / mBin[i]
      Probs <- dDBinom(j, n = mBin[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (mBin[i] * theta[i]))
      weight1 <- weightF1(R, c1)
      weight1[which(R == 0)] <- 1
      weight2 <- weightF2(R, c2)
      weight2[which(R == 0)] <- 1
      
      UMu <- (mBin[i] * j * theta[i]) / mu[i] + (theta[i] * mBin[i] * (1 - j)) / (mu[i] - 1)
      UTheta <- 1 / (2 * theta[i]) + (mBin[i] * j) * log(mu[i] / j) + mBin[i] * (1 - j) * log((1 - mu[i]) / (1 - j))
      
      if(j[1] <= 1e-12){
        UMu[1] <- mBin[i]*theta[i]/(mu[i]-1)
        UTheta[1] <- 1 / (2 * theta[i]) + mBin[i] * log(1 - mu[i])
      }
      if(j[length(j)] >= 1 - 1e-12){
        UMu[length(j)] <- mBin[i]*theta[i]/mu[i]
        UTheta[length(j)] <- 1 / (2 * theta[i]) + mBin[i] * log(mu[i])
      }
      xxt = tcrossprod(fit$B1[i,])
      zzt = tcrossprod(fit$B2[i,])
      xzt = tcrossprod(fit$B1[i,],fit$B2[i,])
      zxt = tcrossprod(fit$B2[i,],fit$B2[i,])
      
      BnumMat = BnumMat + sum(UMu^2 * Probs)*der[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*der[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*der[i]*t(fit$B1[i,])
      OBdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*der[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * der[i] * fit$B1[i,]
      
      GnumMat = GnumMat + sum(UTheta^2 * Probs)*theta[i]^2*zzt
      GdenMat1 = GdenMat1 + sum(weight2 * UTheta^2 * Probs)*weights.xz2[i]*theta[i]^2*zzt
      GdenMat2 = GdenMat2 + sum(UTheta * Probs)*theta[i]*t(fit$B2[i,])
      GdenMat3 = GdenMat3 + sum(weight2^2 * UTheta^2 * Probs)*weights.xz2[i]^2*theta[i]^2*zzt
      a2 = a2 + sum(weight2 * UTheta * Probs) * weights.xz2[i] * theta[i] * fit$B2[i,]
      
      BGnumMat12 = BGnumMat12 + sum(UMu * UTheta * Probs)*der[i]*theta[i]*xzt
      M12 <- M12 + sum(weight1 * UMu * UTheta * Probs)*weights.xz1[i]*der[i]*theta[i]*xzt
      M21 <- M21 + sum(weight2 * UTheta * UMu * Probs)*weights.xz2[i]*der[i]*theta[i]*zxt
      Q12 <- Q12 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*der[i]*theta[i]*xzt
      Q21 <- Q21 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*der[i]*theta[i]*zxt
    }
    K = matrix(0,p1+p2,p1+p2); Q = matrix(0,p1+p2,p1+p2); M = matrix(0,p1+p2,p1+p2)
    K[1:p1,1:p1]=BnumMat; K[1:p1,(p1+1):(p1+p2)]=BGnumMat12; K[(p1+1):(p1+p2),1:p1]=t(BGnumMat12); K[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GnumMat
    Q[1:p1,1:p1]=BdenMat3-1/n*tcrossprod(a1); Q[(p1+1):(p1+p2),1:p1]=Q21-1/n*tcrossprod(a2,a1); Q[1:p1,(p1+1):(p1+p2)]=Q12-1/n*tcrossprod(a1,a2); Q[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p1,1:p1]=BdenMat1-1/n*a1%*%BdenMat2; M[(p1+1):(p1+p2),1:p1]=M21-1/n*a2%*%BdenMat2; M[1:p1,(p1+1):(p1+p2)]=M12-1/n*a1%*%GdenMat2; M[(p1+1):(p1+p2),(p1+1):(p1+p2)]=GdenMat1-1/n*a2%*%GdenMat2
    Binv = solve(1/n*BdenMat1-1/n^2*a1%*%BdenMat2,diag(rep(1,p1)))
    Bden = sum(diag(Binv%*%(1/n*BdenMat3-1/n^2*tcrossprod(a1))%*%t(Binv)))
    Ginv = solve(1/n*GdenMat1-1/n^2*a2%*%GdenMat2,diag(rep(1,p2)))
    Gden = sum(diag(Ginv%*%(1/n*GdenMat3-1/n^2*tcrossprod(a2))%*%t(Ginv)))
    Pinv = solve((1/n*(M[(p1+1):(p1+p2),(p1+1):(p1+p2)] + t(rS)%*%rS) ),diag(rep(1,p2)))
  }else{
    stop("Family not yet implemented.")
  }
  
  RBICGam = t(pseuResp- fit$B2 %*% gamma) %*% weightsAM %*% (pseuResp- fit$B2 %*% gamma) + log(n)/n* sum(diag(Pinv %*% Q[(p1+1):(p1+p2),(p1+1):(p1+p2)]))
  return(list(RBICGam=RBICGam, fit=fit))
}


