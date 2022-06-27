CalcAMSE <- function(y, X, Z,
                           family,
                           m = NULL,
                           weights.xz1, weights.xz2,
                           alpha,
                           cutOff1, weightF1,
                           cutOff2, weightF2
                           )
{
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)

  if(family == "binomial"){
    if(is.null(m)){
      stop('Please provide the variable m.')
    }
  }

  beta <- alpha[1:p]
  gamma <- alpha[(p + 1):(p + q)]

  # Calculate the AMSE
  if (family == "poisson") {
    mu <- as.numeric(exp(X %*% beta))
    theta <- as.numeric(exp(Z %*% gamma))

    # Limit the number of considered observations considerably
    minValue <- pmin(pmax(0, floor(mu - cutOff1 * (1 / theta) * mu)),pmax(0, floor(mu - cutOff2 * (1 / theta) * mu)))
    maxValue <- pmax(ceiling(mu + cutOff1 * (1/theta) * mu),ceiling(mu + cutOff2 * (1/theta) * mu))
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
      weight1 <- weightF1(R, cutOff1)
      weight1[which(R == 0)] <- 1
      weight2 <- weightF2(R, cutOff2)
      weight2[which(R == 0)] <- 1

      UMu <- (j-mu[i])/(mu[i]/theta[i])
      UTheta <- 1 / (2 * theta[i]) - mu[i] + j * log(exp(1) * mu[i] / j)
      if(j[1]==0){
        UTheta[1] <- 1 / (2 * theta[i]) - mu[i]
      }
      xxt = tcrossprod(X[i,])
      zzt = tcrossprod(Z[i,])
      xzt = tcrossprod(X[i,],Z[i,])
      zxt = tcrossprod(Z[i,],X[i,])
      BnumMat = BnumMat + sum(UMu^2 * Probs)*mu[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*mu[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*mu[i]*t(X[i,])
      BdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*mu[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * mu[i] * X[i,]

      GnumMat = GnumMat + sum(UTheta^2 * Probs)*theta[i]^2*zzt
      GdenMat1 = GdenMat1 + sum(weight2 * UTheta^2 * Probs)*weights.xz2[i]*theta[i]^2*zzt
      GdenMat2 = GdenMat2 + sum(UTheta * Probs)*theta[i]*t(Z[i,])
      GdenMat3 = GdenMat3 + sum(weight2^2 * UTheta^2 * Probs)*weights.xz2[i]^2*theta[i]^2*zzt
      a2 = a2 + sum(weight2 * UTheta * Probs) * weights.xz2[i] * theta[i] * Z[i,]

      BGnumMat12 = BGnumMat12 + sum(UMu * UTheta * Probs)*mu[i]*theta[i]*xzt
      M12 <- M12 + sum(weight1 * UMu * UTheta * Probs)*weights.xz1[i]*mu[i]*theta[i]*xzt
      M21 <- M21 + sum(weight2 * UTheta * UMu * Probs)*weights.xz2[i]*mu[i]*theta[i]*zxt
      Q12 <- Q12 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*mu[i]*theta[i]*xzt
      Q21 <- Q21 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*mu[i]*theta[i]*zxt
    }
    K = matrix(0,p+q,p+q); Q = matrix(0,p+q,p+q); M = matrix(0,p+q,p+q)
    K[1:p,1:p]=BnumMat; K[1:p,(p+1):(p+q)]=BGnumMat12; K[(p+1):(p+q),1:p]=t(BGnumMat12); K[(p+1):(p+q),(p+1):(p+q)]=GnumMat
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1); Q[(p+1):(p+q),1:p]=Q21-1/n*tcrossprod(a2,a1); Q[1:p,(p+1):(p+q)]=Q12-1/n*tcrossprod(a1,a2); Q[(p+1):(p+q),(p+1):(p+q)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2; M[(p+1):(p+q),1:p]=M21-1/n*a2%*%BdenMat2; M[1:p,(p+1):(p+q)]=M12-1/n*a1%*%GdenMat2; M[(p+1):(p+q),(p+1):(p+q)]=GdenMat1-1/n*a2%*%GdenMat2
    Binv = solve(1/n*BdenMat1-1/n^2*a1%*%BdenMat2,diag(rep(1,p)))
    Bden = sum(diag(Binv%*%(1/n*BdenMat3-1/n^2*tcrossprod(a1))%*%t(Binv)))
    Ginv = solve(1/n*GdenMat1-1/n^2*a2%*%GdenMat2,diag(rep(1,q)))
    Gden = sum(diag(Ginv%*%(1/n*GdenMat3-1/n^2*tcrossprod(a2))%*%t(Ginv)))
    Minv = solve(1/n*M,diag(rep(1,p+q)))
    AMSEBeta = n*sum(diag(solve(BnumMat,diag(rep(1,p)))))/Bden
    AMSEGamma = n*sum(diag(solve(GnumMat,diag(rep(1,q)))))/Gden
    AMSEBetaGamma = sum(diag(solve(1/n*K,diag(rep(1,p+q)))))/sum(diag(Minv %*% (1/n*Q)%*%t(Minv)))
  }else if(family == "binomial"){
    linpred <- as.numeric(X %*% beta)
    mu <- as.numeric(1 / (1 + exp(-linpred)))
    theta <- as.numeric(exp(Z %*% gamma))
    der <- exp(-linpred)/(1+exp(-linpred))^2
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = m, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)

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
      j <- minValue[i]:maxValue[i] / m[i]
      Probs <- dDBinom(j, n = m[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m[i] * theta[i]))
      weight1 <- weightF1(R, cutOff1)
      weight1[which(R == 0)] <- 1
      weight2 <- weightF2(R, cutOff2)
      weight2[which(R == 0)] <- 1

      UMu <- (m[i] * j * theta[i]) / mu[i] + (theta[i] * m[i] * (1 - j)) / (mu[i] - 1)
      UTheta <- 1 / (2 * theta[i]) + (m[i] * j) * log(mu[i] / j) + m[i] * (1 - j) * log((1 - mu[i]) / (1 - j))

      if(j[1] <= 1e-12){
        UMu[1] <- m[i]*theta[i]/(mu[i]-1)
        UTheta[1] <- 1 / (2 * theta[i]) + m[i] * log(1 - mu[i])
      }
      if(j[length(j)] >= 1 - 1e-12){
        UMu[length(j)] <- m[i]*theta[i]/mu[i]
        UTheta[length(j)] <- 1 / (2 * theta[i]) + m[i] * log(mu[i])
      }
      xxt = tcrossprod(X[i,])
      zzt = tcrossprod(Z[i,])
      xzt = tcrossprod(X[i,],Z[i,])
      zxt = tcrossprod(Z[i,],X[i,])

      BnumMat = BnumMat + sum(UMu^2 * Probs)*der[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*der[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*der[i]*t(X[i,])
      BdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*der[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * der[i] * X[i,]

      GnumMat = GnumMat + sum(UTheta^2 * Probs)*theta[i]^2*zzt
      GdenMat1 = GdenMat1 + sum(weight2 * UTheta^2 * Probs)*weights.xz2[i]*theta[i]^2*zzt
      GdenMat2 = GdenMat2 + sum(UTheta * Probs)*theta[i]*t(Z[i,])
      GdenMat3 = GdenMat3 + sum(weight2^2 * UTheta^2 * Probs)*weights.xz2[i]^2*theta[i]^2*zzt
      a2 = a2 + sum(weight2 * UTheta * Probs) * weights.xz2[i] * theta[i] * Z[i,]

      BGnumMat12 = BGnumMat12 + sum(UMu * UTheta * Probs)*der[i]*theta[i]*xzt
      M12 <- M12 + sum(weight1 * UMu * UTheta * Probs)*weights.xz1[i]*der[i]*theta[i]*xzt
      M21 <- M21 + sum(weight2 * UTheta * UMu * Probs)*weights.xz2[i]*der[i]*theta[i]*zxt
      Q12 <- Q12 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*der[i]*theta[i]*xzt
      Q21 <- Q21 + sum(weight1 * weight2 * UMu * UTheta * Probs)*weights.xz1[i]*weights.xz2[i]*der[i]*theta[i]*zxt
    }
    K = matrix(0,p+q,p+q); Q = matrix(0,p+q,p+q); M = matrix(0,p+q,p+q)
    K[1:p,1:p]=BnumMat; K[1:p,(p+1):(p+q)]=BGnumMat12; K[(p+1):(p+q),1:p]=t(BGnumMat12); K[(p+1):(p+q),(p+1):(p+q)]=GnumMat
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1); Q[(p+1):(p+q),1:p]=Q21-1/n*tcrossprod(a2,a1); Q[1:p,(p+1):(p+q)]=Q12-1/n*tcrossprod(a1,a2); Q[(p+1):(p+q),(p+1):(p+q)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2; M[(p+1):(p+q),1:p]=M21-1/n*a2%*%BdenMat2; M[1:p,(p+1):(p+q)]=M12-1/n*a1%*%GdenMat2; M[(p+1):(p+q),(p+1):(p+q)]=GdenMat1-1/n*a2%*%GdenMat2
    Binv = solve(1/n*BdenMat1-1/n^2*a1%*%BdenMat2,diag(rep(1,p)))
    Bden = sum(diag(Binv%*%(1/n*BdenMat3-1/n^2*tcrossprod(a1))%*%t(Binv)))
    Ginv = solve(1/n*GdenMat1-1/n^2*a2%*%GdenMat2,diag(rep(1,q)))
    Gden = sum(diag(Ginv%*%(1/n*GdenMat3-1/n^2*tcrossprod(a2))%*%t(Ginv)))
    Minv = solve(1/n*M,diag(rep(1,p+q)))
    AMSEBeta = n*sum(diag(solve(BnumMat,diag(rep(1,p)))))/Bden
    AMSEGamma = n*sum(diag(solve(GnumMat,diag(rep(1,q)))))/Gden
    AMSEBetaGamma = sum(diag(solve(1/n*K,diag(rep(1,p+q)))))/sum(diag(Minv %*% (1/n*Q)%*%t(Minv)))
  }else{
    stop("Family not yet implemented.")
  }
  return(c(AMSEBeta,AMSEGamma,AMSEBetaGamma ))
}
