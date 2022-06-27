DispersionTest <- function(y, X, Z,
                           m = NULL,
                           alphaHAlt,
                           alphaHNull,
                           nzeroBeta, #q1
                           nzeroGamma, #q2
                           family,
                           weightFunction1,
                           weightFunction2 = NULL,
                           weights.on.xz1,
                           weights.on.xz2,
                           optionList,
                           M, Q) {
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)

  if(family == "binomial"){
    if(is.null(m)){
      stop('Please provide the variable m.')
    }
  }

  betaHAlt <- alphaHAlt[1:p]
  gammaHAlt <- alphaHAlt[-c(1:p)]
  betaHNull <- alphaHNull[1:p]
  gammaHNull <- alphaHNull[-c(1:p)]


  # Define the weight functions to diminish the effect of large residuals.
  if (weightFunction1 == "Huber") {
    cutOff1 <- optionList$huberC
    weightF1 <- function(Arg1) (HuberResid(Arg1, cutOff1) / Arg1)
  } else if (weightFunction1 == "Tukey") {
    cutOff1 <- optionList$tukeyC
    weightF1 <- function(Arg1) TukeyResid(Arg1, cutOff1) / Arg1
  }
  if(is.null(weightFunction2)){
      cutOff2 <- cutOff1
      weightF2 <- weightF1
      weightFunction2 = weightFunction1
  }else if (weightFunction2 == "Huber") {
    cutOff2 <- optionList$huberC
    weightF2 <- function(Arg2) (HuberResid(Arg2, cutOff2) / Arg2)
  } else if (weightFunction2 == "Tukey") {
    cutOff2 <- optionList$tukeyC
    weightF2 <- function(Arg2) TukeyResid(Arg2, cutOff2) / Arg2
  }

  # Determine the mean and the dispersion under the null hypothesis and the alternative hypothesis.
  if (family == "poisson" ) {
    muHAlt <- as.numeric(exp(X %*% betaHAlt))
    muHNull <- as.numeric(exp(X %*% betaHNull))
    thetaHAlt <- as.numeric(exp(Z %*% gammaHAlt))
    thetaHNull <- as.numeric(exp(Z %*% gammaHNull))
  } else if (family == "binomial") {
    linpredHAlt <- as.numeric(X %*% betaHAlt)
    muHAlt <- as.numeric(exp(linpredHAlt) / (1 + exp(linpredHAlt)))
    linpredHNull <- as.numeric(X %*% betaHNull)
    muHNull <- as.numeric(exp(linpredHNull) / (1 + exp(linpredHNull)))
    thetaHAlt <- as.numeric(exp(Z %*% gammaHAlt))
    thetaHNull <- as.numeric(exp(Z %*% gammaHNull))
  } else{
    stop("Family not yet implemented.")
  }



  # Determine the test value.
  LambdaQM_beta <- 0
  LambdaQM_gamma <- 0
  for (i in 1:n) {
    try(
      Contrib_mu_1 <- integrate(integrand_mu_1, y = y[i], theta = thetaHAlt[i], family = family,
                            weights.on.xz1 = weights.on.xz1[i], weightF = weightF1, weightFunction = weightFunction1,
                            m = m[i],
                            lower = muHAlt[i],
                            upper = muHNull[i],
                            subdivisions = 1000),
      silent = TRUE)
    try(
      Contrib_mu_2 <- integrate(integrand_mu_2, y = y[i], theta = thetaHAlt[i], family = family, weights.on.xz1 = weights.on.xz1[i],
                            weightF = weightF1, weightFunction = weightFunction1,
                            m = m[i],cutOff = cutOff1,
                            lower = muHAlt[i],
                            upper = muHNull[i],
                            subdivisions = 1000),
      silent = TRUE)
    try(
      Contrib1 <- integrate(integrand1, y = y[i], mu = muHAlt[i], family = family,
                            weights.on.xz2 = weights.on.xz2[i], weightF = weightF2, weightFunction = weightFunction2,
                            m = m[i],
                            lower = thetaHAlt[i],
                            upper = thetaHNull[i],
                            subdivisions = 1000),
      silent = TRUE)
    try(
      Contrib2 <- integrate(integrand2, y = y[i], mu = muHAlt[i], family = family, weights.on.xz2 = weights.on.xz2[i],
                            weightF = weightF2, weightFunction = weightFunction2,
                            m = m[i],cutOff = cutOff2,
                          lower = thetaHAlt[i],
                          upper = thetaHNull[i],
                          subdivisions = 1000),
      silent = TRUE)

    if (is.numeric(Contrib1$value) & is.numeric(Contrib2$value)) {
      LambdaQM_gamma <- LambdaQM_gamma + (Contrib1$value - Contrib2$value)
      LambdaQM_beta <- LambdaQM_beta +  (Contrib_mu_1$value - Contrib_mu_2$value)
    }
    else{
      LambdaQM <- -1
      break()
    }
  }
  LambdaQM <- -2*(LambdaQM_beta + LambdaQM_gamma)

  M_beta <- M[1:p, 1:p]
  Q_beta <- Q[1:p, 1:p]
  M_gamma <- M[(p+1):(p +q) , (p+1):(p +q)]
  Q_gamma <- Q[(p+1):(p +q) , (p+1):(p +q)]

  eigValue <- c()

  # Find the eigenvalues for the beta part.
  if (nzeroBeta > 0) {
    M_beta_Invers <- solve(M_beta)
    M_beta_Reduc_Invers <- matrix(0, nrow = p, ncol = p)
    if(p > 1){
      M_beta_Reduc_Invers[1:(p-nzeroBeta),1:(p-nzeroBeta)] <- as.matrix(solve(M[1:(p-nzeroBeta),1:(p-nzeroBeta), drop = FALSE]))
    }
    ModelSelecMatrix <- Q_beta %*% (M_beta_Invers - M_beta_Reduc_Invers)
    ModelSelecMatrix[ModelSelecMatrix <= 10 ^ -5] <- 0
    eigValue <- c(eigValue, Re(eigen(ModelSelecMatrix)$value[1:nzeroBeta]))
  }

  # Find the eigenvalues for the gamma part.
  if (nzeroGamma > 0) {
    M_gamma_Invers <- solve(M_gamma)
    M_gamma_Reduc_Invers <- matrix(0, nrow = q, ncol = q)
    if(q > 1){
      M_gamma_Reduc_Invers[1:(q-nzeroGamma),1:(q-nzeroGamma)] <- as.matrix(solve(M[1:(q-nzeroGamma),1:(q-nzeroGamma), drop = FALSE]))
    }
    ModelSelecMatrix <- Q_gamma %*% (M_gamma_Invers - M_gamma_Reduc_Invers)
    ModelSelecMatrix[ModelSelecMatrix <= 10 ^ -5] <- 0
    eigValue <- c(eigValue, Re(eigen(ModelSelecMatrix)$value[1:nzeroGamma]))
  }

  # Determine the p-value.

  if (is.numeric(eigValue)) {
    pValue <- davies(q = LambdaQM, lambda =  eigValue)$Qq
  } else {
    pValue = NA
  }

  return(list(DispTestValue = LambdaQM,
              DispTestEigValue = eigValue,
              DispTestpValue = pValue))
}


# Contribution to the integral of the main term.
integrand_mu_1 <- function(mu, y, theta, family, weights.on.xz1, weightF, weightFunction, m = NULL){
  Response <- rep(0.0, length(mu))
  if (family == "poisson"){
    for (i in 1:length(mu)) {
      R <- (y - mu[i]) / (sqrt(mu[i]/theta))
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      LLTermBeta <- (y - mu[i]) / (mu[i] / theta)
      Response[i] <- (weight * LLTermBeta * weights.on.xz1)
    }
    return(Response)
  }else if(family == "binomial"){
    for (i in 1:length(mu)) {
      R <- (y - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m * theta ))
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      LLTermBeta <- (y - mu[i]) / (mu[i]*(1-mu[i])/(theta*m))
      Response[i] <- (weight * LLTermBeta * weights.on.xz1)
    }
    return(Response)
  }
  stop("Family not yet implemented.")
}

# Contribution to the integral of the expectancy term.
integrand_mu_2 <- function(mu, y, theta, family, weights.on.xz1, weightF, weightFunction, m = NULL, cutOff){
  Response <- rep(0.0, length(mu))
  if (family == "poisson") {
    for (i in 1:length(mu)) {
      minValue <- max(0, floor(mu[i] - cutOff * (1/theta) * mu[i]))
      maxValue <- ceiling(mu[i] + cutOff * (1/theta) * mu[i])
      tol <- 1e-12
      minCut <- qpois(p = tol, lambda = mu[i])
      maxCut <- qpois(p = 1 - tol , lambda = mu[i])
      minValue <- max(minValue, minCut)
      maxValue <- min(maxValue, maxCut)
      j <- minValue:maxValue
      Probs <- dDPois(j, mu = mu[i], theta = max(theta,0.0001), correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta))
      LLTerm <- R / (sqrt(mu[i]/theta))
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      ExpTermBeta <- sum(weight * LLTerm * Probs) * weights.on.xz1
      Response[i] <- as.numeric(ExpTermBeta)
    }
    return(Response)
  }else if(family == "binomial"){
    for (i in 1:length(mu)) {
      tol <- 1e-12
      minValue <- qbinom(p = tol, size = m, prob = mu[i])
      maxValue <- qbinom(p = 1 - tol, size = m, prob = mu[i])
      j <- minValue:maxValue / m
      Probs <- dDBinom(j, n = m, mu = mu[i], theta = max(theta,0.0001), correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m * theta ))
      LLTerm <- (j - mu[i]) / (mu[i]*(1-mu[i])/(theta*m))
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      ExpTermBeta <- sum(weight * LLTerm * Probs) * weights.on.xz1
      Response[i] <- as.numeric(ExpTermBeta)
    }
    return(Response)
  }
  stop("Family not yet implemented.")
}

# Contribution to the integral of the main term.
integrand1 <- function(theta, y, mu, family, weights.on.xz2, weightF, weightFunction, m = NULL){
  Response <- rep(0.0, length(theta))
  if (family == "poisson"){
    for (i in 1:length(theta)) {
      R <- (y - mu) / (sqrt(mu/theta[i]))
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      LLTermGamma <- (1/(2*theta[i]) - mu + y * log(exp(1)*mu/y))
      if (y == 0) LLTermGamma <- (1/(2*theta[i]) - mu)
      Response[i] <- (weight * LLTermGamma * weights.on.xz2)
    }
    return(Response)
  } else if (family == "binomial"){
    for (i in 1:length(theta)) {
      dLogTheta <- 1 / (2 * theta[i]) + (m * y) * log(mu / y) + m * (1 - y) * log((1 - mu) / (1 - y))
      #Adapt for boundary values
      if (y <= 1e-12) {
        dLogTheta <- 1/(2*theta[i]) + m * log(1 - mu)
      }
      if (y >= 1 - 1e-12) {
        dLogTheta <- 1/(2*theta[i]) + m * log(mu)
      }
      PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta[i] ))
      weight <- weightF(PearsonResid)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(PearsonResid == 0)] <- 1
      Response[i] <- (weight * dLogTheta * weights.on.xz2)
    }
    return(Response)
  }
  stop("Family not yet implemented.")
}

# Contribution to the integral of the expectancy term.
integrand2 <- function(theta, y, mu, family, weights.on.xz2, weightF, weightFunction, m = NULL, cutOff){
  Response <- rep(0.0, length(theta))
  if (family == "poisson") {
    for (i in 1:length(theta)) {
      minValue <- max(0, floor(mu - cutOff * (1/theta[i]) * mu))
      maxValue <- ceiling(mu + cutOff * (1/theta[i]) * mu)
      tol <- 1e-12
      minCut <- qpois(p = tol, lambda = mu)
      maxCut <- qpois(p = 1 - tol , lambda = mu)
      minValue <- max(minValue, minCut)
      maxValue <- min(maxValue, maxCut)
      j <- minValue:maxValue
      Probs <- dDPois(j, mu = mu, theta = max(theta[i],0.0001), correction = FALSE)
      R <- (j - mu) / (sqrt(mu / theta[i]))
      LLTerm <- 1/(2*theta[i]) - mu + j * log(exp(1)*mu/j)
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      if (j[1] == 0) LLTerm[1] <- 1/(2*theta[i]) - mu
      ExpTermGamma <- sum(weight * LLTerm * Probs) * weights.on.xz2
      Response[i] <- as.numeric(ExpTermGamma)
    }
    return(Response)
  } else if(family == "binomial"){
    for (i in 1:length(theta)) {
      tol <- 1e-12
      minValue <- qbinom(p = tol, size = m, prob = mu)
      maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)
      j <- minValue:maxValue / m
      Probs <- dDBinom(j, n = m, mu = mu, theta = max(theta[i],0.0001), correction = FALSE)
      R <- (j - mu) / sqrt(mu * (1 - mu) / (m * theta[i] ))
      LLTerm <- 1/(2*theta[i]) + (m * j) * log(mu/j) + m * (1 - j) * log((1 - mu) / (1 - j))
      if (j[1] <= tol) LLTerm[1] <- 1/(2*theta[i]) + m * log(1 - mu)
      if (j[length(j)] >= 1 - tol) LLTerm[length(j)] <- 1/(2*theta[i]) + m * log(mu)
      weight <- weightF(R)
      if (weightFunction == "Huber") {
        weight <- weight^2
      }
      weight[which(R == 0)] <- 1
      ExpTermGamma = sum(weight * LLTerm * Probs) * weights.on.xz2
      Response[i] <- as.numeric(ExpTermGamma)
    }
    return(Response)
  }
  stop("Family not yet implemented.")
}
