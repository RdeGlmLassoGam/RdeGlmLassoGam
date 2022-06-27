normConstant_Binom <- function(n, mu, theta){

  tol <- 1e-15

  minValue <- qbinom(p = tol, size = n, prob = mu)
  maxValue <- qbinom(p = 1 - tol, size = n, prob = mu)

  j <- minValue:maxValue

  densResult <- rep(0.0, length(j))
  Ind <- which(abs(j) <= tol)
  if (length(Ind) > 0) densResult[Ind] <- theta ^ (1/2) * ((1 - mu) ^ n) ^ theta
  Ind <- which(abs(j)/n >= 1 - tol)
  if (length(Ind) > 0) densResult[Ind] <- theta ^ (1/2) * (mu ^ n) ^ theta

  dBinom <- dbinom(x = j, size = n, prob = mu, log = TRUE) #avoids extra step of taking log of a big vector
  Ind <- which(dBinom > log(tol)) # Since the density is in log format
  TInd <- which(j[Ind] <= tol | j[Ind]/n >= 1 - tol)
  if (length(TInd) > 0) Ind <- Ind[-TInd]
  dBinom <- dBinom[Ind]
  X <- j[Ind] / n

  densResult[Ind] <- exp(1/2*log(theta) +
                           dBinom + # already in log format
                           (n * X) * theta * log(mu / X) +
                           (n * X) * log(X / mu) +
                           theta * n * (1 - X) * log((1 - mu) / (1 - X)) +
                           n * (1 - X) * log((1 - X) / (1 - mu))
                         )

  return(sum(densResult))

}

dDBinom <- function(x, n, mu, theta, correction = TRUE){
  # print(x)
  if (!is.numeric(x)) stop("x must be a numeric vector.")
  if (sum(x < 0 | x > 1)) stop("x must be between zero and one.")
  if (!is.numeric(mu)) stop("mu must be strictly positive.")
  if (mu <= 0||mu>1) stop("mu must be between zero and one.")
  if(theta==Inf){
    return(1*(x==mu))
  }
  if (!is.numeric(theta)) stop("theta must be strictly positive number.")
  if (theta <= 0) stop("theta must be strictly positive number.")
  if (!is.numeric(n)) stop("n must be strictly positive integer.")
  if (as.integer(n) != n) stop("n must be strictly positive integer.")
  if (n < 0) stop("n must be strictly positive integer.")

  tol <- 1e-15

  XInd <- 1:length(x)
  densResult <- rep(0.0, length(x))

  #Treat the null values in x seperatly
  TInd <- abs(x) <= tol
  if (sum(TInd) > 0) {
    densResult[TInd] <- theta^(1/2) * ((1 - mu)^n)^theta
   }
  TInd2 <- abs(x) >= 1 - tol
  if (sum(TInd2) > 0) {
    densResult[TInd2] <- theta ^ (1/2) * (mu ^ n) ^ theta
  }
  XInd <- XInd[!(TInd + TInd2)]

  #Only consider those for which dbinom is not zero,
  #set density of other observations to zero.
  minValue <- qbinom(p = tol, size = n, prob = mu)
  maxValue <- qbinom(p = 1 - tol, size = n, prob = mu)
  Temp <- x[XInd] * n
  XInd <- XInd[(minValue <= Temp) & (Temp <= maxValue)]

  #Calculate the density if the remaining points
  dBinom <- dbinom(x = x[XInd] * n, size = n, prob = mu)
  Ind <- which(dBinom > tol)
  dBinom <- dBinom[Ind]
  XInd <- XInd[Ind]
  y <- x[XInd]

  densResult[XInd] <- exp(1/2*log(theta) +
                            log(dBinom) +
                            (n * y) * theta * log(mu / y) +
                            (n * y) * log(y / mu) +
                            theta * n * (1 - y) * log((1 - mu) / (1 - y)) +
                            n * (1 - y) * log((1 - y) / (1 - mu))
                          )
  if (correction) {
    NormConst <- normConstant_Binom(n, mu, theta)
    densResult <- densResult/NormConst
  }

  return(densResult)

}

pDBinom <- function(q, n, mu, theta){

  if (!is.numeric(q)) stop("q must be a numeric vector.")
  if (sum(q < 0 | q > 1)) stop("q must be between zero and one.")
  if (!is.numeric(mu)) stop("mu must be strictly positive.")
  if (mu <= 0||mu>1) stop("mu must be between zero and one.")
  if (!is.numeric(theta)) stop("theta must be strictly positive number.")
  if (theta <= 0) stop("theta must be strictly positive number.")
  if (!is.numeric(n)) stop("n must be strictly positive integer.")
  if (as.integer(n) != n) stop("n must be strictly positive integer.")
  if (n < 0) stop("n must be strictly positive integer.")

  tol <- 1e-15

  minValue <- qbinom(p = tol, size = n, prob = mu)
  maxValue <- qbinom(p = 1 - tol, size = n, prob = mu)

  rangeOfValues <- minValue:maxValue
  densityOfValues <- dDBinom(x = rangeOfValues / n, n = n, mu = mu, theta = theta)
  cumDensOfValues <- cumsum(densityOfValues)

  Result <- rep(0.0, length(q))
  Temp <- match(rangeOfValues / n, q)
  Ind <- which(!is.na(Temp))
  Result[Temp[Ind]] <- cumDensOfValues[Ind]


  Result[q > (maxValue / n)] <- 1

  return(Result)

}

qDBinom <- function(p, n, mu, theta){

  if (!is.numeric(p)) stop("p must be a numeric vector.")
  if (!is.numeric(mu)) stop("mu must be strictly positive.")
  if (mu <= 0||mu>1) stop("mu must be between zero and one.")
  if (!is.numeric(theta)) stop("theta must be strictly positive number.")
  if (theta <= 0) stop("theta must be strictly positive number.")
  if (!is.numeric(n)) stop("n must be strictly positive integer.")
  if (as.integer(n) != n) stop("n must be strictly positive integer.")
  if (n < 0) stop("n must be strictly positive integer.")

  minValue = 0
  maxValue = n
  rangeOfValues <- minValue:maxValue
  cdf <- pDBinom(rangeOfValues / n, n, mu, theta)

  Quant <- rep(0.0, length(p))
  for (i in 1:length(p)) { Quant[i] <- rangeOfValues[min(which(cdf >= p[i]))]}
  if (rangeOfValues[1] == 0) { Quant[p < cdf[1]] <- 0 }
  Quant <- Quant / n

  return(Quant)
}

rDBinom <- function(m, n, mu, theta){
  if (mu < 0) stop("mu must be greater than 0.")
  if (theta <= 0) stop("theta must be greater than 0.")
  if (n <= 0) stop("n must be a positive integer.")
  if (m <= 0) stop("m must be a positive integer.")

  n <- ceiling(n)
  m <- ceiling(m)
  p <- runif(m)
  randObs <- qDBinom(p, n = n, mu = mu, theta = theta)

  return(randObs)

}
