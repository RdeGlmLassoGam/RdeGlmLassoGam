normConstant <- function(mu, theta){

  tol <- 1e-15

  minValue <- qpois(p = tol, lambda = mu)
  maxValue <- qpois(p = 1-tol , lambda = mu)

  j <- minValue:maxValue

  densResult <- rep(0.0, length(j))
  Ind <- which(abs(j)<=tol)
  if(length(Ind)>0) densResult[Ind] <- theta^(1/2) * exp(-theta * mu)

  dPoisson <- dpois(x = j, lambda = mu, log = TRUE) #avoids extra step of taking log of a big vector
  Ind <- which(dPoisson>log(tol)) # Since the density is in log format
  TInd <- which(j[Ind]<=tol)
  if(length(TInd)>0) Ind <- Ind[-TInd]
  dPoisson <- dPoisson[Ind]
  X <- j[Ind]

  densResult[Ind] <- exp(1/2*log(theta) +
                           dPoisson + #Normally it is a log here
                           (theta-1)*(X-mu) +
                           (theta-1)*(X*(log(mu)-log(X)))
  )

  return(sum(densResult))

}

dDPois <- function(x, mu, theta, correction = TRUE){

  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(mu)) stop("mu must be strictly positive integer.")
  if(mu<=0) stop("mu must be strictly positive integer.")
  if(theta==Inf){
    return(1*(x==mu))
  }
  if(!is.numeric(theta)) stop("theta must be strictly positive number.")
  if(theta<=0) stop("theta must be strictly positive number.")

  tol <- 1e-15

  XInd <- 1:length(x)
  densResult <- rep(0.0, length(x))

  #Treat the null values in x seperatly
  TInd <- abs(x)<=tol
  if(sum(TInd)>0){
    densResult[TInd] <- theta^(1/2) * exp(-theta * mu)
    XInd <- XInd[!TInd]
  }

  #Only consider those for which dpois is not zero,
  #set density of other observations to zero.
  minValue <- qpois(p = tol, lambda = mu)
  maxValue <- qpois(p = 1-tol , lambda = mu)
  Temp <- x[XInd]
  XInd <- XInd[(minValue<=Temp)&(Temp<=maxValue)]

  #Calculate the density if the remaining points
  dPoisson <- dpois(x = x[XInd], lambda = mu)
  TInd <- which(dPoisson>tol)
  dPoisson <- dPoisson[TInd]
  XInd <- XInd[TInd]

  densResult[XInd] <- exp(1/2*log(theta) +
                          log(dPoisson) +
                          (theta-1)*(x[XInd]-mu) +
                          (theta-1)*(x[XInd]*(log(mu)-log(x[XInd])))
                          )
  if(correction){
    NormConst <- normConstant(mu, theta)
    densResult <- densResult/NormConst
  }

  return(densResult)

}

pDPois <- function(q, mu, theta){

  if(!is.numeric(q)) stop("q must be a numeric vector.")
  if(!is.numeric(mu)) stop("mu must be strictly positive integer.")
  # if(as.integer(mu)!=mu) stop("mu must be strictly positive integer.")
  if(mu<=0) stop("mu must be strictly positive integer.")
  if(!is.numeric(theta)) stop("theta must be strictly positive number.")
  if(theta<=0) stop("theta must be strictly positive number.")

  tol <- 1e-15

  minValue <- qpois(p = tol, lambda = mu)
  maxValue <- qpois(p = 1-tol , lambda = mu)

  rangeOfValues <- minValue:maxValue
  densityOfValues <- dDPois(x = rangeOfValues, mu = mu, theta = theta)
  cumDensOfValues <- cumsum(densityOfValues)

  Result <- rep(0.0, length(q))
  Temp <- match(rangeOfValues,q)
  Ind <- which(!is.na(Temp))
  Result[Temp[Ind]] <- cumDensOfValues[Ind]


  Result[q>maxValue] <- 1

  return(Result)

}

qDPois <- function(p, mu, theta){

  if(!is.numeric(p)) stop("p must be a numeric vector.")
  if(!is.numeric(mu)) stop("mu must be strictly positive integer.")
  if(mu<=0) stop("mu must be strictly positive integer.")
  if(!is.numeric(theta)) stop("theta must be strictly positive number.")
  if(theta<=0) stop("theta must be strictly positive number.")

  factor = 3
  minValue = max(0,floor(mu - factor * sqrt((1/theta)*mu)))
  maxValue = ceiling(mu + factor * sqrt((1/theta)*mu))
  rangeOfValues <- minValue:maxValue
  Expand <- TRUE
  while(Expand){
    cdf <- pDPois(rangeOfValues, mu, theta)
    if(((min(p)<cdf[1])&(rangeOfValues[1]!=0)) | (max(p)>cdf[length(rangeOfValues)])){
      factor = factor + 1
      minValue = max(0,floor(mu - factor * sqrt((1/theta)*mu)))
      maxValue = ceiling(mu + factor * sqrt((1/theta)*mu))
      rangeOfValues <- minValue:maxValue
    }
    else{
      Expand = FALSE
    }
  }

  Quant <- rep(0.0, length(p))
  for(i in 1:length(p)){ Quant[i] <- rangeOfValues[min(which(cdf>=p[i]))]}
  if(rangeOfValues[1]==0) { Quant[p<cdf[1]] <- 0 }

  return(Quant)
}

rDPois <- function(n, mu, theta){
  if (mu <= 0) stop("mu must be greater than 0.")
  if (theta <= 0) stop("theta must be greater than 0.")
  if (n <= 0)stop("n must be a positive integer.")

  n <- ceiling(n)
  p <- runif(n)
  randObs <- qDPois(p, mu = mu, theta = theta)

  return(randObs)

}
