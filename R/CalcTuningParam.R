CalcTuningParam <- function(cValues, AMSECrit,
                               y, X, Z, m, family,
                               weights.on.xz1 = "none", weights.on.xz2=NULL,
                               weightFunction1, weightFunction2 = NULL){
  functionforapply <- function(cValue,y,Xmod,Z, m, family,
                               weights.on.xz1, weights.on.xz2,
                               weightFunction1, weightFunction2){
    model = RDE(y=y, X=Xmod, Z=Z, mBin = m,
                           family=family,
                           weights.on.xz1 = weights.on.xz1,
                           weights.on.xz2 = weights.on.xz2,
                           weightFunction1 = weightFunction1,
                           weightFunction2 = weightFunction2,
                           optionList = list(huberC = cValue,
                                             tukeyC = cValue,
                                             tol = 1e-4,
                                             maxIt = 100
                           ))
    AMSEVal = CalcAMSE(y=y, X=Xmod, Z=Z, family=family,
                       m = m,
                       weights.xz1 = model$weights.xz1, weights.xz2 = model$weights.xz2,
                       alpha = c(model$betaEstimate, model$gammaEstimate),
                       cutOff1 = model$c1, weightF1 = model$v_fun1,
                       cutOff2 = model$c2, weightF2 = model$v_fun2)

    return(c(cValue,AMSEVal[1], AMSEVal[2],AMSEVal[3]))
  }
  output = sapply(X=cValues,FUN=functionforapply,y=y,Xmod=X,Z=Z, m=m, family=family,weights.on.xz1 = weights.on.xz1, weights.on.xz2 = weights.on.xz2,
                  weightFunction1 = weightFunction1, weightFunction2 = weightFunction2)
  cBeta <- cValues[which(output[2,]>AMSECrit)[1]]
  cGamma <- cValues[which(output[3,]>AMSECrit)[1]]
  cBetaGamma <- cValues[which(output[4,]>AMSECrit)[1]]
  return(list(tunPar = c(cBeta,cGamma,cBetaGamma),AMSEBeta = output[2,],AMSEGamma = output[3,], AMSEAlpha = output[4,]))
}


