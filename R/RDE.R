RDE <- function(y, X, Z, mBin = NULL, family,
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
                )
)
{
  #################
  # General setup #
  #################
  n <- length(y)
  
  X <- matrix(X, nrow = n)
  nxs <- ncol(X)
  
  Z <- matrix(Z, nrow = n)
  nzs <- ncol(Z)
  
  if( (n < nxs) && (muEstim != "pen") ){
    stop('Please choose a penalized method for the mean model.')
  }
  
  if((n < nzs) && (thetaEstim != "pen")){
    stop('Please choose a penalized method for the dispersion model.')
  }
  
  if( (lambdaBeta<=0) && (muEstim == "pen") ){
    stop('Please choose a strictly positive penalty parameter for the mean model.')
  }
  
  if( (lambdaGamma<=0) && (thetaEstim == "pen") ){
    stop('Please choose a strictly positive penalty parameter for the dispersion model.')
  }
  
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
  
  if(is.null(p2)){p2 = p1}
  if(is.null(K2)){K2 = K1}
  
  if(muEstim == 'pen'){
    if(intercept){
      if(!all(X[,1]==1)){
        X <- cbind(rep(1,n),X)
        contX = contX +1
      }
    }else{
      if(all(X[,1]==1)){
        X <- X[, -1, drop=FALSE]
        contX = contX -1
      }
    }
  }
  
  if(muEstim == 'glm'){
    if(!all(X[,1]==1)){
      X <- cbind(rep(1,n),X)
      contX = contX +1
    }
  }
  
  if(thetaEstim == 'pen'){
    if(intercept){
      if(!all(Z[,1]==1)){
        Z <- cbind(rep(1,n),Z)
        contZ = contZ +1
      }
    }else{
      if(all(Z[,1]==1)){
        Z <- Z[, -1, drop=FALSE]
        contZ = contZ -1
      }
    }
  }
  
  if(thetaEstim == 'glm'){
    if(!all(Z[,1]==1)){
      Z <- cbind(rep(1,n),Z)
      contZ = contZ +1
    }
  }
  
  if(!standardize && !is.null(contX)){
    if(length(contX)==1){
      MED <- median(X[,contX])
      MAD <- mad(X[,contX])
      
      if(MAD==0){  X[,contX] <- (X[,contX]-MED)/0.01}
      else{  X[,contX] <- (X[,contX]-MED)/MAD}
    }else{
      MED <- apply(X[,contX],2,median)
      MAD <- apply(X[,contX],2,mad)
      
      for (i in 1:length(contX)){
        if(MAD[i]==0){  X[,contX[i]] <- (X[,contX[i]]-MED[i])/0.01}
        else{  X[,contX[i]] <- (X[,contX[i]]-MED[i])/MAD[i]}
      }
    }
  }
  
  if(!standardize && !is.null(contZ)){
    if(length(contZ)==1){
      MED <- median(Z[,contZ])
      MAD <- mad(Z[,contZ])
      
      if(MAD==0){  Z[,contZ] <- (Z[,contZ]-MED)/0.01}
      else{  Z[,contZ] <- (Z[,contZ]-MED)/MAD}
    }else{
      MED <- apply(Z[,contZ],2,median)
      MAD <- apply(Z[,contZ],2,mad)
      
      for (i in 1:length(contZ)){
        if(MAD[i]==0){  Z[,contZ[i]] <- (Z[,contZ[i]]-MED[i])/0.01}
        else{  Z[,contZ[i]] <- (Z[,contZ[i]]-MED[i])/MAD[i]}
      }
    }
  }
  
  Xmod <- if(intercept){X[, -1, drop=FALSE]}else{X}
  Zmod <- if(intercept){Z[, -1, drop=FALSE]}else{Z}
  
  p = ncol(X)
  q = ncol(Z)
  
  Zstart = Zmod
  Xstart = Xmod
  ystart = y
  ###########
  # Weights #
  ###########
  if(is.numeric(weights.on.xz1)){
    if(length(weights.on.xz1) == n & !(any(weights.on.xz1 < 0)) & !(any(weights.on.xz1 > 1))){
      weights.xz1 <- weights.on.xz1
    }else{
      stop("All 'weights.on.xz1' must be non-negative and smaller or equal to 1")
    }
  }else{
    if(nxs<=n){
      if(weights.on.xz1 == "none"){
        weights.xz1 = rep(1,n)
      }else{
        if(weights.on.xz1 == "covMcd"){
          if(!("alpha" %in% names(optionList))){
            optionList[["alpha"]] <- 0.75
          }
          if(!is.null(contX)){
            if(!is.null(contZ)){
              weights.xz1 <-  covMcd(X[,contX], alpha = optionList$alpha)$mcd.wt*covMcd(Z[,contZ], alpha = optionList$alpha)$mcd.wt
            }else{
              weights.xz1 <-  covMcd(X[,contX], alpha = optionList$alpha)$mcd.wt
            }
          }else{
            if(!is.null(contZ)){
              weights.xz1 <- covMcd(Z[,contZ], alpha = optionList$alpha)$mcd.wt
            }else{
              stop("When 'weights.on.xz1' is equal to 'covMcd', 'contX' or 'contZ' should be specified.")
            }
          }
        }else{
          stop("This value for 'weights.on.xz1' is undefined.")
        }
      }
    }else{#(nxs>n)
      weights.xz1 = rep(1,n)
      if(!is.null(contX)){
        DetDev = DDC(X[,contX],DDCpars=list(fastDDC=FALSE))
        if(rowLev==TRUE){
          X[,contX] = DetDev$Xest
          Xmod <- if(intercept){X[, -1, drop=FALSE]}else{X}
          IndexLev = DetDev$indrows
          if(length(IndexLev)==0){
            Xstart = Xmod
            ystart = y
          }else{
            Xstart = Xmod[-IndexLev,]
            ystart = y[-IndexLev]
            weights.xz1[IndexLev] = 0
          }
        }else{
          IndexLev = unique((DDC(X[,contX],DDCpars=list(fastDDC=FALSE))$indcells)%%n)
          IndexLev = replace(IndexLev, IndexLev==0, n)
          if(length(IndexLev)==0){
            Xstart = Xmod
            ystart = y
          }else{
            Xstart = Xmod[-IndexLev,]
            ystart = y[-IndexLev]
            weights.xz1[IndexLev] = 0
          }
        }
      }
    }
  }
  
  if(is.null(weights.on.xz2)){
    weights.xz2 = weights.xz1
  }else{
    if(is.numeric(weights.on.xz2)){
      if(length(weights.on.xz2) == n & !(any(weights.on.xz2 < 0)) & !(any(weights.on.xz2 > 1))){
        weights.xz2 <- weights.on.xz2
      }else{
        stop("All 'weights.on.xz2' must be non-negative and smaller or equal to 1")
      }
    }else{
      if(weights.on.xz2 == "none"){
        weights.xz2 = rep(1,n)
      }else if(nzs<=n){
        if(weights.on.xz2== "covMcd"){
          if(!("alpha" %in% names(optionList))){
            optionList[["alpha"]] <- 0.75
          }
          if(!is.null(contX)){
            if(!is.null(contZ)){
              weights.xz2 <-  covMcd(X[,contX], alpha = optionList$alpha)$mcd.wt*covMcd(Z[,contZ], alpha = optionList$alpha)$mcd.wt
            }else{
              weights.xz2 <-  covMcd(X[,contX], alpha = optionList$alpha)$mcd.wt
            }
          }else{
            if(!is.null(contZ)){
              weights.xz2 <- covMcd(Z[,contZ], alpha = optionList$alpha)$mcd.wt
            }else{
              stop("When 'weights.on.xz2' is equal to 'covMcd', 'contX' or 'contZ' should be specified.")
            }
          }
        }else{
          stop("This value for 'weights.on.xz2' is undefined.")
        }
      }else{#(nzs>n)
        weights.xz2 = rep(1,n)
        if(!is.null(contZ)){
          DetDev = DDC(Z[,contZ],DDCpars=list(fastDDC=FALSE))
          if(rowLev==TRUE){
            Z[,contZ] = DetDev$Xest
            Zmod <- if(intercept){Z[, -1, drop=FALSE]}else{Z}
            IndexLev = DetDev$indrows
            if(length(IndexLev)==0){
              Zstart = Zmod
              ystart = y
            }else{
              Zstart = Zmod[-IndexLev,]
              ystart = y[-IndexLev]
              weights.xz2[IndexLev] = 0
            }
          }else{
            IndexLev = unique((DDC(Z[,contZ],DDCpars=list(fastDDC=FALSE))$indcells)%%n)
            IndexLev = replace(IndexLev, IndexLev==0, n)
            if(length(IndexLev)==0){
              Zstart = Zmod
              ystart = y
            }else{
              Zstart = Zmod[-IndexLev,]
              ystart = y[-IndexLev]
              weights.xz2[IndexLev] = 0
            }
          }
        }
      }
    }
  }
  
  Tukey <- function(x, c){
    response <- ((x/c)^2-1)^2*x
    Ind <- which(abs(x)>c)
    if(length(Ind)>0){ response[Ind] <- 0 }
    return( response )
  }
  Huber <- function(x, c){
    x[x <= -c] <- -c
    x[x >= c] <- c
    return(x)
  }
  
  if (weightFunction1 == "Huber") {
    c1 <- optionList$huberC
    v_fun1<- function(resid, c){
      ind = which(resid==0)
      fun = (Huber(resid,c)/resid)^2
      #fun = (Huber(resid,c)/resid)
      fun[ind] = 1
      return(fun)
    }
  } else if (weightFunction1 == "Tukey") {
    c1 <- optionList$tukeyC
    v_fun1<- function(resid, c){
      ind = which(resid==0)
      fun = Tukey(resid,c)/resid
      fun[ind] = 1
      return(fun)
    }
  }
  if(is.null(weightFunction2)){
    c2 <- c1
    v_fun2 <- v_fun1
  }else if (weightFunction2 == "Huber") {
    c2 <- optionList$huberC
    v_fun2<- function(resid, c){
      ind = which(resid==0)
      fun = (Huber(resid,c)/resid)^2
      #fun = (Huber(resid,c)/resid)
      fun[ind] = 1
      return(fun)
    }
  } else if (weightFunction2 == "Tukey") {
    c2 <- optionList$tukeyC
    v_fun2<- function(resid, c){
      ind = which(resid==0)
      fun = Tukey(resid,c)/resid
      fun[ind] = 1
      return(fun)
    }
  }
  
  B1=NULL; B2=NULL; sD1=NULL; sD2=NULL; basis1=NULL;basis2=NULL;
  ##################################################### Basis functions #####################################################
  if(muEstim == 'gam'){
    if (length(sp1)==1){
      sp1 <- rep(sp1,nxs)
    }
    if ((length(sp1)!=nxs)||(!prod(sp1>=0))){
      stop('Please specify smoothing parameter for X!\n')
    }
    
    data1 <- data.frame(data.frame(X),data.frame(y))
    # select splines
    basis1 <- list()
    for (j in (1:nxs)){
      if (smooth.basis1=="tr"){
        stemp1 <- s(X[,j], bs="tr", m=p1, k=p1+K1+1)
        stemp1$term <- names(data1)[j]; stemp1$label <- paste(c("s(",names(data1)[j],")"),collapse="")
        basis1[[j]]  <- smooth.construct.tr.smooth.spec(stemp1,data1,NULL)
      } else if (smooth.basis1=="tp") {
        stemp1 <- s(X[,j], bs="tp", m=p1)
        stemp1$term <- names(data1)[j]; stemp1$label <- paste(c("s(",names(data1)[j],")"),collapse="")
        basis1[[j]]  <- smooth.construct.tp.smooth.spec(stemp1,data1,NULL)
      } else if (smooth.basis1=="cr") {
        stemp1 <- s(X[,j], bs="cr", k=K1)
        stemp1$term <- names(data1)[j]; stemp1$label <- paste(c("s(",names(data1)[j],")"),collapse="")
        basis1[[j]]  <- smooth.construct.cr.smooth.spec(stemp1,data1,NULL)
      } else if (smooth.basis1=="ps") {
        stemp1 <- s(X[,j], bs="ps", m=c(2,2), k=2+K1+2)
        stemp1$term <- names(data1)[j]; stemp1$label <- paste(c("s(",names(data1)[j],")"),collapse="")
        basis1[[j]]  <- smooth.construct.ps.smooth.spec(stemp1,data1,NULL)
        basis1[[j]]$df <- dim(basis1[[j]]$X)[2]
      }
    }
    
    # impose the identifiablity constraint
    dfs1 <- sapply(basis1,function(b){b$df})
    B1 <- basis1[[1]]$X
    sD1 <- matrix(0, nrow=dfs1[1]+sum(dfs1[-1]-1), ncol=dfs1[1]+sum(dfs1[-1]-1))
    sD1[1:dfs1[1],1:dfs1[1]] <- sp1[1]*basis1[[1]]$S[[1]]
    G1 <- list()
    G1[[1]] <- NULL
    if (nxs>1){
      tempindex1 <- dfs1[1]
      for (j in (2:nxs)){
        G1[[j]] <- qr.Q(qr(t(basis1[[j]]$X)%*%rep(1,n)), TRUE)[,-1]
        B1 <- cbind(B1,basis1[[j]]$X%*%G1[[j]])
        sD1[(tempindex1+1):(tempindex1+dfs1[j]-1),(tempindex1+1):(tempindex1+dfs1[j]-1)] <- sp1[j]*(t(G1[[j]])%*%basis1[[j]]$S[[1]]%*%G1[[j]])
        tempindex1 <- tempindex1 + dfs1[j] - 1
      }
    }
    rS1 <- mat.sqrt(sD1)
  }
  
  if(thetaEstim == 'gam'){
    if(is.null(sp2)){sp2 = sp1}
    
    if (length(sp2)==1){
      sp2 <- rep(sp2,nzs)
    }
    if ((length(sp2)!=nzs)||(!prod(sp2>=0))){
      stop('Please specify smoothing parameter for Z!\n')
    }
    
    data2 <- data.frame(data.frame(Z),data.frame(y))
    # select splines
    basis2 <- list()
    for (j in (1:nzs)){
      if (smooth.basis2=="tr"){
        stemp2 <- s(Z[,j], bs="tr", m=p2, k=p2+K2+1)
        stemp2$term <- names(data2)[j]; stemp2$label <- paste(c("s(",names(data2)[j],")"),collapse="")
        basis2[[j]]  <- smooth.construct.tr.smooth.spec(stemp2,data2,NULL)
      } else if (smooth.basis2=="tp") {
        stemp2 <- s(Z[,j], bs="tp", m=p2)
        stemp2$term <- names(data2)[j]; stemp2$label <- paste(c("s(",names(data2)[j],")"),collapse="")
        basis2[[j]]  <- smooth.construct.tp.smooth.spec(stemp2,data2,NULL)
      } else if (smooth.basis2=="cr") {
        stemp2 <- s(Z[,j], bs="cr", k=K2)
        stemp2$term <- names(data2)[j]; stemp2$label <- paste(c("s(",names(data2)[j],")"),collapse="")
        basis2[[j]]  <- smooth.construct.cr.smooth.spec(stemp2,data2,NULL)
      } else if (smooth.basis2=="ps") {
        stemp2 <- s(Z[,j], bs="ps", m=c(2,2), k=2+K2+2)
        stemp2$term <- names(data2)[j]; stemp2$label <- paste(c("s(",names(data2)[j],")"),collapse="")
        basis2[[j]]  <- smooth.construct.ps.smooth.spec(stemp2,data2,NULL)
        basis2[[j]]$df <- dim(basis2[[j]]$X)[2]
      }
    }
    
    # impose the identifiablity constraint
    dfs2 <- sapply(basis2,function(b){b$df})
    B2 <- basis2[[1]]$X
    sD2 <- matrix(0, nrow=dfs2[1]+sum(dfs2[-1]-1), ncol=dfs2[1]+sum(dfs2[-1]-1))
    sD2[1:dfs2[1],1:dfs2[1]] <- sp2[1]*basis2[[1]]$S[[1]]
    G2 <- list()
    G2[[1]] <- NULL
    if (nzs>1){
      tempindex2 <- dfs2[1]
      for (j in (2:nxs)){
        G2[[j]] <- qr.Q(qr(t(basis2[[j]]$X)%*%rep(1,n)), TRUE)[,-1]
        B2 <- cbind(B2,basis2[[j]]$X%*%G2[[j]])
        sD2[(tempindex2+1):(tempindex2+dfs2[j]-1),(tempindex2+1):(tempindex2+dfs2[j]-1)] <- sp2[j]*(t(G2[[j]])%*%basis2[[j]]$S[[1]]%*%G2[[j]])
        tempindex2 <- tempindex2 + dfs2[j] - 1
      }
    }
    rS2 <- mat.sqrt(sD2)
  }
  
  ##################################################### Expectations #####################################################
  
  # choose the fisher consistency correction
  if (family[[1]]=="poisson"){
    nu_fun2 <- nu_fun_pois2
  } else if (family[[1]]=="binomial"){
    nu_fun2 <- nu_fun_bin2
  }
  
  nu_fun1 <- function(v_fun, resid, c, sqrtVar){
    return(v_fun(resid, c) * resid / sqrtVar)
  }
  ##################################################### initial estimate #####################################################
  if(muEstim == 'gam'){
    if(thetaEstim == 'gam'){
      deltaOld <- rep(0.001, ncol(B1) + ncol(B2))
      muInd = 1:ncol(B1)
      thetaInd = (ncol(B1) + 1):(ncol(B1) + ncol(B2))
      if(is.null(gamma.ini)){
        theta.initial = as.vector(exp(B2%*%deltaOld[(1+ncol(B1)):(ncol(B1) + ncol(B2))]))
      }else{
        theta.initial = as.vector(exp(B2%*%gamma.ini))
        deltaOld[(1+ncol(B1)):(ncol(B1) + ncol(B2))] = gamma.ini
      }
      theta = theta.initial
    }else if(thetaEstim == 'glm'){
      deltaOld <- rep(0.001, ncol(B1) + q)
      muInd = 1:ncol(B1)
      thetaInd = (ncol(B1) + 1):(ncol(B1) + q)
      if(!is.null(gamma.ini)){
        deltaOld[(1+ncol(B1)):(ncol(B1) + q)] <- gamma.ini
      }
      theta = as.vector(exp(Z%*%deltaOld[(1+ncol(B1)):(ncol(B1) + q)]))
      theta.initial = theta
    }else if(thetaEstim == 'pen'){
      deltaOld <- rep(0.001, ncol(B1) + q)
      muInd = 1:ncol(B1)
      thetaInd = (ncol(B1) + 1):(ncol(B1) + q)
      if(is.null(gamma.ini)){
        theta.initial = as.vector(exp(Z%*%deltaOld[(1+ncol(B1)):(ncol(B1) + q)]))
      }else{
        theta.initial = as.vector(exp(Z%*%gamma.ini))
        deltaOld[(1+ncol(B1)):(ncol(B1) + q)] = gamma.ini
      }
      theta = theta.initial
    }
    
    if(is.null(beta.ini)){
      fit <- fit.gam.sp1(y, B1, rS1, family)
      deltaOld[1:ncol(B1)] <- as.vector(fit$coefficients)
      m.initial <- as.vector(fit$fitted.values)
    }else{
      deltaOld[1:ncol(B1)] <- beta.ini
      m.initial <- as.vector(family$linkinv(B1 %*% beta.ini))
    }
  }else if((muEstim == 'glm') || (muEstim == 'pen')){
    if(thetaEstim == 'gam'){
      deltaOld <- rep(0.001, p + ncol(B2))
      muInd = 1:p
      thetaInd = (p + 1):(p + ncol(B2))
      if(is.null(gamma.ini)){
        theta.initial = as.vector(exp(B2%*%deltaOld[(1+p):(p + ncol(B2))]))
      }else{
        theta.initial = as.vector(exp(B2%*%gamma.ini))
        deltaOld[(1+p):(p + ncol(B2))] = gamma.ini
      }
      theta = theta.initial
    }else if(thetaEstim == 'glm'){
      deltaOld <- rep(0.001, p + q)
      muInd = 1:p
      thetaInd = (p + 1):(p + q)
      
      if(!is.null(gamma.ini)){
        deltaOld[(p+1):(p+q)] <- gamma.ini
      }
      theta = as.vector(exp(Z%*%deltaOld[(1+p):(p + q)]))
      theta.initial = theta
    }else if(thetaEstim == 'pen'){
      deltaOld <- rep(0.001, p + q)
      muInd = 1:p
      thetaInd = (p + 1):(p + q)
      if(is.null(gamma.ini)){
        theta.initial = as.vector(exp(Z%*%deltaOld[(1+p):(p + q)]))
      }else{
        theta.initial = as.vector(exp(Z%*%gamma.ini))
        deltaOld[(1+p):(p + q)] = gamma.ini
      }
      theta = theta.initial
    }
    
    if (is.null(beta.ini)){
      if(nxs <= n){
        if(family[[1]]=="binomial"){
          yStart <- cbind(mBin*y, mBin-mBin*y)
        }else{
          yStart <- y
        }
        
        if(p==1){
          suppressWarnings(glmrob.ini <- glmrob(yStart ~ 1, method = "Mqle", family = family[[1]],
                               control = glmrobMqle.control(acc = 0.0001 * optionList$tol,
                                                            maxit = 500*optionList$maxIt,
                                                            tcc = optionList$huberC),
                               weights.on.x = weights.xz1,
                               data = data.frame(y, X)))
        } else {
          suppressWarnings(glmrob.ini <- glmrob(yStart~X[,-1],
                               method = "Mqle",
                               family = family[[1]],
                               weights.on.x =weights.xz1 ,
                               control = glmrobMqle.control(acc = 0.0001 * optionList$tol,
                                                            maxit = 500*optionList$maxIt,
                                                            tcc = optionList$huberC)))
          
        }
        betaOld <- as.numeric(glmrob.ini$coefficients)
        deltaOld[1:p] <- betaOld
        m.initial <- as.vector(family$linkinv(X %*% betaOld))
      }else{ #nxs > n
        med = apply(Xstart,2,median)
        Qn = apply(Xstart,2,Qn)
        XstartSc = (Xstart -  matrix(med,nrow(Xstart),ncol(Xstart),byrow = TRUE)) %*% diag(1/Qn)
        YSc = (ystart-median(ystart))/Qn(ystart)
        
        Enorm = sqrt(rowSums(XstartSc^2)+YSc^2)
        
        GoodIndex = order(Enorm)[1:floor(n/2)]
        
        XmodGood = Xstart[GoodIndex,]
        YGood = ystart[GoodIndex]
        
        model <- glmnet(x=XmodGood, y=YGood, family=family[[1]], intercept=intercept)
        deltaOld[1:p] <- tail(as.vector(predict(model, x=XmodGood,  family=family[[1]],y=YGood, intercept=intercept, s = lambdaBeta/n, type="coefficients", exact = T)),nxs)
        m.initial <- as.vector(family$linkinv(X %*% deltaOld[1:p]))
      }
    }else {
      deltaOld[1:p] <- beta.ini
      m.initial <- as.vector(family$linkinv(X %*% beta.ini))
    }
  }
  mu = as.vector(m.initial)
  deltaTemp <- deltaOld
  deltaStart <- deltaOld
  
  #####################
  # Estimation        #
  #####################
  
  i <- 1
  while(i <= optionList$maxIt){
    deltaOld <- deltaTemp
    conv_G <- FALSE
    conv_B <- FALSE
    
    #Update beta
    theta = as.vector(theta)
    mu = as.vector(mu)
    if(muEstim == 'glm'){
      maxIt.beta <- 4 * optionList$maxIt
      for(jBeta in 1:maxIt.beta){
        betaUpdate <- try(
          CalculateBetaUpdate(y = y, X = X, Z = Z,
                              family = family[[1]],
                              m = mBin,
                              mu, theta,
                              c = c1, weightF = v_fun1,
                              weights.xz = weights.xz1,
                              optionList = optionList),
          silent = TRUE)
        
        if(!is.matrix(betaUpdate)){
          break()
        }
        if (any(!is.finite(betaUpdate))) {
          break()
        }
        relE_B1 <- sqrt(sum(betaUpdate^2)/max(1e-20, sum(deltaTemp[1:p]^2)))
        relE_B2 <- max(abs(betaUpdate))
        conv_B <- (relE_B1 <= optionList$tol)|(relE_B2 <= optionList$tol)
        deltaTemp[muInd] <- deltaTemp[muInd] + betaUpdate
        beta.fit = deltaTemp[muInd]
        mu = as.vector(family$linkinv(X%*%beta.fit))
        if (conv_B) break()
      }
    }else if(muEstim == 'gam'){
      maxIt.beta <- 4 * optionList$maxIt
      
      
      for(jBeta in 1:maxIt.beta){
        betaNew <- try(
          CalculateBetaUpdateGAM(y = y, B = B1, rS = rS1, family = family,
                                 mBin = mBin,
                                 mu, theta,
                                 c = c1, v_fun = v_fun1,
                                 weights.xz = weights.xz1),
          silent = TRUE)
        if(!is.matrix(betaNew)){
          break()
        }
        if (any(!is.finite(betaNew))) {
          break()
        }
        
        relE_B1 <- sqrt(sum((betaNew-deltaTemp[muInd])^2)/max(1e-20, sum(deltaTemp[muInd]^2)))
        relE_B2 <- max(abs(betaNew-deltaTemp[muInd]))
        conv_B <- (relE_B1 <= optionList$tol)|(relE_B2 <= optionList$tol)
        deltaTemp[muInd] <- betaNew
        beta.fit = deltaTemp[muInd]
        mu = as.vector(family$linkinv(B1%*%beta.fit))
        if (conv_B) break()
      }
    }else if(muEstim == 'pen'){
      maxIt.beta <- 2*optionList$maxIt
      for(jBeta in 1:maxIt.beta){
        betaNew <- try(
          CalculateBetaUpdateLasso(y = y, X = X, Z = Z,
                                   family = family,
                                   mBin = mBin,
                                   mu, theta, lambdaBeta,
                                   c = c1, v_fun = v_fun1,
                                   weights.xz = weights.xz1,
                                   optionList = optionList, intercept),
          silent = TRUE)
        if(!is.vector(betaNew)){
          break()
        }
        if (any(!is.finite(betaNew))) {
          break()
        }
        
        relE_B1 <- sqrt(sum((betaNew-deltaTemp[muInd])^2)/max(1e-20, sum(deltaTemp[muInd]^2)))
        relE_B2 <- max(abs(betaNew-deltaTemp[muInd]))
        conv_B <- (relE_B1 <= optionList$tol)|(relE_B2 <= optionList$tol)
        deltaTemp[muInd] <- betaNew
        beta.fit = deltaTemp[muInd]
        mu = as.vector(family$linkinv(X%*%beta.fit))
        if (conv_B) break()
      }
    }
    # print(Sys.time())
    # print(deltaTemp[muInd])
    theta = as.vector(theta)
    mu = as.vector(mu)
    maxIt.gamma <- optionList$maxIt
    if (i == 1) maxIt.gamma <- 3 * maxIt.gamma
    if(thetaEstim == 'glm'){
      for (jGamma in 1:maxIt.gamma){
        gammaUpdate <- try(
          CalculateGammaUpdate(y = y, X = X, Z = Z,
                               family = family[[1]],
                               m = mBin,
                               mu, theta,
                               c = c2, weightF = v_fun2,
                               weights.xz = weights.xz2,
                               optionList = optionList),
          silent = TRUE)
        if (!is.matrix(gammaUpdate)){
          break()
        }
        if (any(!is.finite(gammaUpdate))) {
          break()
        }
        relE_G1 <- sqrt(sum(gammaUpdate^2)/max(1e-20, sum(deltaTemp[(p+1):(p+q)]^2)))
        relE_G2 <- max(abs(gammaUpdate))
        conv_G <- (relE_G1 <= optionList$tol*1e-3)|(relE_G2 <= optionList$tol*1e-3)
        deltaTemp[thetaInd] <- deltaTemp[thetaInd] + gammaUpdate
        gamma.fit = deltaTemp[thetaInd]
        theta = as.vector(exp(Z%*%gamma.fit))
        if (conv_G) break()
      }
    }else if(thetaEstim == 'gam'){
      for(jGamma in 1:maxIt.gamma){
        gammaNew <- try(
          CalculateGammaUpdateGAM(y = y, B = B2, rS = rS2, family = family,
                                  mBin = mBin,
                                  mu, theta,
                                  c = c2, v_fun = v_fun2, nu_fun = nu_fun2,
                                  weights.xz = weights.xz2),
          silent = TRUE)
        if(!is.matrix(gammaNew)){
          break()
        }
        if (any(!is.finite(gammaNew))) {
          break()
        }
        relE_G1 <- sqrt(sum((gammaNew-deltaTemp[thetaInd])^2)/max(1e-20, sum(deltaTemp[thetaInd]^2)))
        relE_G2 <- max(abs(gammaNew-deltaTemp[thetaInd]))
        conv_G <- (relE_G1 <= optionList$tol)|(relE_G2 <= optionList$tol)
        deltaTemp[thetaInd] <- gammaNew
        gamma.fit = deltaTemp[thetaInd]
        theta = as.vector(exp(B2%*%gamma.fit))
        if (conv_G) break()
      }
    }else if(thetaEstim == 'pen'){
      for(jGamma in 1:maxIt.gamma){
        gammaNew <- try(
          CalculateGammaUpdateLasso(y = y, X = X, Z = Z,
                                    family = family,
                                    mBin = mBin,
                                    mu, theta, lambdaGamma,
                                    c = c2, v_fun = v_fun2, nu_fun = nu_fun2,
                                    weights.xz = weights.xz2,
                                    optionList = optionList, intercept),
          silent = TRUE)
        if(!is.vector(gammaNew)){
          break()
        }
        if (any(!is.finite(gammaNew))) {
         break()
        }
        
        relE_G1 <- sqrt(sum((gammaNew-deltaTemp[thetaInd])^2)/max(1e-20, sum(deltaTemp[thetaInd]^2)))
        relE_G2 <- max(abs(gammaNew-deltaTemp[thetaInd]))
        conv_G <- (relE_G1 <= optionList$tol)|(relE_G2 <= optionList$tol)
        deltaTemp[thetaInd] <- gammaNew
        gamma.fit = deltaTemp[thetaInd]
        theta = as.vector(exp(Z%*%gamma.fit))
        if (conv_G) break()
      }
    }
    # print(Sys.time())
    # print(deltaTemp[thetaInd])
    relE <- sqrt(sum((deltaOld - deltaTemp)^2/max(deltaOld^2, 10^-20)))
    conv <- relE <= optionList$tol
    if (conv){
      i <- optionList$maxIt + 1
    }
    i <- i + 1
  }
  
  beta <- deltaTemp[muInd]
  gamma <- deltaTemp[thetaInd]
  EstTarget <- deltaTemp
  
  # reconstructing the beta in original basis representation
  if (nxs>1 && muEstim == 'gam'){
    tempindex1 <- dfs1[1]
    for (j in (2:nxs)){
      beta <- c(beta,G1[[j]]%*%beta.fit[(tempindex1+1):(tempindex1+dfs1[j]-1)])
      tempindex1 <- tempindex1 + dfs1[j]
    }
  }
  
  # reconstructing the gamma in original basis representation
  if (nzs>1 && thetaEstim == 'gam'){
    tempindex2 <- dfs2[1]
    for (j in (2:nzs)){
      gamma <- c(gamma,G2[[j]]%*%gamma.fit[(tempindex2+1):(tempindex2+dfs2[j]-1)])
      tempindex2 <- tempindex2 + dfs2[j]
    }
  }
  if((muEstim == 'glm' ) && thetaEstim =='glm'){
    ASInfo <- CalculateAsymptoticInfo(y = y, X = X, Z=Z, m = mBin, mu, theta,
                                      family = family[[1]],
                                      c1 = c1,
                                      weightF1 = v_fun1,
                                      c2 = c2,
                                      weightF2 = v_fun2,
                                      weights.xz1 = weights.xz1,
                                      weights.xz2 = weights.xz2,
                                      optionList = optionList)
    
    Result <- list(betaEstimate = beta, gammaEstimate = gamma,
                   thetaStart = deltaStart,
                   fitted.mu.values=mu,fitted.theta.values=theta,
                   mBin = mBin,
                   c1 = c1, c2 = c2,
                   v_fun1 = v_fun1, v_fun2 = v_fun2,
                   nu_fun1 = nu_fun1, nu_fun2 = nu_fun2, 
                   weights.xz1 = weights.xz1, weights.xz2 = weights.xz2,
                   M = ASInfo$M, Q = ASInfo$Q,
                   ASVar = ASInfo$ASVar
    )
  }else{
    Result <- list(betaEstimate = beta, gammaEstimate = gamma,
                   thetaStart = deltaStart, 
                   fitted.mu.values=mu,fitted.theta.values=theta,
                   mBin = mBin,
                   c1 = c1, c2 = c2,
                   v_fun1 = v_fun1, v_fun2 = v_fun2,
                   nu_fun1 = nu_fun1, nu_fun2 = nu_fun2, 
                   weights.xz1 = weights.xz1, weights.xz2 = weights.xz2,
                   B1=B1,B2=B2,sD1=sD1,sD2=sD2
    )
  }  
  
  return(Result)
}

# This function calculates asymptotic information: M matrix, Q matrix and asymptotic variance matrix.
CalculateAsymptoticInfo <- function(y, X, Z, m, mu, theta,
                                    family,
                                    c1, weightF1,
                                    c2, weightF2,
                                    weights.xz1, weights.xz2,
                                    optionList){
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)
  
  if(family == "binomial"){
    if(is.null(m)){
      stop('Please provide the variable m.')
    }
  }
  
  # Calculate the AMSE
  if (family == "poisson") {
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
      weight1 <- weightF1(R,c1)
      weight2 <- weightF2(R,c2)
      
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
    Q = matrix(0,p+q,p+q); M = matrix(0,p+q,p+q)
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1); Q[(p+1):(p+q),1:p]=Q21-1/n*tcrossprod(a2,a1); Q[1:p,(p+1):(p+q)]=Q12-1/n*tcrossprod(a1,a2); Q[(p+1):(p+q),(p+1):(p+q)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2; M[(p+1):(p+q),1:p]=M21-1/n*a2%*%BdenMat2; M[1:p,(p+1):(p+q)]=M12-1/n*a1%*%GdenMat2; M[(p+1):(p+q),(p+1):(p+q)]=GdenMat1-1/n*a2%*%GdenMat2
    Q = Q/n
    M = M/n
    Minv = solve(M)
    ASVar = Minv %*% Q %*% t(Minv)
  }else if(family == "binomial"){
    linpred = log(mu/(1-mu))
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
      weight1 <- weightF1(R, c1)
      weight2 <- weightF2(R,c2)
      
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
    Q = matrix(0,p+q,p+q); M = matrix(0,p+q,p+q)
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1); Q[(p+1):(p+q),1:p]=Q21-1/n*tcrossprod(a2,a1); Q[1:p,(p+1):(p+q)]=Q12-1/n*tcrossprod(a1,a2); Q[(p+1):(p+q),(p+1):(p+q)]=GdenMat3-1/n*tcrossprod(a2)
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2; M[(p+1):(p+q),1:p]=M21-1/n*a2%*%BdenMat2; M[1:p,(p+1):(p+q)]=M12-1/n*a1%*%GdenMat2; M[(p+1):(p+q),(p+1):(p+q)]=GdenMat1-1/n*a2%*%GdenMat2
    Q = Q/n
    M = M/n
    
    Minv = solve(M)
    ASVar = Minv %*% Q %*% t(Minv)
  }else{
    stop("Family not yet implemented.")
  }
  return(list(M = M,
              Q = Q,
              ASVar = ASVar))
}

# This function is used to update the beta parameter in GLM setting.
CalculateBetaUpdate <- function(y, X, Z,
                                family,
                                m,
                                mu, theta,
                                c, weightF,
                                weights.xz,
                                optionList){
  # Data matrix X rows correspond to observations
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)
  
  if (family == "poisson") {
    dLogMu <- (y - mu) / (mu / theta)
    dMuBeta <- X * mu # contains at row i the vector mu * X[i,]
    dLogTheta <- 1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
    tInd <- which(y == 0)
    if (length(tInd) > 0) dLogTheta[tInd] <- 1 / (2 * theta[tInd]) - mu[tInd]
    dThetaGamma <- Z * theta # contains at row i the vector theta * Z[i,]
    PearsonResid <- (y - mu) / sqrt(mu / theta)
  } else if (family == "binomial") {
    linpred = log(mu/(1-mu))
    dLogMu <- (m * y * theta) / mu + (theta * m * (1 - y)) / (mu - 1)
    dMuBeta <- X * (mu / (1 + exp(linpred))) # contains at row i the vector vec[i] * X[i,]
    dLogTheta <- 1 / (2 * theta) + (m * y) * log(mu / y) + m * (1 - y) * log((1 - mu) / (1 - y))
    dThetaGamma <- Z * theta # contains at row i the vector theta * Z[i,]
    
    #Adapt for boundary values
    tInd <- which(y <= 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(1 - mu[tInd])
    }
    tInd <- which(y >= 1 - 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(mu[tInd])
    }
    
    PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta))
  } else {
    stop("Requested family not yet implemented.")
  }
  
  # Calculate the expectation
  ExpTermBeta <- rep(0.0, p)
  B11 <- rep(0.0, n)
  if (family == "poisson") {
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
      weight <- weightF(R,c)
      
      ExpTermBeta <- ExpTermBeta + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dMuBeta[i,]))
      B11[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * mu[i] ^ 2
    }
    
    ExpTermBeta = (1 / n) * ExpTermBeta
    B11 = Diagonal(x=B11)
  } else if (family == "binomial") {
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = m, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)
    for (i in 1:n) {
      j <- minValue[i]:maxValue[i] / m[i]
      Probs <- dDBinom(j, n = m[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m[i] * theta[i]))
      LLTerm <- (m[i] * j * theta[i]) / mu[i] + (theta[i] * m[i] * (1 - j)) / (mu[i] - 1)
      if(j[1] <= tol) LLTerm[1] <- m[i] * theta[i] / (mu[i] - 1)
      if(j[length(j)] >= 1-tol) LLTerm[length(j)] <- m[i] * theta[i] / mu[i]
      weight <- weightF(R,c)
      
      ExpTermBeta <- ExpTermBeta + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dMuBeta[i,]))
      B11[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * mu[i] ^ 2
    }
    ExpTermBeta = (1 / n) * ExpTermBeta
    B11 = Diagonal(x=B11)
  } else {
    stop("Requested family not yet implemented.")
  }
  
  weight <- weightF(PearsonResid,c)
  
  # Calculate the LogLikelihoodValue
  LLBeta <- rep(0.0, p)
  # for (i in 1:n) {
  #   LLBeta <- LLBeta + weight[i] * dLogMu[i] * weights.xz1[i] * as.matrix(dMuBeta[i,])
  # }
  LLBeta <- as.matrix(colSums(dMuBeta * (weight * dLogMu * weights.xz))) # matrix expression for loop above.
  LLBeta <- LLBeta - n * ExpTermBeta
  # Calculate the PsiMatrix
  PsiDeriv <- crossprod(X, B11) %*% X
  PsiDeriv <- as.matrix(PsiDeriv)
  return(solve(PsiDeriv,diag(p)) %*% LLBeta)
}

# This function is used to update the beta parameter in GLM setting with lasso penalty.
CalculateBetaUpdateLasso <- function(y, X, Z, family, mBin, mu, theta, lambdaBeta,
                                     c, v_fun,  weights.xz,
                                     optionList, intercept){
  # Data matrix X rows correspond to observations
  
  n <- length(y)
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
  resp = eta +(nu - expNu)/(expNuUMu*family$mu.eta(eta))
  weightsAM = diag(c(expNuUMu * weights.xz * family$mu.eta(eta)^2))
  
  if(intercept){
    if(ncol(X)<3){
      beta <- lassoshooting::lassoshooting(y=sqrt(weightsAM)%*%resp,X=sqrt(weightsAM)%*%X,lambda = lambdaBeta,nopenalize=0)$coefficients
    }else{
      fitModel=glmnet(X[,-1],resp, standardize = F, weights = diag(weightsAM))
      beta=as.vector(predict(fitModel, x=X[,-1], y=resp, weights = diag(weightsAM), s = (lambdaBeta/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T))  
    }
  }else{
    if(ncol(X)<2){
      beta <- lassoshooting::lassoshooting(y=sqrt(weightsAM)%*%resp,X=sqrt(weightsAM)%*%X,lambda = lambdaBeta)$coefficients
    }else{
      fitModel=glmnet(X,resp, standardize = F, weights = diag(weightsAM), intercept = F)
      beta=as.vector(predict(fitModel, x=X, y=resp, weights = diag(weightsAM), s = (lambdaBeta/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T, intercept = F))[-1]  
    }
  }
  #Perform adaptive lasso on non-zero elements
  nonZeroBeta <- which(beta!=0)
  
  if(intercept & length(nonZeroBeta)>2){
    fitModel=glmnet(X[,nonZeroBeta[-1]],resp, standardize = F, weights = diag(weightsAM), penalty.factor = abs(1/beta)[nonZeroBeta[-1]])
    betaAdd=as.vector(predict(fitModel, x=X[,nonZeroBeta[-1]], y=resp, weights = diag(weightsAM), penalty.factor = abs(1/beta)[nonZeroBeta[-1]], s = (lambdaBeta/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T))
    beta[nonZeroBeta] = betaAdd
  }else if(!intercept & length(nonZeroBeta)>=2){
    fitModel=glmnet(X[,nonZeroBeta],resp, standardize = F, weights = diag(weightsAM), penalty.factor = abs(1/beta)[nonZeroBeta], intercept = F)
    betaAdd=as.vector(predict(fitModel, x=X[,nonZeroBeta], y=resp, weights = diag(weightsAM), penalty.factor = abs(1/beta)[nonZeroBeta], s = (lambdaBeta/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T, intercept = F))[-1]
    beta[nonZeroBeta] = betaAdd
  }
  return(beta)
}

# This function is used to update the beta parameter in GAM setting.
CalculateBetaUpdateGAM <- function(y, B, rS,
                                   family,
                                   mBin,
                                   mu, theta,
                                   c, v_fun,
                                   weights.xz){
  # Data matrix X rows correspond to observations
  n <- length(y)
  
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
    stop("Requested family not yet implemented.")
  }
  resp = eta +(nu - expNu)/(expNuUMu*family$mu.eta(eta))
  weightsAM = diag(c(expNuUMu * weights.xz * family$mu.eta(eta)^2, rep(1,ncol(B))))
  
  #add penalty part
  resp = c(resp, rep(0, ncol(B)))
  X1 = rbind(B, rS)
  
  return(solve(t(X1)%*%weightsAM%*%X1) %*% t(X1)%*%weightsAM%*%resp)
}

# This function is used to update the gamma parameter in GLM setting.
CalculateGammaUpdate <- function(y, X, Z,
                                 family,
                                 m,
                                 mu, theta,
                                 c, weightF,
                                 weights.xz,
                                 optionList){
  
  # Data matrix X rows correspond to observations
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)
  
  if (family == "poisson") {
    dLogMu <- (y - mu) / (mu / theta)
    dMuBeta <- X * mu # contains at row i the vector mu * X[i,]
    dLogTheta <- 1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
    tInd <- which(y == 0)
    if (length(tInd) > 0) dLogTheta[tInd] <- 1 / (2 * theta[tInd]) - mu[tInd]
    dThetaGamma <- Z * theta # contains at row i the vector theta * Z[i,]
    PearsonResid <- (y - mu) / sqrt(mu / theta)
  } else if (family == "binomial") {
    linpred = log(mu/(1-mu))
    dLogMu <- (m * y * theta) / mu + (theta * m * (1 - y)) / (mu - 1)
    dMuBeta <- X * (mu / (1 + exp(linpred))) # contains at row i the vector vec[i] * X[i,]
    dLogTheta <- 1 / (2 * theta) + (m * y) * log(mu / y) + m * (1 - y) * log((1 - mu) / (1 - y))
    dThetaGamma <- Z * theta # contains at row i the vector theta * Z[i,]
    
    #Adapt for boundary values
    tInd <- which(y <= 1e-12)
    if(length(tInd)>0) {
      dLogMu[tInd] <- m[tInd] * theta[tInd] / (mu[tInd] - 1)
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(1 - mu[tInd])
    }
    tInd <- which(y >= 1 - 1e-12)
    if(length(tInd)>0) {
      dLogMu[tInd] <- m[tInd] * theta[tInd] / mu[tInd]
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(mu[tInd])
    }
    PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta))
  } else {
    stop("Requested family not yet implemented.")
  }
  
  # Calculate the expectation
  ExpTermGamma <- rep(0.0, q)
  B22 <- rep(0.0, n)
  if (family == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmax(0, floor(mu - c * (1 / theta) * mu))
    maxValue <- ceiling(mu + c * (1/theta) * mu)
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)
    
    for(i in 1:n) {
      j <- minValue[i]:maxValue[i]
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i]/theta[i]))
      LLTerm <- 1/(2*theta[i]) - mu[i] + j * log(exp(1)*mu[i]/j)
      if (j[1] == 0) LLTerm[1] <- 1/(2*theta[i]) - mu[i]
      weight <- weightF(R, c)
      
      ExpTermGamma = ExpTermGamma + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dThetaGamma[i,]))
      B22[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * theta[i] ^ 2
    }
    ExpTermGamma = (1 / n) * ExpTermGamma
    B22 <- Diagonal(x=B22)
  } else if (family == "binomial"){
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = m, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)
    for(i in 1:n) {
      j <- minValue[i]:maxValue[i] / m[i]
      Probs <- dDBinom(j, n = m[i], mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (theta[i] * m[i]))
      LLTerm <- 1/(2*theta[i]) + (m[i] * j) * log(mu[i]/j) + m[i] * (1 - j) * log((1 - mu[i]) / (1 - j))
      if(j[1] <= tol) LLTerm[1] <- 1/(2*theta[i]) + m[i] * log(1 - mu[i])
      if(j[length(j)] >= 1-tol) LLTerm[length(j)] <- 1/(2*theta[i]) + m[i] * log(mu[i])
      weight <- weightF(R, c)
      
      ExpTermGamma = ExpTermGamma + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dThetaGamma[i,]))
      B22[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * theta[i] ^ 2
    }
    ExpTermGamma = (1 / n) * ExpTermGamma
    B22 <- Diagonal(x=B22)
  } else {
    stop("Requested family not yet implemented.")
  }
  
  weight <- weightF(PearsonResid, c)
  
  # Calculate the LogLikelihoodValue
  LLGamma <- rep(0.0, q)
  LLGamma <- as.matrix(colSums(dThetaGamma * (weight * dLogTheta * weights.xz))) # matrix expression for loop above
  LLGamma <- LLGamma - n * ExpTermGamma
  
  # Calculate the PsiMatrix
  PsiDeriv <- crossprod(Z, B22) %*% Z
  PsiDeriv <- as.matrix(PsiDeriv)
  
  
  return(solve(PsiDeriv, diag(q)) %*% LLGamma)
}

# This function is used to update the gamma parameterin GLM setting with lasso penalty.
CalculateGammaUpdateLasso <- function(y, X, Z, family,
                                      mBin,
                                      mu, theta, lambdaGamma,
                                      c, v_fun, nu_fun,
                                      weights.xz,
                                      optionList, intercept){
  # Data matrix X rows correspond to observations
  n <- length(y)
  
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
  
  nu = v_fun(PearsonResid, c)*UTheta
  expNu = rep(NA, n)
  expNuUTheta = rep(NA, n)
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
      jZero <- which(j == 0)
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      weight <- v_fun(R, c)
      
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
      
      weight <- v_fun(R, c)
      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUTheta[i] <- sum(weight * LLTerm^2* Probs)
    }
  } else {
    stop("Requested family not yet implemented.")
  }
  
  resp =  log(theta) + (nu - expNu)/(expNuUTheta * theta)
  weightsAM =  diag(c(expNuUTheta * weights.xz * theta^2))

  if(intercept){
    if(ncol(Z)<3){
      gamma <- lassoshooting::lassoshooting(y=sqrt(weightsAM)%*%resp,X=sqrt(weightsAM)%*%Z,lambda = lambdaGamma,nopenalize=0)$coefficients
    }else{
      fitModel=glmnet(Z[,-1],resp, standardize = F, weights = diag(weightsAM))
      gamma=as.vector(predict(fitModel, x=Z[,-1], y=resp, weights = diag(weightsAM), s = (lambdaGamma/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T))
    }
  }else{
    if(ncol(Z)<2){
      gamma <- lassoshooting::lassoshooting(y=sqrt(weightsAM)%*%resp,X=sqrt(weightsAM)%*%Z,lambda = lambdaGamma)$coefficients
    }else{
      fitModel=glmnet(Z,resp, standardize = F, weights = diag(weightsAM), intercept = F)
      gamma=as.vector(predict(fitModel, x=Z, y=resp, weights = diag(weightsAM), s = (lambdaGamma/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T, intercept = F))[-1]
    }
  }
  
  #Perform adaptive lasso on non-zero elements
  nonZeroGamma <- which(gamma!=0)
  if(intercept & length(nonZeroGamma)>2){
    fitModel=glmnet(Z[,nonZeroGamma[-1]],resp, standardize = F, weights = diag(weightsAM), penalty.factor = abs(1/gamma)[nonZeroGamma[-1]])
    gammaAdd=as.vector(predict(fitModel, x=Z[,nonZeroGamma[-1]], y=resp, weights = diag(weightsAM), penalty.factor = abs(1/gamma)[nonZeroGamma[-1]], s = (lambdaGamma/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T))
    gamma[nonZeroGamma] = gammaAdd
  }else if(!intercept & length(nonZeroGamma)>=2){
    fitModel=glmnet(Z[,nonZeroGamma],resp, standardize = F, weights = diag(weightsAM), penalty.factor = abs(1/gamma)[nonZeroGamma], intercept = F)
    gammaAdd=as.vector(predict(fitModel, x=Z[,nonZeroGamma], y=resp, weights = diag(weightsAM), penalty.factor = abs(1/gamma)[nonZeroGamma], s = (lambdaGamma/n)/sum(diag(weightsAM))*n, type="coefficients",  exact = T, intercept = F))[-1]
    gamma[nonZeroGamma] = gammaAdd
  }
  return(gamma)
}

# This function is used to update the gamma parameter in GAM setting.
CalculateGammaUpdateGAM <- function(y, B, rS,
                                    family,
                                    mBin,
                                    mu, theta,
                                    c, v_fun, nu_fun,
                                    weights.xz){
  n <- length(y)
  
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
  
  nu = v_fun(PearsonResid, c)*UTheta
  expNu = rep(NA, n)
  expNuUTheta = rep(NA, n)
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
      jZero <- which(j == 0)
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      weight <- v_fun(R, c)
      
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
      LLTerm[jNotZero] = LLTerm[jNotZero] + mBin[i]*y[jNotZero]*log(mu[i]/y[jNotZero])
      LLTerm[jNotOne] = LLTerm[jNotOne] + mBin[i]*(1-y[jNotOne])*log((1-mu[i])/(1-y[jNotOne]))
      
      weight <- v_fun(R, c)
      expNu[i] <- sum(weight * LLTerm * Probs)
      expNuUTheta[i] <- sum(weight * LLTerm^2* Probs)
    }
  } else {
    stop("Requested family not yet implemented.")
  }
  
  resp =  log(theta) + (nu - expNu)/(expNuUTheta * theta)
  weightsAM =  diag(c(expNuUTheta * weights.xz * theta^2,rep(1,ncol(B))))
  #add penalty part
  resp = c(resp, rep(0, ncol(B)))
  X1 = rbind(B, rS)
  
  return(solve(t(X1)%*%weightsAM%*%X1) %*% t(X1)%*%weightsAM%*%resp)
}


HuberResid <- function(x, c){
  x[x <= -c] <- -c
  x[x >= c] <- c
  return(x)
}

TukeyResid <- function(x, c){
  response <- ((x/c)^2-1)^2*x
  Ind <- which(abs(x)>c)
  if(length(Ind)>0){ response[Ind] <- 0 }
  return( response )
}


nu_fun_pois2 <- function(v_fun, y, mu, theta, c, mBin){
  if((length(theta)>1)){
    resid = (y-mu)/sqrt(mu/theta)
    if(length(y)>1){
      dLogTheta =  1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
      tInd <- which(y == 0)
      if (length(tInd) > 0) dLogTheta[tInd] <- 1 / (2 * theta[tInd]) - mu[tInd]
      return(v_fun(resid, c) * dLogTheta)
    }else{
      if(y==0){
        dLogTheta <- 1 / (2 * theta) - mu
        return(v_fun(resid, c) * dLogTheta)
      }else{
        dLogTheta =  1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
        return(v_fun(resid, c) * dLogTheta)
      }
    }
  }else{
    resid = (y-mu)/sqrt(mu/theta)
    dLogTheta =  1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
    tInd <- which(y == 0)
    if (length(tInd)>0) dLogTheta[tInd] <- 1 / (2 * theta) - mu
    return(v_fun(resid, c) * dLogTheta)
  }
}

nu_fun_bin2 <- function(v_fun, y, mu, theta, c, mBin){
  if((length(theta)>1)){
    resid = (y - mu) / sqrt(mu * (1 - mu) / (mBin * theta))
    if(length(y)>1){
      dLogTheta =  1 / (2 * theta) + (mBin * y) * log(mu / y) + mBin * (1 - y) * log((1 - mu) / (1 - y))
      
      tInd <- which(y <= 1e-12)
      if(length(tInd)>0) {
        dLogTheta[tInd] <- 1/(2*theta[tInd]) + mBin[tInd] * log(1 - mu[tInd])
      }
      tInd <- which(y >= 1 - 1e-12)
      if(length(tInd)>0) {
        dLogTheta[tInd] <- 1/(2*theta[tInd]) + mBin[tInd] * log(mu[tInd])
      }
      
      return(v_fun(resid, c) * dLogTheta)
    }else{
      if(y<= 1e-12){
        dLogTheta <- 1/(2*theta) + mBin * log(1 - mu)
        return(v_fun(resid, c) * dLogTheta)
      }else if(y>= 1 - 1e-12){
        dLogTheta =  1/(2*theta) + mBin * log(mu)
        return(v_fun(resid, c) * dLogTheta)
      }else{
        dLogTheta =  1 / (2 * theta) + (mBin * y) * log(mu / y) + mBin * (1 - y) * log((1 - mu) / (1 - y))
        return(v_fun(resid, c) * dLogTheta)
      }
    }
  }else{
    resid = (y - mu) / sqrt(mu * (1 - mu) / (mBin * theta))
    dLogTheta =  1 / (2 * theta) + (mBin * y) * log(mu / y) + mBin * (1 - y) * log((1 - mu) / (1 - y))
    tInd <- which(y <= 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta) + mBin * log(1 - mu)
    }
    tInd <- which(y >= 1 - 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta) + mBin * log(mu)
    }
    return(v_fun(resid, c) * dLogTheta)
  }
}

var <- function(m, theta, varfun, mBin){
  vari = varfun(m)/mBin
  return(vari/theta)
}
