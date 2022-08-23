# glmrob is a function of the package robustbase, that is copied due to a small 
# change in their function glmrobMqle such that it also accepts a vector value 
# for the parameter weights.on.x

# Description (from package robustbase)
# glmrob is used to fit generalized linear models by robust methods. The models 
# are specified by giving a symbolic description of the linear predictor and a 
# description of the error distribution. Currently, robust methods are 
# implemented for family = binomial, = poisson, = Gamma and = gaussian.
#

# Arguments (from package robustbase)
##' @param formula: a formula, i.e., a symbolic description of the model to be fit 
##' @param family: a description of the error distribution and link function to be used 
#         in the model. This can be a character string naming a family function, 
#         a family function or the result of a call to a family function. 
##' @param data: an optional data frame containing the variables in the model. If not 
#       found in data, the variables are taken from environment(formula), 
#       typically the environment from which glmrob is called.
##' @param weights: an optional vector of weights to be used in the fitting process.
##' @param subset: an optional vector specifying a subset of observations to be used in 
#         the fitting process.
##' @param na.action: a function which indicates what should happen when the data contain 
#            NAs. The default is set by the na.action setting in options. The 
#            "factory-fresh" default is na.omit.
##' @param start:    starting values for the parameters in the linear predictor. Note 
#           that specifying start has somewhat different meaning for the 
#           different methods. Notably, for "MT", this skips the expensive 
#           computation of initial estimates via sub samples, but needs to be 
#           robust itself.
##' @param offset: this can be used to specify an a priori known component to be included
#         in the linear predictor during fitting.
##' @param method: a character string specifying the robust fitting method. The details 
#         of method specification are given below.
##' @param weights.on.x: a character string (can be abbreviated), a function or list 
#               (see below), or a numeric vector of length n, specifying how 
#               points (potential outliers) in x-space are downweighted. If 
#               "hat", weights on the design of the form sqrt{1-h_{ii}} are used, 
#               where h_{ii} are the diagonal elements of the hat matrix. 
#               If "robCov", weights based on the robust Mahalanobis distance of 
#               the design matrix (intercept excluded) are used where the 
#               covariance matrix and the centre is estimated by cov.rob from 
#               the package MASS. Similarly, if "covMcd", robust weights are 
#               computed using covMcd. The default is "none".
#               If weights.on.x is a function, it is called with arguments 
#               (X, intercept) and must return an n-vector of non-negative 
#               weights.
#               If it is a list, it must be of length one, and as element 
#               contain a function much like covMcd() or cov.rob() 
#               (package MASS), which computes multivariate location and 
#               "scatter" of a data matrix X.
##' @param control: a list of parameters for controlling the fitting process. 
#          See the documentation for glmrobMqle.control for details.
##' @param model: a logical value indicating whether model frame should be included as a 
#        component of the returned value.
##' @param x, y: logical values indicating whether the response vector and model matrix 
#       used in the fitting process should be returned as components of the 
#       returned value.
##' @param contrasts: an optional list. See the contrasts.arg of model.matrix.default.
##' @param trace.lev: logical (or integer) indicating if intermediate results should be 
#            printed; defaults to 0 (the same as FALSE).
##' @param ... arguments passed to glmrobMqle.control when control is NULL (as per default).glmrob <-
function (formula, family, data, weights, subset,
	  na.action, start = NULL, offset,
          method = c("Mqle", "BY", "WBY", "MT"),
	  weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL,
	  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, trace.lev = 0,
	  ...)
{
    call <- match.call()
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    fami <- family$family
    if(is.null(fami))
	stop(gettextf("'%s' is not a valid family (see ?family)",
		      as.character(call[["family"]])), domain=NA)

    if (!(fami %in% c("binomial", "poisson", "Gamma", "gaussian"))) {
	stop(gettextf("Robust GLM fitting not yet implemented for family %s",
			  fami), domain=NA)
    }
    if (missing(data))
	data <- environment(formula)
    ##
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if(identical(method, "model.frame")) return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")# "numeric" or "factor"
    if (length(dim(Y)) == 1) {
	nm <- rownames(Y)
	dim(Y) <- NULL
	if (!is.null(nm))
	    names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
	model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0))
	stop("'weights' must be non-negative")
    if (!is.null(offset) && length(offset) != NROW(Y))
	stop(gettextf("Number of offsets is %d, should rather equal %d (number of observations)",
		      length(offset), NROW(Y)), domain=NA)
    method <- match.arg(method)
    meth. <- if(method == "WBY") "BY" else method
### FIXME: the whole 'control' should be changed to "copy"  lmrob() and lmrob.control()
## ------- --> *one* exported glmrob.control() function with 'method' and switch() inside...
## see >>> ./lmrob.MM.R

    if(is.null(control)) # -> use e.g., glmrobMqle.control()
	control <- get(paste0("glmrob", meth., ".control"))(...)
    if(missing(weights.on.x) || is.character(weights.on.x))
        weights.on.x <- match.arg(weights.on.x)
    else if(!(is.function(weights.on.x) || is.list(weights.on.x) ||
              (is.numeric(weights.on.x) && length(weights.on.x) == NROW(Y))))
        stop("'weights.on.x' must be a string, function, list or numeric n-vector")
    if(!is.null(start) && !is.numeric(start)) {
	## initialization methods
	if(!is.character(start))
	    stop("'start' must be a numeric vector, NULL, or a character string")
	start <-
	    switch(start,
		   "lmrob" =, "lmrobMM" = {
		       if(!is.null(weights))
			   warnings("weights are not yet used in computing start estimate")
		       lmrob.fit(x = X, y = family$linkinv(Y),
				 control=lmrob.control())$coefficients
		   },
		   stop("invalid 'start' string"))
    }
    fit <- switch(method,
		  "cubif" = stop("For method 'cubif', use glmRob() from package 'robust'")
		  ,
		  "Mqle" = ## --> ./glmrobMqle.R
		  glmrobMqle(X = X, y = Y, weights = weights, start = start,
			     offset = offset, family = family,
			     weights.on.x = weights.on.x, control = control,
			     intercept = attr(mt, "intercept") > 0, trace=trace.lev),
#                   "BY" =, "WBY" = {
#                       if(fami != "binomial")
#                           stop(gettextf(
# 			"method='%s' is only applicable for binomial family, but family=\"\"",
#                               method, fami), domain=NA)
#                       ### FIXME: use glmrobBY(..) with these arguments, including 'weights'
#                       glmrobBY(X=X, y=Y, weights=weights, start=start,
#                                method=method, ## == "BY" / "WBY"
#                                weights.on.x = weights.on.x, control = control,
#                                intercept = attr(mt, "intercept") > 0,
#                                trace.lev=trace.lev)
#                   },
#                   "MT" = {
#                       glmrobMT(x=X,y=Y, weights=weights, start=start, offset = offset,
# 			       family=family, weights.on.x=weights.on.x, control=control,
#                                intercept = attr(mt, "intercept") > 0, trace.lev=trace.lev)
#                   },
		  stop("invalid 'method': ", method))
    ##-	    if (any(offset) && attr(mt, "intercept") > 0) {
    ##-		fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
    ##-		    y = Y, weights = weights, offset = offset,
    ##-		    control = control, intercept = TRUE)$deviance
    ##-	    }
    fit$na.action <- attr(mf, "na.action")
    if (model)
	fit$model <- mf
    if (x)
	fit$x <- X
    if (!y) ## fit$y <- NULL
	warning("setting 'y = FALSE' has no longer any effect")
    fit <- c(fit,
	     list(call = call, formula = formula, terms = mt, data = data,
		  offset = offset, control = control, method = method,
		  prior.weights = if(is.null(weights)) rep.int(1, nrow(X))
		  else weights,
		  contrasts = attr(X, "contrasts"),
		  xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("glmrob", "glm")
    fit
}


summary.glmrob <- function(object, correlation=FALSE, symbolic.cor=FALSE, ...)
{
    dispersion <- object$dispersion
    if(is.null(dispersion)) dispersion <- 1
    coefs <- object$coefficients
    aliased <- is.na(coefs)# needs care; also used in print method
    if(any(aliased))
	coefs <- coefs[!aliased]
    covmat <- object$cov
    s.err <- sqrt(diag(covmat))
    zvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind("Estimate" = coefs, "Std. Error" = s.err,
			"z value" = zvalue, "Pr(>|z|)" = pvalue)

    ans <- c(object[c("call", "terms", "family", "iter", "control", "method",
		      "residuals", "fitted.values", "w.r", "w.x")],
	     ## MM: should rather keep more from 'object' ?
	     ##	    currently, cannot even print the asympt.efficiency!
	     list(deviance=NULL, df.residual=NULL, null.deviance=NULL,
		  df.null= NULL, df= NULL, ## (because of 0 weights; hmm,...)
		  aliased = aliased,
		  coefficients = coef.table, dispersion = dispersion,
		  cov.scaled = covmat))
    if (correlation) {
	ans$correlation <- cov2cor(covmat)
	ans$symbolic.cor <- symbolic.cor
    }
    structure(ans, class = "summary.glmrob")
}

## almost a copy of vcov.glm() [if that didn't have summmary.glm() explicitly]
vcov.glmrob <- function (object, ...)
{
    so <- summary(object, corr = FALSE, ...)
    ## so$dispersion * so$cov.unscaled
    ## changed from cov.unscaled to cov.scaled
    so$cov.scaled
}


print.glmrob <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(coef(x))) {
	cat("Coefficients")
	if (is.character(co <- x$contrasts))
	    cat("  [contrasts: ", apply(cbind(names(co), co),
					1, paste, collapse = "="), "]")
	cat(":\n")
	print.default(format(x$coefficients, digits = digits),
		      print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nNumber of observations:", length(x$residuals),
	"\nFitted by method ", sQuote(x$method), "\n")
    invisible(x)
}

print.summary.glmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(cf <- coef(x))) {
	if(nsingular <- sum(x$aliased)) # glm has   df[3] - df[1]
	    cat("\nCoefficients: (", nsingular,
		" not defined because of singularities)\n", sep = "")
	else cat("\nCoefficients:\n")
	printCoefmat(cf, digits = digits, signif.stars = signif.stars,
		     na.print = "NA", ...)

	summarizeRobWeights(x$w.r * x$w.x, digits = digits,
			    header = "Robustness weights w.r * w.x:", ...)
    }
    else cat("No coefficients\n\n")

    n <- length(x$residuals)
    cat("\nNumber of observations:", n,
	"\nFitted by method", sQuote(x$method)," (in", x$iter, "iterations)\n")

    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",sep = "")
    if(any(!is.null(unlist(x[c("null.deviance", "deviance")]))))
	cat(apply(cbind(paste(format(c("Null", "Residual"), justify="right"),
			      "deviance:"),
			format(unlist(x[c("null.deviance", "deviance")]),
			       digits=max(5, digits + 1)), " on",
			format(unlist(x[c("df.null", "df.residual")])),
			" degrees of freedom\n"),
		  1L, paste, collapse=" "), "\n", sep = "")
    else
	cat("No deviance values available \n")
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if (isTRUE(symbolic.cor)) {
		print(symnum(correl, abbr.colnames=NULL))
	    }
	    else {
		correl <- format(round(correl, 2), nsmall=2, digits=digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote=FALSE)
	    }
	}
    }

    printControl(x$control, digits = digits)
    cat("\n")
    invisible(x)
}

weights.glmrob <- function(object, type = c("prior", "robustness"), ...) {
    type <- match.arg(type)
    w <- if (type == "prior") {
	## Issue warning only if called from toplevel. Otherwise the warning pop
	## up at quite unexpected places, e.g., case.names().
	if (is.null(object[["weights"]]) && identical(parent.frame(), .GlobalEnv))
	    warning("No weights defined for this object. Use type=\"robustness\" argument to get robustness weights.")
	object[["weights"]]
    } else object$w.r * object$w.x ## those also used summarizeRobWeights(x$w.r * x$w.x, ..)
    if (is.null(object$na.action)) w else naresid(object$na.action, w)
}

## Stems from a copy of residuals.glm() in
## ~/R/D/r-devel/R/src/library/stats/R/glm.R
residuals.glmrob <-
    function(object,
	     type = c("deviance", "pearson", "working", "response",
             "partial"),
	     ...)
{
    type <- match.arg(type)
    y <- object$y
    r <- object$residuals
    mu	<- object$fitted.values
    wts <- object$prior.weights # ok
    p <- length(object$coefficients)
    switch(type,
           deviance=, pearson=, response=
           if(is.null(y)) {
               mu.eta <- object$family$mu.eta
               eta <- object$linear.predictors
               ## we cannot use 'r <- ...$residuals' __ FIXME __
               stop("need non-robust working residuals for this model type")
               y <-  mu + r * mu.eta(eta)
           })
    res <- switch(type,
##		  deviance = if(object$df.residual > 0) {
		  deviance = if((nobs(object) - p) > 0) {
		      d.res <- sqrt(pmax.int((object$family$dev.resids)(y, mu, wts), 0))
		      ifelse(y > mu, d.res, -d.res)
		  } else rep.int(0, length(mu)),
		  pearson = (y-mu)*sqrt(wts)/sqrt(object$family$variance(mu)),
		  working = r,
		  response = y - mu,
		  partial = r
		  )
    if(!is.null(object$na.action))
        res <- naresid(object$na.action, res)
    if (type == "partial") ## need to avoid doing naresid() twice.
        res <- res+predict(object, type="terms")
    res
}
