#' A function to estimate the taxon-wise response functions non-parametrically, and return their fitted values
#'
#' @param sampleScore a vector of length n with environmental scores
#' @param muMarg the offset matrix
#' @param X the n-by-p data matrix
#' @param ncols an integer, the number of columns of X
#' @param thetas a vector of length p with dispersion parameters
#' @param n an integer, the number of samples
#' @param coefInit a 2-by-p matrix with current taxon-wise parameter estimates
#' @param coefInitOverall a vector of length 2 with current overall parameters
#' @param dfSpline a scalar, the degrees of freedom for the smoothing spline.
#' @param vgamMaxit Maximal number of iterations in the fitting of the GAM model
#' @param colWeights a vector of length p with column weights
#' @param verbose a boolean, should number of failed fits be reported
#' @param ... further arguments, passed on to the VGAM:::vgam() function
#'
#' The negative binomial likelihood is still maximized, but now the response function is a non-parametric one. To avoid a perfect fit and overly flexible functions, we enforce smoothness restrictions. In practice we use a generalized additive model (GAM), i.e. with splines.  The same fitting procedure is carried out ignoring species labels. We do not normalize the parameters related to the splines: the psis can be calculated afterwards.
#'
#' @return A list with components
#' \item{taxonWiseFitted}{A n-by-p matrix of fitted valuesof the response curves per taxon at the observed values of the environmental scores}
#' \item{taxonCoef}{The fitted coefficients of the sample-wise response curves}
#' \item{overallFitted}{The fitted values of all taxa combined}
#' \item{taxonWiseFits}{A list of length p of normalized response curves of all the taxa. This may be useful to investigate the shape of the response function through plots.}
#' \item{psi}{The importance parameter of the dimension}
#' \item{taxonWise}{The taxonwise response function}
#'
#' @importFrom MASS negative.binomial
estNPresp = function(sampleScore, muMarg, X, psi, ncols, thetas, n, coefInit, coefInitOverall, dfSpline, vgamMaxit, colWeights, verbose,...){
  logMu = log(muMarg)
    taxonWise = lapply(seq_len(ncols), function(i){
    df = data.frame(x = X[,i], sampleScore = sampleScore, logMu = log(muMarg[,i])) #Going through a dataframe slows things down, so ideally we should appeal directly to the vgam.fit function
      tmp = try(suppressWarnings(vgam.edit(data = df,x ~ s(sampleScore, df = dfSpline), psi = psi, offset = logMu, family = negbinomial.size(lmu = "loge", size = thetas[i]), coefstart = coefInit[[i]], maxit = vgamMaxit,...)), silent = TRUE)
    # }
    if(class(tmp)=="try-error") { #If this fails turn to parametric fit
      warning("GAM would not fit, turned to cubic parametric fit ")
      tmp = try(nleqslv(fn =  dNBllcolNP, x = rep(1,4), X = X[,i], reg = model.matrix(~sampleScore + I(sampleScore^2) + I(sampleScore^3)), theta = thetas[i], muMarg = muMarg[,i], jac = NBjacobianColNP)$x)

       # try(glm(x~sampleScore + I(sampleScore^2) + I(sampleScore^3), offset = logMu, family = negative.binomial(thetas[i]), etastart = logMu, data = df), silent = TRUE)
    }
    if(class(tmp)[[1]]=="try-error") {
      warning("GLM would not fit either, returning independence model! ")
    tmp = list(coef = rep(0,4), fitted = muMarg[,i]) #If nothing will fit, stick to an independence model
    }
    list(fit = tmp, int = getInt(tmp, sampleScore = sampleScore))
  })
  sumFit = sum(sapply(taxonWise, function(x){length(x$fit)==2}))
  if(verbose && sumFit) warning("A total number of",sumFit, "response functions did not converge! \n")
  samRep = rep(sampleScore, ncols)
  overall = vgam(c(X) ~ s(samRep, df = dfSpline), offset = c(logMu), family = negbinomial.size(lmu = "loge", size = rep(thetas, each = n)), coefstart = coefInitOverall, maxit = vgamMaxit, ...)
  #taxonWiseFitted = sapply(taxonWise, function(x){if(class(x$fit)=="vgam") fitted(x$fit) else x$fit$fitted})
  rowMat = sapply(taxonWise, function(x){if(class(x$fit)=="vgam") predict(x$fit, type ="link") else model.matrix(~sampleScore + I(sampleScore^2) + I(sampleScore^3)) %*% x$fit})
  taxonCoef = lapply(taxonWise, function(x){if(class(x$fit)=="vgam") coef(x$fit) else x$fit})
  #psi = sqrt(mean(sapply(taxonWise, function(x){x$int})^2*colWeights))
  names(taxonWise) = colnames(X)
  list(rowMat = rowMat, rowMatOverall = matrix(predict(overall, type = "link"), ncol = ncols), taxonCoef = taxonCoef, overallCoef = coef(overall), taxonWise = taxonWise)
  #overallFitted = matrix(fitted(overall), ncol = ncols),
}
