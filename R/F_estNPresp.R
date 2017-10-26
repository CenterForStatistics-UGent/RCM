#' A function to estimate the taxon-wise response functions non-parametrically, and return their fitted values
#'
#' @param sampleScore a vector of length n with environmental scores
#' @param muMarg the offset matrix
#' @param X the n-by-p data matrix
#' @param ncols an integer, the number of columns of X
#' @param thetas a vector of length p with dispersion parameters
#' @param coefInit a 2-by-p matrix with current taxon-wise parameter estimates
#' @param coefInitOveral a vector of length 2 with current overall parameters
#' @param degreeSpline a scalar, the degrees of freedom for the smoothing spline.
#' @param colWeights a vector of length p with column weights
#' @param ... further arguments, passed on to the VGAM:::vgam() function
#'
#' The negative binomial likelihood is still maximized, but now the response function is a non-parametric one. To avoid a perfect fit and overly flexible functions, we enforce smoothness restrictions. In practice we use a generalized additive model (GAM), i.e. with splines.  The same fitting procedure is carried out ignoring species labels. We do not normalize the parameters related to the splines: the psis can be calculated afterwards. The VGAM package turned ut to be very buggy, and chose the mgcv package instead. Unfortunately we can only use this with a single parameter, so we have to turn to VGAM again for the overall response function. Fitting both models with different packages is of course not desirable.
#'
#' @return A list with components
#' \item{taxonWiseFitted}{A n-by-p matrix of fitted valuesof the response curves per taxon at the observed values of the environmental scores}
#' \item{taxonCoef}{The fitted coefficients of the sample-wise response curves}
#' \item{overallFitted}{The fitted values of all taxa combined}
#' \item{taxonWiseFits}{A list of length p of normalized response curves of all the taxa. This may be useful to investigate the shape of the response function through plots.}
#' \item{psi}{The importance parameter of the dimension}
estNPresp = function(sampleScore, muMarg, X, ncols, thetas, n, coefInit, coefInitOverall, dfSpline, vgamMaxit, colWeights, ...){
  taxonWise = lapply(seq_len(ncols), function(i){#print(i);
    df = data.frame(x = X[,i], logMu = log(muMarg[,i]), sampleScore = sampleScore)
    # tmp = try(vgam2(data = df, x ~ VGAM::s(sampleScore, df = 4), offset = logMu, family = negbinomial.size(lmu = "loge", size = ), coefstart = CoefInit, maxit = vgamMaxit,...),silent=TRUE) #No intercept needed: if the spline equals zero there is no departure from independence. Also we have only one additive term
    # if(class(tmp)=="try-error") {
      tmp = try(vgam2(data = df,x ~ bs(sampleScore, degree = 3), offset = logMu, family = negbinomial.size(lmu = "loge", size = thetas[i]), coefstart = NULL, maxit = vgamMaxit,...))
    # }
    if(class(tmp)=="try-error") { #If still fails turn to parametric fit
      warning("GAM would not fit, turned to cubic parametric fit ")
      tmp = glm.fit(y = X[,i], x = model.matrix(~sampleScore + I(sampleScore^2) + I(sampleScore^3) ), offset = log(muMarg[,i]), family = do.call("negative.binomial", list(theta = thetas[i], link = "log"))) #Fix fitting issues!
    }
    list(fit = tmp, int = getInt(tmp, sampleScore = sampleScore))
    #Try other kinds of splines?
# #
#   # mgcv::gam(data = df, x ~ s(sampleScore)-1, offset = logMu, family = negbin(thetas[i]),start = coefInit[i])
#   mgcv::gam(X[,i] ~ s(sampleScore)-1, offset = log(muMarg[,i]), family = negbin(thetas[i]),start = coefInit[,i],control = list(maxit=vgamMaxit))
  })
  samRep = rep(sampleScore,ncols)
  overall = vgam2(c(X) ~ bs(samRep, degree = 3), offset = c(log(muMarg)), family = negbinomial.size(lmu = "loge", size = rep(thetas, each = n)), coefstart = coefInitOverall, maxit = vgamMaxit,...) # Here we have to use VGAM to allow for different thetas, which is really undesirable
  # overall = mgcv::gam(c(X) ~ s(samRep)- 1, offset = c(log(muMarg)), family = negbin(theta = rep(thetas, each = n)), coefstart = coefInitOverall, ...)
  #TO DO: define a new family function negbin2 that accepts multiple overdispersion parameters.
  taxonWiseFitted = sapply(taxonWise, function(x){fitted(x$fit)})
  taxonCoef = sapply(taxonWise, function(x){coef(x$fit)})
  psi = sqrt(sum(sapply(taxonWise, function(x){x$int})^2*colWeights))
  list(taxonWiseFitted = taxonWiseFitted, taxonCoef = taxonCoef, overallFitted = matrix(fitted(overall), ncol = ncols), overallCoef = coef(overall), psi = psi)
}