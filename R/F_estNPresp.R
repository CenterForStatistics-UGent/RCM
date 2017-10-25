#' A function to estimate the taxon-wise response functions non-parametrically, and return their fitted values
#'
#' @param sampleScore a vector of length n with environmental scores
#' @param muMarg the offset matrix
#' @param X the n-by-p data matrix
#' @param ncols an integer, the number of columns of X
#' @param psi a scalar, the importance parameters
#' @param thetas a vector of length p with dispersion parameters
#' @param coefInit a 2-by-p matrix with current taxon-wise parameter estimates
#' @param coefInitOveral a vector of length 2 with current overall parameters
#' @param dfSpline a scalar, the degrees of freedom for the smoothing spline.
#' @param bs see ?mgcv::s
#' @param ... further arguments, passed on to the VGAM:::vgam() function
#'
#' The negative binomial likelihood is still maximized, but now the response function is a non-parametric one. To avoid a perfect fit and overly flexible functions, we enforce smoothness restrictions. In practice we use a generalized additive model (GAM), i.e. with splines. A cubic smoothing spline is used, which is suited for small sample sizes (see ?gam::s). The same fitting procedure is carried out ignoring species labels. We do not normalize the parameters related to the splines: the psis can be calculated afterwards. The VGAM::vgam had to be cloned and adapted because of a bug.
#'
#' @return A list with components
#' \item{taxonWiseFitted}{A n-by-p matrix of fitted valuesof the response curves per taxon at the observed values of the environmental scores}
#' \item{taxonCoef}{The fitted coefficients of the sample-wise response curves}
#' \item{overallFitted}{The fitted values of all taxa combined}
#' \item{taxonWiseFits}{A list of length p of normalized response curves of all the taxa. This may be useful to investigate the shape of the response function through plots.}
estNPresp = function(sampleScore, muMarg, X, ncols, psi,thetas, n, coefInit, coefInitOverall, dfSpline, vgamMaxit, bs, ...){
  taxonWise = lapply(seq_len(ncols), function(i){#print(i);
    # vgam2(X[,i] ~ spline - 1, offset = log(muMarg[,i]), family = negbinomial.size(lmu = "loge", size = thetas[i]), coefstart = coefInit[i], maxit = vgamMaxit,...) #No intercept needed: if the spline equals zero there is no departure from independence. Also we have only one additive term
# df = data.frame(x = X[,i], logMu = log(muMarg[,i]), sampleScore = sampleScore)
  # mgcv::gam(data = df, x ~ s(sampleScore)-1, offset = logMu, family = negbin(thetas[i]),start = coefInit[i])
  mgcv::gam(X[,i] ~ s(sampleScore)-1, offset = log(muMarg[,i]), family = negbin(thetas[i]),start = coefInit[i],control = list(maxit=vgamMaxit))
  })
  samRep = rep(sampleScore,ncols)
  overall = vgam2(c(X) ~ VGAM::s(samRep, df = dfSpline), offset = c(log(muMarg)), family = negbinomial.size(lmu = "loge", size = rep(thetas, each = n)), coefstart = coefInitOverall, maxit = vgamMaxit,...) # Here we have to use VGAM to allow for different thetas, which is really undesirable
  # overall = gam(c(X) ~ s(samRep)- 1, offset = c(log(muMarg)), family = negbin(theta = rep(thetas, each = n)), coefstart = coefInitOverall, ...)
  #TO DO: define a new family function negbin2 that accepts multiple overdispersion parameters.
  taxonWiseFitted = sapply(taxonWise, fitted)
  taxonCoef = sapply(taxonWise, coef)
  list(taxonWiseFitted = taxonWiseFitted, taxonCoef = taxonCoef, overallFitted = matrix(fitted(overall), ncol = ncols), overallCoef = coef(overall))
}