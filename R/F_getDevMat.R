#' A function to calculate the matrix of deviance residuals for given mean and dispersion matrix
#'
#'@param X the data matrix
#'@param thetaMat the matrix of dispersions
#'@param mu the matrix of means
#'
#'@return The matrix of deviance residuals
getDevMat = function(X, thetaMat, mu){
  tmpMat = suppressWarnings(sqrt(2*(X*log(X/mu)-(X+thetaMat)*log((1+X/thetaMat)/(1+mu/thetaMat))))*sign(X-mu))
  tmpMat[X==0] = -sqrt((2*thetaMat*log(1+mu/thetaMat))[X==0]) #zero observations are always smaller than the mean
  tmpMat
}
