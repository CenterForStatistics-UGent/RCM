#' A function to calculate the matrix of deviance residuals from a RCM object.
#'
#' @param RCM an RCM object
#' @param Dim The dimensions to use
#'
#' For the deviance residuals we use the overdispersions from the reduced model.
#' Standard dimensions used are only first and second, since these are also plotted
#'
#'@return A matrix with deviance residuals of the same size as the original data matrix
getDevianceRes = function(RCM, Dim = c(1,2)){
  mu = extractE(RCM,Dim)
  thetaMat = extractDisp(RCM, mu)
  tmpMat = sqrt(2*(RCM$X*log(RCM$X/mu)-(RCM$X+thetaMat)*log((1+RCM$X/thetaMat)/(1+mu/thetaMat))))*sign(RCM$X-mu)
  tmpMat[RCM$X==0] = -sqrt((2*thetaMat*log(1+mu/thetaMat))[RCM$X==0]) #zero observations are always smaller than the mean
  tmpMat
}