#' A function to calculate the matrix of deviance residuals from a RCM object.
#'
#' @param RCM an RCM object
#' @param Dim The dimensions to use
#'
#' For the deviance residuals we use the overdispersions from the reduced model.
#' Standard dimensions used are only first and second, since these are also plotted
#'
#'@return A matrix with deviance residuals of the same size as the original data matrix
getDevianceRes = function(RCM, Dim = seq_len(RCM$k)){
  mu = extractE(RCM,Dim)
  thetaMat = extractDisp(RCM, mu)
  getDevMat(X = RCM$X, thetaMat = thetaMat, mu = mu)
}
