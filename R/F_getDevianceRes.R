#' A function to calculate the matrix of deviance residuals from a RCM object.
#'
#' @param RCM an RCM object
#' @param Dim The dimensions to use
#'
#' For the deviance residuals we use the overdispersions from the reduced model.
#' Standard dimensions used are only first and second, since these are also plotted
#'
#'@return A matrix with deviance residuals of the same size as the original data matrix
getDevianceRes = function(RCM, Dim = RCM$k){
  mu = extractE(RCM,Dim)
  thetaMat = matrix(byrow = TRUE, nrow = nrow(RCM$X), ncol = ncol(RCM$X),
                    data = RCM$thetas[,switch(as.character(Dim), "0" = "Independence", "0.5" = "Filtered", paste0("Dim", Dim))])
  getDevMat(X = RCM$X, thetaMat = thetaMat, mu = mu)
}
