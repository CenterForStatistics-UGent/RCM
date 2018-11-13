#' The influence function for the psis
#'
#' @param rcm an rcm object
#' @param Dim the required dimensions
#'
#' @return The influence of every single observation on the psi value of this dimension
NBpsiInfl = function(rcm, Dim = 1){
  mu = extractE(rcm, seq_len(Dim)) #Take also lower dimensions into account here
  thetaMat = matrix(byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X),
                    data = rcm$thetas[,switch(as.character(Dim), "0" = "Independence", "0.5" = "Filtered", paste0("Dim", Dim))])
  reg = if(is.null(rcm$covariates)) rcm$rMat[,Dim] %*% rcm$cMat[Dim,] else getRowMat(sampleScore = rcm$alpha[,Dim] %*% rcm$covariates, responseFun = rcm$responseFun, NB_params = rcm$NB_params[,,Dim])
  -((rcm$X-mu)*(thetaMat+mu))/(reg*(thetaMat+rcm$X)*mu)
}
