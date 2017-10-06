#' Calculate the log-likelihoods of the indepence and saturated models and all fitted submodels
#'
#'@param rcm an object of the RCM class
#'
#'Dispersions are re-estimated for every dimension of the model.
#'
#'@return a named vector of length rcm$k+1 containing the log-likelihoods of the independence model and all models with dimension 1 to k
liks = function(rcm){
  arrayOut =array("numeric",c(dim(rcm$X), rcm$k+2))

  arrayOut[,,-(rcm$k+2)] = sapply(c(0,seq_len(rcm$k)), function(k){
    exList = extractE(rcm, k)
   dnbinom(x = rcm$X, mu = exList$E, size = exList$thetaMat, log = TRUE)
  })
  Esat = rcm$X
  Esat[Esat==0] = 1e-300 #For a saturated model, set zero values to 1e-300
  arrayOut[,,rcm$k+2] = dnbinom(x = rcm$X, mu = rcm$X, size = matrix(estDisp(X = rcm$X, muMarg = Esat), byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X)), log = TRUE) #
  dimnames(arrayOut)[[3]] = c("independence", paste0("Dim ", 1:rcm$k),"saturated")
  arrayOut
}