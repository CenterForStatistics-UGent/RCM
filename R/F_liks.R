#' Calculate the log-likelihoods of the indepence models and all fitted submodels
#'
#'@param rcm an object of the RCM class
#'
#'Dispersions are re-estimated for every dimension of the model.
#'
#'@return a named vector of length rcm$k+1 containing the log-likelihoods of the independence model and all models with dimension 1 to k
liks = function(rcm){

  thetaMat = matrix(estDisp(X = rcm$X, muMarg = E), byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
  indLL = sum(dnbinom(x = rcm$X, mu = E, size = thetaMat, log = TRUE)) #Independence model
  rcmLL = sapply(seq_len(rcm$k), function(k){
    exList = extractE(rcm, k)

    sum(dnbinom(x = rcm$X, mu = exList$E, size = exList$thetaMat, log = TRUE))
  })
  names(rcmLL) = paste0("Dim ", 1:rcm$k)
  c("independence" = indLL, rcmLL)
}