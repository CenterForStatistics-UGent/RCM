#' Calculate the log-likelihoods of the indepence and saturated models and all fitted submodels
#'
#'@param rcm an object of the RCM class
#'@param Sum a boolean, should log-likelihoods be summed?
#'
#'Dispersions are re-estimated for every dimension of the model.
#'
#'@return If Sum is FALSE, a named array log-likelihoods of the independence model and all models with dimension 1 to k, including after filtering on confounders. Otherwise a table with log-likelihoods, deviance explained and cumulative deviance explained.
#'@export
liks = function(rcm, Sum = TRUE){
  vec = if(length(rcm$confounders$confMat)) c(0,0.5, seq_len(rcm$k)) else c(0:rcm$k)
  outnames = c("independence", if(length(rcm$confounders$confMat)) "filtered" else NULL, paste0("Dim", 1:rcm$k))#,"saturated")
  if(Sum) {
    tmp = sapply(vec, FUN = function(i){
      sum(getLogLik(rcm, i))
    })
    names(tmp) = outnames
  } else {
  tmp = vapply(vec, FUN.VALUE = matrix(0, nrow(rcm$X), ncol(rcm$X)), FUN = function(i){
getLogLik(rcm, i)
  })
  dimnames(tmp)[[3]] = outnames
  }
  if(Sum){ #Also make cumulative comparisons
    cumDevianceExplained = round((tmp-tmp[1])/(tmp[length(tmp)]-tmp[1]),3)
  out = rbind(logLikelihood = tmp,
              logLikExplained = c(0, diff(cumDevianceExplained)),
              cumLogLikExplained = cumDevianceExplained)
  }
  return(out)
}
