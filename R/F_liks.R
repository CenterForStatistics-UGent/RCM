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
  outnames = c("independence", paste0("Dim ", 1:rcm$k))#,"saturated")
  if(Sum) {
    tmp = sapply(c(0:rcm$k), FUN = function(i){
      sum(getLogLik(rcm, i))
    })
    names(tmp) = outnames
  } else {
  tmp = vapply(c(0:rcm$k), FUN.VALUE = matrix(0, nrow(rcm$X), ncol(rcm$X)), FUN = function(i){
getLogLik(rcm, i)
  })
  dimnames(tmp)[[3]] = outnames
  }
  out = c(tmp[1], filtered = if(Sum &&!is.null(rcm$llFilt)) sum(rcm$llFilt) else rcm$llFilt, tmp[-1])
  if(Sum){ #Also make cumulative comparisons
    cumDevianceExplained = round((out-out[1])/(out[length(out)]-out[1]),3)
  out = rbind(logLikelihood = out,
              devianceExplained = c(0, diff(cumDevianceExplained)),
              cumDevianceExplained = cumDevianceExplained)

  }
  return(out)
}
