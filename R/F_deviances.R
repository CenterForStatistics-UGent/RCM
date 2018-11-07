#' A function to extract deviances for all dimension, including after filtering on confounders
#'
#'@param rcm an object of the RCM class
#'@param squaredSum a boolean, should total deviance be returned?
#'
#'Dispersions are re-estimated for every dimension of the model.
#'
#'@return If Sum is FALSE, a named array of deviance residuals of the independence model and all models with dimension 1 to k, including after filtering on confounders. Otherwise a table with total deviances (the sum of squared deviance residuals), deviance explained and cumulative deviance explained.
#'@export
deviances = function(rcm, squaredSum = TRUE){
  vec = if(length(rcm$confounders)) c(0,0.5, seq_len(rcm$k)) else c(0:rcm$k)
  outnames = c("independence", if(length(rcm$confounders)) "filtered" else NULL, paste0("Dim ", 1:rcm$k))
  if(squaredSum) {
    tmp = sapply(vec, FUN = function(i){
      sum(getDevianceRes(rcm, i)^2)
    })
    names(tmp) = outnames
    out = c(tmp[1], filtered = if(is.null(rcm$devFilt)) NULL else sum(rcm$devFilt^2), tmp[-1])
    #Also make cumulative comparisons
    cumDevianceExplained = round((out-out[1])/(out[length(out)]-out[1]),3)
    out = rbind(deviance = out,
                devianceExplained = c(0, diff(cumDevianceExplained)),
                cumDevianceExplained = cumDevianceExplained)
  } else {
    tmp = lapply(vec, function(i){
      getDevianceRes(rcm, i)
    })
    names(tmp) = outnames
  }
  return(out)
}
