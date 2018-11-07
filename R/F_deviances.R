#' A function to extract deviances for all dimension, including after filtering on confounders
#'
#'@param rcm an object of the RCM class
#'@param squaredSum a boolean, should total deviance be returned?
#'
#'Dispersions are re-estimated for every dimension of the model. Total deviances can be deceptive and not correspond to the differences in log-likelihood. as the dispersion is different for each model. To compare models it is better to compare likelihoods.
#'
#'@return If Sum is FALSE, a named array of deviance residuals of the independence model and all models with dimension 1 to k, including after filtering on confounders. Otherwise a table with total deviances (the sum of squared deviance residuals), deviance explained and cumulative deviance explained.
deviances = function(rcm, squaredSum = FALSE){
  vec = if(length(rcm$confounders)) c(0,0.5, seq_len(rcm$k)) else c(0:rcm$k)
  outnames = c("independence", if(length(rcm$confounders)) "filtered" else NULL, paste0("Dim ", 1:rcm$k))
  if(squaredSum) {
    tmp = sapply(vec, FUN = function(i){
      sum(getDevianceRes(rcm, i)^2)
    })
    names(tmp) = outnames
    #Also make cumulative comparisons
    cumDevianceExplained = round((tmp-tmp[1])/(tmp[length(tmp)]-tmp[1]),3)
    out = rbind(deviance = tmp,
                devianceExplained = c(0, diff(cumDevianceExplained)),
                cumDevianceExplained = cumDevianceExplained)
  } else {
    out = lapply(vec, function(i){
      getDevianceRes(rcm, i)
    })
    names(out) = outnames
  }
  return(out)
}
