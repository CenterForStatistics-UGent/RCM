#' Calculate the log-likelihoods of all possible models
#'
#'@param rcm an object of the RCM class
#'@param Sum a boolean, should log-likelihoods be summed?
#'
#'@return If Sum is FALSE, a named array log-likelihoods
#' of the independence model and all models with dimension 1 to k,
#' including after filtering on confounders.
#' Otherwise a table with log-likelihoods,
#' deviance explained and cumulative deviance explained.
#'@export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, round = TRUE)
#' liks(zellerRCM)
liks = function(rcm, Sum = TRUE){
  vec = if(length(rcm$confounders$confounders)) {c(0,0.5, seq_len(rcm$k), Inf)
    } else c(0:rcm$k, Inf)
  outnames = c("independence",
               if(length(rcm$confounders$confounders)) "filtered" else NULL,
               paste0("Dim", seq_len(rcm$k)),"saturated")
  if(Sum) {
    tmp = vapply(FUN.VALUE = numeric(1), vec, FUN = function(i){
      sum(getLogLik(rcm, i))
    })
    names(tmp) = outnames
  } else {
  tmp = vapply(vec, FUN.VALUE = matrix(0, nrow(rcm$X), ncol(rcm$X)),
               FUN = function(i){
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
