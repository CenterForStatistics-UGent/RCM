#' Trim based on confounders to ensure there are no taxa with only zero counts in one of the subgroups defined by the confounders.
#'
#' @param confounders a nxt confounder matrix
#' @param X the nxp data matrix
#' @param prevCutOff a scalar between 0 and 1, the prevalence cut off
#' @param minFraction a scalar between 0 and 1, each taxon's total abundance should equal at least the number of samples n times minFraction, otherwise it is trimmed
#' @param n: the number of samples
#'
#' Should be called prior to fitting the independence model
#'
#' @return A trimmed data matrix nxp'
trimOnConfounders = function(confounders, X, prevCutOff, minFraction, n){
  trimmingID = apply(X, 2, function(x){ #Over taxa
    any(apply(confounders, 2, function(conf){ #Over confounding variables
      tapply(X = x, INDEX = conf, FUN = function(y){mean(y!=0) <= prevCutOff | sum(y)<(n*minFraction)}) #Any all-zero subgroup?
    }))
  })

  if(all(trimmingID)){
    stop("All taxa would be trimmed, please provide a covariate with less levels! \n")
  }

  X[, !trimmingID] #Return trimmed X
}