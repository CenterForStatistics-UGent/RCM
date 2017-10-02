#' A function to estimate the taxon-wise response functions non-parametrically, and return their fitted values
#'
#' @param sampleScore a vector of length n with environmental scores
#' @param muMarg the offset matrix
#' @param X the n-by-p data matrix
#' @param ncols an integer, the number of columns of X
#' @param psi a scalar, the importance parameters
#' @param ... further arguments, passed on to the locfit.raw() function
#'
#' The density of each species along the environmental scores is fitted non-parametrically, and this curve is normalized to sum to 1 using the normCureve() function. Hereby the ratio of observed values to the offset is modelled with a log-link with the environmental score as explanatory variable. The same fitting procedure is carried out ignoring species labels.
#'
#' @return A list with components
#' \item{taxonWise}{A n-by-p matrix of fitted valuesof the response curves per taxon at the observed values of the environmental scores}
#' \item{overall}{The normalized response curve of all taxa combined}
#' \item{taxonWiseFits}{A list of length p of normalized response curves of all the taxa. This may be useful to investigate the shape of the response function through plots.}
estNPresp = function(sampleScore, muMarg, X, ncols, psi, ...){
  taxonWiseFits = lapply(seq_len(ncols), function(i){
    locfit.raw(x = sampleScore, y = X[,i]/muMarg[,i], link = "log",...)
  })
  overall = matrix(normCurve(locfit.raw(x = rep(sampleScore,ncols), y = c(X/muMarg), link = "log",...)),ncol = ncols)
  taxonWise = sapply(taxonWiseFits, normCurve)
  list(taxonWise = taxonWise, overall = matrix(overall,ncol = ncols), taxonWiseFits = taxonWiseFits)
}