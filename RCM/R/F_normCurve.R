#' A function to normalize the non-parametric curves to total area 1. Remember only the psi parameter can grow in size.
#'
#' @param fit a non-parametric fit, resulting from the locfit:::locfit.raw() function
#' @param newDataLength an integer, number of anchor points used in the normalization
#'
#' This normalization can be seen as the equivalent of the normalization of the beta parameters of the parametric response functions.
#'
#' @return a vector of fitted values at the observed values of the environmental scores, from a normalized response function
normCurve = function(fit, newDataLength = 200){
  x = seq(min(fit$eva$xev), max(fit$eva$xev), length.out = newDataLength)
  auc = flux::auc(x = x, y = predict(fit, newdata = x))
  fitted(fit)/auc
}