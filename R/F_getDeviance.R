#' A function to extract the deviance from an RCM object
#'
#' @param RCM an RCM object
#' @param Dim The dimensions to use
#'
#' Deviance calculation also entails fitting the saturated model, and its overdispersions. For this goal the zero counts are taken to have very tiny means.
#' Standard dimensions used are only first and second, since these are also plotted
#'
#'@return A matrix with deviances of the same size as the original data matrix