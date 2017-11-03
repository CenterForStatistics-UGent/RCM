#' An auxiliary function to get weighted taxon scores based on row scores. Mainly used for PCoA, not an advisable way to make a biplot
#'
#' @param X the data matrix
#' @param rowScores the row scores
#'
#' @return Taxon scores, that are sums of row scores weighted by observed abundances
weightedTaxonScores = function(X, rowScores){
  crossprod(X %*% diag(1/colSums(X)),rowScores)
}