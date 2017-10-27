#' A function that returns the ratio of the contributions of the signal taxa to the non-signal taxa to the separation of the clusters
#'
#' @param rowMatPsi A n-by-k matrik of row coordinates
#' @param colMat a p-by-k matrix of
#' @param groupFactor A factor of length n indicating true group membership
#' @param idTax A matrix indicating which taxa are up and downregulated in the different clusters
#' @param sigFrac a scalar, which fraction of the non-signal taxa should be considered, see Details
#' @param Dim a vector indicating which dimensions to consider
#' @param upDown a boolean, are taxa both up- and downregulated?
#' @param centerFun a function that determines which summary statistic the contribution of the two groups undergoes before taking a ratio
#'
#' We only consider the upper half of non-signal taxa to avoid dividing by very small numbers engendered by non-contributing taxa. Using the median as summary statistic provids robustness to outliers. The geometric mean would be more suitable to summarize ratios, but this fails since some median contributions are negative.
#'
#' @return the mean of the ratios of median contributions of signal to the sigFrac most contributing non-signal taxa
contrTaxa = function(rowMatPsi, colMat, groupFactor, idTax, sigFrac = 0.5, Dim = 1:3, upDown = FALSE, centerFun = median){
  #FIX ME! Hellinger!
  rowMatPsi = rowMatPsi[,Dim]
  colMat = t(colMat)[Dim,]
  if(upDown){
    idTaxNum = matrix(-1,nrow(idTax) , ncol(idTax))
    idTaxNum[idTax=="up"] = 1
  }
  centroids = matrix(unlist(tapply(seq_len(nrow(rowMatPsi)), groupFactor, function(i){
    colMeans(rowMatPsi[i,])
  })), nrow = ncol(rowMatPsi)) #The group centroids
  colnames(centroids) = levels(groupFactor)
  overalCenterSam = colMeans(rowMatPsi) # The overall sample centroid(not always c(0,0,0)!!)
  overalCenterTax = colMeans(colMat)
  contrib = sapply(seq_along(levels(groupFactor)), function(i){
    crossprod(colMat, centroids[, i])
  }) # The inner product of all species' arrows with the vectors locating the centroids
  sigContr = if(upDown) {sapply(seq_along(levels(groupFactor)), function(i){
    centerFun(contrib[match(rownames(idTax), colnames(colMat), nomatch = 0),i]*idTaxNum[,i])
  })} else {
    sapply(seq_along(levels(groupFactor)), function(i){
      centerFun(contrib[match(idTax[,i], colnames(colMat), nomatch = 0),i])
    })
  } # The mean/median contribution of the taxa with signal
  allContr = if(upDown) {sapply(seq_along(levels(groupFactor)), function(i){
    centerFun(
      apply(contrib[!colnames(colMat) %in% rownames(idTax),],2,sort, decreasing = TRUE)[seq_len(round(sigFrac*(ncol(colMat)-length(idTax[[i]])))), ])
  })} else {
    sapply(seq_along(levels(groupFactor)), function(i){
      centerFun( apply(contrib[!colnames(colMat) %in% idTax[,i],],2,sort, decreasing = TRUE)[seq_len(round(sigFrac*(ncol(colMat)-length(idTax[[i]])))), ])
    })
  } # The mean/median contribution of the 50% taxa with highest contribution among non-significant taxa
  ratios = sigContr/allContr # The ratio of these mean contributions
  mean(ratios) # mean ratio
}