#' A function to extract the column vectors. For PCoA and CoDa this uses weighted averaging
#'
#' @param x a list of biplot solutions
#' @param row a list with the corresponding row scores
#'
#' @return a list with corresponding column scores
extractCols = function(x, row){
  tmpList = with(x, list(
    RCM = t(RCM$cMat),
    CApearson = CA$v,
    CAchisq = CA$v,
    CAcontRat = diag(1/sqrt(colSums(RCM$X))) %*% CA$v,
    DCA = DCA$cproj,
    CoDa = CoDa$colScores,
    BC = weightedTaxonScores(RCM$X, row$BC),
    JSD = weightedTaxonScores(RCM$X, row$JSD),
    BCrel = weightedTaxonScores(RCM$X, row$BCrel),
    UniFrac = if(!is.null(UniFrac))  {weightedTaxonScores(RCM$X, row$UniFrac)} ,
    wUniFrac = if(!is.null(wUniFrac)) {weightedTaxonScores(RCM$X, row$wUniFrac)} ,
    DPCOA = if(is.null(DPCOA)) NULL else DPCOA$dls,
    BCrelNMDS = BCrelNMDS$species,
    Hellinger = diag(1/diag((crossprod(sqrt(RCM$X/rowSums(RCM$X))-sqrt(colSums(RCM$X)/sum(RCM$X)), diag(rowSums(RCM$X))) %*% (sqrt(RCM$X/rowSums(RCM$X))-sqrt(colSums(RCM$X)/sum(RCM$X)))))) %*% Hellinger$v
  ))
  tmpList = tmpList[!sapply(tmpList, is.null)]
  rowNames = rownames(tmpList$RCM)
  lapply(tmpList, function(x){
    colnames(x) = NULL
    if(nrow(x)==length(rowNames)) rownames(x) = rowNames
    x
  })
}