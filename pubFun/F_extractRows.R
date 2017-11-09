#' A function to extract the row vectors multiplied by importance elements
#'
#' @param x a list of biplot solutions
#'
#' @return a list of corresponding row scores, matrices of n-by-k
extractRows = function(x){
  tmpList = with(x, list(
    RCM = RCM$rMat %*% diag(RCM$psis),
    CApearson = CA$u %*% diag(CA$d),
    CAcontRat = diag(1/sqrt(rowSums(RCM$X))) %*% CA$u %*% diag(CA$d),
    CAchisq  = diag(1/sqrt(rowSums(RCM$X))) %*% CA$u %*% diag(CA$d),
    DCA = DCA$rproj,
    CoDa = CoDa$rowScores,
    BC = BC$vectors %*% diag(BC$values$Eigenvalues[seq_len(ncol(BC$vectors))]),
    JSD = JSD$vectors %*% diag(JSD$values$Eigenvalues[seq_len(ncol(JSD$vectors))]),
    BCrel = BCrel$vectors %*% diag(BCrel$values$Eigenvalues[seq_len(ncol(BCrel$vectors))]),
    BCrelNMDS = BCrelNMDS$points,
    UniFrac = if(is.null(UniFrac)) NULL else UniFrac$vectors %*% diag(UniFrac$values$Eigenvalues[seq_len(ncol(UniFrac$vectors))]),
    wUniFrac = if(is.null(wUniFrac)) NULL else wUniFrac$vectors %*% diag(wUniFrac$values$Eigenvalues[seq_len(ncol(wUniFrac$vectors))]),
    DPCOA = if(is.null(DPCOA)) NULL else DPCOA$li,
    Hellinger = diag(sqrt(1/rowSums(RCM$X))) %*% Hellinger$u %*% diag(Hellinger$d)
  ))
  tmpList = tmpList[!sapply(tmpList, is.null)]
  lapply(tmpList, function(y){
    colnames(y) = NULL;
    rownames(y) = rownames(x$RCM$X);y
  })
}