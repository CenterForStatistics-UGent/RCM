#' A function to calculate correlation of final scores with the starting values
startFun = function(list){
  X = list$RCM$X
  R = rowSums(X)
  C = colSums(X)
  E = outer(R, C)/sum(R)
  cMat = list$RCM$cMat
  rMat = list$RCM$rMat
  trended.dispersion.ind  <- estimateGLMTrendedDisp(y = t(X), design = NULL, method = "bin.loess",offset=t(log(E)))
  thetas = estDisp(X = X, muMarg = E, cMat = cMat, rMat = rMat, psis = rep(0, ncol(rMat)), trended.dispersion = trended.dispersion.ind)
  svd1 = svd(diag(1/sqrt(R)) %*% (X-E) %*% diag(1/sqrt(C)))
  svd2 = svd(diag(1/R) %*% (X-E) %*% diag(1/C))
  svd3 = svd(diag(R^(-3/2)) %*% (X-E) %*% diag(C^(-3/2)))
  svd4 = svd(diag(R^(-3)) %*% (X-E) %*% diag(C^(-3)))
  svd5 = svd(diag(1/R^2) %*% (X-E) %*% diag(1/C^2))
  svd6 =  svd((X-E)/((outer(R,C)+rowMultiply(outer(R^2, C^2), 1/thetas))/sum(R)))
  list(corRows1 = getCor(svd1$u, rMat), corRows2 = getCor(svd2$u, rMat), corRows3 = getCor(svd3$u, rMat), corRows4 = getCor(svd4$u, rMat), corRows5 = getCor(svd5$u, rMat), corRows6 = getCor(svd6$u, rMat) ,corCols1 = getCor(svd1$v, t(cMat)), corCols2 = getCor(svd2$v, t(cMat)), corCols3 = getCor(svd3$v, t(cMat)), corCols4 = getCor(svd4$v, t(cMat)), corCols5 = getCor(svd5$v, t(cMat)), corCols6 = getCor(svd6$v, t(cMat)))
}
getCor = function(a,b){
  sapply(1:ncol(b), function(x){
    cor(a[,x], b[,x])
  })
}