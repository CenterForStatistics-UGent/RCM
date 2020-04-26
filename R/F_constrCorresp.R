#' Constrained correspondence analysis with adapted powers
#'
#' @param X outcome matrix
#' @param Y constraining matrix
#' @param rowExp,colExp see ?RCM_NB
#' @param muMarg mean matrix under independence model
#'
#' @return a list with eigenvalues, aliased variables and environmentam gradients
#' @details the vegan version, adapted for flexible powers rowExp and colExp
#' @importFrom stats weighted.mean
constrCorresp = function(X, Y, rowExp, colExp,
                         muMarg = outer(rowSums(X), colSums(X))/sum(X)){
  X = X/sum(X)
  RW = rowSums(X)
  CW = colSums(X)
  Xinit = diag(1/RW^rowExp) %*% (X - muMarg/sum(muMarg)) %*% diag(1/CW^colExp)
  ZERO <- sqrt(.Machine$double.eps)
  envcentre <- apply(Y, 2, weighted.mean, w = RW)
  Y <- scale(Y, center = envcentre, scale = FALSE)
  Y <- sweep(Y, 1, sqrt(RW), "*")
  Q <- qr(Y)
  rank <- sum(Q$pivot[seq_len(Q$rank)] > 0)
  if (length(Q$pivot) > Q$rank){
    alias <- colnames(Q$qr)[-seq_len(Q$rank)]
  } else{
    alias <- NULL
  }
  kept <- seq_along(Q$pivot) <= Q$rank & Q$pivot > 0
  Yfit <- qr.fitted(Q, Xinit)
  sol <- svd(Yfit)
  lambda <- sol$d^2
  u <- sol$u
  zeroev <- abs(lambda) < max(ZERO, ZERO * lambda[1L])
  if (any(zeroev)) {
    lambda <- lambda[!zeroev]
    u <- u[, !zeroev, drop = FALSE]
  }
  posev <- lambda > 0
  xx <- Y[, Q$pivot[kept], drop = FALSE]
  bp <- (1/sqrt(colSums(xx^2))) * crossprod(xx, u[, posev, drop = FALSE])
  list(eig = lambda, alias = alias, biplot = bp)
}

