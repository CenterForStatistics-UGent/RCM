#' Replace missing entries in X by their expectation to set their
#' contribution to the estimating equations to zero
#' @param X the matrix of counts
#' @param mu the matrix of expectations
#'
#' @return The matrix X with the NA entries replaced by the
#' corresponding entries in mu
#'
#' @note This may seem like a hacky approach, but it avoids having to deal
#'  with NAs in functions like crossprod().
correctXMissingness = function(X, mu, allowMissingness){
    if(allowMissingness){
        naId = is.na(X)
        X[naId] = mu[naId]
    }
    X
}
