#' Replace missing entries in X by their expectation to set their
#' contribution to the estimating equations to zero
#' @param X the matrix of counts
#' @param mu the matrix of expectations
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#' @return The matrix X with the NA entries replaced by the
#' corresponding entries in mu
#'
#' @note This may seem like a hacky approach, but it avoids having to deal
#'  with NAs in functions like crossprod().
correctXMissingness = function(X, mu, allowMissingness, naId){
    if(allowMissingness){
        X[naId] = mu[naId]
    }
    X
}
