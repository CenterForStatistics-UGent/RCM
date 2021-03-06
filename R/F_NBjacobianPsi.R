#' Jacobian for the psi of a given dimension
#'
#' @param beta a scalar, the current estimate
#' @param X the n-by-p count matrix
#' @param muMarg the nxp offset matrix
#' @param reg the regressor matrix,
#' the outer product of current row and column scores
#' @param theta a n-by-p matrix with the dispersion parameters
#' @param preFabMat a prefab matrix, (1+X/thetas)
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#' @return The evaluation of the jacobian function at beta, a 1-by-1 matrix
NBjacobianPsi = function(beta, X, reg, muMarg,
    theta, preFabMat, allowMissingness, naId) {
    mu = muMarg * exp(reg * beta)
    if(allowMissingness){
        preFabMat = 1 + correctXMissingness(X, mu, allowMissingness, naId)/theta
    }
    matrix(-sum(reg^2 * preFabMat * mu/(1 +
        mu/theta)^2), 1, 1)
}
