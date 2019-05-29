#' The score function of the general response function
#'
#' @param betas a vector of regression parameters with length v
#' @param X the nxp data matrix
#' @param reg a matrix of regressors of dimension nxv
#' @param thetasMat A matrix of dispersion parameters
#' @param muMarg offset matrix of dimension nxp
#' @param psi a scalar, the importance parameter
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#' @param ... further arguments passed on to the jacobian
#'
#' @return The evaluation of the score functions (a vector length v)
dNBllcol_constr_noLab = function(betas, X,
    reg, thetasMat, muMarg, psi, allowMissingness, naId, ...) {
    mu = c(exp(reg %*% betas * psi)) * muMarg
    X = correctXMissingness(X, mu, allowMissingness, naId)
    colSums(crossprod((X - mu)/(1 + (mu/thetasMat)),
        reg) * psi)
}
