#' Calculates the Jacobian of the parametric response functions
#'
#' @param betas a vector of length (deg+1)*(p+1) with regression parameters
#'  with deg the degree of the response function and the lagrangian multipliers
#' @param X the nxp data matrix
#' @param reg a vector of regressors with the dimension n-by-v
#' @param thetaMat The n-by-p matrix with dispersion parameters
#' @param muMarg offset matrix of size nxp
#' @param psi a scalar, the importance parameter
#' @param v an integer, one plus the degree of the response function
#' @param p an integer, the number of taxa
#' @param IDmat an logical matrix with indices of non-zero elements
#' @param IndVec  a vector with indices with non-zero elements
#'
#' @return The jacobian, a square matrix of dimension (deg+1)*(p+1)

respFunJacMat = function(betas, X, reg, thetaMat, 
    muMarg, psi, v, p, IDmat, IndVec) {
    NBparams = matrix(betas[seq_len(p * v)], 
        ncol = p)
    mu = exp(reg %*% NBparams * psi) * muMarg
    Jac = matrix(0, (p + 1) * v, (p + 1) * 
        v)
    did = seq_len(p * v)
    didv = seq_len(v)
    # d²Lag/dlambda dBeta
    Jac[didv + p * v, did][IndVec] = 2 * 
        NBparams
    Jac = Jac + t(Jac)  #symmetrize
    
    tmp = (1 + X/thetaMat) * mu/(1 + mu/thetaMat)^2
    tmp2 = vapply(didv, FUN.VALUE = tmp, 
        function(x) {
            reg[, x] * tmp
        })
    # d²Lag/dBeta²
    Jac[did, did][IDmat] = -aperm(tensor(reg, 
        tmp2, 1, 1), c(3, 1, 2)) * psi^2
    # Permute the dimensions to assure
    # correct insertion
    
    diag(Jac)[did] = diag(Jac)[did] + 2 * 
        betas[seq_len(v) + p * v]
    return(Jac)
}
