#' Jacobian for the estimation of the column scores
#'
#' @param beta vector of length p+1+1+(k-1): p row scores, 1 centering,
#'  one normalization
#'  and (k-1) orhtogonality lagrangian multipliers
#' @param X the nxp data matrix
#' @param reg a nx1 regressor matrix: outer product of rowScores and psis
#' @param thetas nxp matrix with the dispersion parameters
#' (converted to matrix for numeric reasons)
#' @param muMarg the nxp offset
#' @param k an integer, the dimension of the RC solution
#' @param p an integer, the number of taxa
#' @param n an integer, the number of samples
#' @param nLambda an integer, the number of restrictions
#' @param colWeights the weights used for the restrictions
#' @param cMatK the lower dimensions of the colScores
#' @param preFabMat a prefab matrix, (1+X/thetas)
#' @param Jac an empty Jacobian matrix
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#' @return A matrix of dimension p+1+1+(k-1) with evaluations of the Jacobian
NBjacobianCol = function(beta, X, reg, thetas,
    muMarg, k, n, p, colWeights, nLambda,
    cMatK, preFabMat, Jac, allowMissingness, naId) {
    cMat = beta[seq_len(p)]

    # Calculate the mean
    mu = exp(reg %*% cMat) * muMarg

    # The symmetric jacobian matrix. The
    # upper part is filled first, then mirror
    # image is taken for lower triangle

    # dLag²/dr_{ik}dlambda_{1k}
    Jac[seq_len(p), p + 2] = Jac[p + 2, seq_len(p)] = colWeights *
        2 * cMat

    if(allowMissingness){
        preFabMat = 1 + correctXMissingness(X, mu, allowMissingness, naId)/thetas
    }

    # dLag²/dr_{ik}²
    diag(Jac)[seq_len(p)] = -crossprod(preFabMat *
        mu/(1 + mu/thetas)^2, reg^2) + 2 *
        beta[p + 2] * colWeights
    Jac
}
