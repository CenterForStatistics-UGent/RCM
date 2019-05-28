#' The jacobian of the response function without taxon labels
#'
#' @param betas a vector of regression parameters with length v
#' @param X the nxp data matrix
#' @param reg a matrix of regressors of dimension nxv
#' @param thetasMat A matrix of dispersion parameters
#' @param muMarg offset matrix of dimension nxp
#' @param preFabMat a prefabricated matrix
#' @param psi a scalar, the importance parameter
#' @param n an integer, number of rows of X
#' @param v an integer, the number of parameters of the response function
#'
#' @return The jacobian (a v-by-v matrix)
JacCol_constr_noLab = function(betas, X,
    reg, thetasMat, muMarg, psi, n, v, preFabMat, allowMissingness) {
    mu = c(exp(reg %*% betas * psi)) * muMarg
    if(allowMissingness){
      preFabMat = 1 + correctXMissingness(X, mu, allowMissingness)/thetasMat
    }
    tmp = preFabMat * mu/(1 + (mu/thetasMat))^2 *
        psi^2  #Don't forget to square psi!
    -crossprod(reg, vapply(seq_len(v), FUN.VALUE = vector("numeric",
        n), function(x) {
        rowSums(reg[, x] * tmp)
    }))
}
