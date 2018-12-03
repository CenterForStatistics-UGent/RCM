#' The influence function for the column scores
#' @param rcm an rcm object
#' @param Dim the required dimension
#'
#' @return A list with components
#' \item{score}{a matrix with components of the score function}
#' \item{InvJac}{A square matrix of dimension p with the components of the
#'  Jacobian related to the column scores}
NBcolInfl = function(rcm, Dim = 1) {
    reg = rcm$psis[Dim] * rcm$rMat[, Dim]
    mu = extractE(rcm, seq_len(Dim))  #Take also lower dimensions into account here
    thetaMat = matrix(byrow = TRUE, nrow = nrow(rcm$X), 
        ncol = ncol(rcm$X), data = rcm$thetas[, 
            switch(as.character(Dim), `0` = "Independence", 
                `0.5` = "Filtered", paste0("Dim", 
                  Dim))])
    lambdaCol = rcm$lambdaCol[seq_k(Dim)]
    cMatK = rcm$cMat[seq_len(Dim - 1), , 
        drop = FALSE]
    tmp = if (Dim > 1) 
        lambdaCol[-c(1, 2)] %*% cMatK else 0
    
    score = (reg * (rcm$X - mu)/(1 + mu/thetaMat)) + 
        rcm$colWeights * (lambdaCol[1] + 
            lambdaCol[2] * 2 * rcm$cMat + 
            tmp)
    
    JacobianInv = solve(NBjacobianCol(beta = c(rcm$cMat[Dim, 
        ], lambdaCol), X = rcm$X, reg = reg, 
        thetas = thetaMat, muMarg = mu, k = Dim, 
        p = nrow(rcm$X), n = ncol(rcm$X), 
        colWeights = rcm$colWeights, nLambda = length(lambdaCol), 
        cMatK = cMatK))
    # Inverse Jacobian
    
    # After a long thought: The X's do not
    # affect the estimation of the lambda
    # parameters! Matrix becomes too large:
    # return score and inverse jacobian
    return(list(score = score, InvJac = JacobianInv[seq_len(ncol(rcm$X)), 
        seq_len(ncol(rcm$X))]))
}
