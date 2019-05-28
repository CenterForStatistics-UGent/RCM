#' The influence function for the row scores
#'
#' @param rcm an rcm object
#' @param Dim the required dimension
#'
#' @return A list with components
#' \item{score}{a matrix with components of the score function}
#' \item{InvJac}{A square matrix of dimension n with the components
#'  of the Jacobian related to the row scores}
NBrowInfl = function(rcm, Dim = 1) {
    reg = rcm$psis[Dim] * rcm$cMat[Dim, ]
    mu = extractE(rcm, seq_len(Dim))
    rcm$X = correctXMissingness(rcm$X, mu, rcm$allowMissingness)
    #Take also lower dimensions into account here
    thetaMat = matrix(byrow = TRUE, nrow = nrow(rcm$X),
        ncol = ncol(rcm$X), data = rcm$thetas[,
            switch(as.character(Dim), `0` = "Independence",
                `0.5` = "Filtered", paste0("Dim",
                Dim))])
    lambdaRow = rcm$lambdaRow[seq_k(Dim)]
    rMatK = rcm$rMat[, seq_len(Dim - 1),
        drop = FALSE]
    tmp = if (Dim > 1)
        rcm$lambdaRow[-c(1, 2)] %*% rMatK else 0

    score = reg * (rcm$X - mu)/(1 + mu/thetaMat) +
        c(rcm$rowWeights * (lambdaRow[1] +
            lambdaRow[2] * 2 * rcm$rMat[,
                Dim] + tmp))

    JacobianInv = solve(NBjacobianRow(beta = c(rcm$rMat[,
        Dim], lambdaRow), X = rcm$X, reg = reg,
        thetas = thetaMat, muMarg = mu, k = Dim,
        p = ncol(rcm$X), n = nrow(rcm$X),
        rowWeights = rcm$rowWeights, nLambda = Dim +
            1, rMatK = rMatK))
    # Inverse Jacobian

    # After a long thought: The X's do not
    # affect the estimation of the lambda
    # parameters! Matrix of all influences
    # becomes too large: return score and
    # inverse jacobian
    return(list(score = score, InvJac = JacobianInv[seq_len(nrow(rcm$X)),
        seq_len(nrow(rcm$X))]))
}
