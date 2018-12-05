#' An auxiliary R function to 'array' multiply an array with a vector,
#'  kindly provided by Joris Meys
#'
#' @param x a axbxc array
#' @param y a vector of length c
#'
#' @return a axb matrix. The ij-th element equals sum(x[i,j,]*y)
arrayprod <- function(x, y) {
    
    xdim <- dim(x)
    outdim <- xdim[c(1, 2)]
    outn <- prod(outdim)
    
    yexpand <- rep(y, each = outn)
    
    tmp <- x * yexpand
    
    dim(tmp) <- c(outn, xdim[3])
    out <- rowSums(tmp)
    
    dim(out) <- outdim
    
    out
}
