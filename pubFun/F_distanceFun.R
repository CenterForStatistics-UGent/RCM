#' A function to calculate the within cluster and overall distances and the geometric mean of the ratio between them
#'
#' @param coord an nxk matrix of coordinates
#' @param clusters a factor of length n with cluster memberships
#' @param Dim The dimensions to look at

#' @return a list with
#' \item{withinDist} {within distances}
#' \item{overalDist} {the overal distances}
#' \item{ratio}{The geometric mean of the ratio between the two}
distanceFun = function(coord, clusters, Dim = 1:3){

  coord = coord[,Dim, drop=FALSE]

  withinDist = tapply(seq_len(nrow(coord)), clusters, function(x){
    coordsX = coord[x,]
    centerX = colMeans(coordsX)
    dists = apply(coordsX, 1, function(y){
      dist(rbind(y, centerX))
    })
    return(dists)
  })

  overalCenter = colMeans(coord)
  overalDist = apply(coord, 1, function(y){
    dist(rbind(y, overalCenter))
  })

  ratio = exp(mean(log(unlist(sapply(unique(clusters), function(x){overalDist[clusters==x]/withinDist[[x]]}))))) #The geometric mean
  list(withinDist = withinDist, overalDist =  overalDist, ratio = ratio)
}