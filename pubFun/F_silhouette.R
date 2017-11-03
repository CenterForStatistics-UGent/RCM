#' A function to calculate the silhouettes for the different clusters
#'
#' @param coord: an nxk matrix of coordinates
#' @param clusters: a factor of length n with cluster memberships
#' @param Dim: The dimensions to look at
#'
#'The silhouette was defined by Rousseuw in "Silhouettes: A graphical aid to the interpretation and validation of cluster analysis" in 1987
#'
#' @return a vector of n silhouettes
silhouette = function(coord, clusters, Dim = 1:3){
  coord = coord[,Dim, drop=FALSE]

  distMat = as.matrix(dist(coord)) # all distances
  avDists = apply(distMat, 1 , function(x){
    tapply(x[x!=0], clusters[x!=0], mean) #Exclude observation itself
  }) # Average distances within clusters
  ais = sapply(seq_along(clusters), function(i){
    avDists[clusters[i], i]
  }) #Average distance to own cluster
  bis = sapply(seq_along(clusters), function(i){
    min(avDists[clusters[i]!=levels(clusters), i])
  }) #Average distance to nearest cluster

  silhouette = mapply(ais, bis, FUN = function(a,b){(b-a)/max(a,b)}) #The silhouette

  return(silhouette)
}