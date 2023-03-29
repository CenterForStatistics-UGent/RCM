#' A function to calculate the pseudo F-statistic
#'
#' @param coord an nxk matrix of coordinates
#' @param clusters a factor of length n with cluster memberships

#' @return The pseudo F-statistic
pseudoF = function(coord, clusters){

  N = nrow(coord)
  a = length(unique(clusters))
  overalDist = sum(dist(coord)^2)/N
  withinDist = sum(unlist(tapply(seq_len(nrow(coord)), clusters, function(x){
    dist(coord[x,])^2/length(x)
  })))

  Fratio = (overalDist-withinDist)/withinDist * (a-1)/(N-a)
  return(Fratio)
}