#' Enforces squareness off a ggplot
#'
#' @param plt a ggplot object
#'
#' @return a ggplot object, squared
squarePlot <- function(plt){
  return(plt+coord_fixed()+
           expand_limits(x=ggplot_build(plt)$panel$ranges[[1]]$y.range,
                         y=ggplot_build(plt)$panel$ranges[[1]]$x.range))
}