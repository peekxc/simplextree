#' Cech complex 
#' @description Computes the Cech complex up to some radius and dimension. 
#' @details The Cech complex at \code{radius} is defined as the abstract simplicial complex which contains a k-simplex
#' for every (k+1) non-empty intersection of balls with size \code{radius}. The Cech complex is equivalently defined
#' as the nerve complex of the cover formed by placing \code{radius}-balls at every point in \code{x}.
#' @export
cech <- function(x, radius = "default", dim = 1L, dx = NULL){
  stopifnot(is.matrix(x))
  if (missing(dx) || is.null(dx)){ dx <- dist(x) }
  if (missing(radius) || radius == "default"){
    radius <- enclosing_radius(dx)/2.0
  }
  stopifnot(is.numeric(radius), "dist" %in% class(dx))
  R <- rips(dx, eps = 2*radius, dim = 1L)
  miniball_test <- function(idx){
    valid_ball <- seb::seb(x[idx,,drop=FALSE])$radius <= radius
    return(valid_ball)
  }
  ## Compute the Cech complex
  cech_complex(stx = R$as_XPtr(), k = dim, seb_f = miniball_test)
  return(R)
}