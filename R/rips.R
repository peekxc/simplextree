#' rips
#' @description Constructs the Vietoris-Rips complex. 
#' @param d a numeric 'dist' vector. 
#' @param eps diameter parameter. 
#' @param dim maximum dimension to construct up to. Defaults to 1 (edges only).
#' @param filtered whether to construct the filtration. Defaults to false. See details. 
#' @details This function constructs a rips complex from a 'dist' object up to scale \code{diameter}. 
#' By default, if \code{diameter} is not provided, a rips complex is constructed up to the enclosing radius 
#' of the space, at which the complex has trivial homology. 
#' @export
rips <- function(d, diameter = enclosing_radius(d), dim = 1L, filtered = FALSE, ...){
  params <- list(...)
  if ("eps" %in% names(params)){
    diameter <- params[["eps"]]
  }
	stopifnot(is.numeric(d) || 'dist' %in% class(d))
  stopifnot(is.numeric(diameter))
	n <- inverse.choose(length(d), k = 2)
	ind_to_insert <- which(d <= diameter)
	st <- simplex_tree() %>% 
			insert(as.list(seq(n))) %>% 
			insert(nat_to_sub(ind_to_insert, n, 2)) %>% 
			expand(k = dim)
	if (filtered){ st <- st %>% flag(d[ind_to_insert]) }
	return(st)
}

#' enclosing_radius
#' @description Computes the enclosing radius of a set of distances. 
#' @param d a \code{\link[stats]{dist}} object. 
#' @details The enclosing radius is useful as an upper bound of the scale parameter 
#' for the rips filtration. Scales at or above the enclosing radius render the Vietoris–Rips
#' complex as a simplicial cone, beyond which the homology is trivial. 
#' @export 
enclosing_radius <- function(d){
	stopifnot(is.numeric(d) || 'dist' %in% class(d))
	n <- inverse.choose(length(d), k = 2)
	ii <- seq(n)
	radii <- sapply(ii, function(i){ max(d[sub_to_nat(rbind(i, setdiff(ii, i)), n)]) })
	return(min(radii)/2.0)
}