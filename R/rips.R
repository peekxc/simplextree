#' @title Vietoris–Rips complex
#' @description Constructs the Vietoris–Rips complex. 
#' @param d a numeric 'dist' vector. 
#' @param eps diameter parameter. 
#' @param dim maximum dimension to construct up to. Defaults to 1 (edges only).
#' @param filtered whether to construct the filtration. Defaults to false. See details. 
#' @family simplicial complex constructors
#' @export
rips <- function(d, eps = enclosing_radius(d), dim = 1L, filtered = FALSE){
	stopifnot(is.numeric(d) || 'dist' %in% class(d))
  stopifnot(is.numeric(eps))
	n <- inverse.choose(length(d), k = 2)
	ind_to_insert <- which(d <= eps)
	st <- simplex_tree() %>% 
			insert(as.list(seq(n))) %>% 
			insert(nat_to_sub(ind_to_insert, n, 2)) %>% 
			expand(k = dim)
	if (filtered){ st <- st %>% flag(d[ind_to_insert]) }
	return(st)
}

#' @md
#' @title Enclosing radius of a set of distances
#' @description Computes the enclosing radius of a set of distances. 
#' @param d a [stats::dist()] object. 
#' @details The enclosing radius is useful as an upper bound of the scale parameter 
#' for the rips filtration. Scales above the enclosing radius render the Vietoris–Rips
#' complex as a simplicial cone, beyond which the homology is trivial. 
#' @export 
enclosing_radius <- function(d){
	stopifnot(is.numeric(d) || 'dist' %in% class(d))
	n <- inverse.choose(length(d), k = 2)
	ii <- seq(n)
	radii <- sapply(ii, function(i){ max(d[sub_to_nat(rbind(i, setdiff(ii, i)), n)]) })
	return(min(radii))
}
