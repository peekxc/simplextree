#' @title Vietoris–Rips complex
#' @description Constructs the Vietoris–Rips complex. 
#' @param d a numeric 'dist' vector. 
#' @param eps diameter parameter. 
#' @param dim maximum dimension to construct up to. Defaults to 1 (edges only).
#' @param filtered whether to construct the filtration. Defaults to false. See details. 
#' @return a simplicial complex (object of class `"Rcpp_SimplexTree"`).
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

#' @title Enclosing radius of a set of distances
#' @description Computes the enclosing radius of a set of distances. 
#' @param d a [stats::dist()] object. 
#' @return a numeric scalar.
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

#' flag
#' @description Creates a filtration of flag complexes
#' @param st a simplex tree. See details. 
#' @param d a vector of edge weights, or a 'dist' object. 
#' @details A flag complex is a simplicial complex whose k-simplices for k >= 2 are completely determined 
#' by edges/graph of the complex. This function creates filtered simplicial complex using the supplied edge 
#' weights. The resulting complex is a simplex tree object endowed with additional structure; see. 
#' Vertices have their weights set to 0, and k-simplices w/ k >= 2 have their weights set to the maximum
#' weight of any of its edges. 
#' @export
flag <- function(st, d){
  stopifnot(is.numeric(d), class(st) %in% .st_classes)
  fi <- new(Filtration)
  fi$init_tree(st$as_XPtr())
  fi$flag_filtration(d)
  return(fi)
}