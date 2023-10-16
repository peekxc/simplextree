#' @title Flag complexes
#' @description Creates a filtration of flag complexes
#' @param st a simplex tree. See details. 
#' @param d a vector of edge weights, or a 'dist' object. 
#' @return a simplicial filtration (object of class `"Rcpp_Filtration"`).
#' @details A flag complex is a simplicial complex whose \eqn{k}-simplices for \eqn{k >= 2} are completely determined 
#' by edges/graph of the complex. This function creates filtered simplicial complex using the supplied edge 
#' weights. The resulting complex is a simplex tree object endowed with additional structure; see. 
#' Vertices have their weights set to 0, and \eqn{k}-simplices w/ \eqn{k >= 2} have their weights set to the maximum
#' weight of any of its edges. 
#' @family simplicial complex constructors
#' @export
flag <- function(st, d){
  stopifnot(is.numeric(d), class(st) %in% .st_classes)
  fi <- new(Filtration)
  fi$init_tree(st$as_XPtr())
  fi$flag_filtration(d)
  return(fi)
}
