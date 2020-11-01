
# ---- clear ----

#' @name clear
#' @title Clear a simplex tree
#' @description Removes all simplices from the simplex tree, except the root
#'   node.
#' @param st a simplex tree object.
#' @return the simplex tree `st` with simplices removed, invisibly.
#' @examples
#' st <- simplex_tree()
#' st %>% insert(1:3)
#' print(st) ## Simplex Tree with (3, 3, 1) (0, 1, 2)-simplices
#' st %>% clear()
#' print(st) ## < empty simplex tree >
#' @family complex-level operations
#' @export
clear <- function(st){
  stopifnot(class(st) %in% .st_classes)
  st$clear()
  return(invisible(st))
}

# ---- expand ----

#' @name expand
#' @title \eqn{k}-expansion
#' @description Performs a \eqn{k}-expansion on the 1-skeleton of the complex, adding \eqn{k}-simplices 
#' if all combinations of edges are included. Because this operation uses the edges alone to infer 
#' the existence of higher order simplices, the expansion assumes the underlying complex
#' is a flag complex. 
#' @param st a simplex tree. 
#' @param k the target dimension of the expansion.
#' @return the expanded simplex tree `st`, invisibly.
#' @family complex-level operations
#' @export
expand <- function(st, k=2){
  stopifnot(is.numeric(k))
  st$expand(k)
  return(invisible(st))
}

# ---- collapse ----

#' @name collapse
#' @title Elementary collapse
#' @description Performs an elementary collapse. 
#' @param st a simplex tree.
#' @param pair list of simplices to collapse. 
#' @param w vertex to collapse to, if performing a vertex collapse. 
#' @details This function provides two types of _elementary collapses_.
#' 
#' The first type of collapse is in the sense described by (1), which is 
#' summarized here. A simplex \eqn{\sigma} is said to be collapsible through one of its faces \eqn{\tau} if 
#' \eqn{\sigma} is the only coface of \eqn{\tau} (excluding \eqn{\tau} itself). This function checks whether its possible to collapse \eqn{\sigma} through \eqn{\tau}, 
#' (if \eqn{\tau} has \eqn{\sigma} as its only coface), and if so, both simplices are removed. 
#' `tau` and `sigma` are sorted before comparison.
#' To perform this kind of elementary collapse, call `collapse` with two simplices as arguments, i.e. `tau` before `sigma`.
#' 
#' Alternatively, this method supports another type of elementary collapse, also called a _vertex collapse_, as described 
#' in (2). This type of collapse maps a pair of vertices into a single vertex. To use this collapse, specify three vertex ids, the first 
#' two representing the free pair, and the last representing the target vertex to collapse to. 
#' 
#' @return boolean indicating whether the collapse was performed. 
#' @template ref-boissonnat2014
#' @template ref-dey2014
#' @examples 
#' st <- simplextree::simplex_tree(1:3)
#' st %>% print_simplices()
#' # 1, 2, 3, 1 2, 1 3, 2 3, 1 2 3
#' st %>% collapse(list(1:2, 1:3))
#' # 1, 2, 3, 1 3, 2 3=
#' 
#' st %>% insert(list(1:3, 2:5))
#' st %>% print_simplices("column")
#' # 1 2 3 4 5 1 1 2 2 2 3 3 4 1 2 2 2 3 2
#' #           2 3 3 4 5 4 5 5 2 3 3 4 4 3
#' #                           3 4 5 5 5 4
#' #                                     5
#' 
#' st %>% collapse(list(2:4, 2:5))
#' st %>% print_simplices("column") 
#' # 1 2 3 4 5 1 1 2 2 2 3 3 4 1 2 2 3
#' #           2 3 3 4 5 4 5 5 2 3 4 4
#' #                           3 5 5 5
#' @export
collapse <- function(st, pair, w=NULL){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.list(pair))
  tau <- pair[[1]]
  sigma <- pair[[2]]
  if (missing(w)){
    stopifnot(all(tau %in% sigma))
    st$collapse(tau, sigma)
  } else {
    stopifnot(all(sapply(list(tau,sigma,w), length) == 1))
    st$vertex_collapse(tau, sigma, w)
  }
  return(invisible(st))
}

#' @name threshold
#' @title Filtered complex thresholding
#' @description Thresholds a given filtered simplicial complex.
#' @param st simplex tree.
#' @param index integer index to threshold to.
#' @param value numeric index to threshold filtration. 
#' @return the thresholded simplex tree `st`, invisibly.
#' @family complex-level operations
#' @export
threshold <- function(st, index = NULL, value = NULL){
  stopifnot(class(st) %in% .st_classes)
  if ("Rcpp_Filtration" %in% class(st)){
    if (missing(index) && !missing(value)){
      st$threshold_value(value)
      return(invisible(st))
    } else if (!missing(index) && missing(value)){
      st$threshold_index(index)
      return(invisible(st))
    } else {
      stop("Must supply either an integer index or numeric value.")
    }
  } else { return(invisible(st)) }
}

# ---- contract ----

#' @name contract
#' @title Edge contraction
#' @description Performs an edge contraction. 
#' @param st a simplex tree.
#' @param edge an edge to contract, as a 2-length vector. 
#' @return the contracted simplex tree `st`, invisibly.
#' @details This function performs an _edge contraction_ in the sense described by (1), which is 
#' summarized here. Given an edge \eqn{(v,w)}, \eqn{w} is contracted to \eqn{v} if \eqn{w} is 
#' removed from the complex and the link of \eqn{v} is augmented with the link of \eqn{w}. This may be thought as 
#' applying the mapping:
#' 
#' \deqn{f(u) = v}
#' 
#' if \eqn{u = w}
#' and identity otherwise, to all simplices in the complex.
#' 
#' `edge` is **not** sorted prior to contraction: the second vertex of the edge is always contracted to the first. 
#' Note that edge contraction is not symmetric.
#' @template ref-boissonnat2014
#' @examples 
#' st <- simplex_tree(1:3) 
#' st %>% print_simplices()
#' # 1, 2, 3, 1 2, 1 3, 2 3, 1 2 3
#' st %>% contract(c(1, 3)) %>% print_simplices()
#' # 1, 2, 1 2
#' @family complex-level operations
#' @export
contract <- function(st, edge){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.numeric(edge) && length(edge) == 2)
  st$contract(edge)
  return(invisible(st))
}
