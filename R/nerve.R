#' @title Nerve of a cover
#' @description Compute the nerve of a given cover.
#' @param st a simplex tree.
#' @param cover list of integers indicating set membership. See details. 
#' @param k max simplex dimension to consider. 
#' @param threshold the number of elements in common for `k` sets to be considered intersecting. Defaults to 1. 
#' @param neighborhood which combinations of sets to check. See details. 
#' @details This computes the nerve of a given cover, adding a \eqn{k}-simplex for each combination of \eqn{k+1} sets 
#' in the given `cover` that have at least `threshold` elements in their common intersection.
#' 
#' If `neighborhood` is supplied, it can be either (1) a matrix, (2) a list, or (3) a function. Each 
#' type parameterizes which sets in the cover need be checked for to see if they have at least `threshold`
#' elements in their common intersection. If a matrix is supplied, the columns should indicate the indices 
#' of the cover to check (e.g if `neighborhood = matrix(c(1,2), nrow = 2)`, then only the first two sets of `cover`
#' are tested.). Similarly, if a list is supplied, each element in the list should give the indices to test.
#' 
#' The most flexible option is supplying a function to `neighborhood`. If a function is passed, it's assumed to 
#' accept an integer vector of \eqn{k} indices (of the cover) and return a boolean indicating whether or not to 
#' _test_ if they have at least `threshold` elements in their common intersection. This can be used
#' to filter out subsets of the cover the user knows are  The indices are 
#' generated using the same code that performs [expand()].
#' 
#' @family simplicial complex constructors
#' @export
nerve <- function(st, cover, k = st$dimension, threshold = 1L, neighborhood=NULL){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(all(sapply(cover, is.numeric)))
  
  ## Get the cover names to set as the vertices, if not done already
  if (!is.null(names(cover))){
    stopifnot(all(!is.na(as.integer(names(cover)))))
    cover_ids <- as.vector(as.integer(names(cover)))
  } else {
    if (length(st$n_simplices) == 0){ 
      cover_ids <- seq(length(cover)) 
    } else {
      cover_ids <- st$vertices 
    }
  }
  stopifnot(length(cover) == length(cover_ids))
  
  ## Perform a conditional k-expansion based on the cover intersections
  if (missing(neighborhood) || is.null(neighborhood)){
    nerve_expand(st$as_XPtr(), cover_ids, cover, k, threshold)
  } else if (is.matrix(neighborhood)){
    stopifnot(all(neighborhood %in% cover_ids))
    k <- nrow(neighborhood)
    to_include <- apply(neighborhood, 2, function(ids){ 
      nfold_intersection(cover[match(ids, cover_ids)], threshold) 
    })
    st %>% insert(neighborhood[,to_include,drop=FALSE])
  } else if (is.list(neighborhood)){
    stopifnot(all(sapply(neighborhood, function(x){
      all(is.numeric(x) & (x %in% cover_ids))
    })))
    to_include <- sapply(neighborhood, function(ids){ 
      nfold_intersection(cover[match(ids, cover_ids)], threshold)
    })
    st %>% insert(neighborhood[to_include])
  } else if (is.function(neighborhood)){
    ## Should take as a function a vector of ids and return a true or false indicating whether to include it
    nerve_expand_f(st$as_XPtr(), ids = cover_ids, include_f = neighborhood, k = k)
  }
  
  ## Return the complex invisibly
  return(invisible(st))
}

# nfold_intersection()
