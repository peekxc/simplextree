#' nerve 
#' @description Compute the nerve of a given cover.
#' @param st a simplex tree.
#' @param cover list of integers indicating set membership. See details. 
#' @param k max simplex dimension to consider. 
#' @param threshold the number of elements in common for \code{k} sets to be considered intersecting. Defaults to 1. 
#' @param neighborhood which combinations of sets to check. See details. 
#' @details This computes the nerve of a given cover, adding a \emph{k}-simplex for each combination of \emph{k+1} sets 
#' in the given \code{cover} that have at least \code{threshold} elements in their common intersection. 
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
