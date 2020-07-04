#' UnionFind
#' @aliases union_find 
#' @description Union find structure exposed as an Rcpp Module.
#' @docType class 
#' @section Methods: 
#' \describe{
#'     \item{$\code{print.simplextree}}{ S3 method to print a basic summary of the simplex tree. }
#' }
#' @author Matt Piekenbrock
#' @import methods
#' @param n Number of elements in the set. 
#' @return A disjoint set, as a \code{Rcpp_UnionFind} object (Rcpp module). 
union_find <- function(n = 0L){
  stopifnot((is.numeric(n) && n >= 0) || is.list(n))
  uf <- new(UnionFind, n)
  return(uf)
}