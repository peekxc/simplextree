#' @export
traverse <- function(st, simplex, f, type=c("dfs", "bfs", "link", "star", "cofaces", "k-skeleton"), ...){
  stopifnot(is(st, "Rcpp_SimplexTree"), is.vector(simplex), is.function(f))
  extra <- list(...)
  
}