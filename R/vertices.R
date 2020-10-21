
# ---- degree ----

#' @name degree
#' @title Vertex degree
#' @description Returns the number of edges (degree) for each given vertex id. 
#' @param st a simplex tree. 
#' @param vertices the vertex ids to check the degree of. 
degree <- function(st, vertices){
  stopifnot(is.vector(vertices) && is.numeric(vertices))
  return(st$degree(vertices))
}

# ---- adjacent ----

#' @name adjacent
#' @title Adjacent vertices
#' @description Returns a vector of vertex ids that are immediately adjacent to a given vertex.
#' @param st a simplex tree.
#' @param vertices vertex ids. 
#' @examples
#' st <- simplex_tree(1:3)
#' st %>% adjacent(2) 
#' # 1 3
#' @export
adjacent <- function(st, vertices){
  stopifnot(is.vector(vertices))
  if (length(vertices) == 1){ return(st$adjacent(vertices)) }
  else {
    return(lapply(vertices, st$adjacent))
  }
}
