# ---- degree ----

#' @name degree
#' @title Vertex degree
#' @description Returns a list of vertex ids that are immediately adjacent to a given vertex. If a given vertex
#' does not have any adjacencies, a vector of length 0 is returned. 
#' @template params-vertices
#' @return an integer vector of degrees of `vertices` in `st` (taken to be `0` for vertices not in `st`).
#' @family vertex-level operations
#' @export
degree <- function(st, vertices=NULL){
  if (missing(vertices) || is.null(vertices)){ vertices <- st$vertices }
  stopifnot(is.vector(vertices) && is.numeric(vertices))
  return(st$degree(vertices))
}

# ---- adjacent ----

#' @name adjacent
#' @title Adjacent vertices
#' @description Returns a vector of vertex ids that are immediately adjacent to a given vertex.
#' @template params-vertices
#' @return a list of double vectors of vertices adjacent to each of `vertices` in `st` (or `numeric(0)` for vertices not in `st`), unlisted to a single vector if `length(vertices == 1`.
#' @examples
#' st <- simplex_tree(1:3)
#' st %>% adjacent(2) 
#' # 1 3
#' @family vertex-level operations
#' @export
adjacent <- function(st, vertices=NULL){
  if (missing(vertices) || is.null(vertices)){ vertices <- st$vertices }
  stopifnot(is.vector(vertices) && is.numeric(vertices))
  return(st$adjacent(vertices))
}
