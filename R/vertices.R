
# ---- reindex ----

#' @name reindex 
#' @title Reindex vertex ids
#' @description This function allows one to 'reorder' or 'reindex' vertex ids.  
#' @param st a simplex tree. 
#' @param ids vector of new vertex ids. See details. 
#' @return the simplex tree `st` with vertices reindexed by `ids`, invisibly.
#' @details The `ids` parameter must be a sorted integer vector of new ids with length matching the 
#' number of vertices. The simplex tree is modified to replace the vertex label at index `i` with 
#' `ids[i]`. See examples.
#' @examples 
#' st <- simplex_tree()
#' st %>% insert(1:3) %>% print_simplices("tree")
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0):
#' st %>% reindex(4:6) %>% print_simplices("tree")
#' # 4 (h = 2): .( 5 6 )..( 6 )
#' # 5 (h = 1): .( 6 )
#' # 6 (h = 0):
#' @family vertex-level operations
#' @export
reindex <- function(st, ids){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.numeric(ids) && all(order(ids) == seq(length(ids))))
  st$reindex(ids)
  return(invisible(st))
}

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
