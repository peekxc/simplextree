
# ---- generate_ids ----

#' @name generate_ids
#' @aliases id_policy
#' @title Generate vertex ids
#' @description Generates vertex ids representing 0-simplices not in the tree.
#' @param st a simplex tree. 
#' @param n the number of ids to generate. 
#' @details This function generates new vertex ids for use in situations which involve generating new 
#' new 0-simplices, e.g. insertions, contractions, collapses, etc. There are two 'policies' which designate 
#' the generating mechanism of these ids: 'compressed' and 'unique'. 'compressed' generates vertex ids 
#' sequentially, starting at 0. 'unique' tracks an incremental internal counter, which is updated on every 
#' call to \code{generate_ids}. The new ids under the 'unique' policy generates the first sequential \code{n} 
#' ids that are strictly greater  \code{max}(\emph{counter}, \emph{max vertex id}). \cr
#' \cr
#' 
#' @examples 
#' st <- simplex_tree()
#' print(st$id_policy)
#' ## "compressed"
#' st %>% generate_ids(3) 
#' ## 0 1 2
#' st %>% generate_ids(3) 
#' ## 0 1 2
#' st %>% insert(list(1,2,3))
#' print(st$vertices) 
#' ## 1 2 3
#' st %>% insert(as.list(st %>% generate_ids(2)))
#' st %>% print_simplices() 
#' # 0, 1, 2, 3, 4
#' st %>% remove(4)
#' st %>% generate_ids(1) 
#' # 4
#' @export
generate_ids <- function(st, n){
  stopifnot(is.numeric(n) && length(n) == 1)
  return(st$generate_ids(as.integer(n)))
}

# ---- insert ----

#' @name insert
#' @title Insert simplices
#' @description Inserts simplices into the simplex tree. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param st a simplex tree.  
#' @param simplices simplices to insert, either as a vector, a list of vectors, or a column-matrix. See details. 
#' @details This function allows insertion of arbitrary order simplices. If the simplex already exists in the tree, 
#' no insertion is made, and the tree is not modified. \code{simplex} is sorted before traversing the trie. 
#' Faces of \code{simplex} not in the simplex tree are inserted as needed. \cr
#' \cr
#' If \code{simplices} is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @seealso find remove
#' @examples 
#' st <- simplex_tree()
#' st %>% insert(1:3) ## inserts the 2-simplex { 1, 2, 3 }
#' st %>% insert(list(4:5, 6)) ## inserts a 1-simplex { 4, 5 } and a 0-simplex { 6 }.
#' st %>% insert(combn(5,3)) ## inserts all the 2-faces of a 4-simplex
#' @export
insert <- function(st, simplices){
  stopifnot(class(st) %in% .st_classes)
  st$insert(simplices)
  return(invisible(st))
}

# ---- remove ----

#' @name remove
#' @title Remove simplices
#' @description Removes simplices from the simplex tree. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param st a simplex tree.  
#' @param simplices simplices to insert, either as a vector, a list of vectors, or a column-matrix. See details. 
#' @details This function allows removal of a arbitrary order simplices. If \code{simplex} already exists in the tree, 
#' it is removed, otherwise the tree is not modified. \code{simplex} is sorted before traversing the trie.
#' Cofaces of \code{simplex} are also removed. \cr
#' \cr
#' If \code{simplices} is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @seealso find remove
#' @export
remove <- function(st, simplices){
  stopifnot(class(st) %in% .st_classes)
  st$remove(simplices)
  return(invisible(st))
}

# ---- find ----

#' @name find
#' @title Find simplices
#' @description Returns whether supplied simplices exist in the complex.  
#' @param st a simplex tree.  
#' @param simplices simplices to insert, either as a vector, a list of vectors, or a column-matrix. See details. 
#' @section Usage:
#' st %>% find(simplices)
#' @details Traverses the simplex tree looking for \code{simplex}, returning whether or not it exists.
#' \code{simplex} can be specified as vector to represent a single simplex, and a list to represent a set of simplices. 
#' Each \code{simplex} is sorted before traversing the trie. \cr
#' \cr
#' If \code{simplices} is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @return boolean indicating whether or not \code{simplex} exists in the tree. 
#' @seealso insert remove
#' @export
find <- function(st, simplices){
  stopifnot(class(st) %in% .st_classes)
  st$find(simplices)
}

# ---- reindex ----

#' @name reindex 
#' @title Reindex vertex ids
#' @description This function allows one to 'reorder' or 'reindex' vertex ids.  
#' @param st a simplex tree. 
#' @param ids vector of new vertex ids. See details. 
#' @details The \code{ids} parameter must be a sorted integer vector of new ids with length matching the 
#' number of vertices. The simplex tree is modified to replace the vertex label at index \code{i} with 
#' \code{ids}[i]. See examples.
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
#' @export
reindex <- function(st, ids){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.numeric(ids) && all(order(ids) == seq(length(ids))))
  st$reindex(ids)
  return(invisible(st))
}
