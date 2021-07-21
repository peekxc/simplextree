#' @name simplex_tree
#' @aliases SimplexTree Rcpp_SimplexTree
#' @title Simplex Tree
#' @description Simplex tree class exposed as an Rcpp Module.
#' @docType class 
#' @details A simplex tree is an ordered trie-like structure specialized for storing and doing general computation 
#' simplicial complexes. Here is figure of a simplex tree, taken from the original paper (see Boissonnat and Maria, 2014):
#' 
#' \if{html}{\figure{simplextree.png}{options: width="80\%" alt="Figure: simplextree.png"}}
#' \if{latex}{\figure{simplextree.pdf}{options: width=12cm}}
#' 
#' The current implementation provides a subset of the functionality described in the paper.
#' 
#' @field n_simplices A vector, where each index \eqn{k} denotes the number \eqn{(k-1)}-simplices.
#' @field dimension The dimension of the simplicial complex.
#' @section Properties:
#' Properties are actively bound shortcuts to various methods of the simplex tree that may be thought of as fields. 
#' Unlike fields, however, properties are not explicitly stored: they are generated on access. 
#' \describe{
#'     \item{`$id_policy`}{ The policy used to generate new vertex ids. May be assigned `"compressed"` or `"unique"`. See [generate_ids()]. }
#'     \item{`$vertices`}{ The 0-simplices of the simplicial complex, as a matrix. }
#'     \item{`$edges`}{ The 1-simplices of the simplicial complex, as a matrix. }
#'     \item{`$triangles`}{ The 2-simplices of the simplicial complex, as a matrix. }
#'     \item{`$quads`}{ The 3-simplices of the simplicial complex, as a matrix. }
#'     \item{`$connected_components`}{ The connected components of the simplicial complex. }
#' } 
#' @section Methods: 
#' \describe{
#'     \item{`$as_XPtr()`}{ Creates an external pointer. }
#'     \item{`$`[generate_ids()]}{ Generates new vertex ids according to the set policy. }
#'     \item{`$`[insert()]}{ Inserts a simplex into the trie. }
#'     \item{`$`[remove()]}{ Removes a simplex from the trie. }
#'     \item{`$`[find()]}{ Returns whether a simplex exists in the trie. }
#'     \item{`$`[degree()]}{ Returns the degree of each given vertex. }
#'     \item{`$`[adjacent()]}{ Returns vertices adjacent to a given vertex. }
#'     \item{`$`[clear()]}{ Clears the simplex tree. }
#'     \item{`$`[expand()]}{ Performs an \eqn{k}-expansion. }
#'     \item{`$`[collapse()]}{ Performs an elementary collapse. }
#'     \item{`$`[contract()]}{ Performs an edge contraction. }
#'     \item{`$`[traverse()]}{ Traverses a subset of the simplex tree, applying a function to each simplex. }
#'     \item{`$`[ltraverse()]}{ Traverses a subset of the simplex tree, applying a function to each simplex and returning the result as a list. }
#'     \item{`$`[is_face()]}{ Checks for faces. }
#'     \item{`$`[is_tree()]}{ Checks if the simplicial complex is a tree. }
#'     \item{`$as_list()`}{ Converts the simplicial complex to a list. }
#'     \item{`$as_adjacency_matrix()`}{ Converts the 1-skeleton to an adjacency matrix. }
#'     \item{`$as_adjacency_list()`}{ Converts the 1-skeleton to an adjacency list. }
#'     \item{`$as_edgelist()`}{ Converts the 1-skeleton to an edgelist. }
#' }
#' @author Matt Piekenbrock
#' @import methods
#' @param simplices optional simplices to initialize the simplex tree with. See [insert()].
#' @return A queryable simplex tree, as an object (Rcpp module) of class `"Rcpp_SimplexTree"`.
#' @template ref-boissonnat2014
#' @examples
#' ## Recreating simplex tree from figure. 
#' st <- simplex_tree()
#' st %>% insert(list(1:3, 2:5, c(6, 7, 9), 7:8, 10))
#' plot(st)
#' 
#' ## Example insertion
#' st <- simplex_tree(list(1:3, 4:5, 6)) ## Inserts one 2-simplex, one 1-simplex, and one 0-simplex
#' print(st) 
#' # Simplex Tree with (6, 4, 1) (0, 1, 2)-simplices
#' 
#' ## More detailed look at structure
#' print_simplices(st, "tree")
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0): 
#' # 4 (h = 1): .( 5 )
#' # 5 (h = 0): 
#' # 6 (h = 0): 
#' ## Print the set of simplices making up the star of the simplex '2'
#' print_simplices(st %>% cofaces(2))
#' # 2, 2 3, 1 2, 1 2 3
#' 
#' ## Retrieves list of all simplices in DFS order, starting with the empty face 
#' dfs_list <- ltraverse(st %>% preorder(empty_face), identity)
#' 
#' ## Get connected components 
#' print(st$connected_components)
#' # [1] 1 1 1 4 4 5
#' 
#' ## Use clone() to make copies of the complex (don't use the assignment `<-`)
#' new_st <- st %>% clone()
#' 
#' ## Other more internal methods available via `$` 
#' list_of_simplices <- st$as_list()
#' adj_matrix <- st$as_adjacency_matrix()
#' # ... see also as_adjacency_list(), as_edge_list(), etc 
#' @family simplicial complex constructors
#' @export
simplex_tree <- function(simplices = NULL){
  st <- new(SimplexTree)
  # assign('insert', function(x, check=TRUE) { return(st$insert(x, check)) }, envir = st)
  if (!missing(simplices)){ st %>% insert(simplices) }
  return(st)
}

# ---- generate_ids ----

#' @name generate_ids
#' @aliases id_policy
#' @title Generate vertex ids
#' @description Generates vertex ids representing 0-simplices not in the tree.
#' @param st a simplex tree. 
#' @param n the number of ids to generate. 
#' @return a double vector of the `n` smallest natural numbers (starting at `0`) that are not vertex ids of `st`.
#' @details This function generates new vertex ids for use in situations which involve generating new 
#' new 0-simplices, e.g. insertions, contractions, collapses, etc. There are two 'policies' which designate 
#' the generating mechanism of these ids: 'compressed' and 'unique'. 'compressed' generates vertex ids 
#' sequentially, starting at 0. 'unique' tracks an incremental internal counter, which is updated on every 
#' call to `generate_ids()`. The new ids under the 'unique' policy generates the first sequential `n` 
#' ids that are strictly greater `max(`\emph{counter}`, `\emph{max vertex id}`)`.
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
#' @family simplex-level operations
#' @export
generate_ids <- function(st, n){
  stopifnot(is.numeric(n) && length(n) == 1)
  return(st$generate_ids(as.integer(n)))
}


# ---- clone ----
#' @name clone 
#' @title Clones the given simplex tree.
#' @param st a simplex tree.
#' @description Performs a deep-copy on the supplied simplicial complex. 
#' @export
clone <- function(st){
  stopifnot(class(st) %in% .st_classes)
  new_st <- deserialize(st %>% serialize())
  return(new_st)
}

# ---- insert ----

#' @name insert
#' @title Insert simplices
#' @description Inserts simplices into the simplex tree. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param st a simplex tree.  
#' @param simplices simplices to insert, either as a vector, a list of vectors, or a column-matrix. See details. 
#' @return the simplex tree `st` with the simplices `simplices` inserted, invisibly.
#' @details This function allows insertion of arbitrary order simplices. If the simplex already exists in the tree, 
#' no insertion is made, and the tree is not modified. `simplices` is sorted before traversing the trie. 
#' Faces of `simplices` not in the simplex tree are inserted as needed.
#' 
#' If `simplices` is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @family simplex-level operations
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
#' @return the simplex tree `st` with the simplices `simplices` removed, invisibly.
#' @details This function allows removal of a arbitrary order simplices. If `simplices` already exists in the tree, 
#' it is removed, otherwise the tree is not modified. `simplices` is sorted before traversing the trie.
#' Cofaces of `simplices` are also removed.
#' 
#' If `simplices` is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @family simplex-level operations
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
#' @return a boolean vector indicating whether each simplex in `simplices` exists in `st`.
#' @section Usage:
#' st %>% find(simplices)
#' @details Traverses the simplex tree looking for `simplices`, returning whether or not it exists.
#' `simplices` can be specified as vector to represent a single simplex, and a list to represent a set of simplices. 
#' Each `simplices` is sorted before traversing the trie.
#' 
#' If `simplices` is a vector, it's assumed to be a simplex. If a list, its assumed each element in the list 
#' represents a simplex (as vectors). If the simplices to insert are all of the same dimension, you can also 
#' optionally use a matrix, where each column is assumed to be a simplex. 
#' @family simplex-level operations
#' @export
find <- function(st, simplices){
  stopifnot(class(st) %in% .st_classes)
  st$find(simplices)
}
