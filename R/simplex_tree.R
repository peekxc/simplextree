#' Simplex Tree
#' @aliases simplex_tree simplextree SimplexTree
#' @description Simplex tree class exposed as an Rcpp Module.
#' @docType class 
#' @details A simplex tree is an ordered trie-like structure specialized for storing and doing general computation 
#' simplicial complexes. Here is figure of a simplex tree, taken from the original paper (see 1): \cr  
#' \if{html}{\figure{simplextree.png}{options: width="80\%" alt="Figure: simplextree.png"}}
#' \if{latex}{\figure{simplextree.pdf}{options: width=12cm}}
#' \cr 
#' The current implementation provides a subset of the functionality described in the paper.
#' 
#' @field n_simplices A vector, where each index k denotes the number (k-1)-simplices.
#' @field dimension The dimension of the simplicial complex.
#' @section Properties:
#' Properties are actively bound shortcuts to various methods of the simplex tree that may be thought of as fields. 
#' Unlike fields, however, properties are not explicitly stored: they are generated on access. 
#' \describe{
#'     \item{$\code{id_policy}}{ The policy used to generate new vertex ids. May be assigned "compressed" or "unique". See \code{\link{generate_ids}}. }
#'     \item{$\code{vertices}}{ The 0-simplices of the simplicial complex, as a matrix. }
#'     \item{$\code{edges}}{ The 1-simplices of the simplicial complex, as a matrix. }
#'     \item{$\code{triangles}}{ The 2-simplices of the simplicial complex, as a matrix. }
#'     \item{$\code{quads}}{ The 3-simplices of the simplicial complex, as a matrix. }
#'     \item{$\code{connected_components}}{ The connected components of the simplicial complex. }
#' } 
#' @section Methods: 
#' \describe{
#'     \item{$\code{print.simplextree}}{ S3 method to print a basic summary of the simplex tree. }
#'     \item{$\code{\link{plot.simplextree}}}{ S3 method to plot the simplicial complex. }
#'     \item{$\code{\link{print_tree}}}{ Prints the simplex tree structure. }
#'     \item{$\code{\link{as_XPtr}}}{ Creates an external pointer. }
#'     \item{$\code{\link{clear}}}{ Clears the simplex tree. }
#'     \item{$\code{\link{generate_ids}}}{ Generates new vertex ids according to the set policy. }
#'     \item{$\code{\link{degree}}}{ Returns the degree of each given vertex. }
#'     \item{$\code{\link{adjacent}}}{ Returns vertices adjacent to a given vertex. }
#'     \item{$\code{\link{insert}}}{ Inserts a simplex into the trie. }
#'     \item{$\code{\link{remove}}}{ Removes a simplex from the trie. }
#'     \item{$\code{\link{find}}}{ Returns whether a simplex exists in the trie. }
#'     \item{$\code{\link{collapse}}}{ Performs an elementary collapse. }
#'     \item{$\code{\link{contract}}}{ Performs an edge contraction. }
#'     \item{$\code{\link{expand}}}{ Performs an k-expansion. }
#'     \item{$\code{\link{traverse}}}{ Traverses a subset of the simplex tree, applying a function to each simplex. }
#'     \item{$\code{\link{ltraverse}}}{ Traverses a subset of the simplex tree, applying a function to each simplex and returning the result as a list. }
#'     \item{$\code{\link{is_face}}}{ Checks for faces. }
#'     \item{$\code{\link{is_tree}}}{ Checks if the simplicial complex is a tree. }
#'     \item{$\code{\link{serialize}}}{ Serializes the simplex tree. }
#'     \item{$\code{\link{deserialize}}}{ Unserializes a stored simplex tree. }
#'     \item{$\code{\link{save}}}{ Saves the simplex tree to a file. }
#'     \item{$\code{\link{load}}}{ Loads a simplex tree from a file. }
#'     \item{$\code{as_list}}{ Converts the simplicial complex to a list. }
#'     \item{$\code{as_adjacency_matrix}}{ Converts the 1-skeleton to an adjacency matrix. }
#'     \item{$\code{as_adjacency_list}}{ Converts the 1-skeleton to an adjacency list. }
#'     \item{$\code{as_edgelist}}{ Converts the 1-skeleton to an edgelist. }
#' }
#' @author Matt Piekenbrock
#' @import methods
#' @return A queryable simplex tree, as a \code{Rcpp_SimplexTree} object (Rcpp module). 
#' @references Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples
#' ## Recreating simplex tree from figure. 
#' st <- simplex_tree()
#' st$insert(list(1:3, 2:5, c(6, 7, 9), 7:8, 10))
#' plot(st)
#' 
#' ## Example insertion
#' st$insert(list(1:3, 4:5, 6)) ## Inserts one 2-simplex, one 1-simplex, and one 0-simplex
#' @export
simplex_tree <- function(){
  return(new(SimplexTree))
}

# ---- empty_face ----
#' empty_face 
#' @description Simple alias to the NULL value, used to indicate the empty face. 
#' @seealso traverse
#' @export
empty_face <- NULL

# ---- print_tree ----
#' @name print_tree
#' @title Prints the simplex tree
#' @description Prints the simplicial complex to standard out. 
#' By default, this is set to R's buffered output, which is shown in the R console. 
#' The printed format is: \cr 
#' \cr
#' [vertex] (h = [subtree height]): [subtree depth]([subtree]) \cr 
#' \cr
#' Where each lists the top node (\emph{vertex}) and its corresponding subtree. The 
#' \emph{subtree height} displays the highest order k-simplex in that subtree. Each 
#' level in the subtree tree is a set of sibling k-simplices whose order is given  
#' by the number of dots ('.') proceeding the print level.
NULL

# ---- print.Rcpp_SimplexTree ----
setClass("Rcpp_SimplexTree")
.print_simplex_tree <- setMethod("show", "Rcpp_SimplexTree", function (object) {
  max_k <- length(object$n_simplices)
  if (max_k == 0){ cat("< empty simplex tree >\n") }
  else {
    cat(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
  }
})

# ---- as_XPtr ----
#' @name as_XPtr
#' @title Convert to external pointer
#' @description Exports the simplex tree as an \code{externalptr} (i.e. \code{Rcpp::Xptr}) for passing to and from C++ and R. 
#' This method does not register a finalizer. An example is given below using the Rcpp \emph{depends} attribute.
#' @examples 
#' 
#' ## Below is an example of casting an XPtr created in R to a SimplexTree type in C++. 
#' \dontrun{
#' // my_source.cpp
#' #include "Rcpp.h"
#' // [[Rcpp::depends(simplextree)]]
#' #include "simplextree.h"
#' 
#' [[Rcpp::export]]
#' void print_tree(SEXP stree){
#'  Rcpp::XPtr<SimplexTree> stree_ptr(stree);
#'  stree_ptr->print_tree();
#' }
#' }
#' ## Pass to Rcpp as follows
#' st <- simplextree::simplex_tree()
#' print(st$as_XPtr())
NULL

# ---- clear ----
#' @name clear
#' @title Clears the simplex tree
#' @description Removes all simplices from the simplex tree, except the root node.
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3)
#' print(st) ## Simplex Tree with (3, 3, 1) (0, 1, 2)-simplices
#' st$clear()
#' print(st) ## < empty simplex tree >
NULL

# ---- generate_ids ----
#' @name generate_ids
#' @aliases id_policy
#' @title Generates vertex ids.
#' @param n the number of ids to generate. 
#' @description Generates vertex ids representing 0-simplices not in the tree.
#' @details This function generates new vertex ids for use in situations involving, e.g. insertions, 
#' contractions, collapses, etc. There are two 'policies' which designate the generating mechansim  
#' these ids: 'compressed' and 'unique'. 'compressed' generates vertex ids sequentially, starting at 0.
#' 'unique' tracks an incremental internal counter, which is updated on every call to \code{generate_ids}. 
#' The new ids under the 'unique' policy generates the first sequential \code{n} ids that are strictly greater 
#' \code{max}(\emph{counter}, \emph{max vertex id}). \cr
#' \cr
#' 
#' @examples 
#' st <- simplex_tree()
#' st$generate_ids(3) ## 0 1 2
#' st$insert(st$generate_ids(3))
#' print(st$vertices) ## 0 1 2
#' st$insert(st$generate_ids(2))
#' st$print_tree() 
#' st$remove(4)
#' st$generate_ids(1) # 4
NULL

# ---- degree ----
#' @name degree
#' @title The vertex degree.
#' @param ids the vertex ids to check the degree of. 
#' @description Returns the number of edges (degree) for each given vertex id. 
NULL

# ---- traverse ----
#' @name traverse
#' @aliases ltraverse straverse
#' @title traverse
#' @param sigma The simplex to initialize the traversal. See details.  
#' @param f An arbitrary function which accepts as input a simplex. See details. 
#' @param type One of "dfs", "bfs", "cofaces", "star", "link", "skeleton", or "maximal-skeleton".
#' @description Traverses subsets of a simplicial complex.
#' @details \code{\link{traverse}} allows for traversing ordered subsets of the simplex tree. 
#' The simplices within each subset are determined by two aspects: the initial simplex \code{sigma} 
#' and the traversal \code{type}. Given an initial simplex \code{sigma}, \code{traverse} generates 
#' an ordered set of simplices based on the traversal \code{type}, which are iteratively passed to 
#' the supplied function \code{f} as the first argument to \code{f}. \cr
#' \cr
#' \code{sigma} can either be omitted, a simplex, or the \code{empty_face} (which is an alias to NULL).
#' @return NULL; for list or vector-valued returns, use \code{ltraverse} and \code{straverse} respectively.
#' @examples
#' ## Starter example complex 
#' st <- simplex_tree()
#' st$insert(list(1:3, 2:5))
#' 
#' ## Print out complex using depth-first traversal. 
#' ## 'empty_face' implies that the DFS will start at the root. 
#' st$traverse(empty_face, print, "dfs")
#' st$traverse(print, "dfs") ## overload available that assumes start is the empty_face 
#' 
#' ## Print of subtree rooted at vertex 1 using depth-first traversal. 
#' st$traverse(1L, print, "dfs")
#' 
#' ## Print simplices in the star of the edge [4, 5]
#' st$traverse(c(4, 5), print, "star")
#' 
#' ## Traversals can be chained. Here's an example that prints the link of each vertex.
#' st$traverse(function(simplex){
#'   if (length(simplex) == 1){
#'     print(sprintf("Link of %d:", simplex))    
#'     st$traverse(simplex, print, "link")
#'   }
#' }, "bfs")
#' 
#' ## To see the cofaces of a given simplex 
#' st <- simplex_tree()
#' st$insert(c(1, 2, 3))
#' st$traverse(1L, print, "cofaces")
#' 
#' ## Alternatively, collect results into a list 
#' three_cofaces <- st$ltraverse(3L, identity, "cofaces")
NULL

# ---- adjacent ----
#' @name adjacent
#' @title Adjacent vertices.
#' @description Returns a vector of vertex ids that are immediately adjacent to a given vertex.
#' @examples
#' st <- simplex_tree()
#' st$insert(1:3)
#' st$adjacent(2) ## 1 3
NULL

# ---- insert ----
#' @name insert
#' @title Insert simplices
#' @description Inserts simplices into the simplex tree. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @param simplices a list of simplices.
#' @section Usage:
#' st$insert(simplex)
#' st$insert(simplices)
#' @details This function allows insertion of arbitrary order simplices. If the simplex already exists in the tree, 
#' no insertion is made, and the tree is not modified. \code{simplex} is sorted before traversing the trie. 
#' Faces of \code{simplex} not in the simplex tree are inserted as needed.
#' @seealso find remove
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3) ## inserts the 2-simplex { 1, 2, 3 }
#' st$insert(list(4:5, 6)) ## inserts a 1-simplex { 4, 5 } and a 0-simplex { 6 }.
NULL

# ---- remove ----
#' @name remove
#' @title Remove simplices
#' @description Removes simplices from the simplex tree. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex.
#' @section Usage: 
#' st$remove(simplex)
#' @details This function allows removal of a arbitrary order simplices. If \code{simplex} already exists in the tree, 
#' it is removed, otherwise the tree is not modified. \code{simplex} is sorted before traversing the trie.
#' Cofaces of \code{simplex} are also removed.
#' @seealso find remove
NULL

# ---- find ----
#' @name find
#' @title Find simplices
#' @description Finds whether simplices exist the simplex tree.  
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex. Individual simplices are specified as vectors, and a set of simplices as a list of vectors. 
#' @param simplices a list of simplices. 
#' @section Usage:
#' st$find(simplex)
#' st$find(simplices)
#' @details Traverses the simplex tree looking for \code{simplex}, returning whether or not it exists.
#' \code{simplex} can be specified as vector to represent a single simplex, and a list to represent a set of simplices. 
#' \code{simplex} is sorted before traversing the trie.
#' @return boolean indicating whether or not \code{simplex} exists in the tree. 
#' @seealso find remove
NULL

# ---- is_face ----
#' @name is_face
#' @title Is face 
#' @description Checks whether a simplex is a face of another simplex.
#' @param tau a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @param sigma a l-length vector of vertex ids representing a (l-1)-simplex. 
#' @details A simplex \eqn{\tau} is a face of \eqn{\sigma} if \eqn{\tau \subset \sigma}. This function 
#' checks whether that is true. \code{tau} and \code{sigma} are sorted before comparison.
#' @seealso \href{https://en.cppreference.com/w/cpp/algorithm/includes}{std::includes}
#' @return boolean indicating whether \code{tau} is a face of \code{sigma}. 
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3)
#' st$is_face(2:3, 1:3)
#' st$is_face(1:3, 2:3)
NULL

# ---- collapse ----
#' @name collapse
#' @title Elementary collapse
#' @description Performs an elementary collapse. 
#' @param tau a k-length vector of vertex ids representing a (k-1)-simplex. Must be a face of \code{sigma}.
#' @param sigma a n-length vector of vertex ids representing a (n-1)-simplex. Must be the only coface of \code{tau}.
#' @param u a vertex id representing one of the vertices in the free pair.
#' @param v a vertex id representing one of the vertices in the free pair. 
#' @param w a vertex id representing the target of the collapse.
#' @section Usage:
#' st$collapse(tau, sigma)
#' st$collapse(u, v, w)
#' @details This function provides two types of \emph{elementary collapses}. \cr 
#' \cr 
#' The first type of collapse is in the sense described by (1), which is 
#' summarized here. A simplex \eqn{\sigma} is said to be collapsible through one of its faces \eqn{\tau} if 
#' \eqn{\sigma} is the only coface of \eqn{\tau} (excluding \eqn{\tau} itself). This function checks whether its possible to collapse \eqn{\sigma} through \eqn{\tau}, 
#' (if \eqn{\tau} has \eqn{\sigma} as its only coface), and if so, both simplices are removed. 
#' \code{tau} and \code{sigma} are sorted before comparison.
#' To perform this kind of elementary collapse, call \code{collapse} with two simplices as arguments, i.e. \code{tau} before \code{sigma}.
#' 
#' Alternatively, this method supports another type of elementary collapse, also called a \emph{vertex collapse}, as described 
#' in (2). This type of collapse maps a pair of vertices into a single vertex. To use this collapse, specify three vertex ids, the first 
#' two representing the free pair, and the last representing the target vertex to collapse to. 
#' 
#' @return boolean indicating whether the collapse was performed. 
#' @references 
#' 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' 
#' 2. Dey, Tamal K., Fengtao Fan, and Yusu Wang. "Computing topological persistence for simplicial maps." Proceedings of the thirtieth annual symposium on Computational geometry. ACM, 2014.
#' 
#' @examples 
#' st <- simplextree::simplex_tree()
#' st$insert(1:3)
#' st$collapse(1:2, 1:3)
#' st$print_tree()
#' # 1 (h = 1): .( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0):
#' 
#' st$insert(list(1:3, 2:5))
#' st$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 3): .( 3 4 5 )..( 4 5 5 )...( 5 )
#' # 3 (h = 2): .( 4 5 )..( 5 )
#' # 4 (h = 1): .( 5 )
#' # 5 (h = 0): 
#' st$collapse(2:4, 2:5)
#' st$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 2): .( 3 4 5 )..( 5 5 )
#' # 3 (h = 2): .( 4 5 )..( 5 )
#' # 4 (h = 1): .( 5 )
#' # 5 (h = 0): 
NULL

# ---- contract ----
#' @name contract
#' @title Edge contraction
#' @description Performs an edge contraction. 
#' @param edge an edge to contract, as a 2-length vector. 
#' @section Usage: 
#' st$contract(edge)
#' @details This function performs an \emph{edge contraction} in the sense described by (1), which is 
#' summarized here. Given an edge \eqn{ {va, vb}}, \eqn{vb} is contracted to \eqn{va} if \eqn{vb} is 
#' removed from the complex and the link of \eqn{va} is augmented with the link of \eqn{vb}. This may be thought as 
#' applying the mapping: \cr
#' \deqn{f(u) = va}
#' if \eqn{u = vb}
#' and identity otherwise, to all simplices in the complex. \cr 
#' \code{edge} is \strong{not} sorted prior to contraction: the second vertex of the edge is always contracted to the first. 
#' Note that edge contraction is not symmetric.
#' @references 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples 
#' st <- simplextree::simplex_tree()
#' st$insert(1:3)
#' st$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0): 
#' st$contract(c(1, 3))
#' st$print_tree()
#' # 1 (h = 1): .( 2 )
#' # 2 (h = 0): 
NULL


# ---- serialize / deserialize ----
#' @name serialize
#' @aliases deserialize
#' @title Serializes/Deserializes the simplex tree. 
#' @description Provides basic deserialization interface for the simplex tree.
#' @param x a list of simplices to insert into the tree.
#' @details Saves the simplex tree as a compressed RDS file with \code{\link{saveRDS}}. Only the (generally higher order) 
#' simplices which have themselves as a unique coface are saved. 
#' @examples 
#' st <- simplex_tree()
#' st$insert(c(1, 2, 3))
#' tmp <- st$serialize()
#' print(tmp)
#' # [[1]]
#' # [1] 1 2 3
#' st$clear()
#' st$deserialize(tmp)
#' st$print_tree()
NULL

# ---- save ----
#' @name save 
#' @aliases load
#' @title Saves/loads the simplex tree to a file
#' @param filename the filename to save/load the simplex tree to/from. 
#' @description Allows saving/loading the simplex tree with the help of \code{\link{readRDS}}. 
#' @details Both saving and loading requires a filename to save the tree to. Loading 
#' @seealso serialize deserialize
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3)
#' tf <- tempfile()
#' st$print_tree()
#' st$save(tf)
#' st$clear()
#' print(st)
#' # < empty simplex tree >
#' st$load(tf)
#' st$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0):
NULL

# ---- reindex ----
#' @name reindex 
#' @title reindexes vertex ids
#' @param ids Either a vector of new ids, or a named list mapping olds ids to new ids. See details. 
#' @description This function allows one to 'reorder' or 'reindex' vertex ids.  
#' @details The \code{ids} parameter can either be an integer vector or a list. If it's an integer
#' vector, it must be the same length as the number of vertices. The simplex tree is then modified 
#' to replace the vertex label at index \code{i} with \code{ids}[i]. See examples. \cr
#' \cr
#' Alternatively, a named list whose names yield the vertex ids to map from and whose values yield the vertex
#' ids to map to. 
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3)
#' st$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0):
#' st$reindex(4:6)
#' # 4 (h = 2): .( 5 6 )..( 6 )
#' # 5 (h = 1): .( 6 )
#' # 6 (h = 0):
#' st$reindex(list("5"=7))
#' # 4 (h = 2): .( 6 7 )..( 7 )
#' # 6 (h = 1): .( 7 )
#' # 7 (h = 0): 
NULL


# ---- is_tree ----
#' @name is_tree 
#' @title Checks if the simplicial complex is a tree.
#' @description This function performs a breadth-first search on the simplicial complex, checking if the complex is acyclic.
#' @examples 
#' st <- simplex_tree()
#' st$insert(list(1:2, 2:3))
#' st$is_tree() # true
#' st$insert(c(1, 3))
#' st$is_tree() # false
NULL

# ---- plot.Rcpp_SimplexTree ----
#' @name plot.simplextree
#' @title Plots the simplex tree
#' @param x a simplex tree.
#' @param coords Optional (n x 2) matrix of coordinates, where n is the number of 0-simplices. 
#' @param vertex_opt Optional parameters to modify default vertex plotting options. Passed to \code{\link[graphics]{points}}.
#' @param text_opt Optional parameters to modify default vertex text plotting options. Passed to \code{\link[graphics]{text}}.
#' @param edge_opt Optional parameters to modify default edge plotting options. Passed to \code{\link[graphics]{segments}}.
#' @param polygon_opt Optional parameters to modify default k-simplex plotting options for k > 1. Passed to \code{\link[graphics]{polygon}}.
#' @param color_pal Optional vector of colors. See details.
#' @param maximal Whether to draw only the maximal faces of the complex. Defaults to true. 
#' @param by_dim Whether to apply (and recycle or truncate) the color palette to the dimensions rather than to the individual simplices. Defaults to true.
#' @param add Whether to add to the plot or redraw. Defaults to false. See details.
#' @param ... unused
#' @details This function allows generic plotting of simplicial complexes using base \code{\link[graphics:graphics-package]{graphics}}.\cr
#' \cr
#' All parameters passed via list to \code{vertex_opt}, \code{text_opt}, \code{edge_opt}, \code{polygon_opt} 
#' override default parameters and are passed to \code{\link[graphics]{points}}, \code{\link[graphics]{text}}, \code{\link[graphics]{segments}}, 
#' and \code{\link[graphics]{polygon}}, respectively.\cr
#' \cr
#' If \code{add} is true, the plot is not redrawn. \cr
#' \cr
#' If \code{maximal} is true, only the maximal simplices are drawn. \cr
#' \cr
#' The \code{color_pal} argument controls how the simplicial complex is colored. It can be specified in multiple ways.
#' \enumerate{
#'   \item A vector of colors of length \emph{dim+1}, where \emph{dim}=\code{x$dimension}
#'   \item A vector of colors of length \emph{n}, where \emph{n}=\code{sum(x$n_simplices)}
#'   \item A named list of colors
#' }
#' Option (1) assigns every simplex a color based on its dimension. \cr
#' \cr
#' Option (2) assigns each individual simplex a color. The vector must be specified in level-order 
#' (see \code{\link{ltraverse}} or examples below). \cr
#' \cr
#' Option (3) allows specifying individual simplices to draw. It expects a named list, where the names
#' must correspond to simplices in \code{x} as comma-separated strings and whose values are colors. If 
#' option (3) is specified, this method will \emph{only} draw the simplices given in \code{color_pal}.\cr
#' \cr
#' If \code{length(color_pal)} does not match the dimension or the number of simplices in the complex, 
#' the color palette is recyled and simplices are as such. 
#' @importFrom utils modifyList
#' @examples 
#' ## Simple 3-simplex 
#' st <- simplex_tree()
#' st$insert(list(1:4))
#' 
#' ## Default is categorical colors w/ diminishing opacity
#' plot(st)
#' 
#' ## If supplied colors have alpha defined, use that 
#' vpal <- viridis::viridis(st$dimension + 1)
#' plot(st, color_pal = vpal)
#' 
#' ## If alpha not supplied, decreasing opacity applied
#' plot(st, color_pal = substring(vpal, first=1, last=7))
#' 
#' ## Bigger example; observe only maximal faces (+vertices and edges) are drawn
#' st <- simplex_tree()
#' st$insert(list(1:3, 2:5, 5:9, 7:8, 10))
#' plot(st, color_pal = viridis::viridis(st$dimension + 1))
#' 
#' ## If maximal == FALSE, every simplex is drawn (even on top of each other)
#' vpal <- viridis::viridis(st$dimension + 1)[c(1,2,5,4,3)]
#' pal_alpha <- c(1, 1, 0.2, 0.35, 0.35)
#' vpal <- sapply(seq_along(vpal), function(i) adjustcolor(vpal[i], alpha.f = pal_alpha[i]))
#' plot(st, color_pal = vpal, maximal = FALSE)
#' 
#' ## You can also color each simplex individually by supplying a vector 
#' ## of the same length as the number of simplices. 
#' plot(st, color_pal = sample(rainbow(sum(st$n_simplices))))
#' 
#' ## The order is assumed to follow the level order traversal (first 0-simplices, 1-, etc.)
#' ## This example colors simplices on a rainbow gradient based on the sum of their labels
#' si_sum <- unlist(st$ltraverse(sum, "bfs"))[-1] # -1 to remove empty_face
#' rbw_pal <- rev(rainbow(50, start=0,end=4/6))
#' plot(st, color_pal=rbw_pal[cut(si_sum, breaks=50, labels = FALSE)])
#' 
#' ## This also makes highlighting simplicial operations fairly trivial 
#' four_cofaces <- st$ltraverse(4, identity, "cofaces")
#' coface_pal <- st$ltraverse(function(simplex){
#'   ifelse(list(simplex) %in% four_cofaces, "orange", "blue")
#' }, "bfs")
#' plot(st, color_pal=unlist(coface_pal)[-1])
#' 
#' ## You can also give a named list to draw individual simplices. 
#' ## **Only the maximal simplices in the list are drawn** 
#' blue_vertices <- structure(as.list(rep("blue", 5)), names=as.character(seq(5, 9)))
#' plot(st, color_pal=append(blue_vertices, list("5,6,7,8,9"="red")))
#' 
#' ## You can specify add=TRUE to add to the current plot
#' ## This makes e.g. animations easier. 
#' \dontrun{
#'   ## Fix coordinates
#'   g <- igraph::graph_from_adjacency_matrix(st$as_adjacency_matrix())
#'   coords <- igraph::layout.auto(g)
#'   
#'   ## Create rainbow colors 
#'   si_to_str <- function(simplex) { paste0(simplex, collapse=",") }
#'   si_names <- sapply(st$ltraverse(identity, "bfs"), si_to_str)[-1]
#'   si_colors <- structure(as.list(rainbow(sum(st$n_simplices))), names=si_names)
#'   
#'   ## Make animation
#'   animation::saveGIF({
#'     for (i in seq(sum(st$n_simplices))){
#'       plot(st, coords=coords, color_pal=si_colors[1L:i])
#'     }
#'   }, movie.name = "si_animation.gif", interval=0.2)
#' }
#' @export
plot.Rcpp_SimplexTree <- function (x, coords = NULL, vertex_opt=NULL, text_opt=NULL, edge_opt=NULL, polygon_opt=NULL, color_pal=NULL, maximal=TRUE, by_dim=TRUE, add=FALSE, ...) {
  stopifnot(methods::is(x, "Rcpp_SimplexTree"))
  if (sum(x$n_simplices) == 0){ graphics::plot.new(); return() } 

  ## Default color palette; categorical diverging if (# colors) <= 9, o/w rainbow
  if (missing(color_pal) || is.null(color_pal)){  
    n_colors <- if (by_dim) x$dimension+1 else sum(x$n_simplices)
    if (n_colors <= 9){
      color_pal <- .default_st_colors[seq(n_colors)]
    } else {
      color_pal <- substr(rev(grDevices::rainbow(n_colors, start=0, end=4/6)), start=1,stop=7)
    }
  }
  
  ## Regardless of type of palette given, the result is parsed into a vector of hexadecimal colors 
  ## or length (# simplices) in breadth-first order, not including the empty face. 
  simplex_colors <- NULL # placeholder
  draw_simplex <- rep(TRUE, sum(x$n_simplices))
  dim_idx <- (unlist(x$ltraverse(length, "bfs"))[-1]-1L)
  is_char_vec <- all(is.character(color_pal))
  is_in <- function(lst){ 
    return(function(element) { !is.null(element) && (list(as.integer(element)) %in% lst) })
  }
  
  ## If the maximal faces are requested, set non-maximal `draw_simplex` indices to FALSE 
  if (maximal){
    maximal_faces <- lapply(x$serialize(), as.integer)
    draw_simplex <- unlist(x$ltraverse(is_in(maximal_faces), "bfs"))[-1]
    draw_simplex[dim_idx %in% c(0L, 1L)] <- TRUE ## always draw points and edges
  }
  
  ## Converts non-hex colors to hex. Additionally, any 7-length hex has 
  ## simplex dimension opacity scaling applied to it
  col_to_hex <- function(cp){
    is_hex <- (substring(cp,first=1,last=1) == "#")
    is_rgb <- (nchar(cp) == 7)
    is_col <- (!is_hex | (is_hex & is_rgb))
    cp[is_col] <- apply(grDevices::col2rgb(cp[is_col]), 2, function(col){ do.call(grDevices::rgb, as.list(col/255)) })
    cp[is_col] <- alpha4sc(cp)[is_col]
    return(cp)
  }
  
  ## Case 1: color_pal is a named list where each name is a comma-separated simplex 
  if (is.list(color_pal)){
    stopifnot(is.character(names(color_pal)))
    
    ## Extract simplices in names. Check named labels are ordered + simplices exist.
    simplices <- lapply(lapply(strsplit(names(color_pal), ","), as.integer), sort)
    names(color_pal) <- sapply(simplices, function(simplex){ paste0(simplex, collapse=",") })
    stopifnot(all(x$find(simplices)))
    
    ## Color named simplex w/ color if given, otherwise use default
    si_in <- is_in(simplices)
    si_color <- function(simplex){ ifelse(!is.null(simplex) && si_in(simplex), color_pal[[paste0(simplex, collapse=",")]], NA) }
    draw_simplex <- unlist(x$ltraverse(si_in, "bfs"))[-1]
    simplex_colors <- unlist(x$ltraverse(si_color, "bfs"))[-1]
  } else if (is_char_vec && (length(color_pal) == sum(x$n_simplices))){
    ## Case 2: color_pal is character vector w/ length == # simplices
    simplex_colors <- col_to_hex(color_pal)
  } else if (is_char_vec && ((length(color_pal) == x$dimension+1L) || by_dim)){
    ## Case 3: color_pal is character vector, recycled to dimensions
    color_pal <- rep(color_pal, length.out = x$dimension+1L)
    simplex_colors <- col_to_hex(color_pal)[dim_idx+1L]
  } else if (is_char_vec){
    simplex_colors <- rep(color_pal, length.out=sum(x$n_simplices))
  } else {
    stop("Invalid color palette given. Must be either a character vector or named list. See `?plot.simplextree`.")
  }
  
  ## Get coordinates of vertices
  if (!missing(coords)){ stopifnot(is.matrix(coords) && all(dim(coords) == c(x$n_simplices[1], 2))) }
  else {
    requireNamespace("igraph", quietly = TRUE)
    g <- igraph::graph_from_adjacency_matrix(x$as_adjacency_matrix())
    coords <- igraph::layout_with_fr(g)
  }
  
  ## Create a new plot by default unless specified otherwise 
  if (!add){
    graphics::plot.new()
    graphics::plot.window(xlim=range(coords[,1]), ylim=range(coords[,2]))
  }
  
  # plot polygons for simplices of dimension 2+; omits edges and vertices
  # this just plots the triangles
  v <- x$vertices # cache vertices
  if (x$dimension >= 2L){
    for (d in seq(x$dimension, 2)){
      if (any(draw_simplex[dim_idx == d])){
        polys <- x$ltraverse(empty_face, function(simplex){ 
          # idx <- match(simplex[combn(d+1L, 3L)], v)
          poly <- coords[match(simplex, v),] 
          rbind(poly[grDevices::chull(poly),], c(NA, NA))
        }, "maximal-skeleton", list(k=d))
        subset <- (draw_simplex & (dim_idx==d))
        d_subset <- subset[dim_idx==d]
        params <- list(x=do.call(rbind, polys[d_subset]), border=NA, col=simplex_colors[subset])
        do.call(graphics::polygon, modifyList(params, as.list(polygon_opt)))
      }
    }
  }
  # plot segments for edges
  if (length(x$n_simplices) >= 2 && any(draw_simplex[dim_idx == 1L])){
    lc <- apply(x$edges, 1, function(e){ t(coords[match(e, x$vertices),,drop=FALSE]) })
    subset <- (draw_simplex & (dim_idx==1L))
    e_subset <- subset[dim_idx==1L]
    params <- list(x0=lc[1,e_subset], y0=lc[2,e_subset], x1=lc[3,e_subset], y1=lc[4,e_subset], lwd=2, col=simplex_colors[subset])
    do.call(graphics::segments, modifyList(params, as.list(edge_opt)))
  }
  # plot vertices
  if (length(x$n_simplices) >= 1 && any(draw_simplex[dim_idx == 0L])){
    subset <- (draw_simplex & (dim_idx==0L))
    v_subset <- subset[dim_idx==0L]
    do.call(graphics::points, modifyList(list(x=coords[v_subset,,drop=FALSE], pch=21, bg=simplex_colors[subset], col=simplex_colors[subset], cex=2), as.list(vertex_opt)))
    do.call(graphics::text, modifyList(list(x=coords[v_subset,,drop=FALSE], labels=as.character(x$vertices)[v_subset], col="white", cex=0.75), as.list(text_opt))) 
  }
}

.default_st_colors <- c("#e41a1c", "#377eb8", "#ffff33", "#984ea3", "#ff7f00", "#4daf4a", "#a65628", "#f781bf", "#999999")

# Adjusts the simplex alphas for each dimension; expects hexadecimal
alpha4sc <- function(col_pal) {
  nc <- length(col_pal)
  if (nc == 0) { return(col_pal) }
  ext <- if (nc > 2){ seq(0.80, 0.45, length.out = nc-2) } else { NULL }
  si_alpha <- c(1, 1, ext) 
  sapply(seq_along(col_pal), function(i){ grDevices::adjustcolor(col_pal[i], alpha.f = si_alpha[i]) })
}

# ids <- apply(utils::combn(d+1L, 2), 2, function(i){ simplex[i] })
# apply(ids, 2, function(c_id){
#   idx <- match(c_id, v)
#   coords[idx,,drop=FALSE]
#   
# })
# #st$ltraverse(empty_face, identity, "bfs")
# stopifnot(is.character(color_pal))
# if (length(color_pal) == 1){ color_pal <- rep(color_pal, x$dimension+1L) }
# 
# ## Final check:
# stopifnot(length(color_pal) == x$dimension+1L)
