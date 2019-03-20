#' Simplex Tree
#' @aliases simplex_tree simplextree SimplexTree
#' @description Simplex tree class exposed as an Rcpp Module.
#' @docType class 
#' @details A simplex tree is an ordered trie-like structure specialized for storing and doing general computation 
#' simplicial complexes. Here is figure of a simplex tree, taken from the original paper (see 1): \cr  
#' \if{html}{\figure{simplextree.png}{options: width="80\%" alt="Figure: simplextree.png"}}
#' \if{latex}{\figure{simplextree.pdf}{options: width=12cm}}
#' \cr 
#' The current implementation provides a limited API and a subset of the functionality described in the paper.
#' @section Methods: 
#' \describe{
#'     \item{$\code{print_tree}}{ Prints the simplex tree. }
#'     \item{$\code{\link[simplextree:apply.simplex_tree]{apply}}()}{ Applies a function to a subset of the simplex tree. }
#'     \item{$\code{\link{insert_simplex}}}{ Inserts a simplex into the trie, if it doesn't exist. }
#'     \item{$\code{\link{remove_simplex}}}{ Removes a simplex from the trie, if it exists. }
#'     \item{$\code{\link{find_simplex}}}{ Searches the trie for a simplex. }
#'     \item{$\code{\link{collapse}}}{ Performs an elementary collapse. }
#'     \item{$\code{\link{contract}}}{ Performs an edge contraction. }
#'     \item{$\code{\link{is_face}}}{ Checks for faces. }
#'     \item{$\code{\link{serialize}}}{ Serializes the simplex tree. }
#'     \item{$\code{\link{unserialize}}}{ Unserializes a stored simplex tree. }
#'     \item{$\code{as_list}}{ Converts the complex to a list. }
#'     \item{$\code{as_adjacency_matrix}}{ Converts the 1-skeleton to an adjacency matrix. }
#'     \item{$\code{as_adjacency_list}}{ Converts the 1-skeleton to an adjacenecy list. }
#'     \item{$\code{as_edgelist}}{ Converts the 1-skeleton to an edgelist. }
#' }
#' @field n_simplexes A vector, where each index k denotes the number (k-1)-simplices.
#' @field max_depth The maximum height of the tree. The root of the tree has height 0, vertices have height 1, etc.
#' @author Matt Piekenbrock
#' @return A queryable simplex tree, as a \code{Rcpp_SimplexTree} object (Rcpp module). 
#' @references 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples
#' ## Recreating simplex tree from figure. 
#' stree <- simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$insert_simplex(c(2, 3, 4, 5))
#' stree$insert_simplex(c(5, 6, 9))
#' stree$insert_simplex(c(7, 8))
#' stree$insert_simplex(10)
#' @export
simplex_tree <- function(){
  return(new(SimplexTree))
}

setClass("Rcpp_SimplexTree")
.print_simplex_tree <- setMethod("show", "Rcpp_SimplexTree", function (object) {
  max_k <- length(object$n_simplexes)
  if (max_k == 0){ cat("< empty simplex tree >\n") }
  else {
    cat(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplexes, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
  }
})

#' @name apply.simplex_tree
#' @title apply
#' @param sigma The simplex to initialize the traversal, or NULL to use the root. See details.  
#' @param f An arbitrary function which accepts as input a simplex. See details. 
#' @param type One of "dfs", "bfs", "cofaces", "star", or "link"
#' @description Apply operations across subsets of a simplicial complex.
#' @details \code{apply} allows for traversing subsets of the simplex tree. 
#' A subset of the simplex tree is represented by a set of simplices. The simplices within each subset is determined by
#' two aspects: the traversal \code{type} and the initial simplex \code{sigma}. Given a simplex \code{sigma}, a subset of 
#' the simplex tree is generated based on \code{type}, and then each simplex is passed as the first argument to \code{f}.
#' See examples for use-cases. 
#' @examples
#' ## Starter example complex 
#' stree <- simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$insert_simplex(c(2, 3, 4, 5))
#' 
#' ## Print out complex using depth-first traversal. NULL implies that the DFS will start at the root. 
#' stree$apply(NULL, print, "dfs")
#' 
#' ## Print of subtree rooted at vertex 1 using depth-first traversal. 
#' stree$apply(1L, print, "dfs")
#' 
#' ## Print simplices in the star of the edge [4, 5]
#' stree$apply(c(4, 5), print, "star")
#' 
#' ## Traversals can be chained. Here's an example that prints the link of each vertex.
#' stree$apply(NULL, function(simplex){
#'   if (length(simplex) == 1){
#'     print(sprintf("Link of %d:", simplex))    
#'     stree$apply(simplex, print, "link")
#'   }
#' }, "bfs")
#' 
#' ## To see the cofaces of a given simplex 
#' stree <- simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$apply(1L, print, "cofaces")
#' stree$apply(2L, print, "cofaces")
#' stree$apply(3L, print, "cofaces")
NULL

#' @name as_XPtr
#' @title Convert to external pointer
#' @description Exports the simplex tree as an Rcpp::Xptr for e.g. passing from C++ --> R --> C++. Does not register a finalizer. 
NULL

#' @name adjacent_vertices
#' @title Adjacent vertices.
#' @description Returns a vector of vertex ids that are immediately adjacent to a given vertex.
NULL

#' @name insert_simplex
#' @title Insert simplex
#' @description Inserts a simplex. 
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @usage $insert(simplex)
#' @details This function allows insertion of arbitrary order simplices. If the simplex already exists in the tree, 
#' no insertion is made, and the tree is not modified. \code{simplex} is sorted before traversing the trie. 
#' Lower order simplices are inserted as needed if they do not already exist.
NULL

#' @name remove_simplex
#' @title Remove simplex
#' @description Removes a simplex. 
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @usage $remove(simplex)
#' @details This function allows removal of a arbitrary order simplices. If \code{simplex} already exists in the tree, 
#' it is removed, otherwise the tree is not modified. \code{simplex} is sorted before traversing the trie.
#' Higher order simplices which are cofaces of \code{simplex} are also removed.
NULL

#' @name find_simplex
#' @title Remove simplex
#' @description Removes a simplex. 
#' @param simplex a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @usage $find_simplex(simplex)
#' @details Traverses the simplex tree looking for \code{simplex}, returning whether or not it exists. 
#' \code{simplex} is sorted before traversing the trie.
#' @return boolean indicating whether or not \code{simplex} exists in the tree. 
NULL

#' @name is_face
#' @title Is face 
#' @description Checks whether a simplex is a face of another simplex.
#' @param tau a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @param sigma a l-length vector of vertex ids representing a (l-1)-simplex. 
#' @usage $is_face(tau, sigma)
#' @details A simplex \eqn{\tau} is a face of \eqn{\sigma} if \eqn{\tau \subset \sigma}. This function 
#' checks whether that is true. \code{tau} and \code{sigma} are sorted before comparison.
#' @seealso \href{https://en.cppreference.com/w/cpp/algorithm/includes}{std::includes}
#' @return boolean indicating whether \code{tau} is a face of \code{sigma}. 
NULL

#' @name collapse
#' @title Elementary collapse
#' @description Performs an elementary collapse. 
#' @param tau a k-length vector of vertex ids representing a (k-1)-simplex. 
#' @param sigma a l-length vector of vertex ids representing a (l-1)-simplex. 
#' @usage $collapse(tau, sigma)
#' @details This function performs an \emph{elementary collapse} in the sense described by (1), which is 
#' summarized here. A simplex \eqn{\sigma} is said to be collapsible through one of its faces \eqn{\tau} if 
#' \eqn{\sigma} is the only coface of \eqn{\tau} (excluding \eqn{\tau} itself). This function checks whether its possible to collapse \eqn{\sigma} through \eqn{\tau}, 
#' (if \eqn{\tau} has \eqn{\sigma} as its only coface), and if so, both simplices are removed. 
#' \code{tau} and \code{sigma} are sorted before comparison.
#' @return boolean indicating whether \code{tau} is a coface of \code{sigma}, and the two were removed. 
#' @references 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples 
#' stree <- simplextree::simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$collapse(c(1, 2), c(1, 2, 3))
#' stree$print_tree()
#' # 1 (h = 1): .( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0):
#' 
#' stree$insert_simplex(1:3)
#' stree$insert_simplex(2:5)
#' stree$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 3): .( 3 4 5 )..( 4 5 5 )...( 5 )
#' # 3 (h = 2): .( 4 5 )..( 5 )
#' # 4 (h = 1): .( 5 )
#' # 5 (h = 0): 
#' stree$collapse(2:4, 2:5)
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 2): .( 3 4 5 )..( 5 5 )
#' # 3 (h = 2): .( 4 5 )..( 5 )
#' # 4 (h = 1): .( 5 )
#' # 5 (h = 0): 
NULL

#' @name contract
#' @title Edge contraction
#' @description Performs an edge contraction. 
#' @param edge an edge to contract, as a 2-length vector. 
#' @usage $contract(edge)
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
#' stree <- simplextree::simplex_tree()
#' stree$insert_simplex(1:3)
#' stree$print_tree()
#' # 1 (h = 2): .( 2 3 )..( 3 )
#' # 2 (h = 1): .( 3 )
#' # 3 (h = 0): 
#' stree$contract(c(1, 3))
#' stree$print_tree()
#' # 1 (h = 1): .( 2 )
#' # 2 (h = 0): 
NULL


#' @name serialize
#' @title Serializes the simplex tree. 
#' @description Provides basic serialization of the simplex tree with the help of \code{\link{saveRDS}}. 
#' @param filename The file to write the simplex tree too (path). 
#' @usage $serialize(filename)
#' @details Saves the simplex tree as a compressed RDS file with \code{\link{saveRDS}}. Only the (generally higher order) 
#' simplices which have themselves as a unique coface are saved. 
#' @references 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples 
#' stree <- simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$serialize("test.rds")
#' readRDS("test.rds")
#' # [[1]]
#' # [1] 1 2 3
NULL

#' @name unserialize
#' @title Unserializes a simplex tree. 
#' @description Provides basic unserialization of the simplex tree with the help of \code{\link{readRDS}}. 
#' @param filename The file to read the simplex tree from (path). 
#' @usage $unserialize(filename)
#' @details Reads the simplex tree stored in the compressed RDS file given by \code{filename} with \code{\link{readRDS}}, 
#' successively reinserting the simplices into the current tree. 
#' @references 1. Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @examples 
#' stree <- simplex_tree()
#' stree$insert_simplex(c(1, 2, 3))
#' stree$serialize("test.rds")
#' stree <- simplex_tree()
#' print(stree)
#' # < empty simplex tree >
#' stree$unserialize("test.rds")
#' print(stree)
#' # Simplex Tree with (3, 3, 1) (0, 1, 2)-simplices
NULL




# plot.Rcpp_SimplexTree <- function (x) {

#   g <- igraph::graph_from_edgelist(x$as_edgelist(), directed = FALSE)
#   coords <- igraph::layout.auto(g, dim = 2L)
# }
# rev(viridis::viridis(100, alpha = 0.8, begin = 2/6, end = 1))
# .default_st_colors <- c("#FDE725CC","#F9E621CC","#F5E61FCC","#F1E51DCC","#ECE51BCC","#E8E419CC","#E4E419CC","#DFE318CC","#DBE319CC","#D7E219CC","#D2E21BCC","#CDE11DCC","#C9E020CC","#C4E022CC","#C0DF25CC","#BBDE28CC","#B7DE2ACC","#B2DD2DCC","#ADDC30CC","#A9DB33CC","#A4DB36CC","#A0DA39CC","#9BD93CCC","#96D83FCC","#92D741CC","#8ED645CC","#8AD547CC","#85D54ACC","#81D34DCC","#7DD250CC","#78D152CC","#75D054CC","#70CF57CC","#6DCD59CC","#68CD5BCC","#65CB5ECC","#61CA60CC","#5DC863CC","#59C864CC","#56C667CC","#53C569CC","#4FC46ACC","#4CC26CCC","#48C16ECC","#45BF70CC","#41BE71CC","#3FBC73CC","#3BBB75CC","#39BA76CC","#37B878CC","#34B679CC","#31B67BCC","#2FB47CCC","#2DB27DCC","#2BB07FCC","#29AF7FCC","#27AD81CC","#25AC82CC","#24AA83CC","#23A983CC","#22A785CC","#21A585CC","#20A486CC","#1FA287CC","#1FA188CC","#1F9F88CC","#1F9E89CC","#1E9C89CC","#1F9A8ACC","#1F998ACC","#1F978BCC","#1F958BCC","#20938CCC","#20928CCC","#21918CCC","#218F8DCC","#228D8DCC","#228C8DCC","#238A8DCC","#23888ECC","#24878ECC","#25858ECC","#25838ECC","#26828ECC","#26818ECC","#277F8ECC","#287D8ECC","#287C8ECC","#297A8ECC","#2A788ECC","#2A768ECC","#2B758ECC","#2C738ECC","#2C718ECC","#2D718ECC","#2E6F8ECC","#2E6D8ECC","#2F6B8ECC","#306A8ECC","#31688ECC")