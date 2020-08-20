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
#'     \item{$\code{as_XPtr}}{ Creates an external pointer. }
#'     \item{$\code{clear}}{ Clears the simplex tree. }
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
#'     \item{$\code{as_list}}{ Converts the simplicial complex to a list. }
#'     \item{$\code{as_adjacency_matrix}}{ Converts the 1-skeleton to an adjacency matrix. }
#'     \item{$\code{as_adjacency_list}}{ Converts the 1-skeleton to an adjacency list. }
#'     \item{$\code{as_edgelist}}{ Converts the 1-skeleton to an edgelist. }
#' }
#' @author Matt Piekenbrock
#' @import methods
#' @param simplices optional simplices to initialize the simplex tree with. See \code{\link{insert}}.
#' @return A queryable simplex tree, as a \code{Rcpp_SimplexTree} object (Rcpp module). 
#' @references Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
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
#' @export
simplex_tree <- function(simplices = NULL){
  st <- new(SimplexTree)
  # assign('insert', function(x, check=TRUE) { return(st$insert(x, check)) }, envir = st)
  if (!missing(simplices)){ st %>% insert(simplices) }
  return(st)
}

#' flag
#' @description Creates a filtration of flag complexes
#' @param st a simplex tree. See details. 
#' @param d a vector of edge weights, or a 'dist' object. 
#' @details A flag complex is a simplicial complex whose k-simplices for k >= 2 are completely determined 
#' by edges/graph of the complex. This function creates filtered simplicial complex using the supplied edge 
#' weights. The resulting complex is a simplex tree object endowed with additional structure; see. 
#' Vertices have their weights set to 0, and k-simplices w/ k >= 2 have their weights set to the maximum
#' weight of any of its edges. 
#' @export
flag <- function(st, d){
  stopifnot(is.numeric(d), class(st) %in% .st_classes)
  fi <- new(Filtration)
  fi$init_tree(st$as_XPtr())
  fi$flag_filtration(d)
  return(fi)
}

# ---- empty_face ----
#' empty_face 
#' @description Alias to the empty integer vector (integer(0L)). Used to indicate the empty face of the tree. 
#' @seealso traverse
#' @export
empty_face <- integer(0L)

.traversal_types = c("Preorder", "Level order", "Face", "Coface", "Coface roots", "K-skeleton", "K-simplices", "Maximal simplex", "Link")

# ---- print.st_traversal ----
#' print.st_traversal
#' @param x traversal object.
#' @param ... unused. 
#' @export
print.st_traversal <- function(x, ...){
  sigma_str <- ifelse(length(x$sigma) == 0 || is.null(x$sigma), "empty face", paste0(x$sigma, collapse = " "))
  tt <- .traversal_types[x$traversal_type+1L]
  writeLines(sprintf("%s traversal @ { %s }", tt, sigma_str))
}

# ---- as.list.st_traversal ----
#' as.list.st_traversal
#' @param x traversal object.
#' @param ... unused. 
#' @export
as.list.st_traversal <- function(x, ...){
  return(ltraverse(x, identity))
}


# ---- traverse ----
#' @name traverse
#' @title traverse
#' @param traversal The type of traversal to use.
#' @param f An arbitrary function to apply to eac simplex of the traversal. See details. 
#' @param ... unused. 
#' @description Traverses specific subsets of a simplicial complex.
#' @details \code{\link{traverse}} allows for traversing ordered subsets of the simplex tree. 
#' The specific subset and order are determined by the choice of \emph{traversal}: examples include 
#' the \code{\link{preorder}} traversal, the \code{\link{cofaces}} traversal, etc. See the links below. 
#' Each simplex in the traversal is passed as the first and only argument to \code{f}, one per simplex in the traversal.
#' \code{\link{traverse}} does nothing with the result; if you want to collect the results of applying \code{f} to each simplex 
#' into a list, use \code{\link{ltraverse}} (or \code{\link{straverse}}), which are meant to be used like \code{\link{lapply}} 
#' and \code{\link{sapply}}, respectively. 
#' @family traversals 
#' @return NULL; for list or vector-valued returns, use \code{ltraverse} and \code{straverse} respectively.
#' @examples
#' ## Starter example complex 
#' st <- simplex_tree()
#' st %>% insert(list(1:3, 2:5))
#' 
#' ## Print out complex using depth-first traversal. 
#' st %>% preorder() %>% traverse(print)
#' 
#' ## Collect the last labels of each simplex in the tree. 
#' last_labels <- st %>% preorder() %>% straverse(function(simplex){ tail(simplex, 1) })
#' @export
traverse <- function(traversal, f, ...){
  stopifnot("st_traversal" %in% class(traversal))
  # if (missing(f)){ return(function(traversal, f){ traverse_R(traversal, f) }) }
  # if (is.function(traversal)){ traversal(f) }
  traverse_R(traversal, f)
}

#' @name straverse 
#' @param traversal the type of traversal. 
#' @param f the function to apply to each simplex. 
#' @param ... unused. 
#' @rdname traverse
#' @export
straverse <- function(traversal, f, ...){
  stopifnot("st_traversal" %in% class(traversal))
  stopifnot(is.function(f))
  # if (missing(f)){ return(function(traversal, f){ straverse_R(traversal, f) }) }
  # if (is.function(traversal)){ traversal(f) }
  return(straverse_R(traversal, f))
}

#' ltraverse 
#' @param traversal the type of traversal. 
#' @param f the function to apply to each simplex. 
#' @param ... unused. 
#' @rdname traverse
#' @export
ltraverse <- function(traversal, f, ...){
  stopifnot("st_traversal" %in% class(traversal))
  stopifnot(is.function(f))
  # if (missing(f)){ return(function(traversal, f){ ltraverse_R(traversal, f) }) }
  # if (is.function(traversal)){ traversal(f) }
  return(ltraverse_R(traversal, f))
}

# ---- preorder ----- 
#' @name preorder 
#' @title Generates a preorder traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
preorder <- function(st, sigma = NULL){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "preorder", NULL)
}

# ---- level_order ----- 
#' @name level_order 
#' @title Generates a level order traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
level_order <- function(st, sigma = NULL){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "level_order", NULL)
}

# ---- faces ----- 
#' @name faces 
#' @title Generates a face traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
faces <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "faces", NULL)
}

# ---- cofaces ----- 
#' @name cofaces 
#' @title Generates a coface traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
cofaces <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "cofaces", NULL)
}

# ---- k_skeleton ----- 
#' @name k_skeleton 
#' @title Generates a k-skeleton traversal on the simplex tree.
#' @param st the simplex tree to traverse.
#' @param k the dimension of the skeleton to include.
#' @param sigma simplex to start the traversal at. 
#' @export
k_skeleton <- function(st, k, sigma = NULL){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "k_skeleton", list(k=k))
}

# ---- coface_roots ----- 
#' @name coface_roots 
#' @title Generates a coface roots traversal on the simplex tree.
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @description The coface roots of a given simplex \code{sigma} are the roots of subtrees 
#' in the trie whose descendents (including the roots themselves) are cofaces of \code{sigma}.
#' This traversal is more useful when used in conjunction with other traversals, e.g. a \emph{preorder} 
#' or \emph{level_order} traversal at the roots enumerates the cofaces of \code{sigma}.
#' @export
coface_roots <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "coface_roots", NULL)
}


# ---- maximal ----- 
#' @name maximal 
#' @title Generates a traversal on the maximal of the simplex tree.
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
maximal <- function(st, sigma = NULL){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "maximal", NULL)
}

# ---- k_simplices ----- 
#' @name k_simplices 
#' @title Generates a traversal on the k-simplices of the simplex tree.
#' @param st the simplex tree to traverse.
#' @param k the dimension of the skeleton to include.
#' @param sigma simplex to start the traversal at. 
#' @export
k_simplices <- function(st, k, sigma = NULL){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "k_simplices", list(k=k))
}

# ---- link ----- 
#' @name link 
#' @title Generates a traversal on the link of a given simplex in the simplex tree.
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
link <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "link", NULL)
}


# setClass("Rcpp_SimplexTree")
# .format_simplex_tree <- setMethod("format", "Rcpp_SimplexTree", function (object) {
#   max_k <- length(object$n_simplices)
#   if (max_k == 0){ return("< empty simplex tree >") }
#   else {
#     return(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
#   }
# })

# format.Rcpp_SimplexTree <- function(x){
#   max_k <- length(x$n_simplices)
#   if (max_k == 0){ return("< empty simplex tree >") }
#   else {
#     return(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(x$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
#   }
# }
# 
# format.Rcpp_Filtration<- function(x){
#   paste0(format.Rcpp_SimplexTree(x), sprintf("Current filtration index: %d", x$current_index), collapse = "\n")
# }

## One printer to rule them all
# ----- print method ------
#' @name print_simplices
#' @title Print simplices to the console
#' @description Prints simplices in a formatted way 
#' @param st a simplex tree. 
#' @param format the choice of how to format the printing. See details.   
#' @description Prints a traversal, a simplex tree, or a list of simplices to the R console, with 
#' options to customize how the simplices are printed. The \code{format} must be one of 
#' "summary", "tree", "cousins", "short", "column", or "row", with the default being "short".
#' In general, the "tree" and "cousins" format give more details on the structure of the trie, 
#' whereas the other formats just change how the given set of simplices are formatted.
#' \cr
#' The "tree" method prints the nodes grouped by the same last label and indexed by depth.
#' The printed format is: \cr 
#' \cr
#' [vertex] (h = [subtree height]): [subtree depth]([subtree]) \cr 
#' \cr
#' Where each lists the top node (\emph{vertex}) and its corresponding subtree. The 
#' \emph{subtree height} displays the highest order k-simplex in that subtree. Each 
#' level in the subtree tree is a set of sibling k-simplices whose order is given  
#' by the number of dots ('.') proceeding the print level.\cr 
#' \cr
#' The "cousin" format prints the simplex relations used by various algorithms to speed 
#' up finding adjacencies in the complex. The cousins are grouped by label and depth. \cr 
#' The format looks like: 
#' \cr
#' (last=[label], depth=[depth of label]): [simplex] \cr
#' \cr
#' This function is useful for understanding how the simplex tree is stored, and for debugging purposes. 
#' @export
print_simplices <- function (st, format=c("summary", "tree", "cousins", "short", "column", "row")){
  if (missing(format)){ format <- "short" } 
  if (format == "summary" && (class(st) %in% .st_classes)){ show(st) }
  else if (format == "tree" && (class(st) %in% .st_classes)){ st$print_tree() }
  else if (format == "cousins" && (class(st) %in% .st_classes)){ st$print_cousins()}
  else {
    if (is.list(st)){
      stopifnot(all(sapply(st, is.numeric)))
      simplex_str <- lapply(st, as.character)
    } else if ("st_traversal" %in% class(st)){
      simplex_str <- straverse(st, as.character)
    } else if (class(st) %in% .st_classes){
      simplex_str <- straverse(level_order(st), as.character)
    } else {
      stop("Unknown type of 'st' passed in.")
    }
    
    if (format == "short"){
      format_simplex <- function(sigma){ paste0(sigma, collapse = " ") }
      writeLines(paste0(sapply(simplex_str, format_simplex), collapse = ", "))
    } else if (format == "column"){
      d <- max(sapply(simplex_str, length))
      simplices_str <- sapply(seq(d), function(i){
        paste0(sapply(simplex_str, function(labels){ 
          width <- max(sapply(labels, nchar))
          ifelse(length(labels) < i, paste0(rep(" ", width), collapse=""), sprintf(paste0("%", width, "d"), as.integer(labels[i])))
        }), collapse = " ")
      })
      writeLines(simplices_str)
    } else if (format == "row"){
      writeLines(sapply(simplex_str, function(sigma){ paste0(sigma, collapse = " ") }))
    } else {
      stop("Unknown format specified.")
    }
  } 
}

# ---- print.Rcpp_SimplexTree ----
# nocov start
setClass("Rcpp_SimplexTree")
.print_simplex_tree <- setMethod("show", "Rcpp_SimplexTree", function (object) { 
  max_k <- length(object$n_simplices)
  if (max_k == 0){ cat("< empty simplex tree >\n") }
  else {
    cat(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
  }
}) 
# nocov end

# ---- print.Rcpp_Filtration ----
# nocov start
setClass("Rcpp_Filtration")
.print_filtration <- setMethod("show", "Rcpp_Filtration", function (object) {
  # cat(format(object))
  max_k <- length(object$n_simplices)
  if (max_k == 0){ cat("< empty filtration >\n") }
  else {
    writeLines(c(
      sprintf("Simplex Tree with (%s) (%s)-simplices", paste0(object$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")),
      sprintf("Current filtration index: %d", object$current_index)
    ))
  }
})
# nocov end

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
#' @param st a simplex tree object. 
#' @description Removes all simplices from the simplex tree, except the root node.
#' @examples 
#' st <- simplex_tree()
#' st$insert(1:3)
#' print(st) ## Simplex Tree with (3, 3, 1) (0, 1, 2)-simplices
#' st$clear()
#' print(st) ## < empty simplex tree >
clear <- function(st){
  stopifnot(class(st) %in% .st_classes)
  st$clear()
  return(invisible(st))
}

# ---- generate_ids ----
#' @name generate_ids
#' @aliases id_policy
#' @title Generates vertex ids.
#' @param st a simplex tree. 
#' @param n the number of ids to generate. 
#' @description Generates vertex ids representing 0-simplices not in the tree.
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
#' st %>% insert(as.list(1:3))
#' print(st$vertices) ## 0 1 2
#' st %>% insert(as.list(st %>% generate_ids(2)))
#' st %>% print_simplices() 
#' # 0, 1, 2, 3, 4, 0 4
#' st %>% remove(4)
#' st %>% generate_ids(1) 
#' # 4
#' @export
generate_ids <- function(st, n){
  stopifnot(is.numeric(n) && length(n) == 1)
  return(st$generate_ids(as.integer(n)))
}

# ---- degree ----
#' @name degree
#' @title The vertex degree.
#' @param st a simplex tree. 
#' @param vertices the vertex ids to check the degree of. 
#' @description Returns the number of edges (degree) for each given vertex id. 
degree <- function(st, vertices){
  stopifnot(is.vector(vertices) && is.numeric(vertices))
  return(st$degree(vertices))
}

# ---- expand ----
#' @name expand
#' @title k-expansion.
#' @param st a simplex tree. 
#' @param k the target dimension of the expansion.
#' @description Performs a k-expansion on the 1-skeleton of the complex, adding k-simplices 
#' if all combinations of edges are included. Because this operation uses the edges alone to infer 
#' the existence of higher order simplices, the expansion assumes the underlying complex
#' is a flag complex. 
#' @export
expand <- function(st, k=2){
  stopifnot(is.numeric(k))
  st$expand(k)
  return(invisible(st))
}


# ---- adjacent ----
#' @name adjacent
#' @title Adjacent vertices.
#' @param st a simplex tree.
#' @param vertices vertex ids. 
#' @description Returns a vector of vertex ids that are immediately adjacent to a given vertex.
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


# ---- is_face ----
#' @name is_face
#' @title Is face 
#' @description Checks whether a simplex is a face of another simplex and is in the complex.
#' @param st a simplex tree.  
#' @param tau a simplex which may contain \code{sigma} as a coface. 
#' @param sigma a simplex which may contain \code{tau} as a face. 
#' @details A simplex \eqn{\tau} is a face of \eqn{\sigma} if \eqn{\tau \subset \sigma}. This function 
#' checks whether that is true. \code{tau} and \code{sigma} are sorted before comparison.
#' @seealso \href{https://en.cppreference.com/w/cpp/algorithm/includes}{std::includes}
#' @return boolean indicating whether \code{tau} is a face of \code{sigma}. 
#' @examples 
#' st <- simplex_tree()
#' st %>% insert(1:3)
#' st %>% is_face(2:3, 1:3)
#' st %>% is_face(1:3, 2:3)
#' @export
is_face <- function(st, tau, sigma){
  tau_exists <- find(st, tau)
  sigma_exists <- find(st, sigma)
  return(tau_exists && sigma_exists && all(tau %in% sigma))
}


# ---- collapse ----
#' @name collapse
#' @title Elementary collapse
#' @description Performs an elementary collapse. 
#' @param st a simplex tree.
#' @param pair list of simplices to collapse. 
#' @param w vertex to collapse to, if performing a vertex collapse. 
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
#' st <- simplextree::simplex_tree(1:3)
#' st %>% print_simplices()
#' # 1, 2, 3, 1 2, 1 3, 2 3, 1 2 3
#' st %>% collapse(list(1:2, 1:3))
#' # 1, 2, 3, 1 3, 2 3=
#' 
#' st %>% insert(list(1:3, 2:5))
#' st %>% print_simplices("column")
#' # 1 2 3 4 5 1 1 2 2 2 3 3 4 1 2 2 2 3 2
#' #           2 3 3 4 5 4 5 5 2 3 3 4 4 3
#' #                           3 4 5 5 5 4
#' #                                     5
#' 
#' st %>% collapse(list(2:4, 2:5))
#' st %>% print_simplices("column") 
#' # 1 2 3 4 5 1 1 2 2 2 3 3 4 1 2 2 3
#' #           2 3 3 4 5 4 5 5 2 3 4 4
#' #                           3 5 5 5
#' @export
collapse <- function(st, pair, w=NULL){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.list(pair))
  tau <- pair[[1]]
  sigma <- pair[[2]]
  if (missing(w)){
    stopifnot(all(tau %in% sigma))
    st$collapse(tau, sigma)
  } else {
    stopifnot(all(sapply(list(tau,sigma,w), length) == 1))
    st$vertex_collapse(tau, sigma, w)
  }
  return(invisible(st))
}


#' threshold
#' @description Thresholds a given filtered simplicial complex.
#' @param st simplex tree.
#' @param index integer index to threshold to.
#' @param value numeric index to threshold filtration. 
#' @export
threshold <- function(st, index = NULL, value = NULL){
  stopifnot(class(st) %in% .st_classes)
  if ("Rcpp_Filtration" %in% class(st)){
    if (missing(index) && !missing(value)){
      st$threshold_value(value)
      return(invisible(st))
    } else if (!missing(index) && missing(value)){
      st$threshold_index(index)
      return(invisible(st))
    } else {
      stop("Must supply either an integer index or numeric value.")
    }
  } else { return(invisible(st)) }
}

# ---- contract ----
#' @name contract
#' @title Edge contraction
#' @description Performs an edge contraction. 
#' @param st a simplex tree.
#' @param edge an edge to contract, as a 2-length vector. 
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
#' st <- simplex_tree(1:3) 
#' st %>% print_simplices()
#' # 1, 2, 3, 1 2, 1 3, 2 3, 1 2 3
#' st %>% contract(c(1, 3)) %>% print_simplices()
#' # 1, 2, 1 2
#' @export
contract <- function(st, edge){
  stopifnot(class(st) %in% .st_classes)
  stopifnot(is.numeric(edge) && length(edge) == 2)
  st$contract(edge)
  return(invisible(st))
}


# ---- serialize / deserialize ----
#' @name serialize
#' @title Serializes the simplex tree. 
#' @description Provides a compressed serialization interface for the simplex tree.
#' @param st a simplex tree.
#' @family serialization
#' @details The serialize/deserialize commands can be used to compress/uncompress the complex into 
#' smaller form amenable for e.g. storing on disk (see \code{saveRDS}) or saving for later use. 
#' The serialization.
#' @examples 
#' st <- simplex_tree(list(1:5, 7:9))
#' st2 <- deserialize(serialize(st))
#' all.equal(as.list(preorder(st)), as.list(preorder(st2)))
#' # TRUE 
#' 
#' set.seed(1234)
#' R <- rips(dist(replicate(2, rnorm(150))), eps = pnorm(0.20), dim = 2)
#' print(R$n_simplices)
#' # 50 137 229
#' 
#' ## Approx. size of the full complex 
#' print(utils::object.size(as.list(preorder(R))), units = "Kb")
#' # 345.7 Kb
#' 
#' ## Approx. size of serialized version 
#' print(utils::object.size(serialize(R)), units = "Kb")
#' # 14.4 Kb
#' 
#' ## You can save these to disk via e.g. saveRDS(serialize(R), ...)
#' @export
serialize <- function(st){
  stopifnot(class(st) %in% .st_classes)
  n <- st$n_simplices[1]
  complex <- local({
    ids <- st$vertices
    minimal <- straverse(maximal(st), function(simplex){ 
      c(length(simplex), sub_to_nat(match(simplex, ids), n)) 
    })
    minimal <- minimal[,order(minimal[1,])]
    ids <- structure(rle(diff(ids)), head=ids[1])
    list(ids = ids, dims = rle(minimal[1,]), maps = minimal[2,])
  })
  return(complex)
}

#' @name deserialize 
#' @title Deserializes the simplex tree. 
#' @description Provides a compressed serialization interface for the simplex tree.
#' @param complex The result of \code{\link{serialize}}.
#' @param st optionally, the simplex tree to insert into. Otherwise a new one is created. 
#' @family serialization
#' @details The serialize/deserialize commands can be used to compress/uncompress the complex into 
#' smaller form amenable for e.g. storing on disk (see \code{saveRDS}) or saving for later use. 
#' @export
deserialize <- function(complex, st = NULL){
  if (is.null(complex)){ return(simplex_tree()) }
  stopifnot(all(c("ids", "dims", "maps") %in% names(complex)))
  if (missing(st) || is.null(st)){ st <- simplex_tree() } 
  else { stopifnot(st %in% .st_classes) }
  with(complex, {
    ids <- c(attr(ids, "head"), cumsum(inverse.rle(ids))+attr(ids, "head"))
    st %>% insert(as.list(ids))
    d <- inverse.rle(dims)
    n <- st$n_simplices[1]
    for (di in unique(d)){
      d_simplices <- nat_to_sub(maps[d == di], n, k = di)
      st %>% insert(matrix(ids[d_simplices], ncol = ncol(d_simplices), nrow = nrow(d_simplices), byrow = FALSE))
    }
  })
  return(st)
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

# ---- reindex ----
#' @name reindex 
#' @title reindexes vertex ids
#' @param st a simplex tree. 
#' @param ids vector of new vertex ids. See details. 
#' @description This function allows one to 'reorder' or 'reindex' vertex ids.  
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


# ---- is_tree ----
#' @name is_tree 
#' @title Checks if the simplicial complex is a tree.
#' @param st a simplex tree. 
#' @description This function performs a breadth-first search on the simplicial complex, checking if the complex is acyclic.
#' @examples 
#' st <- simplex_tree()
#' st %>% insert(list(1:2, 2:3))
#' st %>% is_tree() # true
#' st %>% insert(c(1, 3))
#' st %>% is_tree() # false
#' @export
is_tree <- function(st){
  stopifnot(class(st) %in% .st_classes)
  return(st$is_tree())
}

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
## @param clip_polygons Whether to clip the polygons. Useful when visualizing large complexes. See details. 
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
#' st %>% insert(list(1:4))
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
#' st %>% insert(list(1:3, 2:5, 5:9, 7:8, 10))
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
#' si_sum <- straverse(st %>% level_order, sum) 
#' rbw_pal <- rev(rainbow(50, start=0,end=4/6))
#' plot(st, color_pal=rbw_pal[cut(si_sum, breaks=50, labels = FALSE)])
#' 
#' ## This also makes highlighting simplicial operations fairly trivial 
#' four_cofaces <- as.list(cofaces(st, 4))
#' coface_pal <- straverse(level_order(st), function(simplex){ 
#'     ifelse(list(simplex) %in% four_cofaces, "orange", "blue") 
#' })
#' plot(st, color_pal=unlist(coface_pal))
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
plot.Rcpp_SimplexTree <- function(x, coords = NULL, vertex_opt=NULL, text_opt=NULL, edge_opt=NULL, polygon_opt=NULL, color_pal=NULL, maximal=TRUE, by_dim=TRUE, add=FALSE,...) { # nocov start
  stopifnot(class(x) %in% .st_classes)
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
  dim_idx <- straverse(level_order(x), length)-1L
  is_char_vec <- all(is.character(color_pal))
  is_in <- function(lst){ 
    return(function(element) { any(sapply(lst, function(x) (all.equal(x, element) == TRUE))) })
  }
  
  ## If the maximal faces are requested, set non-maximal `draw_simplex` indices to FALSE 
  if (maximal){
    all_simplices <- as.list(level_order(x))
    max_idx <- match(as.list(maximal(x)), all_simplices)
    draw_simplex <- vector("logical", length=sum(x$n_simplices))
    draw_simplex[max_idx] <- TRUE
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
    draw_simplex <- straverse(level_order(x), si_in)
    simplex_colors <- straverse(level_order(x), si_color)
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
    params <- list(...)
    rel_params <- intersect(c("xlim", "ylim", "log", "asp", "xaxs", "yaxs", "lab"), names(params))
    graphics::plot.new()
    if (length(rel_params) > 0){
      default_p <- list(xlim=range(coords[,1]), ylim=range(coords[,2]))
      do.call(graphics::plot.window, modifyList(default_p, params[rel_params]))
    } else {
      graphics::plot.window(xlim=range(coords[,1]), ylim=range(coords[,2])) 
    }
  }
  
  # plot polygons for simplices of dimension 2+; omits edges and vertices
  # this just plots the triangles
  v <- x$vertices # cache vertices
  if (x$dimension >= 2L){
    for (d in seq(x$dimension, 2)){
      if (any(draw_simplex[dim_idx == d])){
        
        # safe_to_clip <- is_char_vec && ((length(color_pal) == x$dimension+1L) || by_dim)
        # if (clip_polygons && safe_to_clip){
        #   polys <- ltraverse(k_simplices(x, k=d), function(simplex){ 
        #     poly <- coords[match(simplex, v),] 
        #     poly <- poly[grDevices::chull(poly),,drop=FALSE]
        #     list(x=poly[,1], y=poly[,2])
        #   })
        #   subset <- (draw_simplex & (dim_idx==d))
        #   polys <- polys[subset[dim_idx==d]]
        #   clipped_polys <- Reduce(f = function(A, B){ polyclip::polyclip(A, B, op = "union") }, 
        #                           x = polys[-1], init = polys[[1]])
        #   clipped_polys <- do.call(rbind, lapply(clipped_polys, function(p){ rbind(cbind(p$x, p$y), c(NA,NA))}))
        #   params <- list(x=clipped_polys, border=NA, col=simplex_colors[head(which(dim_idx==d),1)])
        #   do.call(graphics::polygon, modifyList(params, as.list(polygon_opt)))
        # } else {
          polys <- ltraverse(k_simplices(x, k=d), function(simplex){ 
            poly <- coords[match(simplex, v),] 
            rbind(poly[grDevices::chull(poly),,drop=FALSE], c(NA, NA))
          })
          subset <- (draw_simplex & (dim_idx==d))
          polys_to_draw <- polys[subset[dim_idx==d]]
          params <- list(x=do.call(rbind, polys_to_draw), border=NA, col=simplex_colors[dim_idx==d])
          do.call(graphics::polygon, modifyList(params, as.list(polygon_opt)))
        # }
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
} # nocov end

#' plot.Rcpp_Filtration
#' @param ... passed to \code{\link{plot.Rcpp_SimplexTree}}
#' @describeIn plot_simplextree family of plotting methods. 
#' @export
plot.Rcpp_Filtration <- function(...){
  plot.Rcpp_SimplexTree(...)
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
