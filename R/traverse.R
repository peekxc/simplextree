
.traversal_types = c("Preorder", "Level order", "Face", "Coface", "Coface roots", "K-skeleton", "K-simplices", "Maximal simplex", "Link")

#' @name traversals
#' @title Methods for traversal objects
NULL

# ---- print.st_traversal ----

#' @rdname traversals
#' @param x traversal object.
#' @param ... unused. 
#' @export
print.st_traversal <- function(x, ...){
  sigma_str <- ifelse(length(x$sigma) == 0 || is.null(x$sigma), "empty face", paste0(x$sigma, collapse = " "))
  tt <- .traversal_types[x$traversal_type+1L]
  writeLines(sprintf("%s traversal @ { %s }", tt, sigma_str))
}

# ---- as.list.st_traversal ----

#' @rdname traversals
#' @param x traversal object.
#' @param ... unused. 
#' @export
as.list.st_traversal <- function(x, ...){
  return(ltraverse(x, identity))
}


# ---- traverse ----

#' @name traverse
#' @title Apply a function along a traversal
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

#' @rdname traverse
#' @export
straverse <- function(traversal, f, ...){
  stopifnot("st_traversal" %in% class(traversal))
  stopifnot(is.function(f))
  # if (missing(f)){ return(function(traversal, f){ straverse_R(traversal, f) }) }
  # if (is.function(traversal)){ traversal(f) }
  return(straverse_R(traversal, f))
}

#' @rdname traverse
#' @export
ltraverse <- function(traversal, f, ...){
  stopifnot("st_traversal" %in% class(traversal))
  stopifnot(is.function(f))
  # if (missing(f)){ return(function(traversal, f){ ltraverse_R(traversal, f) }) }
  # if (is.function(traversal)){ traversal(f) }
  return(ltraverse_R(traversal, f))
}

# ---- empty_face ----

#' @name empty_face
#' @title Empty faces
#' @description Alias to the empty integer vector (integer(0L)). Used to indicate the empty face of the tree. 
#' @seealso traverse
#' @export
empty_face <- integer(0L)

# ---- preorder ----- 

#' @name preorder 
#' @title Preorder traversal
#' @description Generate a preorder traversal on the simplex tree. 
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
#' @title Level order traversal
#' @description Generates a level order traversal on the simplex tree. 
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
#' @title Face traversal
#' @description Generates a face traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
faces <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "faces", NULL)
}

# ---- cofaces ----- 

#' @name cofaces 
#' @title Coface traversal
#' @description Generates a coface traversal on the simplex tree. 
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at. 
#' @export
cofaces <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "cofaces", NULL)
}

# ---- k_skeleton ----- 

#' @name k_skeleton 
#' @title k-Skeleton traversal
#' @description Generates a k-skeleton traversal on the simplex tree.
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
#' @title Coface roots traversal
#' @description Generates a coface roots traversal on the simplex tree.
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at.
#' @description The coface roots of a given simplex \code{sigma} are the roots
#'   of subtrees in the trie whose descendents (including the roots themselves)
#'   are cofaces of \code{sigma}. This traversal is more useful when used in
#'   conjunction with other traversals, e.g. a \emph{preorder} or
#'   \emph{level_order} traversal at the roots enumerates the cofaces of
#'   \code{sigma}.
#' @export
coface_roots <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  parameterize_R(st$as_XPtr(), sigma, "coface_roots", NULL)
}

# ---- maximal ----- 

#' @name maximal 
#' @title Maximal traversal
#' @description Generates a traversal on the maximal of the simplex tree.
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
#' @title k-Simplex traversal
#' @description Generates a traversal on the k-simplices of the simplex tree.
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
#' @title Link traversal
#' @description Generates a traversal on the link of a given simplex in the
#'   simplex tree.
#' @param st the simplex tree to traverse.
#' @param sigma simplex to start the traversal at.
#' @export
link <- function(st, sigma){
  stopifnot(class(st) %in% .st_classes)
  if (is.null(sigma)){ sigma <- empty_face }
  parameterize_R(st$as_XPtr(), sigma, "link", NULL)
}
