
# ---- is_face ----

#' @md
#' @name is_face
#' @title Face test
#' @description Checks whether a simplex is a face of another simplex and is in the complex.
#' @param st a simplex tree.  
#' @param tau a simplex which may contain `sigma` as a coface. 
#' @param sigma a simplex which may contain `tau` as a face. 
#' @details A simplex \eqn{\tau} is a face of \eqn{\sigma} if the vertices of \eqn{\tau} are vertices of \eqn{\sigma}. This function 
#' checks whether that is true. `tau` and `sigma` are sorted before comparison.
#' @seealso \href{https://en.cppreference.com/w/cpp/algorithm/includes}{std::includes}
#' @return boolean indicating whether `tau` is a face of `sigma`. 
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

# ---- is_tree ----

#' @md
#' @name is_tree 
#' @title Tree (acyclicity) test
#' @description This function performs a breadth-first search on the simplicial complex, checking if the complex is acyclic.
#' @param st a simplex tree. 
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
