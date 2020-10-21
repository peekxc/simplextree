
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

#' @name print_simplices
#' @title Print simplices to the console
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
