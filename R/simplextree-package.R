#' @name simplextree-package
#' @aliases simplextree
#' @title 'simplextree' package
#' @description Provides an R/Rcpp implementation of a Simplex Tree data structure and its related tools. 
#' @author Matt Piekenbrock
#' @useDynLib simplextree, .registration = TRUE
#' @import methods Rcpp
#' @details This package provides a lightweight implementation of a Simplex Tree data structure, exported as an Rcpp Module.
#' The current implementation provides a limited API and a subset of the functionality described in the paper.
#' @docType package
NULL

# @importFrom Rcpp evalCpp Module cpp_object_initializer

## Simplex Tree module (can be loaded anywhere)
Rcpp::loadModule("simplex_tree_module", TRUE)
Rcpp::loadModule("filtration_module", TRUE)
Rcpp::loadModule("union_find_module", TRUE) 

.st_classes <- c("Rcpp_SimplexTree", "Rcpp_Filtration")