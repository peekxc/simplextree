#' @title nat_to_sub
#' @description Computes the \code{x}^th (\code{n} choose \code{2}) combination.
#' @details The mapping is done via an lexicographically-ordered combinadic mapping. \cr
#' In general, this function is \emph{not} intended to be used to \emph{generate} all (n choose \code{k})  combinations in the combinadic mapping. 
#' @param x non-negative integers in the range \code{c(1, choose(n, k))}
#' @param n numerator of the binomial coefficient 
#' @param k denominator of the binomial coefficient
#' @return integer matrix whose columns give the combinadics of \code{x}.
#' @references McCaffrey, J. D. "Generating the mth lexicographical element of a mathematical combination." MSDN Library (2004).
#' @examples 
#' library(simplextree)
#' all(nat_to_sub(seq(choose(100,2)), n = 100, k = 2) == combn(100,2))
#' @export  
nat_to_sub  <- function(x, n, k){
	if (length(x) == 0){ return(matrix(integer(0), ncol = k)) }
	stopifnot(all(x >= 1 && x <= choose(n,k)))
	return(to_subscript_R(as.integer(x)-1L, n, k)+1L)
}

#' @title sub_to_nat
#' @description Given a combination \code{x}, computes its position out of all lexicographically-ordered (\code{n} choose \code{2}) combinations.
#' @details The mapping is done via an lexicographically-ordered combinadic mapping. 
#' @param x matrix whose columns represent \code{k}-combinations.
#' @param n numerator of the binomial coefficient 
#' @references McCaffrey, J. D. "Generating the mth lexicographical element of a mathematical combination." MSDN Library (2004).
#' @return integer vector of the positions of the given combinations. 
#' @export
sub_to_nat  <- function(x, n){
	if (length(x) == 0){ return(integer(0)) }
	if (is.null(dim(x))){ x <- matrix(x, ncol = 1L) }
	return(to_natural_R(x-1L, n)+1L)
}