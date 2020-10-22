#' @name combinadic
#' @title Lexicographic mappings of \emph{k}-combinations
#' @description Map between \emph{k}-combinations of \code{n} and their
#'   lexicographic positions, and recover binomial coefficient numerators.
#' @details \code{nat_to_sub} computes the \code{x}^th (\code{n} choose
#'   \code{k}) combination in the lexicographic order. \code{sub_to_nat}
#'   computes the position of a combination \code{x} out of all
#'   lexicographically-ordered (\code{n} choose \code{k}) combinations. The
#'   values are calculated via a lexicographically-ordered combinadic mapping.
#'
#'   In general, \code{nat_to_sub} is \emph{not} intended to be used to
#'   \emph{generate} all 2-combinations in the combinadic mapping.
#'
#'   \code{inverse.choose} inverts the binomial coefficient for general
#'   (\emph{n}, \emph{k}). That is, given the denominator \code{k} and \code{x =
#'   choose(n, k)}, find \code{n}.
#' @param n integer numerator of the binomial coefficient.
#' @param k integer denominator of the binomial coefficient.
#' @param i vector of integers in the range \code{c(1L, choose(n, k))}.
#' @param s matrix whose columns represent \code{k}-combinations.
#' @param x a binomial coefficient (\emph{n} choose \emph{k}).
#' @return an integer matrix whose columns give the combinadics of \code{x}
#'   (\code{nat_to_sub}), an integer vector of the positions of the given
#'   combinations (\code{sub_to_nat}), or the integer numerator of the binomial
#'   coefficient (\code{inverse.choose}).
#' @references McCaffrey, J. D. "Generating the mth lexicographical element of a
#'   mathematical combination." MSDN Library (2004).
#' @examples
#' library(simplextree)
#' all(nat_to_sub(seq(choose(100,2)), n = 100, k = 2) == combn(100,2))
#'
#' ## Generating pairwise combinadics is particularly fast
#' ## Below: test to generate ~ 45k combinadics
#' ## (note: better to use microbenchmark)
#' system.time({
#'   x <- seq(choose(300,2))
#'   nat_to_sub(x, n = 300, k = 2L)
#' })
#'
#' ## Compare with generating raw combinations
#' system.time(combn(300,2))
#'
#' 100 == inverse.choose(choose(100,2), k = 2)
#' # TRUE
#' 12345 == inverse.choose(choose(12345, 5), k = 5)
#' # TRUE
NULL

#' @rdname combinadic
#' @export  
nat_to_sub  <- function(i, n, k){
	if (length(i) == 0){ return(matrix(integer(0), nrow = k)) }
	stopifnot(all(i >= 1 & i <= choose(n,k)))
	if (choose(n, k) > ((2^64) - 2)){ stop("(n,k) combination too big; combinadics limited to 64-bit integer arithmetic.") }
	return(to_subscript_R(as.integer(i)-1L, n, k)+1L)
}

#' @rdname combinadic
#' @export
sub_to_nat  <- function(s, n){
	if (length(s) == 0){ return(integer(0)) }
	if (is.null(dim(s))){ s <- matrix(s, ncol = 1L) }
  if (choose(n, nrow(s)) > ((2^64) - 2)){ stop("(n,k) combination too big; combinadics limited to 64-bit integer arithmetic.") }
	return(to_natural_R(s-1L, n)+1L)
}

#' @rdname combinadic
#' @export 
inverse.choose <- function(x, k) {
  stopifnot(k >= 1)
  if (k == 1){ final_n <- x }
  else if (k == 2){
    rng <- floor(sqrt(2*x)):ceiling(sqrt(2*x)+2)
	  final_n <- rng[choose(rng, 2) == x]
  } else {
    # From: https://math.stackexchange.com/a/2381576/214799
    if (x < 10^7){
      lb <- (factorial(k)*x)^(1/k)
      potential_n <- seq(floor(lb), ceiling(lb+k))
      final_n <- potential_n[choose(potential_n, k) == x]
    } else {
      lb <- floor((4^k)/(2*k + 1))
      C <- factorial(k)*x
      n <- 1
      while(n^k < C){ n <- n*2 }
      m <- utils::head(which(seq(1,n)^k >= C),1)
      potential_n <- seq(max(c(m, 2*k)), m+k)
      final_n <- potential_n[choose(potential_n, k) == x]
    }
  }
  if (length(final_n) == 0){ stop("Unable to invert binomial coefficient.") }
  return(final_n)
}
