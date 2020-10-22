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
#' 
#' ## Generating pairwise combinadics is particularly fast
#' ## Below: test to generate ~ 45k combinadics (note: better to use microbenchmark) 
#' system.time({
#'   x <- seq(choose(300,2))
#'   nat_to_sub(x, n = 300, k = 2L)
#' })
#' 
#' ## Compare with generating raw combinations
#' system.time(combn(300,2))
#' @export  
nat_to_sub  <- function(x, n, k){
	if (length(x) == 0){ return(matrix(integer(0), nrow = k)) }
	stopifnot(all(x >= 1 & x <= choose(n,k)))
	if (choose(n, k) > ((2^64) - 2)){ stop("(n,k) combination too big; combinadics limited to 64-bit integer arithmetic.") }
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
  if (choose(n, nrow(x)) > ((2^64) - 2)){ stop("(n,k) combination too big; combinadics limited to 64-bit integer arithmetic.") }
	return(to_natural_R(x-1L, n)+1L)
}

#' inverse.choose
#' @description Inverts the binomial coefficient for general (n,k).
#' @details Given a quantity x = choose(n, k) with fixed k, finds n. 
#' @param x the binomial coefficient. 
#' @param k the denominator of the binomial coefficient \code{x}. 
#' @examples 
#' 100 == inverse.choose(choose(100,2), k = 2)
#' # TRUE 
#' 12345 == inverse.choose(choose(12345, 5), k = 5)
#' # TRUE
#' @return the numerator of the binomial coefficient, if the Otherwise  
#' @export 
inverse.choose <- function(x, k) {
  stopifnot(k >= 1)
  if (k == 1){ final_n <- x }
  else if (k == 2){
    rng <- floor(sqrt(2*x)):ceiling(sqrt(2*x)+2)
	  final_n <- rng[choose(rng, 2) == x]
  } else {
    # From: https://math.stackexchange.com/questions/103377/how-to-reverse-the-n-choose-k-formula
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



