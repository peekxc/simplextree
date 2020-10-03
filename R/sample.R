#' @name sample
#' @title Sample random simplicial complexes
#' @description Generate random simplicial complexes following the model of
#'   Costa and Farber (2016).
#' @param n an integer number of starting vertices.
#' @param prob a numeric vector of simplex insertion probabilities. The
#'   dimension of the resulting simplicial complex will be at most
#'   \code{length(prob) - 1L}.

#' @details The random simplicial complex model of Costa and Farber (2016)
#'   begins with a finite number of vertices \eqn{n} (\code{n}) and proceeds as
#'   follows, based on the \eqn{d+1}-dimensional vector of probabilities
#'   \eqn{p_0,\ldots,p_d} (\code{prob}):

#'   \itemize{
#'     \item{Delete each vertex
#'           with probability \eqn{1-p_0}.}
#'     \item{Insert an edge
#'           on each pair of vertices
#'           with probability \eqn{p_1}.}
#'     \item{Insert a \eqn{2}-simplex
#'           on each triangle
#'           with probability \eqn{p_2}.}
#'     \item{For \eqn{k=3,\ldots,d}, insert a \eqn{k}-simplex
#'           on each subcomplex that forms a \eqn{k}-simplex boundary
#'           with probability \eqn{p_k}.}
#'   }

#' @references Costa A., Farber M. (2016) Random Simplicial Complexes. In:
#'   Callegaro F., Cohen F., De Concini C., Feichtner E., Gaiffi G., Salvetti M.
#'   (eds) Configuration Spaces. Springer INdAM Series, vol 14. Springer, Cham.
#'   https://doi.org/10.1007/978-3-319-31580-5_6
#' @examples
#' set.seed(1)
#' ## Generate Costa-Farber random simplicial complexes
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' ## Construct a complete complex of a given size and dimension
#' sample_costa_farber(n = 6L, prob = rep(1, 5L))
#' ## Construct the clique complex of a random 1-skeleton
#' plot(sample_costa_farber(n = 10L, prob = c(.7, .4, rep(1, 11L))))
#' ## Generate Linial-Meshulam random simplicial complexes
#' sample_costa_farber(n = 6L, prob = c(rep(1, 1L), .5))
#' sample_costa_farber(n = 6L, prob = c(rep(1, 2L), .5))
#' sample_costa_farber(n = 6L, prob = c(rep(1, 3L), .5))

#' @rdname sample
#' @export
sample_costa_farber <- function(n, prob) {
  stopifnot(
    n >= 0L,
    inherits(prob, "numeric")
  )
  
  ## Create an empty simplicial complex
  st <- simplex_tree()
  
  ## Retain vertices independently with probability p_0
  vs <- which(as.logical(stats::rbinom(n, 1L, prob[[1L]])))
  st %>% insert(as.list(vs))
  ## Return the complex if done
  if (length(prob) == 1L) return(st)
  
  ## Insert edges independently with probability p_1
  m1 <- choose(st$n_simplices[[1L]], 2L)
  n1 <- stats::rbinom(n = 1L, size = m1, prob = prob[[2L]])
  ex <- sort(sample.int(n = m1, size = n1))
  es <- nat_to_sub(ex, n = st$n_simplices[[1L]], k = 2L)
  es[] <- st$vertices[es]
  st$insert_lex(es)
  ## Return the complex if done
  if (length(prob) == 2L) return(st)
  
  ## Iteravely sample higher-dimensional simplices
  for (k in seq(3L, length(prob))) {
    expand_f_bernoulli(st$as_XPtr(), k = k - 1L, p = prob[[k]])
    if (st$dimension < k - 1L) break
  }
  ## Return the complex
  return(st)
}
