#' @name sample
#' @title Sample random simplicial complexes
#' @description Generate random simplicial complexes following the models of
#'   Meshulam and Wallach (2009) and of Costa and Farber (2016).
#' @param n an integer number of starting vertices.
#' @param dimension an integer dimension at which to randomly insert simplices.
#' @param prob a numeric simplex insertion probability (Linial-Meshulam-Wallach)
#'   or a vector of probabilities for all dimensions (Costa-Farber). The
#'   dimension of a Costa-Farber random simplicial complex will be at most
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

#'   The model of Meshulam and Wallach (2009), generalized from that of Linial
#'   and Meshulam (2006) is a special case in which \eqn{p_k=1} for \eqn{0\le
#'   k\le d-1}; the only parameters are \eqn{n} (\code{n}) and \eqn{p_d}
#'   (\code{prob}).

#' @references Linial N. and Meshulam R. (2006) Homological Connectivity of
#'   Random 2-Complexes. Combinatorica 26, 4, 475–487.
#'   DOI:https://doi.org/10.1007/s00493-006-0027-9
#' @references Meshulam, R. and Wallach, N. (2009) Homological Connectivity of
#'   Random k‐Dimensional Complexes. Random Struct. Alg., 34: 408-417.
#'   doi:10.1002/rsa.20238
#' @references Costa A. and Farber M. (2016) Random Simplicial Complexes. In:
#'   Callegaro F., Cohen F., De Concini C., Feichtner E., Gaiffi G., Salvetti M.
#'   (eds) Configuration Spaces. Springer INdAM Series, vol 14. Springer, Cham.
#'   https://doi.org/10.1007/978-3-319-31580-5_6
#' @examples
#' set.seed(1)
#' ## Generate Linial-Meshulam random simplicial complexes
#' sample_linial_meshulam_wallach(n = 6L, dimension = 0L, prob = .5)
#' sample_linial_meshulam_wallach(n = 6L, dimension = 1L, prob = .5)
#' sample_linial_meshulam_wallach(n = 6L, dimension = 2L, prob = .5)
#' sample_linial_meshulam_wallach(n = 6L, dimension = 3L, prob = .5)
#' ## Generate Costa-Farber random simplicial complexes
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' plot(sample_costa_farber(n = 12L, prob = c(.5, .5, .5)))
#' ## Construct a complete complex of a given size and dimension
#' sample_linial_meshulam_wallach(n = 6L, dimension = 4L, prob = 0)
#' sample_costa_farber(n = 6L, prob = rep(1, 4L))
#' ## Construct the clique complex of a random 1-skeleton
#' plot(sample_costa_farber(n = 10L, prob = c(.7, .4, rep(1, 11L))))

#' @rdname sample
#' @export
sample_linial_meshulam_wallach <- function(n, dimension, prob) {
  stopifnot(
    n >= 0L,
    dimension >= 0L,
    inherits(prob, "numeric"),
    length(prob) == 1L
  )
  
  ## Create an empty simplicial complex
  st <- simplex_tree()
  
  if (dimension == 0L) {
    ## Retain vertices independently with probability p
    vs <- which(as.logical(stats::rbinom(n, 1L, prob)))
    insert(st, as.list(vs))
  } else if (dimension == 1L) {
    ## Insert edges independently with probability p
    m1 <- choose(n, 2L)
    n1 <- stats::rbinom(n = 1L, size = m1, prob = prob)
    ex <- sort(sample.int(n = m1, size = n1))
    es <- nat_to_sub(ex, n = n, k = 2L)
    st$insert_lex(es)
  } else {
    ## Make `st` (d-1)-complete
    insert(st, simplices = utils::combn(n, dimension))
    ## Insert d-simplices with probability p 
    expand_f_bernoulli(st$as_XPtr(), k = dimension, p = prob)
  }
  
  ## Return the complex
  return(st)
}

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
  for (d in seq(2L, length(prob) - 1L)) {
    expand_f_bernoulli(st$as_XPtr(), k = d, p = prob[[d + 1L]])
    if (st$dimension < d) break
  }
  
  ## Return the complex
  return(st)
}

sample_complex <- function(n, prob, dim=NULL, method = c("costa farber", "linial meshulam wallach")){
  
}
