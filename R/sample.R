#' @name sample
#' @title Sample random simplicial complexes
#' @description Generate random simplicial complexes following the models of
#'   Meshulam and Wallach (2009), Kahle (2009), and Costa and Farber (2016).
#' @param n an integer number of starting vertices.
#' @param prob a numeric simplex insertion probability (Linial-Meshulam-Wallach,
#'   Kahle) or a vector of probabilities for all dimensions (Costa-Farber). The
#'   dimension of a Costa-Farber random simplicial complex will be at most
#'   \code{length(prob) - 1L}.
#' @param dimension an integer dimension at which to randomly insert simplices.
#' @param method a character string indicating the model to use; matched to
#'   \code{"erdos_renyi"}, \code{"kahle"}, \code{"linial_meshulam_wallach"}, and
#'   \code{"costa_farber"}, allowing for spaces in place of underscores.

#' @details The random graph model \eqn{G(n,p)} of Erdős and Rényi (1959) powers
#'   parts of other models and is exported for convenience.
#'
#'   The random clique complex model of Kahle (2009) samples an Erdős-Rényi
#'   random graph, then uses \code{\link{expand}} to insert all complete
#'   subgraphs.
#'
#'   The random simplicial complex model of Costa and Farber (2016) begins with
#'   a finite number of vertices \eqn{n} (\code{n}) and proceeds as follows,
#'   based on the \eqn{d+1}-dimensional vector of probabilities
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
#'   and Meshulam (2006), is a special case in which \eqn{p_k=1} for \eqn{0\le
#'   k\le d-1}; the only parameters are \eqn{n} (\code{n}) and \eqn{p_d}
#'   (\code{prob}).

#' @references Erdős P. and Rényi A. (1959) On Random Graphs I. Publicationes
#'   Mathematicae 6: 290–297.
#' @references Linial N. and Meshulam R. (2006) Homological Connectivity of
#'   Random 2-Complexes. Combinatorica 26, 4, 475–487.
#'   doi:10.1007/s00493-006-0027-9
#' @references Meshulam, R. and Wallach, N. (2009) Homological Connectivity of
#'   Random k‐Dimensional Complexes. Random Struct. Alg., 34: 408–417.
#'   doi:10.1002/rsa.20238
#' @references Kahle, M. (2009) Topology of Random Clique Complexes. Discrete
#'   Math., 309(6): 1658–1671. doi:10.1016/j.disc.2008.02.037
#' @references Costa A. and Farber M. (2016) Random Simplicial Complexes. In:
#'   Callegaro F., Cohen F., De Concini C., Feichtner E., Gaiffi G., Salvetti M.
#'   (eds) Configuration Spaces. Springer INdAM Series, vol 14. Springer, Cham.
#'   doi:10.1007/978-3-319-31580-5_6
#' @examples
#' set.seed(1)
#' ## Generate Erdos-Renyi random graphs
#' plot(sample_abstract(n = 12L, prob = .2, method = "erdos_renyi"))
#' plot(sample_abstract(n = 12L, prob = .5, method = "erdos_renyi"))
#' plot(sample_abstract(n = 12L, prob = .8, method = "erdos_renyi"))
#' ## Generate Kahle random clique complexes
#' sample_abstract(n = 6L, prob = .2, method = "kahle")
#' sample_abstract(n = 6L, prob = .5, method = "kahle")
#' sample_abstract(n = 6L, prob = .8, method = "kahle")
#' ## Generate Linial-Meshulam random simplicial complexes
#' sample_abstract(n = 6L, dimension = 0L, prob = .6,
#'                 method = "linial_meshulam_wallach")
#' sample_abstract(n = 6L, dimension = 1L, prob = .6,
#'                 method = "linial_meshulam_wallach")
#' sample_abstract(n = 6L, dimension = 2L, prob = .6,
#'                 method = "linial_meshulam_wallach")
#' sample_abstract(n = 6L, dimension = 3L, prob = .6,
#'                 method = "linial_meshulam_wallach")
#' ## Generate Costa-Farber random simplicial complexes
#' plot(sample_abstract(n = 12L, prob = c(.5, .5, .5), method = "costa_farber"))
#' plot(sample_abstract(n = 12L, prob = c(.5, .5, .5), method = "costa_farber"))
#' plot(sample_abstract(n = 12L, prob = c(.5, .5, .5), method = "costa_farber"))
#' ## Construct a complete complex of a given size and dimension
#' sample_abstract(n = 6L, dimension = 4L, prob = 0,
#'                 method = "linial_meshulam_wallach")
#' sample_abstract(n = 6L, prob = rep(1, 4L), method = "costa_farber")
#' ## Construct the clique complex of a random 1-skeleton
#' plot(sample_abstract(n = 10L, prob = c(.7, .6, rep(1, 11L)),
#'                      method = "costa_farber"))

#' @rdname sample
#' @export
sample_abstract <- function(
  n, prob, dimension = NULL,
  method = c("costa_farber", "linial_meshulam_wallach", "kahle", "erdos_renyi")
) {
  
  ## Match (unique) method, allowing for spaces
  method <- match.arg(
    gsub(" ", "_", method),
    c("erdos_renyi", "kahle", "linial_meshulam_wallach", "costa_farber")
  )
  
  ## Check validity of parameters
  stopifnot(
    n >= 0L,
    inherits(prob, "numeric")
  )
  if (method != "costa_farber") stopifnot(length(prob) == 1L)
  if (method == "linial_meshulam_wallach") 
    stopifnot(dimension >= 0L) else
      if (! is.null(dimension))
        warning("`dimension` argument unused by `", method, "` method.")
  
  ## Transform parameters for the Costa-Farber model
  if (method == "erdos_renyi" || method == "kahle") {
    prob <- c(1, prob)
  } else if (method == "linial_meshulam_wallach") {
    prob <- c(rep(1, dimension), prob)
  }
  
  ## Execute Costa-Farber model on transformed parameters
  
  if (length(prob) == 0L) {
    ## Create an empty simplicial complex
    st <- simplex_tree()
    ## Return the complex if done
    return(st)
  }
  
  ## Retain vertices independently with probability p
  vs <- which(as.logical(stats::rbinom(n, 1L, prob[[1L]])))
  
  ## Create and return the complex if done
  if (length(prob) == 1L || length(vs) == 0L) {
    ## Create an empty simplicial complex
    st <- simplex_tree()
    insert(st, as.list(vs))
    return(st)
  }
  
  ## Create an empty simplicial complex
  st <- simplex_tree()
  
  ## Select edges independently with probability p
  m1 <- choose(length(vs), 2L)
  n1 <- stats::rbinom(n = 1L, size = m1, prob = prob[[2L]])
  ex <- sort(sample.int(n = m1, size = n1))
  es <- nat_to_sub(ex, n = length(vs), k = 2L)
  
  ## Populate simplicial complex with vertices and edges
  insert(st, as.list(seq(length(vs))))
  st$insert_lex(es)
  
  ## Reindex vertices
  reindex(st, vs)
  ## Return the complex if done
  if (length(prob) == 2L || st$dimension < 1L) return(st)
  
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
