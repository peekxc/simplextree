#' @name sample-geometric
#' @title Sample random geometric simplicial complexes
#' @description Generate Vietoris–Rips complexes on random point clouds.
#' @param d a distance matrix (\code{"dist"} class) object, or a numeric matrix
#'   of (row) coordinates of points (which will be transformed into a distance
#'   matrix).
#' @param radius a numeric distance within which subsets of points will form
#'   simplices.
#' @param dimension an integer maximum dimension of simplices to form.
#' @param method a character string indicating the model to use; matched only to
#'   \code{"vietoris_rips"}, allowing for spaces in place of underscores,
#'   anticipating future additional methods like \code{"cech"}.
#' @param coords a logical instruction to retain the coordinates from a numeric
#'   matrix \code{d} as an attribute of the simplicial complex.
#' @param n an integer number of starting points.
#' @param torus a logical instruction to identify opposite faces of the sampling
#'   region.
#' @param ... additional parameters passed to \code{\link{rips}()}.

#' @details The geometric random graph model (see Penrose, 2003) begins with a
#'   random sample of points from a distribution on a manifold (usually
#'   Euclidean space), which are taken to be vertices, and introduces edges
#'   between vertices within a fixed distance of each other.
#'
#'   The geometric random simplicial complex model extends this model by
#'   constructing a Vietoris–Rips or Čech complex on the sample. See Kahle
#'   (2011) and Bobrowski and Weinberger (2017) for key results and Kahle (2017)
#'   for a review.

#' @references Penrose M. (2003) Random Geometric Graphs. Oxford University
#'   Press. doi:10.1093/acprof:oso/9780198506263.001.0001/acprof-9780198506263
#' @references Kahle M. (2011) Random Geometric Complexes. Discrete Comput.
#'   Geom. 45, 553–573. doi:10.1007/s00454-010-9319-3
#' @references Bobrowski O. and Weinberger S. (2017) On the vanishing of
#'   homology in random Čech complexes. Random Struct. Alg., 51: 14–51.
#'   doi:10.1002/rsa.20697
#' @references Kahle M. (2017) In: J.E. Goodman, J. O'Rourke, and C.D. Tóth
#'   (eds) Handbook of Discrete and Computational Geometry, 3rd edition. CRC
#'   Press, Boca Raton, FL.
#' @examples
#' set.seed(1)
#' ## Construct geometric simplicial complexes from a sample point cloud
#' theta <- stats::runif(n = 24L, min = 0, max = 2*pi)
#' x <- cbind(x = cos(theta), y = sin(theta))
#' plot(x)
#' make_geometric(x, radius = .03, dimension = 2L)
#' make_geometric(x, radius = .3, dimension = 2L)
#' ## Check distance ranges for square and toroidal samples
#' sqrt(2)
#' range(sample_unit(n = 1e3L))
#' sqrt(2)/2
#' range(sample_unit(n = 1e3L, torus = TRUE))
#' ## Construct random geometric simplicial complexes, square and toroidal
#' plot(sample_geometric(24L, radius = .1, dimension = 1L))
#' plot(sample_geometric(24L, radius = .1, dimension = 1L, torus = TRUE))
#' plot(sample_geometric(24L, radius = .1, dimension = 2L))
#' plot(sample_geometric(24L, radius = .1, dimension = 2L, torus = TRUE))
#' @return A \code{\link{simplex_tree}} (\code{*_geometric()}) or a
#'   \code{"dist"} object or coordinate matrix (\code{sample_unit()}).

#' @rdname sample-geometric
#' @export
make_geometric <- function(
  d, radius = NULL, dimension = 1L,
  method = c("vietoris_rips"), coords = FALSE, ...
) {
  
  ## Check that `coords` are available
  if (coords && ! "dist" %in% class(d)) {
    warning("Coordinates cannot be recovered from a 'dist' object.")
  } else {
    x <- as.matrix(d)
  }
  ## Ensure that `d` is a 'dist' object
  if (! "dist" %in% class(d)) d <- stats::dist(d)
  ## Take a `NULL` radius to be the enclosing radius of `d`
  if (is.null(radius)) radius <- enclosing_radius(d)
  ## Check validity of parameters
  stopifnot(
    inherits(radius, "numeric"),
    radius >= 0,
    dimension > 0L
  )
  ## Match (unique) method, allowing for spaces
  method <- match.arg(
    gsub(" ", "_", method),
    c("vietoris_rips", "rips_vietoris")
  )
  method <- switch(
    method,
    vietoris_rips = "vietoris_rips",
    rips_vietoris = "vietoris_rips"
  )
  
  ## Construct the geometric simplicial complex
  st <- rips(d, eps = 2 * radius, dim = dimension, ...)
  ## Attribute coordinates if requested
  if (coords) attr(st, "coords") <- x
  ## Return geometric complex
  st
}

#' @rdname sample-geometric
#' @export
sample_unit <- function(n, torus = FALSE, coords = FALSE) {
  
  ## Sample points from the unit square
  x <- cbind(x = stats::runif(n = n), y = stats::runif(n = n))
  
  ## Return coordinates if feasible and desired
  if (coords) {
    if (torus) {
      warning("Coordinates cannot be returned for a toroidal sample.")
    } else {
      return(x)
    }
  }
  
  ## Calculate Euclidean distances
  d <- stats::dist(x)
  
  ## Calculate additional flat toroidal distances and take pairwise minima
  if (torus) {
    ## shift only horizontal coordinate
    x[, 1L] <- (x[, 1L] + .5) %% 1
    d <- pmin(d, stats::dist(x))
    ## shift also vertical coordinate
    x[, 2L] <- (x[, 2L] + .5) %% 1
    d <- pmin(d, stats::dist(x))
    ## restore horizontal coordinate
    x[, 1L] <- (x[, 1L] - .5) %% 1
    d <- pmin(d, stats::dist(x))
  }
  
  ## Return distance matrix
  d
}

#' @rdname sample-geometric
#' @export
sample_geometric <- function(
  n, torus = FALSE, radius = NULL, dimension = 1L,
  method = c("vietoris_rips"), coords = FALSE, ...
) {
  
  ## Generate sample
  d <- sample_unit(n = n, torus = torus, coords = coords)
  
  ## Construct geometric simplicial complex
  make_geometric(
    d, radius = radius, dimension = dimension,
    method = method, coords = coords, ...
  )
}
