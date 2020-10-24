#' @name serialize
#' @title Serialize the simplex tree
#' @description Provides a compressed serialization interface for the simplex tree.
#' @param st a simplex tree.
#' @family serialization methods
#' @details The serialize/deserialize commands can be used to compress/uncompress the complex into 
#' smaller form amenable for e.g. storing on disk (see \code{saveRDS}) or saving for later use. 
#' The serialization.
#' @examples 
#' st <- simplex_tree(list(1:5, 7:9))
#' st2 <- deserialize(serialize(st))
#' all.equal(as.list(preorder(st)), as.list(preorder(st2)))
#' # TRUE 
#' 
#' set.seed(1234)
#' R <- rips(dist(replicate(2, rnorm(100))), eps = pnorm(0.10), dim = 2)
#' print(R$n_simplices)
#' # 100 384 851
#' 
#' ## Approx. size of the full complex 
#' print(utils::object.size(as.list(preorder(R))), units = "Kb")
#' # 106.4 Kb
#' 
#' ## Approx. size of serialized version 
#' print(utils::object.size(serialize(R)), units = "Kb")
#' # 5.4 Kb
#' ## You can save these to disk via e.g. saveRDS(serialize(R), ...)
#' @export
serialize <- function(st){
  stopifnot(class(st) %in% .st_classes)
  n <- st$n_simplices[1]
  complex <- local({
    ids <- st$vertices
    minimal <- matrix(straverse(maximal(st), function(simplex){ 
      c(length(simplex), sub_to_nat(match(simplex, ids), n)) 
    }), nrow = 2)
    minimal <- minimal[,order(minimal[1,]),drop=FALSE]
    ids <- structure(rle(diff(ids)), head=ids[1])
    list(ids = ids, dims = rle(minimal[1,]), maps = minimal[2,])
  })
  return(complex)
}

#' @name deserialize 
#' @title Deserialize the simplex tree
#' @description Provides a compressed serialization interface for the simplex tree.
#' @param complex The result of \code{\link{serialize}}.
#' @param st optionally, the simplex tree to insert into. Otherwise a new one is created. 
#' @family serialization methods
#' @details The serialize/deserialize commands can be used to compress/uncompress the complex into 
#' smaller form amenable for e.g. storing on disk (see \code{saveRDS}) or saving for later use. 
#' @export
deserialize <- function(complex, st = NULL){
  if (is.null(complex)){ return(simplex_tree()) }
  stopifnot(all(c("ids", "dims", "maps") %in% names(complex)))
  if (missing(st) || is.null(st)){ st <- simplex_tree() } 
  else { stopifnot(st %in% .st_classes) }
  with(complex, {
    ids <- c(attr(ids, "head"), cumsum(inverse.rle(ids))+attr(ids, "head"))
    st %>% insert(as.list(ids))
    d <- inverse.rle(dims)
    n <- st$n_simplices[1]
    for (di in unique(d)){
      d_simplices <- nat_to_sub(maps[d == di], n, k = di)
      st %>% insert(matrix(ids[d_simplices], ncol = ncol(d_simplices), nrow = nrow(d_simplices), byrow = FALSE))
    }
  })
  return(st)
}

# ---- clone ----

#' @name clone
#' @title Clone a simplex tree
#' @description Performs a deep-copy on the supplied simplicial complex.
#' @param st a simplex tree.
#' @return a new simplex tree copied from `st`.
#' @family serialization methods
#' @details A clone is produced by serializing the input simplex tree \code{st}
#'   and deserializing it into a new simplex tree that is then returned.
#' @export
clone <- function(st){
  stopifnot(class(st) %in% .st_classes)
  new_st <- deserialize(st %>% serialize())
  return(new_st)
}
