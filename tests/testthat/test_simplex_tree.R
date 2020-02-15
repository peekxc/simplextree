## test_simplex_tree.R
library("simplextree")
library("testthat")
testthat::context("Testing Simplex Tree")

test_that("Can construct and deconstruct a Simplex Tree object", {
  st <- simplextree::simplex_tree()
  expect_is(st, "Rcpp_SimplexTree")
  expect_is(st$as_XPtr(), "externalptr")
  expect_equal(st$n_simplices, numeric(0))
  expect_silent({ rm(st); invisible(gc(verbose = FALSE)) })
})

test_that("Can add and remove vertices", {
  st <- simplextree::simplex_tree()
  invisible(sapply(seq(5L), function(i) st$insert(i)))
  expect_equal(head(st$n_simplices, 1), 5L)
  invisible(sapply(seq(5L), function(i) st$remove(i)))
  expect_equal(head(st$n_simplices, 1), numeric(0))
  rm(st)
})

test_that("Can add and remove edges", {
  stree <- simplextree::simplex_tree()
  n_vertices <- sample(5L:25L, size = 1)
  edges <- t(combn(n_vertices, 2L))
  
  ## Insert vertices
  expect_silent(invisible(sapply(seq(n_vertices), function(v){
    stree$insert(as.integer(v))
  })))
  ## Insert edges
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$insert(as.integer(e))
  })))
  expect_equal(stree$n_simplices, c(n_vertices, choose(n_vertices, 2L)))
  
  ## Remove edges, check each time
  cc <- choose(n_vertices, 2L)
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$remove(as.integer(e))
  })))
  expect_equal(stree$n_simplices, c(n_vertices, numeric(0L)))
  rm(stree)
})

test_that("Export types work", {
  stree <- simplextree::simplex_tree()
  n_vertices <- sample(2:25, size = 1)
  edges <- t(combn(n_vertices, 2L))
  invisible(apply(edges, 1, function(e){
    stree$insert(as.integer(e))
  }))
  
  ## Test can export to adjacency matrix 
  expect_is(stree$as_adjacency_matrix(), class = "matrix")
  am <- stree$as_adjacency_matrix()
  expect_equal(nrow(am), n_vertices)
  expect_equal(sum(am == 1), choose(n_vertices, 2)*2) 
  
  ## Test can export to adjacency list 
  expect_is(stree$as_adjacency_list(), class = "list")
  al <- stree$as_adjacency_list()
  expect_equal(length(al), n_vertices)
  expect_true(all(sapply(names(al), function(key) length(al[[key]]) == n_vertices-1)))
  
  ## Test can export to an edgelist
  expect_is(stree$as_edge_list(), class = "matrix")
  el <- stree$as_edge_list()
  expect_equal(dim(el), c(choose(n_vertices, 2), 2L))
})

test_that("serialization/deserialization works", {
  st <- simplex_tree()
  testthat::expect_silent(st$insert(list(1:3, 4:5, 6)))
  testthat::expect_equal(list(1:3, 4:5, 6), st$serialize())
  st2 <- simplex_tree()
  testthat::expect_silent(st2$deserialize(st$serialize()))
  testthat::expect_equal(st, st2)
})

## Example from simplicial maps paper (after converting letters to numbers)
test_that("vertex collapse works", {
  st <- simplex_tree()
  st$insert(list(c(1,5), c(4,5), 2:4))
  st$collapse(1, 2, 1)
  
  testthat::expect_true(st$find(c(1, 3, 4)))
  testthat::expect_true(st$find(c(1, 5)))
  testthat::expect_true(st$find(c(4, 5)))
  testthat::expect_false(st$find(2))
  testthat::expect_equal(st$n_simplices, c(4, 5, 1))
  
  ## Test that {u,v} -> {w} is the same as {u,w} -> {w}, {v,w} -> {w}
  st1 <- simplex_tree()
  st1$insert(list(1:2, 2:3, 3:5, c(3,5,6)))
  
  st2 <- simplex_tree()
  st2$insert(list(1:2, 2:3, 3:5, c(3,5,6)))
  
  st1$collapse(1, 5, 5) # {u,w} -> {w}
  st1$collapse(6, 5, 5) # {v,w} -> {w}
  st2$collapse(1, 6, 5) # {u,v} -> {w}
  all.equal(st1$as_list(), st2$as_list())
})

testthat::test_that("reindexing work", {
  st <- simplex_tree()
  st$insert(list(c(1,2), c(2,3), c(1,3), c(1,4), c(2,5), c(2,6), c(3,7), c(3,8), c(3,9)))
  st2 <- simplex_tree()
  st2$deserialize(st$serialize())
  st$reindex(c(3,1,2,4:9))  
  st2$reindex(list("1"=3, "2"=1, "3"=2))
  expect_equal(st$serialize(), st2$serialize())
})

## Testing k-expansion 
testthat::test_that("K-expansion works", {
  st$serialize()
  base_complex <- function(){
    st <- simplex_tree()
    sigma1 <- apply(combn(3, 2), 2, function(idx){ (1:3)[idx] })
    sigma2 <- apply(combn(4, 2), 2, function(idx){ (4:7)[idx] })
    apply(cbind(sigma1, sigma2), 2, st$insert) 
    return(st)
  }
  st <- base_complex()
  testthat::expect_true(methods::is(st, "Rcpp_SimplexTree"))
  testthat::expect_equal(st$dimension, 1)
  testthat::expect_silent(st$expand(2))
  testthat::expect_equal(st$dimension, 2)
  testthat::expect_true(all(st$find(list(1:3, 4:6, c(4, 5, 7), c(4, 6, 7)))))
  testthat::expect_false(st$find(4:7))
  testthat::expect_error(st$expand(3))
  st <- base_complex()
  testthat::expect_silent(st$expand(3))
  testthat::expect_equal(st$dimension, 3)
  testthat::expect_true(st$find(4:7))
  
  ## More tests 
  # st <- simplextree::simplex_tree()
  # st$insert(list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),c(3,5),c(4,5),c(5,6),c(5,7),c(5,8),c(6,7,8,9)))
  # st$expand(k=2)
  
})


testthat::test_that("Facet traversal works", {
  st <- simplex_tree()
  st$insert(list(1:3, 2:5))
  testthat::expect_equivalent(st$ltraverse(2:5, identity, "facets"), list(c(2,3,4), c(2,3,5), c(2,4,5), c(3,4,5)))
})

testthat::test_that("Rips filtration works", {
  
  ## Helper function to turn column matrices into edge lists
  matrix_to_list <- function(M){ split(as.vector(M), rep(1:ncol(M), each=nrow(M))) }

  ## Test can build rips filtration 
  xy <- replicate(2, runif(20)) # bench at 20
  X_dist <- as.matrix(dist(xy))

  st <- simplex_tree()
  st$insert(matrix_to_list(combn(nrow(xy),2)))
  st$rips(X_dist[st$edges], 2)
  
  st$threshold_function(0.1)
  
  ## Test the weights match exactly to what one might compute manually
  simplex_weight <- function(sigma){
    if (length(sigma) == 1){ return(0); }
    idx <- matrix(sigma[combn(length(sigma),2)], nrow=2)
    max(apply(idx, 2, function(i){ X_dist[i[1],i[2]] }))
  }
  test_weights <- sapply(st$rips_simplices[-1], simplex_weight)
  all(abs(st$rips_weights[-1] - test_weights) < .Machine$double.eps)
  
  st <- simplex_tree()
  st$insert(list(1:3, 2:5))
  st$faces(c(2,4,5))
  
  st <- simplex_tree()
  st$insert(1:5)
  test_dist <- matrix(runif(5*5),nrow=st$n_simplices[1], ncol=st$n_simplices[1])
  test_dist <- as.matrix(as.dist(test_dist))
  st$filtration(test_dist[lower.tri(test_dist)], 2)
  
  ## Test can build rips filtration 
  g <- igraph::sample_grg(10^4, radius = 0.01, coords = TRUE)
  X_dist <- parallelDist::parDist(cbind(igraph::vertex_attr(g, "x"), igraph::vertex_attr(g,"y")))

  library(simplextree)
  st <- simplex_tree()
  edges <- matrix_to_list(t(igraph::as_edgelist(g)))
  st$insert(edges)
  w_edges <- as.matrix(X_dist)[st$edges]
  microbenchmark::microbenchmark({ st$rips(w_edges, 2) }, times = 1L)
  
  # sapply(unique(sort(test_dist[lower.tri(test_dist)])), function(eps){ 
  #   st$threshold_function(eps) 
  #   st$connected_components
  # })
  
  ## Test the thresholding function works
  testthat::expect_silent(st$threshold_function(Inf))
  testthat::expect_equal(st$n_simplices, as.integer(c(5, 10, 10, 5, 1)))

  testthat::expect_silent(st$threshold_function(0))
  testthat::expect_equal(st$n_simplices, 5L)
  
  testthat::expect_silent(st$threshold_function(-1))
  testthat::expect_equal(st$n_simplices, numeric(0L))
  
  testthat::expect_silent(st$threshold_function(Inf))
  testthat::expect_equal(st$n_simplices, as.integer(c(5, 10, 10, 5, 1)))
  
  ## Check the number of simplices changes at each iteration
  eps_vals <- unique(sort(test_dist[lower.tri(test_dist)]))
  n_simplices_filt <- lapply(eps_vals, function(eps){ 
    st$threshold_function(eps) 
    st$n_simplices
  })
  testthat::expect_equal(n_simplices_filt, unique(n_simplices_filt))
  
  ## Test number of edges at threshold value is always 
  num_edges <- sapply(eps_vals, function(eps){ 
    st$threshold_function(eps) 
    nrow(st$edges) == sum(eps_vals <= eps)
  })
  testthat::expect_true(all(num_edges))
})

# g <- igraph::sample_grg(7, radius = 0.35, coords = TRUE)
#   X <- cbind(igraph::V(g)$x, igraph::V(g)$y)
#   
#   st <- simplextree::simplex_tree()
#   st$insert(as.list(igraph::V(g)))
#   st$insert(unname(as.list(as.data.frame(t(igraph::as_edgelist(g))))))
#   
#   mst <- dbscan:::prims(dist(X), n = nrow(X))
#   e_sorted <- apply(mst[,1:2], 1, sort)
#   
#   get_edge_weight <- function(edge){
#     e_idx <- which(apply(e_sorted, 2, function(sigma){ all(sigma == edge) }))
#     if (length(e_idx) == 1) { mst[e_idx,3] } else { return(Inf) }
#   }
#   X_dist <- as.matrix(dist(X))
#   get_edge_weight <- function(edge){
#     X_dist[edge[1], edge[2]]
#   }
# 
#   weighted_si <- st$ltraverse(empty_face, function(simplex){
#     if (length(simplex) == 0){ return(list(NULL, 0)) }
#     else if (length(simplex) == 1){ return(list(simplex, 0)) }
#     else if (length(simplex) == 2){ 
#       return(list(simplex, get_edge_weight(simplex)))
#     } else {
#       k <- length(simplex)
#       eps <- max(tapply(as.vector(simplex[combn(k, 2)]), rep(1:choose(k, 2), each=2), get_edge_weight))
#       return(list(simplex, eps))
#     }
#   }, "bfs")
  
  # max(mst[,3])+sqrt(.Machine$double.eps)
  # 
  # simplices <- lapply(weighted_si[-1], function(x) x[[1]])
  # weights <- sapply(weighted_si[-1], function(x) x[[2]])
  # 
  # s_simplices <- simplices[order(weights)]
  # s_weights <- weights[order(weights)]
  # 
  # rips_ref <- TDA::ripsFiltration(X, maxdimension = 5, maxscale=max(mst[,3])+sqrt(.Machine$double.eps))
  # 
  # library(simplextree)
  # matrix_to_list <- function(M){ split(as.vector(M), rep(1:ncol(M), each=nrow(M)))  }
  # 
  # eps <- quantile(X_dist, 0.60)
  # rips_ref <- TDA::ripsFiltration(X, maxdimension = 5, maxscale = eps)
  # 
  # valid_idx <- (lower.tri(X_dist) & (X_dist <= eps))
  # row_idx <- row(X_dist)[valid_idx]
  # col_idx <- col(X_dist)[valid_idx]
  # rips_test <- simplex_tree()
  # rips_test$insert(as.list(seq(nrow(X_dist))))
  # rips_test$insert(matrix_to_list(rbind(row_idx, col_idx)))
  # rips_test$filtration(X_dist[valid_idx], 5)
  # 
  # st$insert(list(1:2, 2:3, c(1, 3)))
  # st$expand(2)
  # .fix_order <- order(cbind(rips_ref$values, sapply(rips_ref$cmplx, length)))
  # 
  # s1 <- lapply(rips_ref[.fix_order]$cmplx, sort)
  # s2 <- lapply(rips_test$filtration_simplices[-1], sort)
  # all.equal(s1, s2)
  # 
  # xy <- cmdscale(X_dist)
  # plot(xy, asp = 1)
  # text(xy, labels = seq(nrow(X_dist)), pos=3)
  # 
  # length(wut$cmplx) == length(s_simplices)
  # all(abs(wut$values - s_weights) < sqrt(.Machine$double.eps))
  # all(wut$cmplx == s_simplices)
  # which(!mapply(function(u, v){ all(sort(u) == sort(v)) }, wut$cmplx, s_simplices))
  # wut$cmplx[60:65] == s_simplices[60:65]
  # 
