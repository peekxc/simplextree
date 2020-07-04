## test_simplex_tree.R
library("simplextree")
library("testthat")
testthat::context("Testing Simplex Tree")

## ---- can_construct ----
test_that("Can construct and deconstruct a Simplex Tree object", {
  st <- simplex_tree()
  expect_is(st, "Rcpp_SimplexTree")
  expect_is(st$as_XPtr(), "externalptr")
  expect_equal(st$n_simplices, numeric(0))
  expect_silent({ rm(st); invisible(gc(verbose = FALSE)) })
})

## ---- can_remove ----
test_that("Can add and remove vertices", {
  st <- simplex_tree()
  invisible(sapply(seq(5L), function(i) insert(st, i)))
  expect_equal(head(st$n_simplices, 1), 5L)
  invisible(sapply(seq(5L), function(i) remove(st, i)))
  expect_equal(head(st$n_simplices, 1), numeric(0))
  rm(st)
})

## ---- can_add_remove ----
test_that("Can add and remove edges", {
  st <- simplex_tree()
  n_vertices <- sample(5L:25L, size = 1)
  edges <- t(combn(n_vertices, 2L))
  
  ## Insert vertices
  expect_silent(invisible(sapply(seq(n_vertices), function(v){
    st %>% insert(as.integer(v))
  })))
  ## Insert edges
  expect_silent(invisible(apply(edges, 1, function(e){
    st %>% insert(as.integer(e))
  })))
  expect_equal(st$n_simplices, c(n_vertices, choose(n_vertices, 2L)))
  
  ## Remove edges, check each time
  cc <- choose(n_vertices, 2L)
  expect_silent(invisible(apply(edges, 1, function(e){
    st %>% remove(as.integer(e))
  })))
  expect_equal(st$n_simplices, c(n_vertices, numeric(0L)))
  rm(st)
})

## ---- can_insert_bulk ----
test_that("Bulk matrix insertion and removal works", {
  st <- simplex_tree()
  testthat::expect_silent(st %>% insert(combn(10,4)))
  testthat::expect_equal(st$n_simplices, choose(10,1:4))
  testthat::expect_silent(st %>% remove(combn(10,4)))
  testthat::expect_equal(st$n_simplices, choose(10, seq(3)))
  testthat::expect_silent(st %>% remove(combn(10,3)))
  testthat::expect_equal(st$n_simplices, choose(10, seq(2)))
  testthat::expect_silent(st %>% remove(combn(10,2)))
  testthat::expect_equal(st$n_simplices, choose(10, seq(1)))
  testthat::expect_silent(st %>% remove(combn(10,1)))
  testthat::expect_equal(st$n_simplices, numeric(0L))
})

## ---- cofaces traversal ----
test_that("Cofaces traversal works", {
  st <- simplex_tree(1:5)
  full_test <- straverse(preorder(st), function(simplex){
    cofaces_sigma <- as.list(st %>% cofaces(simplex))
    testthat::expect_true(all(sapply(cofaces_sigma, function(coface){
      all(simplex %in% coface)
    })))
  })
  testthat::expect_true(all(full_test))
})

## ---- cousin map ----
test_that("Cousins works", {
  st <- simplex_tree(1:5)
  testthat::expect_equal(st$n_simplices, c(5,10,10,5,1))
  testthat::expect_silent(st %>% remove(5))
  testthat::expect_equal(st$n_simplices, c(4,6,4,1))
})

## ---- export types ----
test_that("Export types work", {
  st <- simplex_tree()
  n_vertices <- sample(2:25, size = 1)
  edges <- t(combn(n_vertices, 2L))
  invisible(apply(edges, 1, function(e){
    st %>% insert(as.integer(e))
  }))
  
  ## Test can export to adjacency matrix 
  expect_is(st$as_adjacency_matrix(), class = "matrix")
  am <- st$as_adjacency_matrix()
  expect_equal(nrow(am), n_vertices)
  expect_equal(sum(am == 1), choose(n_vertices, 2)*2) 
  
  ## Test can export to adjacency list 
  expect_is(st$as_adjacency_list(), class = "list")
  al <- st$as_adjacency_list()
  expect_equal(length(al), n_vertices)
  expect_true(all(sapply(names(al), function(key) length(al[[key]]) == n_vertices-1)))
  
  ## Test can export to an edgelist
  expect_is(st$as_edge_list(), class = "matrix")
  el <- st$as_edge_list()
  expect_equal(dim(el), c(choose(n_vertices, 2), 2L))
})

## ---- serialization ----
test_that("serialization/deserialization works", {
  st <- simplex_tree()
  testthat::expect_silent(st %>% insert(list(2:5, 7:12, 15:20)))
  testthat::expect_equal(st, deserialize(serialize(st)))
})

## ---- vertex collapse ----
## Example from simplicial maps paper (after converting letters to numbers)
test_that("vertex collapse works", {
  st <- simplex_tree()
  st %>% insert(list(c(1,5), c(4,5), 2:4))
  st %>% collapse(pair = list(1, 2), 1)

  testthat::expect_true(st %>% find(c(1, 3, 4)))
  testthat::expect_true(st %>% find(c(1, 5)))
  testthat::expect_true(st %>% find(c(4, 5)))
  testthat::expect_false(st %>% find(2))
  testthat::expect_equal(st$n_simplices, c(4, 5, 1))

  ## Test that {u,v} -> {w} is the same as {u,w} -> {w}, {v,w} -> {w}
  st1 <- simplex_tree()
  st1 %>% insert(list(1:2, 2:3, 3:5, c(3,5,6)))

  st2 <- simplex_tree()
  st2 %>% insert(list(1:2, 2:3, 3:5, c(3,5,6)))

  st1 %>% collapse(list(1, 5), 5) # {u,w} -> {w}
  st1 %>% collapse(list(6, 5), 5) # {v,w} -> {w}
  st2 %>% collapse(list(1, 6), 5) # {u,v} -> {w}
  all.equal(st1$as_list(), st2$as_list())
})

## ---- edge contractions ----
# testthat::test_that("edge contractions work", {
#   simplextree::traverse()
# })

## ---- reindexing ----
testthat::test_that("reindexing work", {
  st <- simplex_tree(1:5)
  st %>% reindex(6:10)
  expect_equal(st$vertices, 6:10)
})

## ---- clear ----
testthat::test_that("clear works", {
  st <- simplex_tree(combn(5,3))
  testthat::expect_true(all(st$n_simplices == choose(5,1:3)))
  testthat::expect_silent(st$clear())
  testthat::expect_equal(st$dimension, 0)
  testthat::expect_equal(st$n_simplices, numeric(0L))
})

## ---- preorder traversal ----
testthat::test_that("preorder traversal works", {
  st <- simplex_tree(list(2:5))
  testthat::expect_is(preorder(st), "st_traversal")
  testthat::expect_equal(capture.output(preorder(st)), "Preorder traversal @ { empty face }")
  testthat::expect_equal(
    capture.output(print_simplices(preorder(st), "short")), 
    "2, 2 3, 2 3 4, 2 3 4 5, 2 3 5, 2 4, 2 4 5, 2 5, 3, 3 4, 3 4 5, 3 5, 4, 4 5, 5"
  )
  testthat::expect_is(as.list(preorder(st)), "list")
})

## ---- level order traversal ----
testthat::test_that("level order traversal works", {
  st <- simplex_tree(list(2:5))
  testthat::expect_is(level_order(st), "st_traversal")
  testthat::expect_equal(capture.output(level_order(st)), "Level order traversal @ { empty face }")
  testthat::expect_equal(
    capture.output(print_simplices(level_order(st), "short")), 
    "2, 3, 4, 5, 2 3, 2 4, 2 5, 3 4, 3 5, 4 5, 2 3 4, 2 3 5, 2 4 5, 3 4 5, 2 3 4 5"
  )
  testthat::expect_is(as.list(level_order(st)), "list")
})


## ---- maximal traversal ----
testthat::test_that("maximal traversal works", {
  st <- simplex_tree(list(2:5, 8:12))
  testthat::expect_is(maximal(st), "st_traversal")
  testthat::expect_equal(capture.output(maximal(st)), "Maximal simplex traversal @ { empty face }")
  testthat::expect_equal(
    capture.output(print_simplices(maximal(st), "short")), 
    "2 3 4 5, 8 9 10 11 12"
  )
  testthat::expect_is(as.list(maximal(st)), "list")
})


## ---- link traversal ----
testthat::test_that("link traversal works", {
  st <- simplex_tree(list(c(1,2,3), c(1,3,4), c(1,4,5), c(1,5,6), c(1,6,7), c(1,2,7)))
  testthat::expect_is(link(st, 1), "st_traversal")
  testthat::expect_equal(capture.output(link(st, 1)), "Link traversal @ { 1 }")
  testthat::expect_equal(
    capture.output(print_simplices(link(st, 1), "short")), 
    "2, 2 3, 2 7, 3, 3 4, 4, 4 5, 5, 5 6, 6, 6 7, 7"
  )
  testthat::expect_is(as.list(link(st, 1)), "list")
})


## ---- face traversal ----
testthat::test_that("face traversal works", {
  st <- simplex_tree(list(1:5))
  testthat::expect_is(faces(st, 3:5), "st_traversal")
  testthat::expect_equal(capture.output(faces(st, 3:5)), "Face traversal @ { 3 4 5 }")
  testthat::expect_equal(
    capture.output(print_simplices(faces(st, 3:5), "short")), 
    "3, 4, 5, 3 4, 3 5, 4 5, 3 4 5"
  )
  testthat::expect_is(as.list(faces(st, 3:5)), "list")
})


## ---- k simplices traversal ----
testthat::test_that("k simplices traversal works", {
  st <- simplex_tree(list(1:5))
  testthat::expect_is(k_simplices(st, 2), "st_traversal")
  testthat::expect_equal(capture.output(k_simplices(st, 2)), "K-simplices traversal @ { empty face }")
  testthat::expect_equal(
    capture.output(print_simplices(k_simplices(st, 2), "short")), 
    "1 2 3, 1 2 4, 1 2 5, 1 3 4, 1 3 5, 1 4 5, 2 3 4, 2 3 5, 2 4 5, 3 4 5"
  )
  testthat::expect_is(as.list(k_simplices(st, 2)), "list")
})

## ---- k skeleton traversal ----
testthat::test_that("k skeleton traversal works", {
  st <- simplex_tree(list(1:5))
  testthat::expect_is(k_skeleton(st, 2), "st_traversal")
  testthat::expect_equal(capture.output(k_skeleton(st, 2)), "K-skeleton traversal @ { empty face }")
  testthat::expect_equal(
    capture.output(print_simplices(k_skeleton(st, 2), "short")), 
    "1, 1 2, 1 2 3, 1 2 4, 1 2 5, 1 3, 1 3 4, 1 3 5, 1 4, 1 4 5, 1 5, 2, 2 3, 2 3 4, 2 3 5, 2 4, 2 4 5, 2 5, 3, 3 4, 3 4 5, 3 5, 4, 4 5, 5"
  )
  testthat::expect_is(as.list(k_skeleton(st, 2)), "list")
})

## ---- traversal returns ----
testthat::test_that("traversal returns work", {
  st <- simplex_tree(list(1:4))
  testthat::expect_equal(straverse(k_simplices(st, 2), identity), combn(4,3))
  testthat::expect_equal(ltraverse(k_simplices(st, 2), identity), combn(4,3,simplify=FALSE))
})

## Testing k-expansion
testthat::test_that("K-expansion works", {
  base_complex <- function(){
    st <- simplex_tree()
    sigma1 <- apply(combn(3, 2), 2, function(idx){ (1:3)[idx] })
    sigma2 <- apply(combn(4, 2), 2, function(idx){ (4:7)[idx] })
    st %>% insert(cbind(sigma1, sigma2))
    return(st)
  }
  st <- base_complex()
  testthat::expect_true(methods::is(st, "Rcpp_SimplexTree"))
  testthat::expect_equal(st$dimension, 1)
  testthat::expect_silent(st %>% expand(2))
  testthat::expect_equal(st$dimension, 2)
  testthat::expect_true(all(st %>% find(list(1:3, 4:6, c(4, 5, 7), c(4, 6, 7)))))
  testthat::expect_false(st %>% find(4:7))
  st <- base_complex()
  testthat::expect_silent(st %>% expand(3))
  testthat::expect_equal(st$dimension, 3)
  testthat::expect_true(st %>% find(4:7))
})

# testthat::test_that("combinadics work", {
#   testthat::expect_error(nat_to_sub(0, 10, 10))
#   testthat::expect_equal(nat_to_sub(1, 10, 10), matrix(seq(10), ncol=1))
#   testthat::expect_equal(nat_to_sub(1, 10, 1), matrix(seq(10), ncol=1))
#   simplextree::traverse()
# })
# 
# testthat::test_that("flag construction works", {
#   simplextree::traverse()
# })
# 
# testthat::test_that("rips filtration works", {
#   simplextree::traverse()
# })
# 
# 
# testthat::test_that("intersection testing works", {
#   
# })
# testthat::test_that("nerve construction works", {
#   simplextree::traverse()
# })



# 
# 
# testthat::test_that("Facet traversal works", {
#   st <- simplex_tree()
#   st$insert(list(1:3, 2:5))
#   testthat::expect_equivalent(st$ltraverse(2:5, identity, "facets"), list(c(2,3,4), c(2,3,5), c(2,4,5), c(3,4,5)))
# })
# 
# testthat::test_that("Rips filtration works", {
#   
#   ## Helper function to turn column matrices into edge lists
#   matrix_to_list <- function(M){ split(as.vector(M), rep(1:ncol(M), each=nrow(M))) }
# 
#   ## Test can build rips filtration 
#   xy <- replicate(2, runif(20)) # bench at 20
#   X_dist <- as.matrix(dist(xy))
# 
#   st <- simplex_tree()
#   st$insert(matrix_to_list(combn(nrow(xy),2)))
#   st$rips(X_dist[st$edges], 2)
#   
#   st$threshold_function(0.1)
#   
#   ## Test the weights match exactly to what one might compute manually
#   simplex_weight <- function(sigma){
#     if (length(sigma) == 1){ return(0); }
#     idx <- matrix(sigma[combn(length(sigma),2)], nrow=2)
#     max(apply(idx, 2, function(i){ X_dist[i[1],i[2]] }))
#   }
#   test_weights <- sapply(st$rips_simplices[-1], simplex_weight)
#   all(abs(st$rips_weights[-1] - test_weights) < .Machine$double.eps)
#   
#   st <- simplex_tree()
#   st$insert(list(1:3, 2:5))
#   st$faces(c(2,4,5))
#   
#   st <- simplex_tree()
#   st$insert(1:5)
#   test_dist <- matrix(runif(5*5),nrow=st$n_simplices[1], ncol=st$n_simplices[1])
#   test_dist <- as.matrix(as.dist(test_dist))
#   st$filtration(test_dist[lower.tri(test_dist)], 2)
#   
#   ## Test can build rips filtration 
#   g <- igraph::sample_grg(10^4, radius = 0.01, coords = TRUE)
#   X_dist <- parallelDist::parDist(cbind(igraph::vertex_attr(g, "x"), igraph::vertex_attr(g,"y")))
# 
#   library(simplextree)
#   st <- simplex_tree()
#   edges <- matrix_to_list(t(igraph::as_edgelist(g)))
#   st$insert(edges)
#   w_edges <- as.matrix(X_dist)[st$edges]
#   microbenchmark::microbenchmark({ st$rips(w_edges, 2) }, times = 1L)
#   
#   # sapply(unique(sort(test_dist[lower.tri(test_dist)])), function(eps){ 
#   #   st$threshold_function(eps) 
#   #   st$connected_components
#   # })
#   
#   ## Test the thresholding function works
#   testthat::expect_silent(st$threshold_function(Inf))
#   testthat::expect_equal(st$n_simplices, as.integer(c(5, 10, 10, 5, 1)))
# 
#   testthat::expect_silent(st$threshold_function(0))
#   testthat::expect_equal(st$n_simplices, 5L)
#   
#   testthat::expect_silent(st$threshold_function(-1))
#   testthat::expect_equal(st$n_simplices, numeric(0L))
#   
#   testthat::expect_silent(st$threshold_function(Inf))
#   testthat::expect_equal(st$n_simplices, as.integer(c(5, 10, 10, 5, 1)))
#   
#   ## Check the number of simplices changes at each iteration
#   eps_vals <- unique(sort(test_dist[lower.tri(test_dist)]))
#   n_simplices_filt <- lapply(eps_vals, function(eps){ 
#     st$threshold_function(eps) 
#     st$n_simplices
#   })
#   testthat::expect_equal(n_simplices_filt, unique(n_simplices_filt))
#   
#   ## Test number of edges at threshold value is always 
#   num_edges <- sapply(eps_vals, function(eps){ 
#     st$threshold_function(eps) 
#     nrow(st$edges) == sum(eps_vals <= eps)
#   })
#   testthat::expect_true(all(num_edges))
# })

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
