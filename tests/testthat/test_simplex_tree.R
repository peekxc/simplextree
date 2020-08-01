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
test_that("collapse works", {
  st <- simplex_tree(list(c(1,5), c(4,5), 2:4))
  
  ## Testing elementary collapse 
  expect_silent(st %>% collapse(list(c(2,3), c(2,3,4))))
  expect_equal(
    capture.output(print_simplices(st)), 
    "1, 2, 3, 4, 5, 1 5, 2 4, 3 4, 4 5"
  )
  expect_false(st %>% find(c(2,3)))
  expect_false(st %>% find(c(2,3,4)))
  
  ## Testing vertex collapse 
  st <- simplex_tree(list(c(1,5), c(4,5), 2:4))
  expect_silent(st %>% collapse(pair = list(1, 2), 1))
  
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

# ---- is face ----
testthat::test_that("face testing works", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_true(st %>% is_face(1:2, 1:3))  ## test face relation 
  expect_false(st %>% is_face(1:3, 1:2)) ## test non-face relation 
  expect_false(st %>% is_face(1:2, 1:4)) ## test face relation + existence
})

# ---- adjacent ----
testthat::test_that("adjacent", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_equal(st %>% adjacent(1), c(2, 3))  
  expect_equal(st %>% adjacent(2), c(1, 3))
  expect_equal(st %>% adjacent(7), numeric(length = 0))
})

# ---- degree ----
testthat::test_that("degree", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_silent(st %>% degree(st$vertices))  
  expect_equal(st %>% degree(st$vertices), c(2, 2, 2, 1, 1, 0))
  expect_equal(st %>% degree(c(st$vertices, 7)), c(2, 2, 2, 1, 1, 0, 0))
})

# ---- generate ids  ----
testthat::test_that("generating ids", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_equal(st$generate_ids(10), c(0, 7:15))  
  expect_equal(st$generate_ids(10), c(0, 7:15))  
  expect_silent({ st$id_policy <- "unique" })
  expect_equal(st$id_policy, "unique")
})

# ---- clear  ----
testthat::test_that("clearing the complex", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_equal(st$n_simplices, c(6,4,1))
  expect_silent(st %>% clear())
  expect_equal(st$n_simplices, numeric(length = 0L))
})

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

testthat::test_that("combinadics work", {
  # testthat::expect_error(nat_to_sub(0, 10, 10))
  expect_equal(nat_to_sub(numeric(0L), 10, 2), matrix(integer(0), nrow = 2, ncol=0))
  expect_equal(nat_to_sub(1, 10, 10), matrix(seq(10), ncol=1))
  expect_equal(nat_to_sub(seq(10), 10, 1), matrix(seq(10), nrow=1))
  expect_equal(nat_to_sub(seq(choose(10,2)), n = 10, k = 2), combn(10, 2)) 
  expect_equal(nat_to_sub(seq(choose(10,3)), n = 10, k = 3), combn(10, 3)) 
  
  # testthat::expect_error(sub_to_nat(0, 10))
  expect_equal(sub_to_nat(numeric(0L), 10), integer(length = 0L))
  expect_equal(sub_to_nat(seq(10),10), 1)
  expect_equal(sub_to_nat(matrix(seq(10), ncol=1),10), 1)
  expect_equal(sub_to_nat(combn(10,2), n = 10), seq(choose(10,2)))
  expect_equal(sub_to_nat(combn(10,3), n = 10), seq(choose(10,3)))
  
})

testthat::test_that("inverse binomial coefficients work", {
  expect_equal(inverse.choose(choose(10,2), k = 2), 10)
  expect_equal(inverse.choose(choose(10,3), k = 3), 10)
  expect_equal(inverse.choose(choose(120, k = 10), k = 10), 120)
})

testthat::test_that("flag construction works", {
  xy <- cbind(runif(10), runif(10))
	d <- dist(xy, method = "euclidean")
	st <- simplex_tree(as.list(10)) %>% flag(d)
	testthat::expect_true("Rcpp_Filtration" %in% class(st))
	testthat::expect_equal(st$weights, numeric(0))
  testthat::expect_equal(st$simplices, list())
  testthat::expect_equal(st$included, logical(0))
})

testthat::test_that("rips filtration works", {
  xy <- cbind(runif(10), runif(10))
	d <- dist(xy, method = "euclidean")
  expect_silent(rips(d))
  expect_true("Rcpp_SimplexTree" %in% class(rips(d)))
  expect_true("Rcpp_SimplexTree" %in% class(rips(d, eps = 1.0)))
  expect_true("Rcpp_SimplexTree" %in% class(rips(d, dim = 5)))
  expect_true("Rcpp_Filtration" %in% class(rips(d, filtered = TRUE)))
  
  ## Testing rips complex
  R <- rips(d, dim = 3)
  eps <- enclosing_radius(d)
  max_weight <- function(simplex){ max(combn(simplex, 2, function(sigma){ d[sub_to_nat(sigma, 10)] })) }
  expect_true(all(straverse(k_simplices(R, 1), max_weight) <= eps))
  expect_true(all(straverse(k_simplices(R, 2), max_weight) <= eps))
  expect_true(all(straverse(k_simplices(R, 3), max_weight) <= eps))
  expect_equal(min(apply(as.matrix(d), 1, max)), eps)
  expect_true(all(combn(10, 2, function(e){ ifelse(d[sub_to_nat(e,10)] <= eps, R %>% find(e), !(R %>% find(e))) })))
  expect_true("Rcpp_SimplexTree" %in% class(R))
  
  ## Testing rips filtration
  R <- rips(d, dim = 2, filtered = TRUE)
  ns <- R$n_simplices
  expect_silent(R$weights)
  expect_silent(R$simplices)
  s <- R$simplices
  w <- R$weights
  expect_silent(R %>% threshold(index = 0))
  expect_equal(R$n_simplices, numeric(length = 0L))
  expect_silent(R %>% threshold(value = Inf))
  expect_equal(R$simplices, s)
  expect_equal(R$weights, w)
  expect_true(all(R$included))
  N <- length(R$included)
  for (i in sample(seq(N))[seq(10)]){
    expect_silent(R %>% threshold(index = i))
    expect_equal(sum(R$included), i)
    expect_equal(R$weights, w[seq(i)])
    expect_equal(R$simplices, s[seq(i)])
  }
})


testthat::test_that("nerve construction works", {
  alphabet <- seq(50)
  cover <- lapply(seq(15), function(i){
    set_size <- as.integer(runif(n = 1, min = 1, max = 15))
    sample(alphabet, size = set_size, replace = FALSE)
  })
  nfold <- function(nt){ function(cc){ length(Reduce(intersect, cover[cc])) >= nt } }
  
  ## Testing unsorted n-fold intersections
  for (nt in seq(1, 5, by = 1)){
    st <- simplex_tree() %>% nerve(cover, k = 2, threshold = nt)
    e2 <- nat_to_sub(which(combn(15, 2, nfold(nt))), n = 15, k = 2)
    e3 <- nat_to_sub(which(combn(15, 3, nfold(nt))), n = 15, k = 3)
    expect_equal(e2, t(st$edges))
    expect_equal(e3, t(st$triangles))
  }
  
  ## Testing sorted n-fold intersections
  cover <- lapply(cover, sort)
  for (nt in seq(1, 5, by = 1)){
    st <- nerve(simplex_tree(), cover, k = 2, threshold = nt)
    e2 <- nat_to_sub(which(combn(15, 2, nfold(nt))), n = 15, k = 2)
    e3 <- nat_to_sub(which(combn(15, 3, nfold(nt))), n = 15, k = 3)
    expect_equal(e2, t(st$edges))
    expect_equal(e3, t(st$triangles))
  }
  
  ## Testing using neighborhood function
  include_f <- function(ids){ return(all(ids %in% 1:3)) }
  st <- simplex_tree() %>% nerve(cover, k = 2, neighborhood = include_f)
  expect_equal(
    capture.output(print_simplices(st, "short")), 
    "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1 2, 1 3, 2 3, 1 2 3"
  )
  
  ## Test matrix input
  st <- simplex_tree() 
  expect_silent(st %>% nerve(cover, k = 2, threshold = 1, neighborhood = combn(length(cover),2)))
  e2 <- nat_to_sub(which(combn(15, 2, nfold(1))), n = 15, k = 2)
  expect_equal(e2, t(st$edges))
  
  ## Test list input
  cc <- combn(length(cover), 3)
  neighborhood_lst <- lapply(1:ncol(cc), function(i){ cc[,i] })
  expect_silent(st %>% nerve(cover, k = 3, threshold = 1, neighborhood = neighborhood_lst))
  e3 <- nat_to_sub(which(combn(15, 3, nfold(1))), n = 15, k = 3)
  expect_equal(e3, t(st$triangles))
})



## testing union find 
testthat::test_that("union find works", {
  expect_silent(simplextree:::union_find())
  ut <- simplextree:::union_find(10)
  expect_equal(as.character(class(ut)), "Rcpp_UnionFind")
  expect_equal(ut$size, 10)
  expect_equal(ut$connected_components(), seq(0, 9))

  expect_silent(ut$union(0, 1))
  expect_equal(length(unique(ut$connected_components())), 9)
  
  expect_silent(ut$union(1, 2))
  expect_equal(length(unique(ut$connected_components())), 8)
  
  expect_silent(ut$union_all(c(3,4,5,6)))
  expect_equal(length(unique(ut$connected_components())), 5)
  
  expect_silent(ut$add_sets(5))
  expect_equal(ut$size, 15)
  expect_equal(length(unique(ut$connected_components())), 10)
})

testthat::test_that("contraction works", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_equal(
    capture.output(st %>% print_simplices("short")), 
    "1, 2, 3, 4, 5, 6, 1 2, 1 3, 2 3, 4 5, 1 2 3"
  )
  expect_silent(st %>% contract(c(1, 3)))
  expect_equal(
    capture.output(st %>% print_simplices("short")), 
    "1, 2, 4, 5, 6, 1 2, 4 5"
  )
})

testthat::test_that("is tree works", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_false(st$is_tree())
  st <- simplex_tree(list(1:2, 2:3, 3:4, 4:5, c(4,6), c(1,8)))
  expect_true(st$is_tree())
})

testthat::test_that("clone works", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  expect_equal(st, clone(st))
})

testthat::test_that("plotting is silent", {
  st <- simplex_tree(list(1:3, 4:5, 6))
  tfn <- tempfile()
  png(tfn)
  expect_silent(plot(st))
  dev.off()
  file.remove(tfn)
})


