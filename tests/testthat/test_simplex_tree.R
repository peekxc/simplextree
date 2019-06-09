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
})


