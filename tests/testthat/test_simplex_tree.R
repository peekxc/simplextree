## test_simplex_tree.R
library("simplextree")
library("testthat")
testthat::context("Testing Simplex Tree")

test_that("Can construct and deconstruct a Simplex Tree object", {
  stree <- simplextree::simplex_tree()
  expect_is(stree, "Rcpp_SimplexTree")
  expect_equal(stree$n_simplexes, numeric(0))
  expect_silent({ rm(stree); invisible(gc(verbose = FALSE)) })
})

test_that("Can add and remove vertices", {
  stree <- simplextree::simplex_tree()
  invisible(sapply(seq(5L), function(i) stree$insert_simplex(i)))
  expect_equal(head(stree$n_simplexes, 1), 5L)
  invisible(sapply(seq(5L), function(i) stree$remove_simplex(i)))
  expect_equal(head(stree$n_simplexes, 1), numeric(0))
  rm(stree)
})


test_that("Can add and remove edges", {
  stree <- simplextree::simplex_tree()
  n_vertices <- sample(1:25, size = 1)
  edges <- t(combn(n_vertices, 2L))
  
  ## Insert vertices
  expect_silent(invisible(sapply(seq(n_vertices), function(v){
    stree$insert_simplex(as.integer(v))
  })))
  ## Insert edges
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$insert_simplex(as.integer(e))
  })))
  expect_equal(stree$n_simplexes, c(n_vertices, choose(n_vertices, 2L)))
  
  ## Remove edges, check each time
  cc <- choose(n_vertices, 2L)
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$remove_simplex(as.integer(e))
  })))
  expect_equal(stree$n_simplexes, c(n_vertices, numeric(0L)))
  rm(stree)
})

test_that("Export types work", {
  stree <- simplextree::simplex_tree()
  n_vertices <- sample(2:25, size = 1)
  edges <- t(combn(n_vertices, 2L))
  invisible(apply(edges, 1, function(e){
    stree$insert_simplex(as.integer(e))
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






