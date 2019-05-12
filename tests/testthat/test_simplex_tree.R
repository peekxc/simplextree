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



## Example from simplicial maps paper (after converting letters to numbers)
test_that("vertex collapse works", {
  st <- simplex_tree()
  st$insert_simplex(c(1,5))
  st$insert_simplex(c(4,5))
  st$insert_simplex(c(2,3,4))
  st$collapse(1, 2, 1)
  
  testthat::expect_true(st$find_simplex(c(1, 3, 4)))
  testthat::expect_true(st$find_simplex(c(1, 5)))
  testthat::expect_true(st$find_simplex(c(4, 5)))
  testthat::expect_false(st$find_simplex(2))
  testthat::expect_equal(st$n_simplexes, c(4, 5, 1))
  
  ## Test that {u,v} -> {w} is the same as {u,w} -> {w}, {v,w} -> {w}
  
  st1 <- simplex_tree()
  st1$insert_simplex(1:2)
  st1$insert_simplex(2:3)
  st1$insert_simplex(3:5)
  st1$insert_simplex(c(3,5,6))
  
  st2 <- simplex_tree()
  st2$insert_simplex(1:2)
  st2$insert_simplex(2:3)
  st2$insert_simplex(3:5)
  st2$insert_simplex(c(3,5,6))
  
  st1$collapse(1, 5, 5) # {u,w} -> {w}
  st1$collapse(6, 5, 5) # {v,w} -> {w}
  st2$collapse(1, 6, 5) # {u,v} -> {w}
  all.equal(st1$as_list(), st2$as_list())
  
  
})


# st <- simplex_tree()
# st$insert_simplex(1:3)
# st$ltraverse(1, identity, "cofaces")
# st$ltraverse(2, identity, "cofaces")
# st$ltraverse(3, identity, "cofaces")


